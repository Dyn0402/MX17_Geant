#!/usr/bin/env python3
"""
manifest.py — the production run register for S1 kernel products.

design/RESPONSE_SIM_PLAN.md §2: "No un-manifested runs — one CSV manifest row
per production run." Every greens_comb_*.npz already embeds its own provenance
(parameters, git hash of the producing tree, and the sum-rule result) in its
`meta` block, so the manifest is not a separate record that can disagree with
the data — it is read back out of the products themselves.

Reading `meta` is cheap: np.load on an npz is lazy, so this touches only the
small member and never decompresses the multi-GB arrays.

    python3 -m response.solver.manifest ~/x17/response_sim/s1
    python3 -m response.solver.manifest ~/x17/response_sim/s1 --check
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np

# No transport timescale here on purpose. The products' meta carries a
# `sharing.*.tau_1e_s`, but that estimator is window-dependent and withdrawn
# (see response/solver/spread.py); a register people trust must not repeat it.
# Transport numbers come from `python3 -m response.solver.spread`.
COLUMNS = ["file", "rho_s_MOhm_sq", "d_kapton_um", "nx", "ny",
           "sum_rule_err", "tiling_err", "channel_capture_prompt",
           "x_fraction_prompt", "git", "bytes"]

# The closed-form sum rule is the acceptance gate (plan §3). A product that
# fails it must not reach the digitizer.
SUM_RULE_TOL = 2e-2


def row(path):
    with np.load(path) as d:
        m = json.loads(str(d["meta"]))
    sr = m.get("sum_rule") or {}
    return {
        "file": os.path.basename(path),
        "rho_s_MOhm_sq": m["rho_s_ohm_sq"] / 1e6,
        "d_kapton_um": m["d_kapton_m"] * 1e6,
        "nx": m.get("nx"), "ny": m.get("ny"),
        "sum_rule_err": sr.get("err"),
        "tiling_err": sr.get("tiling_err"),
        "channel_capture_prompt": m.get("channel_capture_prompt"),
        "x_fraction_prompt": m.get("x_fraction_prompt"),
        "git": m.get("git", "")[:12],
        "bytes": os.path.getsize(path),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("--out", default=None, help="default <indir>/MANIFEST.csv")
    ap.add_argument("--check", action="store_true",
                    help="exit nonzero if any product fails the sum rule or "
                         "if the rho_s x d_k grid has holes")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.indir, "greens_comb_*.npz")))
    if not files:
        print(f"no greens_comb_*.npz under {a.indir}")
        return 1

    rows = []
    for f in files:
        try:
            rows.append(row(f))
        except Exception as e:
            print(f"  UNREADABLE {os.path.basename(f)}: {e}")
            if a.check:
                return 2
    rows.sort(key=lambda r: (r["d_kapton_um"], r["rho_s_MOhm_sq"]))

    out = a.out or os.path.join(a.indir, "MANIFEST.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} products in {a.indir}\n")
    print(f"{'file':34s} {'rho':>5s} {'d_k':>5s} {'sumrule':>9s} "
          f"{'tiling':>8s} {'capture':>8s} {'X frac':>7s}  git")
    for r in rows:
        print(f"{r['file']:34s} {r['rho_s_MOhm_sq']:5.1f} "
              f"{r['d_kapton_um']:5.0f} {r['sum_rule_err']:9.1e} "
              f"{r['tiling_err']:8.0e} {r['channel_capture_prompt']:8.4f} "
              f"{r['x_fraction_prompt']:7.4f}  {r['git']}")
    print(f"\nwrote {out}  ({sum(r['bytes'] for r in rows)/1e9:.1f} GB total)")

    if not a.check:
        return 0

    bad = [r for r in rows
           if r["sum_rule_err"] is None or r["sum_rule_err"] > SUM_RULE_TOL]
    # Round before comparing: d_kapton round-trips through metres, so 50 µm
    # comes back as 50.000000000000007 and exact-equality reports every
    # existing point as a hole.
    from ..common import constants as C
    key = lambda r, d: (round(float(r), 6), round(float(d), 6))
    have = {key(r["rho_s_MOhm_sq"], r["d_kapton_um"]) for r in rows}
    want = {key(r / 1e6, d) for r in C.RHO_S_SCAN_OHM_SQ
            for d in C.KAPTON_THICK_UM_SCAN}
    holes = sorted(want - have)

    if bad:
        print(f"\nFAIL: {len(bad)} product(s) miss the sum rule "
              f"(tol {SUM_RULE_TOL:.0e}): "
              + ", ".join(r["file"] for r in bad))
    if holes:
        print(f"\nincomplete grid, {len(holes)} point(s) missing: "
              + ", ".join(f"rho={r:g}M/d_k={d:g}um" for r, d in holes))
    if not bad and not holes:
        print("\ngrid COMPLETE and every product passes the sum rule.")
        return 0
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
