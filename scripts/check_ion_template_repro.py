#!/usr/bin/env python3
"""
check_ion_template_repro.py — compare two S3 calibrations' ion templates.

    python3 scripts/check_ion_template_repro.py <ref_calib.json> <new_calib.json> [voltage]

Both arguments are calibration files as written by
`python3 -m response.avalanche.collect` (or the shipped `aval_calib_v*.json`),
NOT raw slice files. Voltage defaults to 490.

WHY THIS EXISTS. The first S3 campaign wrote `i_elec`/`i_ion` arrays that were
identically zero in all 56 slices: `ComponentConstant` has no weighting
potential until `SetWeightingPotential()` is called, and Garfield computes the
induced current from psi rather than the weighting field by default, so every
signal came out zero while the avalanche and ion drift ran perfectly normally —
gain and sigma0 unaffected, which is why it went unnoticed. Commit 42390d1
fixed it. EOS `avalanche/raw/` is still that PRE-FIX campaign; the raw that
actually backs `aval_calib_v2.json` is on AFS in `results_v2/`.

WHAT IT CHECKS. f_ion = Qi/(Qe+Qi), the fraction of induced charge carried by
the ion tail. It is the right observable because it is a RATIO, so it is
insensitive to gain and to slice count — in the 2026-08-08 verification an
independent re-run matched v2's f_ion to 1e-4 while its gain differed by 8.5 %
on a quarter of the statistics. Comparing gains would have shown a "failure"
that was only counting statistics.

It also flags an all-zero template outright, which is the failure mode that
started all of this and which a f_ion comparison alone would report as nan.
"""

from __future__ import annotations

import json
import sys

import numpy as np

TOL = 0.005          # absolute on f_ion


def _point(calib, voltage):
    pts = calib["points"]
    hits = [k for k in pts if f"@{voltage}V" in k]
    if not hits:
        raise SystemExit(f"no {voltage} V point; have: {sorted(pts)}")
    return hits[0], pts[hits[0]]


def _f_ion(p):
    ie = np.asarray(p["i_elec"], dtype=float)
    ii = np.asarray(p["i_ion"], dtype=float)
    dt = float(p["signal_dt_ns"])
    qe, qi = ie.sum() * dt, ii.sum() * dt
    nz = int(np.count_nonzero(ie)) + int(np.count_nonzero(ii))
    tot = qe + qi
    return (qi / tot if tot else float("nan")), qe, qi, nz


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    voltage = argv[2] if len(argv) > 2 else "490"
    rows, zero = [], False
    for label, path in (("reference", argv[0]), ("new", argv[1])):
        key, p = _point(json.load(open(path)), voltage)
        f, qe, qi, nz = _f_ion(p)
        gain = (p.get("polya") or {}).get("gain_mean", float("nan"))
        rows.append((label, path.split("/")[-1], p.get("nev_total"), gain, f))
        if nz == 0:
            print(f"  !! {path}: template is ALL ZERO — this is the "
                  f"SetWeightingPotential bug, not a small disagreement")
            zero = True

    print(f"\n  {'':10s} {'file':28s} {'nev':>6} {'gain':>10} {'f_ion':>10}")
    for label, name, nev, gain, f in rows:
        print(f"  {label:10s} {name:28s} {nev if nev else 0:6d} "
              f"{gain:10.1f} {f:10.6f}")

    d = rows[1][4] - rows[0][4]
    print(f"\n  f_ion difference {d:+.6f}   tol {TOL}")
    print(f"  (gain differs by {100*(rows[1][3]/rows[0][3]-1):+.1f} %, which is "
          f"why the check is on the ratio)")
    ok = (not zero) and abs(d) <= TOL
    print("\n" + ("PASS — the ion template reproduces."
                  if ok else
                  "FAIL — the templates disagree; do not treat them as "
                  "interchangeable."))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
