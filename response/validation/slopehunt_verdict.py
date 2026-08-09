#!/usr/bin/env python3
"""Slope-hunt collection harness — one command, pre-registered decision rule.

WHAT THIS IS FOR
================
The T7 slope hunt (`response/avalanche/run_slopehunt_chain2.sh`, 144 slices:
a Penning rP scan plus a 90/10 isobutane confirmation point) is the last
outstanding item of the 2026-08-09→10 overnight queue. It runs unattended on the
desktop, self-merges to `aval_calib_slopehunt.json` and shuttles to EOS. This
turns its output into the verdict table in one command.

It gates the amplitude thread. From `OVERNIGHT_2026-08-10.md` sections 7 and 12,
the amplitude ledger now has exactly ONE surviving candidate — the avalanche
gain — because every other row is dead, controlled, separated, or (primary
ionisation, checked 2026-08-10) measured and correct. Two independent numbers
point at it:

    demanded gain  ~4e4 against the sim's 24 094 at 490 V   (x1.6)
    HV slope       data 0.449/10 V vs sim 0.296/10 V        (x1.52, ~12 sigma)

A single alpha(E) / Penning-transfer error at the operating point would produce
both, and neither number was derived from the other.

THE DECISION RULE — PRE-REGISTERED, 2026-08-10, BEFORE THE DATA EXISTS
======================================================================
Let rP* be the Penning transfer that best matches the data's gain-vs-HV slope.

  A. rP* fixes the slope AND moves gain at 490 V toward ~4e4 (>= x1.4 of the
     rP=0.40 baseline)
     -> THE LEDGER CLOSES ON A SINGLE DEFECT. The amplitude and slope threads
        are one alpha(E) error. Adopt rP* as the calibration, re-run Stage B/C,
        and the amplitude discrepancy is explained.

  B. rP* fixes the slope but leaves gain essentially unchanged (< x1.2)
     -> NO CANDIDATE SURVIVES. Every ledger row is then eliminated and the
        chain decomposition itself is suspect. This is the more interesting
        outcome and the one to be ready for. Do NOT paper over it by fitting a
        gain factor -- that would be exactly the modelling error the S3 ion
        closeout warned about.

  C. No rP in the scanned range reaches the data slope
     -> Penning is not the knob. The alpha(E) discrepancy is in the cross
        sections or the field map, not the transfer rate, and the gain campaign
        needs a different axis.

  D. rP* fixes the slope and OVERSHOOTS gain (> x2.2)
     -> the two threads are inconsistent with a single defect; report as such
        rather than choosing the rP that splits the difference.

Reporting rule that applies in every case: rP* is FITTED TO THE SLOPE. Once
fitted, slope agreement is not evidence. The gain is the out-of-sample
prediction, and it is the only part of this that can confirm anything. Same
discipline as beta and undershoot in the T14 freeze queue.

USAGE
=====
    python3 -m response.validation.slopehunt_verdict \
        --calib ~/x17/response_sim/avalanche/aval_calib_slopehunt.json \
        --slopes ~/x17/response_sim/hv_slope/slopes.json

Tested against `aval_calib_meshfield_hvscan.json` (which has no rP axis) so the
code path is proven before the real input exists:

    python3 -m response.validation.slopehunt_verdict --self-test
"""
from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict

import numpy as np

V_OP = 490.0
GAIN_SIM_REF = 24094.0        # T7 pooled meshfield 490 V
GAIN_DEMANDED = 4.0e4         # amplitude ledger, saturation-corrected X
DATA_SLOPE_KEY = ("x", "p50_head")


def _fit_slope10(volts, gains):
    """d ln(gain) / dV, expressed per 10 V, by least squares."""
    v = np.asarray(volts, float)
    g = np.asarray(gains, float)
    m = (g > 0) & np.isfinite(g)
    if m.sum() < 3:
        return None, None
    p, cov = np.polyfit(v[m], np.log(g[m]), 1, cov=True)
    return float(p[0] * 10.0), float(np.sqrt(cov[0, 0]) * 10.0)


def _rp_of(key, point):
    """Penning rP for a calib point, from the point or its key. None if absent."""
    for f in ("penning_rP", "penning_rp", "rP", "rp"):
        if f in point:
            return float(point[f])
    pen = point.get("penning")
    if isinstance(pen, dict) and "rP" in pen:
        return float(pen["rP"])
    m = re.search(r"rP0*(\d+)", key)
    if m:
        return float(m.group(1)) / (100.0 if len(m.group(1)) > 1 else 10.0)
    return None


def load(calib_path):
    """-> {(gas, rP): {voltage: gain_mean}}"""
    d = json.load(open(os.path.expanduser(calib_path)))
    pts = d["points"]
    out = defaultdict(dict)
    for key, p in pts.items():
        g = p.get("polya", {}).get("gain_mean")
        v = p.get("voltage_V")
        if g is None or v is None:
            continue
        gas = p.get("gas_file", key.split("@")[0])
        out[(gas, _rp_of(key, p))][float(v)] = float(g)
    return out


def data_slope(slopes_path):
    s = json.load(open(os.path.expanduser(slopes_path)))["data"]
    view, est = DATA_SLOPE_KEY
    r = s[view][est]
    return r["slope10"], r["err10"]


def verdict(arms, d_slope, d_err, baseline_rp=0.40):
    rows = []
    for (gas, rp), gv in sorted(arms.items(), key=lambda kv: (kv[0][0], kv[0][1] if kv[0][1] is not None else -1)):
        sl, err = _fit_slope10(list(gv), list(gv.values()))
        g_op = gv.get(V_OP)
        rows.append(dict(gas=gas, rP=rp, n=len(gv), slope10=sl, slope_err=err,
                         gain_at_490=g_op))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib", default="~/x17/response_sim/avalanche/"
                                       "aval_calib_slopehunt.json")
    ap.add_argument("--slopes", default="~/x17/response_sim/hv_slope/slopes.json")
    ap.add_argument("--self-test", action="store_true",
                    help="run against the 08-08 hvscan calib to prove the path")
    a = ap.parse_args()

    calib = a.calib
    if a.self_test:
        calib = "~/x17/response_sim/avalanche/aval_calib_meshfield_hvscan.json"
        print("SELF-TEST against the 08-08 HV scan (no rP axis; the rP column "
              "will read None and the decision rule will abstain).\n")

    arms = load(calib)
    d_sl, d_err = data_slope(a.slopes)
    rows = verdict(arms, d_sl, d_err)

    print(f"data slope ({DATA_SLOPE_KEY[0]}/{DATA_SLOPE_KEY[1]}): "
          f"{d_sl:.4f} +- {d_err:.4f} per 10 V")
    print(f"sim reference gain at {V_OP:.0f} V: {GAIN_SIM_REF:.0f};  "
          f"ledger demands ~{GAIN_DEMANDED:.0f}\n")
    print(f"{'gas':>34} {'rP':>6} {'n':>3} {'slope10':>9} {'+-':>7} "
          f"{'gain@490':>10} {'g/ref':>7} {'slope/data':>11}")
    for r in rows:
        rp = "  n/a" if r["rP"] is None else f"{r['rP']:5.2f}"
        sl = "      n/a" if r["slope10"] is None else f"{r['slope10']:9.4f}"
        se = "    n/a" if r["slope_err"] is None else f"{r['slope_err']:7.4f}"
        go = "       n/a" if r["gain_at_490"] is None else f"{r['gain_at_490']:10.0f}"
        gr = ("    n/a" if r["gain_at_490"] is None
              else f"{r['gain_at_490']/GAIN_SIM_REF:7.3f}")
        sr = ("        n/a" if r["slope10"] is None
              else f"{r['slope10']/d_sl:11.3f}")
        print(f"{os.path.basename(r['gas']):>34} {rp} {r['n']:>3} {sl} {se} "
              f"{go} {gr} {sr}")

    # ---- apply the pre-registered decision rule -----------------------------
    cand = [r for r in rows if r["rP"] is not None and r["slope10"] is not None]
    print("\n=== PRE-REGISTERED DECISION RULE ===")
    if not cand:
        print("  ABSTAIN — no rP axis in this calib. (Expected for the "
              "self-test; if you see this on the real slope-hunt output, the "
              "merge did not carry the Penning setting and that is a bug.)")
        return 0

    best = min(cand, key=lambda r: abs(r["slope10"] - d_sl))
    floor = d_sl - 2.0 * d_err
    reached = max(r["slope10"] for r in cand) >= floor
    print(f"  best-slope arm: rP = {best['rP']:.2f}, slope "
          f"{best['slope10']:.4f} vs data {d_sl:.4f} "
          f"({(best['slope10'] - d_sl) / d_err:+.1f} sigma)")

    if not reached:
        print("  -> OUTCOME C: no rP in the scanned range reaches the data "
              "slope. Penning is not the knob; the alpha(E) discrepancy is in "
              "the cross sections or the field map. The gain campaign needs a "
              "different axis.")
        return 0

    g = best["gain_at_490"]
    if g is None:
        print("  -> INCOMPLETE: the best-slope arm has no 490 V point, so the "
              "out-of-sample gain prediction cannot be read. Re-run that arm "
              "at 490 V before concluding anything.")
        return 0

    ratio = g / GAIN_SIM_REF
    print(f"  out-of-sample gain prediction: {g:.0f} at {V_OP:.0f} V "
          f"= x{ratio:.2f} of the rP-0.40 reference "
          f"(ledger demands ~x{GAIN_DEMANDED / GAIN_SIM_REF:.2f})")

    if ratio >= 1.4:
        if ratio > 2.2:
            print("  -> OUTCOME D: slope fixed but gain OVERSHOOTS. The two "
                  "threads are not consistent with a single defect. Report "
                  "that; do NOT pick the rP that splits the difference.")
        else:
            print("  -> OUTCOME A: THE LEDGER CLOSES ON A SINGLE DEFECT. The "
                  "amplitude and slope threads are one alpha(E) error. Adopt "
                  f"rP = {best['rP']:.2f}, re-run Stage B/C, and the amplitude "
                  "discrepancy is explained.")
    elif ratio < 1.2:
        print("  -> OUTCOME B: slope fixed at essentially UNCHANGED GAIN. No "
              "ledger candidate survives and the chain decomposition itself is "
              "suspect. Do NOT paper over it with a fitted gain factor -- that "
              "is precisely the modelling error the S3 ion closeout warned "
              "against. This is the interesting outcome.")
    else:
        print(f"  -> BETWEEN A AND B (x{ratio:.2f} sits in the 1.2-1.4 dead "
              "band the rule left deliberately unassigned). Partial closure; "
              "quote it as such and do not round it into either outcome.")

    print("\n  REPORTING RULE: rP was FITTED to the slope, so slope agreement "
          "is not evidence. The gain above is the out-of-sample prediction and "
          "is the only part of this that can confirm anything.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
