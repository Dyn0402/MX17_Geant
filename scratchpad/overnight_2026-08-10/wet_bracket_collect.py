#!/usr/bin/env python3
"""Collect the wet amp-range gain bracket and score it against its
pre-registration.

Pre-registration: design/report/WET_GAIN_BRACKET_PREREG_2026-08-10.md
(committed 287612b, BEFORE condor cluster 16705137 was submitted).

The four predictions and their falsifiers, restated here verbatim in substance
so the scoring cannot drift from what was actually claimed:

  1. Direction: gain falls monotonically with water, A > B > C.
     FALSIFIER: gain rises at either water point.
  2. Magnitude: 1 % H2O moves mean gain by 10-40 %.
     FALSIFIER: outside that band. (This prediction CONTRADICTS the roadmap's
     hoped-for "less than a few %", i.e. I predicted the contaminant axis has
     to STAY in the gain campaign.)
  3. Attachment: survival P(g>0) stays >= 0.95, unchanged within errors.
     FALSIFIER: survival below 0.90 at 1 % H2O.
  4. The 0.5 % point sits between A and C, roughly half way in log-gain.
     FALSIFIER: non-monotonic, or B outside [C, A].

Arms A/B/C are all at manual rP = 0.40 and are the bracket. Arm D (dry at
Penning auto) measures the Penning-MODEL delta and must never be differenced
against a wet arm -- that was the confound this run was designed to avoid.
"""
from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

SRC = os.path.expanduser("~/x17/response_sim/avalanche/wetbracket")
ARMS = {
    "A_dry_rp040": "dry 95/5, rP 0.40",
    "B_h2o05_rp040": "+0.5 % H2O, rP 0.40",
    "C_h2o10_rp040": "+1.0 % H2O, rP 0.40",
    "D_dry_auto": "dry 95/5, Penning auto",
}
BRACKET = ["A_dry_rp040", "B_h2o05_rp040", "C_h2o10_rp040"]


def load():
    acc = defaultdict(lambda: dict(gains=[], attached=0, attempted=0,
                                   batches=0, runtime=0.0, field=None))
    for p in sorted(glob.glob(os.path.join(SRC, "*.json"))):
        d = json.load(open(p))
        a = acc[d["gas_label"]]
        a["gains"].extend(d["gain_raw"])
        a["attached"] += d.get("n_attached", 0)
        a["attempted"] += d["n_events_attempted"]
        a["batches"] += 1
        a["runtime"] += d.get("runtime_s", 0.0)
        a["field"] = d.get("field")
        a["penning"] = (d.get("penning_mode"), d.get("penning_rP"))
    return acc


def stats(a, nboot=4000, seed=20260810):
    g = np.asarray(a["gains"], float)
    surv = float((g > 0).mean())
    gm = g[g > 0]
    rng = np.random.default_rng(seed)
    bs = np.array([gm[rng.integers(0, gm.size, gm.size)].mean()
                   for _ in range(nboot)])
    return dict(n=g.size, survival=surv,
                surv_err=float(np.sqrt(surv * (1 - surv) / g.size)),
                gain_mean=float(gm.mean()), gain_err=float(bs.std()),
                gain_median=float(np.median(gm)),
                rel_var=float(gm.var() / gm.mean() ** 2))


def main():
    acc = load()
    st = {k: stats(v) for k, v in acc.items()}

    print("=== wet amp-range gain bracket, condor 16705137 ===")
    print(f"{'arm':>16} {'batches':>8} {'events':>7} {'gain_mean':>11} "
          f"{'+-':>8} {'median':>9} {'survival':>9} {'attach':>7} {'rel_var':>8}")
    for k in ARMS:
        if k not in st:
            print(f"{k:>16}   MISSING")
            continue
        s, a = st[k], acc[k]
        print(f"{k:>16} {a['batches']:>8d} {s['n']:>7d} {s['gain_mean']:>11.0f} "
              f"{s['gain_err']:>8.0f} {s['gain_median']:>9.0f} "
              f"{s['survival']:>9.4f} {a['attached']:>7d} {s['rel_var']:>8.4f}")
    f = acc[BRACKET[0]]["field"]
    print(f"\n  field {f:.0f} V/cm at 490 V over a 0.015 cm gap; "
          f"total CPU {sum(v['runtime'] for v in acc.values())/3600:.1f} h")

    A, B, C = (st[k] for k in BRACKET)
    print("\n=== the bracket (rP 0.40 arms only) ===")
    for k, s in zip(BRACKET, (A, B, C)):
        rel = s["gain_mean"] / A["gain_mean"]
        nsig = ((s["gain_mean"] - A["gain_mean"])
                / np.hypot(s["gain_err"], A["gain_err"]))
        print(f"  {ARMS[k]:<24} gain {s['gain_mean']:8.0f} +- {s['gain_err']:5.0f}"
              f"   vs dry x{rel:.4f}  ({nsig:+.1f} sigma)")

    print("\n=== Penning-model delta (arm D, kept SEPARATE by design) ===")
    D = st["D_dry_auto"]
    print(f"  dry auto vs dry rP 0.40: {D['gain_mean']:.0f} vs "
          f"{A['gain_mean']:.0f}  = x{D['gain_mean']/A['gain_mean']:.4f}")
    print("  (never differenced against a wet arm -- that is the confound this "
          "run exists to avoid)")

    # ---- score the pre-registration ---------------------------------------
    print("\n=== SCORING AGAINST THE PRE-REGISTRATION ===")
    rB, rC = B["gain_mean"] / A["gain_mean"], C["gain_mean"] / A["gain_mean"]
    dC = abs(1.0 - rC) * 100.0

    p1 = A["gain_mean"] > B["gain_mean"] > C["gain_mean"]
    print(f"  1 direction   gain falls monotonically A > B > C : "
          f"{'CONFIRMED' if p1 else 'FALSIFIED'}  "
          f"({A['gain_mean']:.0f} / {B['gain_mean']:.0f} / {C['gain_mean']:.0f})")

    p2 = 10.0 <= dC <= 40.0
    print(f"  2 magnitude   1 % H2O moves gain 10-40 %          : "
          f"{'CONFIRMED' if p2 else 'FALSIFIED'}  (measured {dC:.1f} %)")

    p3 = C["survival"] >= 0.90
    print(f"  3 attachment  survival >= 0.90 at 1 % H2O         : "
          f"{'CONFIRMED' if p3 else 'FALSIFIED'}  ({C['survival']:.4f})")

    p4 = min(rC, 1.0) <= rB <= max(rC, 1.0)
    print(f"  4 ordering    B between A and C                   : "
          f"{'CONFIRMED' if p4 else 'FALSIFIED'}  "
          f"(B/A {rB:.4f}, C/A {rC:.4f})")

    print("\n=== THE ROADMAP DECISION (§2 step 2) ===")
    if dC < 5.0:
        print(f"  1 % H2O moves gain {dC:.1f} % -- less than a few percent.")
        print("  -> DROP the contaminant axis from the gain campaign. Keep "
              "contaminants only where they are known to matter (drift "
              "transport).")
    else:
        print(f"  1 % H2O moves gain {dC:.1f} % -- NOT negligible.")
        print("  -> The contaminant axis STAYS in the gain campaign, and the "
              "T7 slope hunt's assumption that composition is fixed while "
              "Penning varies needs its own error bar.")

    print("\n  Labelling, unchanged: every water fraction here is "
          "FITTED-TO-DATA. This measures a derivative (would water, if "
          "present, move gain?), not an operating point, and is not evidence "
          "that the bench gas contained water.")
    print("  Scope: one voltage. It says nothing about the gain SLOPE, which "
          "is the live problem. If the effect is large the follow-up is a "
          "two-voltage version, not a gain correction.")

    out = dict(cluster=16705137, field_Vcm=f,
               arms={k: dict(**st[k], label=ARMS[k],
                             penning=list(acc[k]["penning"]),
                             batches=acc[k]["batches"]) for k in st},
               ratios=dict(B_over_A=rB, C_over_A=rC,
                           D_over_A=D["gain_mean"] / A["gain_mean"]),
               prereg=dict(direction=bool(p1), magnitude=bool(p2),
                           attachment=bool(p3), ordering=bool(p4),
                           pct_change_1pct_H2O=dC))
    dst = os.path.expanduser("~/CLionProjects/MX17_Geant/design/report/"
                             "wet_gain_bracket_2026-08-10.json")
    with open(dst, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {dst}")


if __name__ == "__main__":
    main()
