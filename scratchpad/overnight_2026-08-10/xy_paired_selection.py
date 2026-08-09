#!/usr/bin/env python3
"""Does the sim's Y/X amplitude asymmetry survive matched selection?

T14 reads sim Y/X = 0.790 on median peak amplitude against an unselected
detector value of 1.0002, and the strip-anisotropy A/B has already exonerated
the S1 electrostatics (xy_anisotropy_ab.py). But before hunting Stage B/C code,
the selection has to be ruled out: the two views' legs are selected
INDEPENDENTLY (each on its own view's reco quality and theta), and the data's Y
leg is cut ~4x harder than its X leg (detector saturation 0.327 -> 0.110 in Y
vs 0.326 -> 0.260 in X).

Three progressively tighter comparisons, all read-only on the frozen parquets:

  1. as-published -- each view's own leg, the T14 numbers
  2. PAIRED       -- only events present in BOTH views' legs, so the two views
                     see the identical event set and per-view selection cannot
                     contribute at all
  3. paired + unsaturated -- additionally drop events where either view is at
                     or near the 3550 ADC rail, which removes the clipping that
                     biases the data's medians downward

If Y/X survives (2) and (3) on the sim while the data stays near 1, the
asymmetry is genuinely sim-side and downstream of the kernel, and the next
place to look is per-view Stage B/C handling. If it collapses, it was selection
all along.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"
POINT = "t14_compare"
RAIL = 3400.0          # just below the 3550 ADC rail


def load(leg):
    x = pd.read_parquet(BASE / POINT / f"wf_{leg}_x.parquet")
    y = pd.read_parquet(BASE / POINT / f"wf_{leg}_y.parquet")
    return x, y


def med(a):
    return float(np.median(a)) if len(a) else float("nan")


rows = {}
for leg in ("sim", "data"):
    x, y = load(leg)
    p = x.merge(y, on="event_id", suffixes=("_x", "_y"))
    unsat = p[(p.peak_amp_x < RAIL) & (p.peak_amp_y < RAIL)]

    r = {}
    r["published"] = dict(
        n_x=len(x), n_y=len(y),
        X=med(x.peak_amp), Y=med(y.peak_amp),
        YX=med(y.peak_amp) / med(x.peak_amp))
    r["paired"] = dict(
        n=len(p), X=med(p.peak_amp_x), Y=med(p.peak_amp_y),
        YX=med(p.peak_amp_y) / med(p.peak_amp_x),
        YX_per_event=med(p.peak_amp_y / p.peak_amp_x))
    r["paired_unsat"] = dict(
        n=len(unsat), X=med(unsat.peak_amp_x), Y=med(unsat.peak_amp_y),
        YX=med(unsat.peak_amp_y) / med(unsat.peak_amp_x),
        YX_per_event=med(unsat.peak_amp_y / unsat.peak_amp_x))
    r["sat_frac"] = dict(
        x=float((x.peak_amp >= RAIL).mean()),
        y=float((y.peak_amp >= RAIL).mean()))
    rows[leg] = r

print(f"{'leg':>6} {'selection':>16} {'n':>6} {'X med':>9} {'Y med':>9} "
      f"{'Y/X med':>8} {'med(Y/X)':>9}")
print("-" * 70)
for leg in ("sim", "data"):
    for sel in ("published", "paired", "paired_unsat"):
        d = rows[leg][sel]
        n = d.get("n", d.get("n_x"))
        pe = d.get("YX_per_event")
        pes = f"{pe:9.3f}" if pe is not None else " " * 9
        print(f"{leg:>6} {sel:>16} {n:>6d} {d['X']:9.1f} {d['Y']:9.1f} "
              f"{d['YX']:8.3f} {pes}")
    print(f"{'':>6} {'railed frac':>16}  x={rows[leg]['sat_frac']['x']:.3f} "
          f"y={rows[leg]['sat_frac']['y']:.3f}")

print()
for sel in ("published", "paired", "paired_unsat"):
    s, d = rows["sim"][sel]["YX"], rows["data"][sel]["YX"]
    print(f"{sel:>16}:  sim Y/X = {s:.3f}   data Y/X = {d:.3f}   "
          f"sim/data = {s/d:.3f}")

print("\nReference: unselected whole-detector data Y/X = 1.0002 "
      "(2607.0 / 2606.5 ADC, hv_slope audit)")

p = BASE / "t14_ang_trend" / "xy_paired_selection.json"
p.write_text(json.dumps(rows, indent=1))
print(f"\nwrote {p}")
