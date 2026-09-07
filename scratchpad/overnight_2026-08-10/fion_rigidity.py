#!/usr/bin/env python3
"""Does the f_ion dial displace the rise distribution RIGIDLY?

The angled ladder established that the sim-vs-data rise mismatch, measured at
inclined incidence where both legs are compressed to similar widths, is a
RIGID ~65-80 ns offset at every quantile -- not a missing fast population.

A rigid offset is a fingerprint. If f_ion also moves the rise distribution
rigidly, and by a comparable amount, that ties the offset to the ion fraction
specifically -- much stronger evidence than the p5-only demand curve in the S3
closeout, which only ever matched ONE statistic.

If instead f_ion stretches the distribution (moves the tail more than the
core), then f_ion is NOT a pure delay, the ~70 ns rigid offset is not purely
f_ion, and something else fixed-delay-like is unaccounted for -- shaper group
delay, template t0, or seeding latency are the candidates.

Both outcomes are informative, which is why this is worth running.

All legs here are DIAGNOSIS points already on disk (own directories, nothing
frozen touched). They are vertical points, so the sim legs are compared to
EACH OTHER, never to the vertical data leg -- the vertical data leg is known
to be broader than any sim leg in both tails (see OVERNIGHT_2026-08-10.md).
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"
QS = [5, 10, 25, 50, 75, 90]

# label -> (directory, nominal f_ion)
LEGS = [
    ("noions",   "t14_DIAGNOSIS_noions",      0.0000),
    ("fion030",  "t14_DIAGNOSIS_fion030",     0.3000),
    ("fion050",  "t14_DIAGNOSIS_fion050",     0.5000),
    ("fion070",  "t14_DIAGNOSIS_fion070",     0.7000),
    ("default",  "t14_compare",               0.9006),
]


def q(d, view, leg="sim"):
    r = pd.read_parquet(BASE / d / f"wf_{leg}_{view}.parquet")["rise_ns"].to_numpy()
    r = r[np.isfinite(r)]
    return np.percentile(r, QS), r.size


out = {}
for view in ("x", "y"):
    print(f"\n===== {view.upper()} view : sim rise quantiles [ns] vs f_ion =====")
    print(f"{'leg':>9} {'f_ion':>6} {'n':>6} "
          + " ".join(f"{'p'+str(p):>7}" for p in QS))
    tab = {}
    for lab, d, f in LEGS:
        try:
            qq, n = q(d, view)
        except FileNotFoundError:
            print(f"{lab:>9} {f:>6.4f}  MISSING")
            continue
        tab[lab] = (f, qq)
        print(f"{lab:>9} {f:>6.4f} {n:>6d} " + " ".join(f"{v:>7.1f}" for v in qq))

    # rigidity: shift of each leg relative to noions, quantile by quantile
    print(f"\n  shift vs noions [ns] -- rigid means flat across the row")
    print(f"  {'leg':>9} " + " ".join(f"{'p'+str(p):>7}" for p in QS) + "    span")
    base = tab["noions"][1]
    rig = {}
    for lab, (f, qq) in tab.items():
        if lab == "noions":
            continue
        d = qq - base
        span = float(d.max() - d.min())
        rig[lab] = dict(f_ion=f, shift=d.tolist(), span=span,
                        mean=float(d.mean()))
        print(f"  {lab:>9} " + " ".join(f"{v:>7.1f}" for v in d)
              + f"  {span:7.1f}")

    # the production step, which is the one to compare against the ~70 ns
    d_full = tab["default"][1] - base
    print(f"\n  full ion term (f_ion 0 -> 0.9006): mean shift "
          f"{d_full.mean():.1f} ns, span {d_full.max()-d_full.min():.1f} ns")
    print(f"  for contrast, the sim-vs-data offset at 10/20 deg had "
          f"span 11-32 ns on a ~70 ns offset")
    out[view] = dict(quantiles=QS,
                     abs={k: dict(f_ion=v[0], q=v[1].tolist())
                          for k, v in tab.items()},
                     shift_vs_noions=rig)

p = BASE / "t14_ang_trend" / "fion_rigidity.json"
p.write_text(json.dumps(out, indent=1))
print(f"\nwrote {p}")
