#!/usr/bin/env python3
"""Is the sim's rise deficit an angular-RESPONSE failure or a constant OFFSET?

A fast-fraction counter at a fixed threshold conflates the two: if the sim's
whole distribution sits above the threshold, a large shift registers as almost
no change in the count. The threshold-free test is to compare the two legs
quantile by quantile and ask whether sim-minus-data is (a) large and roughly
constant with angle -> a fixed offset, or (b) growing with angle -> the sim's
angular response itself is wrong.

Data legs are angle-matched at +-3 deg (t14_angW_*, the fixed-window set).
Both legs carry the selection caveats of the parent campaign; the comparison
here is sim-vs-data at each angle, so a bias common to all angles cancels out
of the trend even though it does not cancel out of the absolute offset.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"
QS = [5, 10, 25, 50, 75, 90]
LADDER = {
    "x": [("t14_angW_th00", 0.0), ("t14_angW_thx10", 10.0), ("t14_angW_thx20", 20.0)],
    "y": [("t14_angW_th00", 0.0), ("t14_angW_thy10", 10.0), ("t14_angW_thy20", 20.0)],
}


def q(point, leg, view):
    r = pd.read_parquet(BASE / point / f"wf_{leg}_{view}.parquet")["rise_ns"].to_numpy()
    r = r[np.isfinite(r)]
    return np.percentile(r, QS), r.size


out = {}
for view, pts in LADDER.items():
    print(f"\n===== view {view} : rise-time quantiles [ns] =====")
    print(f"{'deg':>4} {'leg':>5} {'n':>6} " + " ".join(f"{'p'+str(p):>7}" for p in QS))
    tab = {}
    for point, deg in pts:
        qs, ns = q(point, "sim", view)
        qd, nd = q(point, "data", view)
        tab[deg] = (qs, qd)
        print(f"{deg:>4.0f} {'sim':>5} {ns:>6d} " + " ".join(f"{v:>7.1f}" for v in qs))
        print(f"{deg:>4.0f} {'data':>5} {nd:>6d} " + " ".join(f"{v:>7.1f}" for v in qd))
        print(f"{'':>4} {'diff':>5} {'':>6} " + " ".join(f"{a-b:>7.1f}" for a, b in zip(qs, qd)))

    print(f"\n  sim - data offset [ns] vs angle (is it constant?)")
    print(f"  {'deg':>4} " + " ".join(f"{'p'+str(p):>7}" for p in QS))
    offs = []
    for deg in sorted(tab):
        qs, qd = tab[deg]
        o = qs - qd
        offs.append(o)
        print(f"  {deg:>4.0f} " + " ".join(f"{v:>7.1f}" for v in o))
    offs = np.array(offs)
    print(f"  {'rng':>4} " + " ".join(f"{v:>7.1f}" for v in offs.max(0) - offs.min(0))
          + "   <- spread of the offset across 0-20 deg")

    print(f"\n  how much each leg SPEEDS UP from 0 to 20 deg [ns, negative = faster]")
    qs0, qd0 = tab[0.0]
    qs20, qd20 = tab[20.0]
    print(f"  {'sim':>5} " + " ".join(f"{v:>7.1f}" for v in qs20 - qs0))
    print(f"  {'data':>5} " + " ".join(f"{v:>7.1f}" for v in qd20 - qd0))
    ratio = (qs20 - qs0) / (qd20 - qd0)
    print(f"  {'ratio':>5} " + " ".join(f"{v:>7.2f}" for v in ratio)
          + "   <- 1.00 = sim's angular response matches data's")
    out[view] = dict(
        quantiles=QS,
        offset_by_deg={str(d): (tab[d][0] - tab[d][1]).tolist() for d in sorted(tab)},
        speedup_sim=(qs20 - qs0).tolist(),
        speedup_data=(qd20 - qd0).tolist(),
        response_ratio=ratio.tolist(),
    )

p = BASE / "t14_ang_trend" / "rise_offset_vs_angle.json"
p.write_text(json.dumps(out, indent=1))
print(f"\nwrote {p}")
