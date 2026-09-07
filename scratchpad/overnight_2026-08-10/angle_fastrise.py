#!/usr/bin/env python3
"""Q2 verdict: does track inclination generate the data's fast-rise population?

The angled ladder ran and t14_ang_trend/report.html lays out the test, but its
"rise floor -- geometry or electronics?" section stops at the figure captions
and never states the answer. This computes it at BOTH thresholds (the report's
<200 ns and the queue's <240 ns) so the number can't be an artifact of one cut.

Reads only the fixed-angle-window (t14_angW_*) parquets, which are the correct
set per ANGLED_LADDER_2026-08-09.md section 5.
"""
import json
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"

# (point, tilted view, |gun angle| deg)
LADDER = [
    ("t14_angW_th00", "x", 0.0),
    ("t14_angW_thx10", "x", 10.0),
    ("t14_angW_thx20", "x", 20.0),
    ("t14_angW_th00", "y", 0.0),
    ("t14_angW_thy10", "y", 10.0),
    ("t14_angW_thy20", "y", 20.0),
]


def stats(df):
    r = df["rise_ns"].to_numpy()
    r = r[np.isfinite(r)]
    return dict(
        n=int(r.size),
        p05=float(np.percentile(r, 5)),
        p25=float(np.percentile(r, 25)),
        p50=float(np.percentile(r, 50)),
        f200=float((r < 200).mean()),
        f240=float((r < 240).mean()),
    )


rows = []
for point, view, ang in LADDER:
    d = {}
    for leg in ("sim", "data"):
        p = BASE / point / f"wf_{leg}_{view}.parquet"
        d[leg] = stats(pd.read_parquet(p))
    rows.append(dict(point=point, view=view, deg=ang, sim=d["sim"], data=d["data"]))

hdr = (f"{'view':>4} {'deg':>5} | {'n_sim':>6} {'n_dat':>6} | "
       f"{'p05_sim':>8} {'p05_dat':>8} {'gap':>7} | "
       f"{'f240_sim':>8} {'f240_dat':>8} {'ratio':>6} | "
       f"{'f200_sim':>8} {'f200_dat':>8}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    s, dd = r["sim"], r["data"]
    ratio = dd["f240"] / s["f240"] if s["f240"] > 0 else float("inf")
    print(f"{r['view']:>4} {r['deg']:>5.0f} | {s['n']:>6d} {dd['n']:>6d} | "
          f"{s['p05']:>8.1f} {dd['p05']:>8.1f} {s['p05']-dd['p05']:>7.1f} | "
          f"{s['f240']:>8.3f} {dd['f240']:>8.3f} {ratio:>6.1f} | "
          f"{s['f200']:>8.3f} {dd['f200']:>8.3f}")

# Angular leverage: how much fast-fraction does each leg BUY per 20 deg of tilt?
print()
print("Angular leverage on the fast fraction (<240 ns), 0 -> 20 deg:")
for view in ("x", "y"):
    v = [r for r in rows if r["view"] == view]
    v.sort(key=lambda r: r["deg"])
    d_sim = v[-1]["sim"]["f240"] - v[0]["sim"]["f240"]
    d_dat = v[-1]["data"]["f240"] - v[0]["data"]["f240"]
    lev = d_dat / d_sim if d_sim != 0 else float("inf")
    print(f"  {view}: sim {v[0]['sim']['f240']:.3f} -> {v[-1]['sim']['f240']:.3f} "
          f"(+{d_sim:.3f} pp/20deg)   "
          f"data {v[0]['data']['f240']:.3f} -> {v[-1]['data']['f240']:.3f} "
          f"(+{d_dat:.3f})   data/sim leverage = x{lev:.1f}")

print()
print("p05 floor: does the sim's floor fall fast enough to reach the data?")
for view in ("x", "y"):
    v = [r for r in rows if r["view"] == view]
    v.sort(key=lambda r: r["deg"])
    print(f"  {view}: sim p05 {v[0]['sim']['p05']:.1f} -> {v[-1]['sim']['p05']:.1f} ns "
          f"({v[-1]['sim']['p05']-v[0]['sim']['p05']:+.1f});  "
          f"data {v[0]['data']['p05']:.1f} -> {v[-1]['data']['p05']:.1f} ns "
          f"({v[-1]['data']['p05']-v[0]['data']['p05']:+.1f});  "
          f"gap {v[0]['sim']['p05']-v[0]['data']['p05']:.1f} -> "
          f"{v[-1]['sim']['p05']-v[-1]['data']['p05']:.1f} ns")

out = BASE / "t14_ang_trend" / "fastrise_verdict.json"
out.write_text(json.dumps(rows, indent=1))
print(f"\nwrote {out}")
