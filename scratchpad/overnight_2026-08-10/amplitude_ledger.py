#!/usr/bin/env python3
"""Amplitude-deficit ledger: size the saturation bias and the demanded gain.

Two quantitative pieces the elimination table needs and nobody had computed:

1. The amplitude deficit measured in CHARGE (q_sum) is f_ion-INDEPENDENT. Peak
   amplitude moves 0.5646 -> 0.6319 across f_ion 0.9006 -> 0, but q_sum sits at
   0.626-0.642 throughout, because f_ion redistributes charge in time and
   conserves it. So the charge deficit cannot be double-counted against the ion
   thread: whatever happens to f_eff, the charge deficit stays.

2. A SIGNED correction for the data leg's saturation depletion. The reco-quality
   cut drops railed waveforms, so the data leg is the bottom (1-f) of the true
   amplitude distribution, where f is the fraction removed from the top. If the
   leg is the bottom (1-f), then the leg's own quantile 0.5/(1-f) is the true
   median. That is computable read-only from the frozen parquet.

f is estimated from the saturation fractions already in the record: the
detector saturates at 0.326 (X) / 0.327 (Y) unselected, the leg at 0.260 /
0.110. Both are quoted in the T14 verdict and the hv_slope audit.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path.home() / "x17/response_sim/stageB_w2"
G_SIM = 24094.0          # T7 pooled meshfield 490 V

# unselected detector vs leg saturation fraction, per view (hv_slope audit)
SAT_DETECTOR = {"x": 0.326, "y": 0.327}
SAT_LEG = {"x": 0.260, "y": 0.110}

s = json.load(open(BASE / "t14_compare" / "t14_summary.json"))["views"]

print("=== 1. the charge deficit is f_ion-independent ===")
print(f"{'leg':>26} {'X peak':>8} {'X q_sum':>8} {'Y peak':>8} {'Y q_sum':>8}")
legs = [("default (f_ion 0.9006)", "t14_compare"),
        ("noions (f_ion 0)", "t14_DIAGNOSIS_noions"),
        ("fion030", "t14_DIAGNOSIS_fion030"),
        ("fion050", "t14_DIAGNOSIS_fion050"),
        ("fion070", "t14_DIAGNOSIS_fion070")]
qx = []
for lab, d in legs:
    v = json.load(open(BASE / d / "t14_summary.json"))["views"]
    qx.append(v["x"]["q_sum_reco_tight3deg"]["ratio"])
    print(f"{lab:>26} {v['x']['peak_amp_med']['ratio']:8.4f} "
          f"{v['x']['q_sum_reco_tight3deg']['ratio']:8.4f} "
          f"{v['y']['peak_amp_med']['ratio']:8.4f} "
          f"{v['y']['q_sum_reco_tight3deg']['ratio']:8.4f}")
print(f"  X q_sum ratio across the whole f_ion range: "
      f"{min(qx):.4f} - {max(qx):.4f}  (spread {100*(max(qx)-min(qx))/np.mean(qx):.1f} %)")

print("\n=== 2. signed saturation-depletion correction ===")
out = {}
for view in ("x", "y"):
    d = pd.read_parquet(BASE / "t14_compare" / f"wf_data_{view}.parquet")
    a = np.sort(d["peak_amp"].to_numpy(dtype=float))
    q = np.sort(d["q_event"].to_numpy(dtype=float))
    # fraction of the TRUE population missing from the top of the leg
    f = (SAT_DETECTOR[view] - SAT_LEG[view]) / (1.0 - SAT_LEG[view])
    qt = 0.5 / (1.0 - f)
    if qt >= 1.0:
        print(f"  {view}: f = {f:.3f} -> quantile {qt:.3f} out of range, skipped")
        continue
    a_obs, a_cor = np.percentile(a, 50), np.percentile(a, 100 * qt)
    q_obs, q_cor = np.percentile(q, 50), np.percentile(q, 100 * qt)
    print(f"  {view}: detector sat {SAT_DETECTOR[view]:.3f}, leg sat "
          f"{SAT_LEG[view]:.3f} -> f = {f:.3f}, true median at leg p{100*qt:.1f}")
    print(f"      peak    {a_obs:8.1f} -> {a_cor:8.1f}  (x{a_cor/a_obs:.3f})")
    print(f"      q_event {q_obs:8.1f} -> {q_cor:8.1f}  (x{q_cor/q_obs:.3f})")
    out[view] = dict(f_missing=float(f), quantile=float(qt),
                     peak_obs=float(a_obs), peak_corr=float(a_cor),
                     peak_factor=float(a_cor / a_obs),
                     q_obs=float(q_obs), q_corr=float(q_cor),
                     q_factor=float(q_cor / q_obs))

print("\n=== 3. the gain the deficit demands ===")
print(f"  T7 pooled sim gain at 490 V: {G_SIM:.0f}")
for view in ("x", "y"):
    r = s[view]["q_sum_reco_tight3deg"]["ratio"]
    g0 = G_SIM / r
    line = f"  {view}: q_sum ratio {r:.4f} -> G_demanded {g0:8.0f}"
    if view in out:
        g1 = g0 * out[view]["q_factor"]
        line += f"   saturation-corrected {g1:8.0f}"
    print(line)
print("\n  literature max stable Ar/iso bulk-MM gain: 3-4e4 "
      "(T14_CAMPAIGN, Saha/Saclay); det3 sparks above 500 V")

p = BASE / "t14_ang_trend" / "amplitude_ledger.json"
p.write_text(json.dumps(dict(sat_correction=out, q_sum_ratios=qx,
                             G_sim=G_SIM), indent=1))
print(f"\nwrote {p}")
