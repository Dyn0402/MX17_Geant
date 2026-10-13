#!/usr/bin/env python3
"""Target run (fiducial box from the M3 position: misses removed) (490 V), threshold-free peak strip vs the TRUE projected angle from M3.
View x inclination = +atan(tan_y), view y = -atan(tan_x) (from the reco correlation, +0.73/-0.71).
Prints: peak median / p25 / sat fraction by |true angle| in the view (other projection |.|<6 deg),
the reco-vs-true relation, and what T14's selections (reco window, quality cut) do to the median."""
import os, numpy as np, pandas as pd
pk = pd.read_parquet(os.path.expanduser('~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet'))
pk = pk[(pk.ref_x > -200) & (pk.ref_x < 120) & (pk.ref_y.abs() < 170)]   # fiducial: beyond ref_x 125 the ray misses det3 (peak < 150 ADC on 90 %)
pk['th'] = np.where(pk.view == 'x', pk.ay, -pk.ax); pk['oth'] = np.where(pk.view == 'x', pk.ax, -pk.ay)
B = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 26]
for v in 'xy':
    d = pk[(pk.view == v) & (pk.oth.abs() < 6)]
    print(f'\n=== view {v}: n {len(d)} (|other| < 6 deg)')
    print(' |true| bin   n   peak med   p25   sat3500  | reco med  reco IQR  qual_ok frac  peak med(qual_ok)')
    for lo, hi in zip(B[:-1], B[1:]):
        m = (d.th.abs() >= lo) & (d.th.abs() < hi); x = d[m]
        if len(x) < 15: continue
        r = (x.reco_theta * np.sign(x.th)).dropna()
        q = x[x.quality_ok == True]
        print(f' {lo:3d}-{hi:<3d} {len(x):5d} {x.peak_amp.median():8.0f} {x.peak_amp.quantile(.25):6.0f} {(x.peak_amp >= 3500).mean():7.3f}   | {r.median():7.2f} {r.quantile(.75) - r.quantile(.25):8.2f} {(x.quality_ok == True).mean():10.3f} {q.peak_amp.median():12.0f}')
    # the T14 selections, on the vertical-true sample and on the T14 reco window
    t = d[d.th.abs() < 3]
    w = pk[(pk.view == v) & (pk.reco_theta.abs() <= 3)]
    print(f' true |th|<3 (other<6): median {t.peak_amp.median():.0f} (n {len(t)}); T14-like reco |th|<=3 + quality: {w[w.quality_ok == True].peak_amp.median():.0f} (n {(w.quality_ok == True).sum()}); reco window no quality: {w.peak_amp.median():.0f}')
    for c in (10, 20):
        t = d[(d.th.abs() > c - 3) & (d.th.abs() < c + 3)]
        w = pk[(pk.view == v) & ((pk.reco_theta.abs() - c).abs() <= 3)]
        print(f' true |th| {c}+-3: median {t.peak_amp.median():.0f} (n {len(t)}); reco window + quality: {w[w.quality_ok == True].peak_amp.median():.0f} (n {(w.quality_ok == True).sum()}); reco window no quality: {w.peak_amp.median():.0f}')
