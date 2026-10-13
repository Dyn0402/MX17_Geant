#!/usr/bin/env python3
"""bG x shaper scan summary (condor 4451431/4451432, plus the tp-180 fullphys twins and noisecm A/B).
Amplitude: ledger_m3.csv (data by TRUE M3 angle, sim all events): p50 and p25 data/sim, sat fraction,
per view at 0/10/20 deg. Shape: sim rise/FWHM medians (T14 summary, sim leg) against the data by TRUE
angle (rise_vs_m3.log; censored rise, as T14): X 10 deg 187/321, 20 deg 169/295; Y 10 deg 188/310,
20 deg 160/286 (rise/FWHM ns). Runs on yoga after ledger_m3.py has covered the new dirs."""
import os, json, glob, numpy as np, pandas as pd
SB = os.path.expanduser('~/x17/response_sim/stageB_w2')
L = pd.read_csv(os.path.expanduser('~/ledger_m3.csv'))
DATA_SHAPE = {('x', 10): (187, 321), ('x', 20): (169, 295), ('y', 10): (188, 310), ('y', 20): (160, 286)}
pts = sorted(set(L.kind))
print('amplitude, data/sim  p50 | p25 | sat d/s   (X0 X10 X20 | Y0 Y10 Y20)')
rows = []
for k in pts:
    cells = []; r50 = []; r25 = []
    for v, a, ang in (('x', 0, 'th00'), ('x', 10, 'thx10'), ('x', 20, 'thx20'), ('y', 0, 'th00'), ('y', 10, 'thy10'), ('y', 20, 'thy20')):
        m = L[(L.kind == k) & (L.view == v) & (L.ang == ang)]
        if not len(m): cells.append('      -       '); continue
        m = m.iloc[0]; r50.append(m.r50); r25.append(m.r25)
        cells.append(f'{m.r50:.2f}|{m.r25:.2f}|{m.sat_d:.2f}/{m.sat_s:.2f}')
    g50 = np.exp(np.mean(np.log(r50))) if r50 else np.nan; g25 = np.exp(np.mean(np.log(r25))) if r25 else np.nan
    sp = np.std(np.log(r50)) if r50 else np.nan
    rows.append(dict(kind=k, gm50=g50, gm25=g25, spread50=sp))
    print(f'{k:16s} ' + ' '.join(f'{c:>20s}' for c in cells) + f'   geo-mean p50 {g50:.3f} p25 {g25:.3f} rms(ln p50) {sp:.3f}')
print('\nshape (inclined), sim - data [ns]: rise | FWHM   (X10 X20 Y10 Y20)')
for k in pts:
    cells = []
    for v, a, ang in (('x', 10, 'thx10'), ('x', 20, 'thx20'), ('y', 10, 'thy10'), ('y', 20, 'thy20')):
        f = f'{SB}/t14_{k}_{ang}/t14_summary.json'
        if not os.path.exists(f): cells.append('    -    '); continue
        s = json.load(open(f))['views']; s = eval(s) if isinstance(s, str) else s
        r, w = s[v]['rise_ns_med']['sim'], s[v]['fwhm_ns_med']['sim']; dr, dw = DATA_SHAPE[(v, a)]
        cells.append(f'{r - dr:+4.0f}|{w - dw:+4.0f}')
    print(f'{k:16s} ' + ' '.join(f'{c:>10s}' for c in cells))
pd.DataFrame(rows).to_csv(os.path.expanduser('~/scan_summary.csv'), index=False)
