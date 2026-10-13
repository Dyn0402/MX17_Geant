#!/usr/bin/env python3
"""Relative width of the peak-amplitude distribution vs voltage, per dataset.
Pure gain scaling keeps ln(pHi/pLo) constant; avalanche feedback with a random number
of secondaries broadens it as bG -> 1. The noise floor inflates low quantiles at low gain
(which fakes broadening), so the noise-only peak max f (median of rays with no signal,
from the lowest-voltage point) is subtracted: W = ln((pHi - f)/(pLo - f)), used only where
pLo >= 3 f and pHi <= 0.7 rail. Spatial gain spread: the same W inside 50x50 mm cells
(>= 60 events), averaged."""
import glob, numpy as np, pandas as pd
C = '/home/dylan/.claude/jobs/913ff27a/tmp/c/'
SETS = [('det0 Ar/iso', 'det0__*Ar_Iso*', (50, 150, -200, 150), 22, [(.5, .9), (.5, .75)], {490}),
        ('det0 He/C2H6', 'det0__*He_HV_Scan_4-2*', (50, 150, -200, 150), 22, [(.75, .9), (.5, .75)], set()),
        ('det3o Ar/iso', 'det3old__*HV_Scan_5-5*', (-200, 100, -200, 50), 38, [(.1, .5), (.25, .5)], {510}),
        ('det3o Ar/CF4', 'det3old__*ArCF4*', (-200, 100, -200, 50), 38, [(.1, .5), (.25, .5)], set())]
for name, pat, (x0, x1, y0, y1), f, qpairs, excl in SETS:
    df = pd.concat([pd.read_parquet(p) for p in glob.glob(C + pat + '.parquet')])
    df = df[(df.ref_x > x0) & (df.ref_x < x1) & (df.ref_y > y0) & (df.ref_y < y1) & ~df.subrun.str.startswith('final')]
    df['cell'] = (np.floor(df.ref_x / 50)).astype(int) * 100 + np.floor(df.ref_y / 50).astype(int)
    for ql, qh in qpairs:
        print(f'{name}: W = ln((p{int(qh*100)}-f)/(p{int(ql*100)}-f)), f={f}   [all | cell mean]')
        for view in 'xy':
            out = []
            for v, g in df[df.view == view].groupby('volt'):
                lo, hi = np.quantile(g.peak_amp, [ql, qh])
                if lo < 3 * f or hi > 0.7 * 3872 or v in excl: continue
                cw = [np.log((np.quantile(c.peak_amp, qh) - f) / (np.quantile(c.peak_amp, ql) - f))
                      for _, c in g.groupby('cell') if len(c) >= 60 and np.quantile(c.peak_amp, ql) > 1.5 * f]
                out.append(f'{v}:{np.log((hi-f)/(lo-f)):.2f}|{np.mean(cw) if cw else np.nan:.2f}')
            print('  ', view, '  '.join(out))
J = pd.read_csv('/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-13/mesh_ladder_det3_2026-06-27.csv')
J = J[J.p02_ok & J.p50_ok]
print('June det3 (T14): W = ln(p50/p02), no floor subtraction (p02 >= 69 ADC)')
for v, g in J.groupby('view'):
    print('  ', v, '  '.join(f"{int(r['volt'])}{r['scan'][-1]}:{np.log(r['p50']/r['p02']):.2f}" for _, r in g.sort_values('volt').iterrows()))
