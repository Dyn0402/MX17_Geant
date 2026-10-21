"""Audit (2026-10-21): a1 and the +-1 sharings per band, Y-pad class. det2 by 64-row band (its field runs in y),
det4 / det3 June-15 / det6 / det7 by 64-column chip (their fields run in x), det3 May control. Output bands.json.
Usage: NB_RAW=~/x17/nb_asym_raw python bands.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
OUT = {}
def mse(v): return (float(v.mean()), float(v.std() / np.sqrt(len(v)))) if len(v) > 20 else (np.nan, np.nan)
for b, axis in (('det2', 'y'), ('det3 June-15', 'x'), ('det4', 'x'), ('det6', 'x'), ('det7', 'x'), ('det3 May', 'x'), ('det3 May', 'y'), ('det3 June-15', 'y')):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    k = (cy if axis == 'y' else cx) // 64
    a1 = E['alx'][:, 5] - E['alx'][:, 3]; sx = E['alx'][:, 5] + E['alx'][:, 3]; sy = E['ysh']
    a1y = E['aly'][:, 5] - E['aly'][:, 3]
    rows = []
    for cls in ('ypad', 'all'):
        cm = E['ypad'] if cls == 'ypad' else np.ones(len(a1), bool)
        for q in range(8):
            m = cm & (k == q)
            rows.append(dict(cls=cls, band=q, n=int(m.sum()), a1=mse(a1[m]), s1x=mse(sx[m]), s1y=mse(sy[m]), a1y=mse(a1y[m]),
                             pm1=mse(E['alx'][m, 3]), pp1=mse(E['alx'][m, 5])))
    OUT[f'{b}|{axis}'] = rows
    print(f'== {b} by {axis} band (Y-pad class): band n a1 | s1X s1Y | a1Y | a(-1) a(+1)')
    for r in rows:
        if r['cls'] == 'ypad':
            print(f"  {r['band']} {r['n']:5d} {r['a1'][0]:+.3f}±{r['a1'][1]:.3f} | {r['s1x'][0]:.3f}±{r['s1x'][1]:.3f} {r['s1y'][0]:.3f} | {r['a1y'][0]:+.3f} | {r['pm1'][0]:.3f} {r['pp1'][0]:.3f}")
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bands.json'), 'w'), indent=1)
