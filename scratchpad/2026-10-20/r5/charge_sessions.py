"""(C3) Bookkeeping across det3 June-15 sessions (2026-10-20): same board and electronics, copy 0.10 -> 0.26.
Per Y-share bin (pad class Y and X): X(n-1)/SY, X(n+1)/SY, Xc/SY, (X(n-1)+X(n+1))/SY, SX/SY, integrals of positive
samples. If the excess on n+1 is moved charge, n-1 falls as the copy rises at fixed Y share; if it is added, n-1 stays
and SX/SY rises. Usage: charge_sessions.py"""
import json, os, numpy as np
import r5lib as R
from probe_charge_bins import ED
SESS = {'6-27/28 1000V F7': lambda r: r['build'] == 'det3 June-15' and r['feu'] in (7, 8) and r['drift'] == '1000',
        '6-22/23 1000V F3': lambda r: r['build'] == 'det3 June-15' and 'det2_det3_overnight_6-22' in r['run'],
        '6-23/24 600V F3': lambda r: r['build'] == 'det3 June-15' and '6-23-26' in r['run'],
        'det3 May': lambda r: r['build'] == 'det3 May', 'det2': lambda r: r['build'] == 'det2'}
OUT = {}
for nm, sel in SESS.items():
    E = R.ev(sel); qx, qy = E['qx'], E['qy']; SY = qy.sum(1); b = np.digitize(E['ysh'], ED) - 1; o = {}
    for cls, cm in (('Y', E['ypad']), ('X', ~E['ypad'])):
        rows = []
        for k in range(len(ED) - 1):
            m = cm & (b == k) & (SY > 0)
            if m.sum() < 40: rows.append(None); continue
            f = lambda v: float(np.mean(v[m] / SY[m]))
            rows.append(dict(n=int(m.sum()), m1=f(qx[:, 3]), p1=f(qx[:, 5]), c=f(qx[:, 4]), s=f(qx[:, 3] + qx[:, 5]), sx=f(qx.sum(1)),
                             a1=float(E['alx'][m, 5].mean() - E['alx'][m, 3].mean())))
        o[cls] = rows
    OUT[nm] = o
    a1 = E['alx'][:, 5] - E['alx'][:, 3]
    print(f'== {nm}: n={len(SY)}  a1 all {a1.mean():+.3f}  Ypad {a1[E["ypad"]].mean():+.3f}  median Ysh {np.median(E["ysh"]):.2f}')
    for cls in ('Y', 'X'):
        print(f'  {cls}pad  n-1   ' + ' '.join('  -  ' if r is None else f'{r["m1"]:.3f}' for r in o[cls]))
        print(f'  {cls}pad  n+1   ' + ' '.join('  -  ' if r is None else f'{r["p1"]:.3f}' for r in o[cls]))
        print(f'  {cls}pad  sum   ' + ' '.join('  -  ' if r is None else f'{r["s"]:.3f}' for r in o[cls]))
        print(f'  {cls}pad  Xc    ' + ' '.join('  -  ' if r is None else f'{r["c"]:.3f}' for r in o[cls]))
        print(f'  {cls}pad  SX/SY ' + ' '.join('  -  ' if r is None else f'{r["sx"]:.3f}' for r in o[cls]))
        print(f'  {cls}pad  n     ' + ' '.join('  -  ' if r is None else f'{r["n"]:5d}' for r in o[cls]))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'charge_sessions.json'), 'w'), indent=1)
