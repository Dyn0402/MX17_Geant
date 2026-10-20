"""Is det2's positive region a 'little det3'? (2026-10-20) det2's X asymmetry changes sign with y (y quarter 1: +0.17,
quarters 0 and 3: -0.08). Split det2 matched events by the Y peak's true row: POS = rows 128-255, NEG = rows 0-127 and
384-511, and compute det3 June-15's fingerprints on each: class asymmetries A_Y / A_X and Y-pad fraction in Y-share
bins (fingerprint.py), the x-position conditional E[alpha(n+1) | alpha(n-1)] (xpos_probe binning; effect against
det2's quiet middle rows 256-383 + det3 May), and the excess timing (timing.py).
Usage: det2_regions.py"""
import numpy as np, os, json
import r5lib as R, fingerprint as FP, timing as T, xpos_probe as XP
E = R.ev(R.BUILDS['det2']); py = E['py'] % 512; row = 64 * (py // 64) + 63 - py % 64
def sub(E, m): return {k: (v[m] if hasattr(v, '__len__') and len(v) == len(m) else v) for k, v in E.items()}
REG = {'det2 POS (rows 128-255)': (row >= 128) & (row < 256), 'det2 NEG (rows 0-127, 384-511)': (row < 128) | (row >= 384),
       'det2 MID (rows 256-383)': (row >= 256) & (row < 384)}
D3 = R.ev(R.BUILDS['det3 June-15']); OUT = {}
for nm, m in list(REG.items()) + [('det3 June-15', None)]:
    S = sub(E, m) if m is not None else D3
    o = FP.fp_arrays(S['Wx'], S['Wy'], S['px'], S['py']); FP.show(nm, o)
    t = T.prof(S['Wx'], S['Wy'], S['ypad']); print('   ' + T.line('timing', t))
    # x-position conditional (pooled Ysh 0.4-0.85)
    u = S['alx'][:, 3]; v = S['alx'][:, 5]; mm = S['ypad'] & (S['ysh'] >= .4) & (S['ysh'] < .85); r = []
    for a, b in zip(XP.UB[:-1], XP.UB[1:]):
        k = mm & (u >= a) & (u < b); r.append(v[k].mean() if k.sum() >= 25 else np.nan)
    print('   E[a+1|u] ' + ' '.join('  -  ' if np.isnan(x) else f'{x:.3f}' for x in r))
    OUT[nm] = dict(marg=o['marg'], A_Y=o['A_Y'].tolist(), A_X=o['A_X'].tolist(), fY=o['fY'].tolist(), xpos=r,
                   timing={k: t[k] for k in ('Yc', 'Xc', 'Xm1', 'exc', 'fit')})
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'det2_regions.json'), 'w'), indent=1, default=float)
