"""Strip beat with the right noise unit (2026-10-20): per-column means (detrended, sigma 30 columns), Fourier amplitude
at f, null = the same with the column means permuted (1000 permutations). Reports p-value at f = 1/40 and the max over
f in 0.015..0.040 (look-elsewhere included via the permutation max). X and Y views; Y is the control.
Usage: beat_col.py"""
import json, os, numpy as np
import r5lib as R
from beat_sig import detrend
F = np.arange(0.015, 0.0401, 0.00025); rng = np.random.default_rng(5); OUT = {}
def colmeans(col, v):
    n = np.bincount(col, minlength=512); s = np.bincount(col, v, minlength=512); ok = n >= 15
    return np.flatnonzero(ok), (s[ok] / n[ok])
def A(c, m, f): m = m - m.mean(); return np.abs(np.exp(2j * np.pi * np.outer(np.atleast_1d(f), c)) @ m) * 2 / len(m)
for b in ('det3 June-15', 'det3 May', 'det2', 'det1', 'det4', 'det6', 'det7'):
    for view in ('X', 'Y'):
        p, s1, d1 = [], [], []
        for r in R.rows():
            if r['build'] != b or r['view'] != view: continue
            d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
            p.append(d['p'][ok] % 512); s1.append((d['al'][1] + d['al'][-1])[ok]); d1.append(d['a1'][ok])
        p, s1, d1 = map(np.concatenate, (p, s1, d1)); col = 64 * (p // 64) + 63 - p % 64
        m = (np.abs(s1) < 3) & (np.abs(d1) < 3); col, s1, d1 = col[m], s1[m], d1[m]
        line = f'{b:13s} {view}'
        for nm, v in (('s1', detrend(col, s1)), ('d1', detrend(col, d1))):
            c, mc = colmeans(col, v)
            a40 = A(c, mc, 1 / 40)[0]; aF = A(c, mc, F); 
            null40 = np.array([A(c, rng.permutation(mc), 1 / 40)[0] for _ in range(1000)])
            nullF = np.array([A(c, rng.permutation(mc), F).max() for _ in range(300)])
            p40 = (null40 >= a40).mean(); pF = (nullF >= aF.max()).mean()
            OUT[f'{b}|{view}|{nm}'] = dict(ncol=int(len(c)), A40=float(a40), p40=float(p40), fbest=float(F[aF.argmax()]), Abest=float(aF.max()), pbest=float(pF))
            line += f' | {nm} ncol {len(c)} A(1/40)={a40:.4f} p={p40:.3f}; best f={F[aF.argmax()]:.4f} A={aF.max():.4f} p={pF:.2f}'
        print(line)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'beat_col.json'), 'w'), indent=1)
