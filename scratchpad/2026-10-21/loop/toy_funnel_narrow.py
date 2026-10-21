"""Saturated funnel with narrow strips (2026-10-21). If the June print had narrow strips and wide charged gaps, a full funnel
puts most of the charge onto narrow lines, the configuration the line-attractor toy needs. Same toy as toy_funnel.py;
strip widths 150-400 um, with and without spreading across the strip. Reports, over a half period of delta, the best
point with |ds1| <= 0.03 and its a1, plus the global a1 max. Usage: toy_funnel_narrow.py"""
import json, os, numpy as np
from toy_funnel import run, respf, P
OUT = []
for sig_ind in (100., 200.):
    Rm = respf(sig_ind)
    for sigma in (300., 450., 600.):
        b0 = run(sigma, P, 0., Rm)
        for w in (150., 250., 400.):
            for eq in (0., 1.0):
                D = np.arange(0, 391, 30.); L = [run(sigma, w, d, Rm, eq) for d in D]
                a = np.array([l[0] for l in L]); s = np.array([l[1] for l in L]) - b0[1]
                ok = np.abs(s) <= 0.03; best = (a[ok].max(), D[ok][np.argmax(a[ok])]) if ok.any() else (np.nan, np.nan)
                OUT.append(dict(sig_ind=sig_ind, sigma=sigma, w=w, eq=eq, delta=D.tolist(), a1=a.tolist(), ds1=s.tolist()))
                print(f'ind {sig_ind:3.0f} sigma {sigma:3.0f} strip {w:3.0f} um eq {eq:.0f}: a1 max {a.max():+.3f} at delta {D[np.argmax(a)]:3.0f} (ds1 {s[np.argmax(a)]:+.3f}); '
                      f'best with |ds1|<=0.03: a1 {best[0]:+.3f} at delta {best[1]}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_funnel_narrow.json'), 'w'), indent=1)
