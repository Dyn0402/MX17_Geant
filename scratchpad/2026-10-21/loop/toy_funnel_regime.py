"""Saturated-gap funnel in the narrow-cloud / wide-induction regime (2026-10-21). toy_local_narrow.py shows a per-electron
pull onto narrow lines reaches a1 0.1-0.2 at |ds1| <= 0.03 when the electron cloud is narrow (sigma <= 200 um) and the
induction wide (~350 um). The physical version: charged gaps (saturated) funnel every gap electron onto the nearest strip
edge; with narrow strips that is a pull onto narrow lines. Scan strip width 100-400 um, cloud 100-200 um, induction 300/350.
Usage: toy_funnel_regime.py"""
import json, os, numpy as np
from toy_funnel import run, respf, P
OUT = []
for sig_ind in (300., 350.):
    Rm = respf(sig_ind)
    for sigma in (100., 150., 200.):
        b0 = run(sigma, P, 0., Rm, n=120)
        for w in (100., 200., 300., 400.):
            D = np.arange(0, 391, 30.); L = [run(sigma, w, d, Rm, 0., n=120) for d in D]
            a = np.array([l[0] for l in L]); s = np.array([l[1] for l in L]) - b0[1]
            ok = np.abs(s) <= 0.03; best = float(a[ok].max()) if ok.any() else float('nan')
            OUT.append(dict(sig_ind=sig_ind, sigma=sigma, w=w, s1_nopull=b0[1], delta=D.tolist(), a1=a.tolist(), ds1=s.tolist(), best=best))
            print(f'ind {sig_ind:3.0f} sigma {sigma:3.0f} (s1 {b0[1]:.2f}) strip {w:3.0f} um: a1 max {a.max():+.3f} at delta {D[np.argmax(a)]:3.0f} (ds1 {s[np.argmax(a)]:+.3f}); best |ds1|<=0.03: {best:+.3f}', flush=True)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_funnel_regime.json'), 'w'), indent=1)
