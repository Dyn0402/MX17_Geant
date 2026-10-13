#!/usr/bin/env python3
"""Cache det0 mean pulses (data) and det0_toy model curves for the slide deck:
det0_curves.npz (t grid, He data, Ar data, toy He at tp 130/140, toy Ar at tp 140 for
(f, T) = (0.9, 340), (0.9, 180), (0.5, 180), (0.2, 340))."""
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import det0_toy as D, he_template_fit as H
def job(c):
    gas, tp, f, T = c
    W = D.sim(gas, tp, f, T, n=9000, seed=11); g, a = H.fine_avg(W); return c, g, a
if __name__ == '__main__':
    He = H.collect('det0__*He_HV_Scan_4-2*', 495, 515, 0, .05); Ar = H.collect('det0__*Ar_Iso*', 440, 480, 0, .05)
    g, dHe = H.fine_avg(He); _, dAr = H.fine_avg(Ar)
    cfg = [('He', 130, .9, 40.), ('He', 140, .9, 40.), ('Ar', 140, .9, 340.), ('Ar', 140, .9, 180.), ('Ar', 140, .5, 180.), ('Ar', 140, .2, 340.)]
    with ProcessPoolExecutor(6) as ex:
        R = list(ex.map(job, cfg))
    out = dict(t=g, He=dHe, Ar=dAr, nHe=len(He), nAr=len(Ar))
    for (gas, tp, f, T), gg, a in R:
        out[f'toy_{gas}_tp{tp}_f{f}_T{int(T)}'] = np.interp(g, gg, a)
    np.savez('det0_curves.npz', **out); print(sorted(out))
