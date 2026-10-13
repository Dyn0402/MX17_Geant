#!/usr/bin/env python3
"""det0_toy.py with more statistics and systematics: He shaper calibration over
tp 110-160 ns x He drift speed 12/15/20 um/ns; then Ar (f, T) at each He-allowed tp."""
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import det0_toy as D, he_template_fit as H
def job(c):
    gas, tp, f, T, v = c
    if v: D.GAS[gas] = dict(D.GAS[gas], v=v)
    W = D.sim(gas, tp, f, T, n=9000, seed=11)
    g, a = H.fine_avg(W); return c, g, a, len(W)
if __name__ == '__main__':
    He = H.collect('det0__*He_HV_Scan_4-2*', 495, 515, 0, .05); Ar = H.collect('det0__*Ar_Iso*', 440, 480, 0, .05)
    g, dHe = H.fine_avg(He); _, dAr = H.fine_avg(Ar); m = (g >= -150) & (g <= 700)
    rms = lambda gg, a, d: float(np.sqrt(np.mean((np.interp(g, gg, a) - d)[m] ** 2)))
    with ProcessPoolExecutor(12) as ex:
        R = list(ex.map(job, [('He', tp, 0.9, 40., v) for tp in (110, 120, 130, 140, 150, 160) for v in (0.012, 0.015, 0.020)]))
    print('He/C2H6: rms by shaper tp (rows) x He drift v (cols 12/15/20 um/ns)')
    tab = {}
    for c, gg, a, n in R: tab[(c[1], c[4])] = rms(gg, a, dHe)
    for tp in (110, 120, 130, 140, 150, 160):
        print(f'   tp {tp}: ' + '  '.join(f'{tab[(tp, v)]:.4f}' for v in (0.012, 0.015, 0.020)))
    cfg = [('Ar', tp, f, T, None) for tp in (120, 130, 140, 150) for f in (0.2, 0.5, 0.9) for T in (100., 140., 180., 220., 340.)]
    with ProcessPoolExecutor(12) as ex:
        R = list(ex.map(job, cfg))
    print('Ar/iso: rms by (tp, f) rows x T cols 100/140/180/220/340')
    tab = {(c[1], c[2], c[3]): rms(gg, a, dAr) for c, gg, a, n in R}
    for tp in (120, 130, 140, 150):
        for f in (0.2, 0.5, 0.9):
            print(f'   tp {tp} f {f:.1f}: ' + '  '.join(f'{tab[(tp, f, T)]:.4f}' for T in (100., 140., 180., 220., 340.)))
