"""Reach test (2026-10-20): over a wide grid of linear channel mixings, the most negative narrow-cluster X-pad asymmetry
any mixing reaches on copy-free templates, with the Y-pad narrow bins held near the data (+0.26 +- 0.06). Usage: mix_reach.py"""
import itertools, numpy as np, mix_fit as F
for tb in ('det2', 'det3 May', 'det1'):
    W, px, py, ysh = F.UP.events(tb); W = W.astype(float); wts = np.ones(len(W)); best = None; best_any = None
    for mm1, mp1, mm2, mp2 in itertools.product(np.arange(0, 1.01, .1), np.arange(-.5, .51, .1), (-.2, 0, .2), (-.2, 0, .2)):
        t = F.analyse(F.mix(W, {-1: mm1, 1: mp1, -2: mm2, 2: mp2}), px, py, ysh, wts)
        xn = np.nanmean([t[2, 0], t[7, 0]]); yn = np.nanmean([t[0, 0], t[5, 0]])
        if best_any is None or xn < best_any[0]: best_any = (xn, yn, mm1, mp1, mm2, mp2)
        if abs(yn - .265) < .06 and (best is None or xn < best[0]): best = (xn, yn, mm1, mp1, mm2, mp2)
    print(f'{tb:9s} lowest narrow X-pad asym with Y-pad narrow ~ data: {best[0]:+.2f} (Ypad {best[1]:+.2f}; m-1 {best[2]:.1f} m+1 {best[3]:+.1f} m-2 {best[4]:+.1f} m+2 {best[5]:+.1f})'
          f' | any mixing: {best_any[0]:+.2f} (Ypad {best_any[1]:+.2f}; {best_any[2]:.1f} {best_any[3]:+.1f} {best_any[4]:+.1f} {best_any[5]:+.1f})   data: -0.35 / +0.265')
