"""Audit (2026-10-21): a1(delta) and s1(delta) - s1(0) curves of the strip-quantisation toy (toy_strips.py) for the
deck: 700 um strips, full equalisation, cloud sigma 250 / 450 um, induction width 50 / 150 / 300 um.
-> toy_curves.json. Usage: python toy_curves.py"""
import json, os, numpy as np
import toy_strips as T
OUT = []
for si in (50., 150., 300.):
    T.RK = np.array([T.R(T.x - k * T.P, si) for k in range(-3, 4)])
    for sg in (250., 450.):
        a0, s0 = T.run(sg, 700., 1., 0.); pts = []
        for d in np.arange(0, 346, 30.):
            a, s = T.run(sg, 700., 1., d); pts.append([float(d), float(a), float(s - s0)])
        OUT.append(dict(sig_ind=si, sigma=sg, w=700., s1_0=float(s0), pts=pts)); print(si, sg, pts[-3:])
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_curves.json'), 'w'), indent=0)
