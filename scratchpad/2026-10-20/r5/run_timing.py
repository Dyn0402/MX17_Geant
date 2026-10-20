"""Run timing.prof on data sessions, controls (det2, det3 May with an injected copy), and sim arms (2026-10-20).
Usage: run_timing.py [extra sim arm patterns...]"""
import os, sys, json, glob, re, numpy as np
import r5lib as R, timing as T
from eval_mix import SESS, arm, SIM, MIX
OUT = {}
def classify(Wx, Wy, px, py):
    cy = Wy[:, 4].astype(float); ny = np.einsum('ij,ij->i', cy, cy)
    ysh = (np.einsum('ij,ij->i', Wy[:, 3], cy) + np.einsum('ij,ij->i', Wy[:, 5], cy)) / ny; par = (px + py) % 2
    return par == int(ysh[par == 1].mean() < ysh[par == 0].mean())
for k in ('6-22/23 1000V', '6-23/24 600V', '6-27/28 1000V', 'det2', 'det3 May'):
    E = R.ev(SESS[k]); o = T.prof(E['Wx'], E['Wy'], E['ypad']); OUT[k] = o; print(T.line(k, o))
    if k in ('det2', 'det3 May'):
        Wx = E['Wx'].astype(float).copy(); Wx[:, 1:] += 0.25 * E['Wx'][:, :-1]
        o = T.prof(Wx, E['Wy'], E['ypad']); OUT[k + ' +copy .25'] = o; print(T.line(k + ' +copy .25', o))
pats = {'sim D0': f'{SIM}/xo0_{{ang}}_s{{s}}_x.npz', 'sim D-200': f'{SIM}/xom200_{{ang}}_s{{s}}_x.npz', 'sim D-300': f'{SIM}/xom300_{{ang}}_s{{s}}_x.npz'}
for a in sys.argv[1:]: pats[a] = a
for nm, pat in pats.items():
    a = arm(pat)
    if a is None: continue
    Wx, Wy, px, py = a; yp = classify(Wx, Wy, px, py); o = T.prof(Wx, Wy, yp); OUT[nm] = o; print(T.line(nm, o))
    if nm == 'sim D0':
        Wc = Wx.astype(float).copy(); Wc[:, 1:] += 0.25 * Wx[:, :-1]; o = T.prof(Wc, Wy, yp); OUT['sim D0 +copy .25'] = o; print(T.line('sim D0 +copy .25', o))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'timing.json'), 'w'))
