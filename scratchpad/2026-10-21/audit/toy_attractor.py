"""Audit toy (2026-10-21): any mechanism that pulls a fraction f of each cloud toward a pad-pitch-periodic attractor
(strip centre) at offset delta from the pad centre makes an ODD a1(delta) and an EVEN change of the +-1 sharing s1.
1-D, pads of pitch P, point-charge pad response = box(P) * gaussian(sig_ind); cloud gaussian sigma; uniform positions;
peak-pad selection. Pulled fraction f is moved to the nearest attractor and smeared by sig_line (strip width).
Prints, per (sigma, f, sig_line), the delta needed for a1 = 0.10 / 0.20 and the accompanying s1 change.
Usage: python toy_attractor.py"""
import json, os, numpy as np
from math import erf
P = 780.; x = np.arange(-6 * P, 6 * P, 5.)
def prf(sig_ind):
    from scipy.special import erf as E
    s = max(sig_ind, 1.)
    return lambda u: 0.5 * (E((u + P / 2) / (np.sqrt(2) * s)) - E((u - P / 2) / (np.sqrt(2) * s)))
def run(sigma, f, sig_line, delta, sig_ind=150., n=400):
    R = prf(sig_ind); a1s, s1s = [], []
    for x0 in np.linspace(-P / 2, P / 2, n, endpoint=False):
        rho = (1 - f) * np.exp(-0.5 * ((x - x0) / sigma) ** 2) / sigma
        xs = np.round((x0 - delta) / P) * P + delta           # nearest attractor
        rho += f * np.exp(-0.5 * ((x - xs) / max(sig_line, 5.)) ** 2) / max(sig_line, 5.)
        q = np.array([(rho * R(x - k * P)).sum() for k in range(-3, 4)])
        c = np.argmax(q); 
        if c in (0, 6): continue
        a1s.append((q[c + 1] - q[c - 1]) / q[c]); s1s.append((q[c + 1] + q[c - 1]) / q[c])
    return np.mean(a1s), np.mean(s1s)
OUT = []
for sigma in (250., 450.):
    for f in (0.2, 0.5, 1.0):
        for sl in (50., 150.):
            a0, s0 = run(sigma, f, sl, 0.)
            ds = np.arange(0, 391, 15.); A = []
            for d in ds: A.append(run(sigma, f, sl, d, n=160))
            A = np.array(A); row = dict(sigma=sigma, f=f, sig_line=sl, s1_0=float(s0), a1_max=float(A[:, 0].max()))
            for tgt in (0.1, 0.2):
                k = np.flatnonzero(A[:, 0] >= tgt)
                if len(k): row[f'd_{tgt}'] = float(ds[k[0]]); row[f'ds1_{tgt}'] = float(A[k[0], 1] - s0)
            OUT.append(row); print(row)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_attractor.json'), 'w'), indent=1)
