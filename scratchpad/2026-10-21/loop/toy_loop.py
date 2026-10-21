"""Attractor loop vs a no-pull baseline (2026-10-21). The even-term exclusion compared s1 against delta = 0, but a1(delta) is
zero at delta = 0 (attractor over the pad centre: narrowest) AND at delta = P/2 (over the pad boundary: widest). Here the full
loop (a1(delta), s1(delta) - s1(f = 0)) over one period, for the audit toy's model (1-D; pad response box(P) * gauss(sig_ind);
gaussian cloud sigma, uniform positions, peak-pad selection; fraction f pulled to the nearest attractor, smeared by sig_line).
Usage: toy_loop.py"""
import json, os, numpy as np
from scipy.special import erf
P = 780.; x = np.arange(-6 * P, 6 * P, 5.); K = np.arange(-3, 4)
def resp(sig_ind):
    u = x[None, :] - K[:, None] * P; s = np.sqrt(2) * sig_ind
    return 0.5 * (erf((u + P / 2) / s) - erf((u - P / 2) / s))           # (7, nx)
def run(sigma, f, sig_line, delta, Rm, n=240):
    x0 = np.linspace(-P / 2, P / 2, n, endpoint=False)
    rho = (1 - f) * np.exp(-0.5 * ((x[None] - x0[:, None]) / sigma) ** 2) / sigma
    xs = np.round((x0 - delta) / P) * P + delta
    rho += f * np.exp(-0.5 * ((x[None] - xs[:, None]) / sig_line) ** 2) / sig_line
    q = rho @ Rm.T                                                       # (n, 7)
    c = np.argmax(q, 1); ok = (c > 0) & (c < 6); i = np.arange(n)[ok]; c = c[ok]
    a1 = (q[i, c + 1] - q[i, c - 1]) / q[i, c]; s1 = (q[i, c + 1] + q[i, c - 1]) / q[i, c]
    return float(a1.mean()), float(s1.mean())
OUT = []
D = np.arange(0, 781, 30.)
for sig_ind in (100., 200.):
    Rm = resp(sig_ind)
    for sigma in (300., 450., 600.):
        _, s0 = run(sigma, 0., 50., 0., Rm)
        for f in (0.2, 0.5, 1.0):
            for sl in (50., 200.):
                L = [run(sigma, f, sl, d, Rm) for d in D]
                a = np.array([l[0] for l in L]); s = np.array([l[1] for l in L]) - s0
                OUT.append(dict(sig_ind=sig_ind, sigma=sigma, f=f, sig_line=sl, s1_nopull=s0, delta=D.tolist(), a1=a.tolist(), ds1=s.tolist()))
                imax = int(np.argmax(a))
                # where |a1| >= 0.1: range of ds1
                big = np.abs(a) >= 0.1
                rng = f'[{s[big].min():+.3f}, {s[big].max():+.3f}]' if big.any() else '   (a1 < 0.1)   '
                print(f'ind {sig_ind:3.0f} sigma {sigma:3.0f} f {f:.1f} line {sl:3.0f}: s1(no pull) {s0:.3f} | ds1 at delta 0 {s[0]:+.3f}, '
                      f'at P/2 {s[13]:+.3f} | a1 max {a.max():+.3f} at delta {D[imax]:.0f} (ds1 {s[imax]:+.3f}) | ds1 where |a1|>=0.1 {rng}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_loop.json'), 'w'))
