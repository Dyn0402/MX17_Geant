"""Transfer shape of the attractor toy (2026-10-21): change of alpha(k), k = -2..2, relative to no pull, at offsets where
a1 is 0.10-0.20 and |ds1| <= 0.03 (the data's operating point: det3 June-15 vs May), per unit a1. Data: det2 per unit a1
odd parts -0.54/+0.47 (+-1), -0.08/+0.04 (+-2); det3 June-15 - May odd ~0.5 (+-1), ~0.11 (+-2).
Usage: toy_shape.py"""
import json, os, numpy as np
from toy_loop import resp, x, P
def prof(sigma, f, sig_line, delta, Rm, n=240):
    x0 = np.linspace(-P / 2, P / 2, n, endpoint=False)
    rho = (1 - f) * np.exp(-0.5 * ((x[None] - x0[:, None]) / sigma) ** 2) / sigma
    xs = np.round((x0 - delta) / P) * P + delta
    rho += f * np.exp(-0.5 * ((x[None] - xs[:, None]) / sig_line) ** 2) / sig_line
    q = rho @ Rm.T; c = np.argmax(q, 1); ok = (c >= 2) & (c <= 4); i = np.arange(n)[ok]; c = c[ok]
    return np.array([(q[i, c + k] / q[i, c]).mean() for k in (-2, -1, 1, 2)])
OUT = []
for sig_ind in (100., 200.):
    Rm = __import__('toy_loop').resp(sig_ind)
    for sigma in (300., 450., 600.):
        p0 = prof(sigma, 0., 50., 0., Rm)
        for f in (0.1, 0.2, 0.3, 0.5):
            for sl in (50., 200.):
                for d in np.arange(0, 391, 15.):
                    p = prof(sigma, f, sl, d, Rm) - p0; a1 = p[2] - p[1]; ds1 = p[2] + p[1]
                    if 0.10 <= a1 <= 0.20 and abs(ds1) <= 0.03:
                        o1 = (p[2] - p[1]) / 2; o2 = (p[3] - p[0]) / 2
                        OUT.append(dict(sig_ind=sig_ind, sigma=sigma, f=f, sig_line=sl, delta=float(d), a1=float(a1), ds1=float(ds1), odd2_over_odd1=float(o2 / o1),
                                        even2=float((p[3] + p[0]) / 2 / a1)))
for r in OUT:
    print(f'ind {r["sig_ind"]:3.0f} sigma {r["sigma"]:3.0f} f {r["f"]:.1f} line {r["sig_line"]:3.0f} delta {r["delta"]:5.0f}: a1 {r["a1"]:+.3f} ds1 {r["ds1"]:+.3f} | '
          f'odd(+-2)/odd(+-1) {r["odd2_over_odd1"]:+.2f}  even(+-2) per a1 {r["even2"]:+.2f}')
v = np.array([r['odd2_over_odd1'] for r in OUT]); print(f'{len(OUT)} operating points: odd2/odd1 median {np.median(v):+.2f}, range {v.min():+.2f} .. {v.max():+.2f}  (data: det2 0.12, det3 ~0.22; linear antisym ~0.4)')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_shape.json'), 'w'), indent=1)
# representative operating point for the deck: per-unit-a1 change of alpha(k), k = -3..3 (k = 0 omitted)
def prof7(sigma, f, sig_line, delta, Rm, n=240):
    x0 = np.linspace(-P / 2, P / 2, n, endpoint=False)
    rho = (1 - f) * np.exp(-0.5 * ((x[None] - x0[:, None]) / sigma) ** 2) / sigma
    xs = np.round((x0 - delta) / P) * P + delta
    rho += f * np.exp(-0.5 * ((x[None] - xs[:, None]) / sig_line) ** 2) / sig_line
    q = rho @ Rm.T; c = np.argmax(q, 1); ok = c == 3; i = np.arange(n)[ok]
    return np.array([(q[i, 3 + k] / q[i, 3]).mean() for k in (-3, -2, -1, 1, 2, 3)])
Rm = resp(200.); REP = dict(sigma=450., f=0.3, sig_line=200., delta=300., sig_ind=200.)
d = prof(REP['sigma'], REP['f'], REP['sig_line'], REP['delta'], Rm) - prof(REP['sigma'], 0., 50., 0., Rm); a1 = d[2] - d[1]
REP.update(a1=float(a1), ds1=float(d[2] + d[1]), per_a1={'-2': d[0] / a1, '-1': d[1] / a1, '1': d[2] / a1, '2': d[3] / a1})
print('representative', REP)
json.dump(REP, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_shape_rep.json'), 'w'), indent=1)
