"""Audit toy 2 (2026-10-21): strip quantisation. Strips of width w at pitch P (= pad pitch), centres offset delta from
the pad centres. The part of the cloud landing on a strip is spread across that strip's width (fraction e equalised,
1-e left where it landed); gap charge stays. Pads read the 1-D charge through a box(P)*gauss(sig_ind) response.
Gives a1(delta) (odd) and s1(delta) - s1(0) (even). Prints the s1 change that comes with a1 = 0.1 and 0.2, and the
ratio ds1/a1^2 to compare with the data's per-cell slope (share_vs_a1.py). Usage: python toy_strips.py"""
import json, os, numpy as np
from scipy.special import erf
P = 780.; dx = 5.; x = np.arange(-6 * P, 6 * P, dx)
def R(u, s=150.): return 0.5 * (erf((u + P / 2) / (np.sqrt(2) * s)) - erf((u - P / 2) / (np.sqrt(2) * s)))
RK = np.array([R(x - k * P) for k in range(-3, 4)])
def run(sigma, w, e, delta, n=120):
    a1s, s1s = [], []
    rel = ((x - delta + P / 2) % P) - P / 2; on = np.abs(rel) < w / 2; sid = np.round((x - delta) / P)
    for x0 in np.linspace(-P / 2, P / 2, n, endpoint=False):
        rho = np.exp(-0.5 * ((x - x0) / sigma) ** 2)
        out = rho * (~on) + (1 - e) * rho * on
        q_on = rho * on
        for k in np.unique(sid[on]):
            m = on & (sid == k); out[m] += e * q_on[m].sum() / m.sum()
        q = RK @ out; c = np.argmax(q)
        if c in (0, 6): continue
        a1s.append((q[c + 1] - q[c - 1]) / q[c]); s1s.append((q[c + 1] + q[c - 1]) / q[c])
    return np.mean(a1s), np.mean(s1s)
OUT = []
for sigma in (250., 450.):
    for w in (540., 700.):
        for e in (0.3, 1.0):
            a0, s0 = run(sigma, w, e, 0.); ds = np.arange(0, 391, 15.)
            A = np.array([run(sigma, w, e, d) for d in ds])
            row = dict(sigma=sigma, w=w, e=e, s1_0=float(s0), a1_max=float(A[:, 0].max()), d_at_max=float(ds[A[:, 0].argmax()]))
            for tgt in (0.05, 0.1, 0.2):
                k = np.flatnonzero(A[:, 0] >= tgt)
                if len(k): row[f'ds1_{tgt}'] = float(A[k[0], 1] - s0); row[f'd_{tgt}'] = float(ds[k[0]])
            OUT.append(row); print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_strips.json'), 'w'), indent=1)
