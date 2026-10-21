"""Upper bound for an electrostatic pull (2026-10-21). Charged insulating gaps saturate at sigma = eps0 E0 (~3.5 nC/cm^2 at
39 kV/cm); then no field line ends on a gap and every avalanche electron over a gap lands on the nearest strip edge (a full
funnel). Optionally the strip charge also equalises across the strip width (resistive spreading, fraction eq). 1-D toy as
toy_loop.py: strips of width w at the pad pitch (780 um) offset by delta; a1 and ds1 vs no strips, over a half period.
Usage: toy_funnel.py"""
import json, os, numpy as np
from toy_loop import resp, x, P
def run(sigma, w, delta, Rm, eq=0.0, n=240):
    g = P - w; x0s = np.linspace(-P / 2, P / 2, n, endpoint=False)
    u = ((x - delta + P / 2) % P) - P / 2                 # position relative to the nearest strip centre
    k = np.round((x - delta) / P)                         # strip index
    xe = np.where(np.abs(u) > w / 2, delta + k * P + np.sign(u) * w / 2, x)   # gap -> nearest edge
    a1s, s1s = [], []
    for x0 in x0s:
        rho = np.exp(-0.5 * ((x - x0) / sigma) ** 2)
        # funnel: move each bin's charge to xe; then a fraction eq of each strip's charge spread uniformly over the strip
        xc = xe
        q = np.zeros(7)
        for j, kk in enumerate(range(-3, 4)):
            q[j] = (rho * Rm[j](xc)).sum() if callable(Rm[j]) else 0
        if eq > 0:
            cen = delta + k * P; q2 = np.zeros(7)
            for j in range(7):
                # uniform over the strip: average of the pad response across the strip width
                ss = np.linspace(-w / 2, w / 2, 21)
                avg = np.mean([Rm[j](cen + s) for s in ss], 0)
                q2[j] = (rho * avg).sum()
            q = (1 - eq) * q + eq * q2
        c = int(np.argmax(q))
        if c in (0, 6): continue
        a1s.append((q[c + 1] - q[c - 1]) / q[c]); s1s.append((q[c + 1] + q[c - 1]) / q[c])
    return float(np.mean(a1s)), float(np.mean(s1s))
from scipy.special import erf
def respf(sig_ind):
    s = np.sqrt(2) * sig_ind
    return [lambda xx, K=K: 0.5 * (erf((xx - K * P + P / 2) / s) - erf((xx - K * P - P / 2) / s)) for K in range(-3, 4)]
OUT = []
for sig_ind in (100., 200.):
    Rm = respf(sig_ind)
    for sigma in (300., 450.):
        b0 = run(sigma, P, 0., Rm)        # w = P: no gaps, no funnel = no-pull baseline
        for w in (540., 700.):
            for eq in (0., 0.5, 1.0):
                L = [run(sigma, w, d, Rm, eq) for d in np.arange(0, 391, 30.)]
                a = np.array([l[0] for l in L]); s = np.array([l[1] for l in L]) - b0[1]
                OUT.append(dict(sig_ind=sig_ind, sigma=sigma, w=w, eq=eq, a1=a.tolist(), ds1=s.tolist()))
                print(f'ind {sig_ind:3.0f} sigma {sigma:3.0f} strip {w:3.0f} um, equalise {eq:.1f}: max a1 {a.max():+.3f}  (ds1 there {s[np.argmax(a)]:+.3f}; range {s.min():+.3f}..{s.max():+.3f})')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_funnel.json'), 'w'), indent=1)
