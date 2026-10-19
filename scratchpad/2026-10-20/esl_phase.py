"""Copy vs ESL stripe phase, with TRUE pad columns (each 64-ch connector reversed), 2026-10-20.
ESL: 550 um stripes, 250 um gaps, 800 um pitch, along y; pads 780 um -> column k sits at stripe phase
phi_k = (k*780 + x0) mod 800, i.e. 20 um per column backwards; 40 columns = one beat. x0 unknown -> scan.
Model: copy = a + b * f(phi), f = fraction of a Gaussian cloud (sigma) over the stripe at the Y pad's x."""
import os, numpy as np
from scipy.special import erf
z = np.load(os.path.expanduser('~/x17/nb_asym_raw/esl_cache.npz'), allow_pickle=True)
col = lambda ch: 64 * (ch % 512 // 64) + 63 - ch % 64
a1 = z['ax'][:, 5] - z['ax'][:, 3]
def frac(xc, sig, x0):
    # fraction of N(xc, sig) on stripes centred at x0 + 800 m (width 550)
    t = 0
    for m in range(-3, 4):
        c = x0 + 800 * m; lo, hi = c - 275, c + 275
        t = t + 0.5 * (erf((hi - xc) / (np.sqrt(2) * sig)) - erf((lo - xc) / (np.sqrt(2) * sig)))
    return t
def fit(y, xc, sig):
    best = None
    for x0 in np.arange(0, 800, 10):
        f = frac(xc % 800, sig, x0); X = np.c_[np.ones_like(f), f - f.mean()]
        c, res, *_ = np.linalg.lstsq(X, y, rcond=None); r = ((y - X @ c) ** 2).sum()
        if best is None or r < best[0]: best = (r, x0, c)
    return best
for bu in ('det3 June-15', 'det4', 'det7', 'det3 May', 'det2', 'det6', 'det1'):
    for cls in ('ypad', 'xpad'):
        m = (z['build'] == bu) & (z['ypad'] if cls == 'ypad' else ~z['ypad'])
        m &= (z['px'] % 64 >= 4) & (z['px'] % 64 <= 59)
        if m.sum() < 200: continue
        k = col(z['px'][m]); ph = k % 40; y = a1[m]
        mu = np.array([y[ph == i].mean() for i in range(40)]); se = np.array([y[ph == i].std() / np.sqrt((ph == i).sum()) for i in range(40)])
        chi2 = (((mu - y.mean()) / se) ** 2).sum()
        # sine fit amplitude in phase
        w = 1 / se ** 2; th = 2 * np.pi * np.arange(40) / 40
        X = np.c_[np.ones(40), np.cos(th), np.sin(th)]; c = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * mu))
        cov = np.linalg.inv(X.T @ (w[:, None] * X)); amp = np.hypot(c[1], c[2])
        ea = np.sqrt((c[1] ** 2 * cov[1, 1] + c[2] ** 2 * cov[2, 2]) / amp ** 2)
        out = f'{bu:13s} {cls}  n={m.sum():5d}  mean={y.mean():+.3f}±{y.std()/np.sqrt(len(y)):.3f}  chi2/39={chi2/39:.2f}  sine amp={amp:.3f}±{ea:.3f}'
        for sig in (150, 300):
            r, x0, cc = fit(y, k * 780., sig)
            r0 = ((y - y.mean()) ** 2).sum(); dchi = (r0 - r) / y.var()
            out += f' | sig{sig}: b={cc[1]:+.3f} x0={x0} dchi2={dchi:.1f}'
        print(out)
