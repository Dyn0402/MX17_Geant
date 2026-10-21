"""Per-electron vs cloud-level pull (2026-10-21). toy_loop.py moves a fraction f of the WHOLE cloud to the attractor nearest the
cloud centre: a non-local, cloud-level pull. Physical pulls (funnel, spreading on a strip, field focusing) act per electron:
each electron, with probability f, moves to the attractor nearest where IT lands. Here both versions side by side, same
parameters; operating points with |ds1| <= 0.03 vs no pull. Usage: toy_local.py"""
import json, os, numpy as np
from toy_loop import resp, x, P, run as run_cloud
def run_local(sigma, f, sig_line, delta, Rm, n=240):
    x0 = np.linspace(-P / 2, P / 2, n, endpoint=False)
    rho = np.exp(-0.5 * ((x[None] - x0[:, None]) / sigma) ** 2) / sigma          # (n, nx)
    kk = np.round((x - delta) / P).astype(int); lines = np.unique(kk)
    pulled = np.zeros_like(rho)
    for k in lines:
        m = kk == k; mass = rho[:, m].sum(1)                                  # cloud mass landing in this attractor's cell
        pulled += mass[:, None] * np.exp(-0.5 * ((x[None] - (k * P + delta)) / sig_line) ** 2) / (sig_line * np.sqrt(2 * np.pi)) * 5.
    rho2 = (1 - f) * rho + f * pulled
    q = rho2 @ Rm.T; c = np.argmax(q, 1); ok = (c > 0) & (c < 6); i = np.arange(n)[ok]; c = c[ok]
    return float(((q[i, c + 1] - q[i, c - 1]) / q[i, c]).mean()), float(((q[i, c + 1] + q[i, c - 1]) / q[i, c]).mean())
OUT = []
for sig_ind in (100., 200.):
    Rm = resp(sig_ind)
    for sigma in (300., 450.):
        _, s0 = run_cloud(sigma, 0., 50., 0., Rm)
        for f in (0.2, 0.3, 0.5, 1.0):
            for sl in (50., 200.):
                D = np.arange(0, 391, 30.)
                for name, fn in (('cloud', run_cloud), ('local', run_local)):
                    L = [fn(sigma, f, sl, d, Rm) for d in D]; a = np.array([l[0] for l in L]); s = np.array([l[1] for l in L]) - s0
                    ok = np.abs(s) <= 0.03; best = float(a[ok].max()) if ok.any() else float('nan')
                    OUT.append(dict(kind=name, sig_ind=sig_ind, sigma=sigma, f=f, sig_line=sl, a1=a.tolist(), ds1=s.tolist(), best=best))
                    print(f'{name:5s} ind {sig_ind:3.0f} sigma {sigma:3.0f} f {f:.1f} line {sl:3.0f}: a1 max {a.max():+.3f}; best a1 with |ds1|<=0.03: {best:+.3f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_local.json'), 'w'), indent=1)
