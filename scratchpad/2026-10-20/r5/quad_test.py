"""Quadrature test of the print-registration model (2026-10-20). With offset phase phi(x, y) = phi0 + kx x + ky y from
the a1 map fit (map_fit.py), the model predicts a1 ~ A sin(phi) and an X-sharing change ~ B cos(phi), the Y view
unaffected. Per 8 x 8 cell: s1X, s1Y (matched events). Fit s1X = c0 + c1 s1Y + B cos(phi) + C sin(phi); report B, C with
errors (cells weighted by their s1X error). B != 0 at fixed phase = the sharing follows the same phase as the asymmetry.
Null: the same fit with phi from a random (kx, ky, phi0) (200 draws) -> the spread of B. Usage: quad_test.py"""
import json, os, numpy as np
import r5lib as R
from map_fit import NC
F = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'map_fit.json')))
rng = np.random.default_rng(3); OUT = {}
for b in ('det3 June-15', 'det2', 'det6', 'det4'):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    sx = E['alx'][:, 5] + E['alx'][:, 3]; sy = E['ysh']; i = cx * NC // 512; j = cy * NC // 512; C = []
    for a in range(NC):
        for c in range(NC):
            m = (i == a) & (j == c)
            if m.sum() >= 25: C.append(((a + .5) * 512 / NC, (c + .5) * 512 / NC, sx[m].mean(), sx[m].std() / np.sqrt(m.sum()), sy[m].mean()))
    C = np.array(C); x, y, v, e, s = C.T; w = 1 / np.maximum(e, .01) ** 2
    def reg(ph):
        X = np.stack([np.ones_like(v), s, np.cos(ph), np.sin(ph)], 1); Xw = X * np.sqrt(w)[:, None]
        c, *_ = np.linalg.lstsq(Xw, v * np.sqrt(w), rcond=None); cov = np.linalg.inv(Xw.T @ Xw)
        r = v - X @ c; s2 = (w * r ** 2).sum() / (len(v) - 4)
        return c, np.sqrt(np.diag(cov) * max(s2, 1)), float((w * r ** 2).sum())
    f = F[b]; ph = f['phi0'] + f['kx'] * x + f['ky'] * y
    c, ce, chi = reg(ph)
    null = []
    for _ in range(200):
        phr = rng.uniform(-np.pi, np.pi) + rng.uniform(-3, 3) * np.pi / 256 * x + rng.uniform(-3, 3) * np.pi / 256 * y
        null.append(np.hypot(*reg(phr)[0][2:]))
    amp = np.hypot(c[2], c[3])
    OUT[b] = dict(B=float(c[2]), Be=float(ce[2]), C=float(c[3]), Ce=float(ce[3]), amp=float(amp), p=float((np.array(null) >= amp).mean()))
    print(f'{b:13s} cells {len(v)}: s1X = {c[0]:+.2f} + {c[1]:.2f} s1Y + ({c[2]:+.3f}±{ce[2]:.3f}) cos phi + ({c[3]:+.3f}±{ce[3]:.3f}) sin phi'
          f' | quadrature amplitude {amp:.3f}, random-phase p = {OUT[b]["p"]:.2f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'quad_test.json'), 'w'), indent=1)
