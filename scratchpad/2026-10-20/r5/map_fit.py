"""Registration-field fit (2026-10-20). Model: a1(x, y) = A sin(phi0 + kx x + ky y) + b (a printed periodic layer with
pitch ~ the pad pitch whose offset from the pad columns varies linearly over the board: scale error kx, rotation ky).
x, y in units of columns/rows (0..511) from the true X column / Y row of matched events; X-only events add x-only info
but need y, so matched events only. 8 x 8 cells; grid search over (kx, ky), linear in (A cos, A sin, b).
Per det3 session as well: is the ridge in the same place (fixed geometry) with only the amplitude changing?
Usage: map_fit.py"""
import json, os, numpy as np
import r5lib as R
from eval_mix import SESS
NC = 8
def cells(E):
    px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64; a1 = E['alx'][:, 5] - E['alx'][:, 3]
    i = cx * NC // 512; j = cy * NC // 512; G = []
    for a in range(NC):
        for b in range(NC):
            m = (i == a) & (j == b)
            if m.sum() >= 25: G.append(((a + .5) * 512 / NC, (b + .5) * 512 / NC, a1[m].mean(), a1[m].std() / np.sqrt(m.sum()), m.sum()))
    return np.array(G)
def fit(G, fixk=None):
    x, y, v, e = G[:, 0], G[:, 1], G[:, 2], np.maximum(G[:, 3], .01); w = 1 / e ** 2; best = None
    ks = np.linspace(-3 * np.pi / 256, 3 * np.pi / 256, 61)  # up to 3 periods over the board
    grid = [fixk] if fixk is not None else [(kx, ky) for kx in ks for ky in ks]
    for kx, ky in grid:
        ph = kx * x + ky * y; X = np.stack([np.sin(ph), np.cos(ph), np.ones_like(ph)], 1)
        c, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], v * np.sqrt(w), rcond=None); chi = float((w * (v - X @ c) ** 2).sum())
        if best is None or chi < best[0]: best = (chi, kx, ky, c)
    chi, kx, ky, c = best; m0 = (v * w).sum() / w.sum(); chi0 = float((w * (v - m0) ** 2).sum())
    return dict(chi=chi, ndf=len(v) - 5, chi_const=chi0, kx=float(kx), ky=float(ky), A=float(np.hypot(c[0], c[1])),
                phi0=float(np.arctan2(c[1], c[0])), b=float(c[2]), n=len(v))
OUT = {}
for b in ('det3 June-15', 'det4', 'det6', 'det7', 'det2', 'det3 May', 'det1'):
    G = cells(R.ev(R.BUILDS[b])); F = fit(G); OUT[b] = F
    print(f'{b:13s} cells {F["n"]:2d}: chi2 {F["chi"]:6.1f}/{F["ndf"]} (const {F["chi_const"]:6.1f}/{F["n"]-1}) | A {F["A"]:.3f} b {F["b"]:+.3f} '
          f'| kx {F["kx"]*512/np.pi:+.2f} pi, ky {F["ky"]*512/np.pi:+.2f} pi per board width | phi0 {np.degrees(F["phi0"]):+.0f}')
ref = OUT['det3 June-15']
for s in ('6-22/23 1000V', '6-23/24 600V', '6-27/28 1000V'):
    G = cells(R.ev(SESS[s])); F = fit(G); Ff = fit(G, fixk=(ref['kx'], ref['ky'])); OUT['sess ' + s] = dict(free=F, fixed=Ff)
    print(f'det3 {s:14s} free: chi2 {F["chi"]:5.1f}/{F["ndf"]} A {F["A"]:.3f} b {F["b"]:+.3f} kx {F["kx"]*512/np.pi:+.2f}pi ky {F["ky"]*512/np.pi:+.2f}pi phi0 {np.degrees(F["phi0"]):+.0f}'
          f' | det3-pooled (kx,ky): chi2 {Ff["chi"]:5.1f}/{Ff["ndf"]+2} A {Ff["A"]:.3f} b {Ff["b"]:+.3f} phi0 {np.degrees(Ff["phi0"]):+.0f} | const {F["chi_const"]:.1f}/{F["n"]-1}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'map_fit.json'), 'w'), indent=1)
