"""Across-board profiles in 16-column blocks (true column order), per board (2026-10-20):
  a1   X asymmetry (alpha(+1) - alpha(-1))       -- strip-offset model: A sin(phi(c))
  s1X  X sharing (alpha(+1) + alpha(-1))         -- strip-offset model: B cos(phi(c)) on top of width/gain trends
  s1Y  Y sharing of the same matched events      -- control for gas/gain-driven width trends
  Edge test: a1 in the 8 columns either side of every connector boundary (a wiring flip would step there).
Fits a1(c) = A sin(2 pi (c - c0)/N) + b over N in 300..3000 columns (grid), reports best N, A, c0, b, chi2,
and the same with N -> infinity (constant). Usage: block_profile.py"""
import json, os, numpy as np
import r5lib as R
OUT = {}
NB = 16; cb = np.arange(0, 512, NB) + NB / 2
def blocks(col, v, nb=NB):
    k = col // nb; m = np.array([v[k == i].mean() if (k == i).sum() >= 30 else np.nan for i in range(512 // nb)])
    e = np.array([v[k == i].std() / np.sqrt((k == i).sum()) if (k == i).sum() >= 30 else np.nan for i in range(512 // nb)])
    return m, e
def fit(m, e):
    ok = np.isfinite(m); x, y, w = cb[ok], m[ok], 1 / np.maximum(e[ok], .005) ** 2; best = None
    for N in np.r_[np.arange(300, 1000, 20), np.arange(1000, 3001, 100)]:
        X = np.stack([np.sin(2 * np.pi * x / N), np.cos(2 * np.pi * x / N), np.ones_like(x)], 1)
        c, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None); chi = float((w * (y - X @ c) ** 2).sum())
        if best is None or chi < best[0]: best = (chi, N, c)
    c0 = float(((y * w).sum() / w.sum())); chic = float((w * (y - c0) ** 2).sum())
    chi, N, c = best; A = float(np.hypot(c[0], c[1])); ph = float(np.arctan2(c[1], c[0]))
    return dict(N=float(N), A=A, phase=ph, b=float(c[2]), chi=chi, ndf=int(ok.sum() - 4), chi_const=chic, ndf_const=int(ok.sum() - 1))
for b in ('det3 June-15', 'det4', 'det6', 'det7', 'det2', 'det3 May', 'det1'):
    E = R.ev(R.BUILDS[b]); col = 64 * ((E['px'] % 512) // 64) + 63 - (E['px'] % 512) % 64
    a1 = E['alx'][:, 5] - E['alx'][:, 3]; sx = E['alx'][:, 5] + E['alx'][:, 3]; sy = E['ysh']
    # X-only (unmatched) a1 for statistics
    p, ax = [], []
    for r in R.rows():
        if r['build'] != b or r['view'] != 'X': continue
        d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59); p.append(d['p'][ok] % 512); ax.append(d['a1'][ok])
    p, ax = map(np.concatenate, (p, ax)); colx = 64 * (p // 64) + 63 - p % 64
    ma, ea = blocks(colx, ax); msx, esx = blocks(col, sx); msy, esy = blocks(col, sy)
    F = fit(ma, ea); OUT[b] = dict(a1=ma.tolist(), a1_e=ea.tolist(), s1X=msx.tolist(), s1Y=msy.tolist(), fit=F)
    f = lambda a: ' '.join('  -  ' if np.isnan(x) else f'{x:+.2f}' for x in a)
    print(f'== {b}: a1 fit N={F["N"]:.0f} cols A={F["A"]:.3f} b={F["b"]:+.3f} chi2 {F["chi"]:.0f}/{F["ndf"]} (const {F["chi_const"]:.0f}/{F["ndf_const"]})')
    print(f'  a1  {f(ma)}'); print(f'  s1X {f(msx)}'); print(f'  s1Y {f(msy)}')
    # connector-edge steps: columns 56..63 vs 64..71 of each boundary
    st = []
    for B in range(64, 512, 64):
        l = (colx >= B - 8) & (colx < B); r_ = (colx >= B) & (colx < B + 8)
        if l.sum() > 30 and r_.sum() > 30: st.append(f'{B}:{ax[r_].mean() - ax[l].mean():+.2f}±{np.hypot(ax[r_].std()/np.sqrt(r_.sum()), ax[l].std()/np.sqrt(l.sum())):.2f}')
    print('  steps at connector boundaries (right - left, 8 cols each):', ' '.join(st))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'block_profile.json'), 'w'), indent=1)
