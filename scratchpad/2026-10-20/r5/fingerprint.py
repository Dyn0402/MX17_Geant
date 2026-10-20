"""One fingerprint for data sessions and sim arms (2026-10-20), so models are scored on everything at once.
Per sample of matched X/Y events (isolated, good pulse, not at a chip edge), in absolute Y-share bins ED:
  A_Y, A_X  : mean a1 = alpha(+1) - alpha(-1) in the Y-pad / X-pad class
  fY        : Y-pad fraction
  m1, p1, c : X(n-1)/SY, X(n+1)/SY, Xc/SY for Y-pad events (integrals; the bookkeeping)
plus marginals: a1 (all), coin = P(A(n+1) > A(n-1)), median Y share.
score(model, data): chi2 over all binned quantities with both errors and a floor per quantity (FLOOR).
Usage: imported (fp_data(sel), fp_arrays(Wx, Wy, px, py), score)."""
import numpy as np
import r5lib as R
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9])
FLOOR = dict(A_Y=.02, A_X=.02, fY=.03, m1=.01, p1=.01, c=.02)
def fp_arrays(Wx, Wy, px, py):
    Wx = Wx.astype(float); Wy = Wy.astype(float)
    cx, cy = Wx[:, 4], Wy[:, 4]; nx = np.einsum('ij,ij->i', cx, cx); ny = np.einsum('ij,ij->i', cy, cy)
    alx = lambda k: np.einsum('ij,ij->i', Wx[:, 4 + k], cx) / nx
    aly = lambda k: np.einsum('ij,ij->i', Wy[:, 4 + k], cy) / ny
    a1 = alx(1) - alx(-1); ysh = aly(1) + aly(-1); par = (px + py) % 2
    yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean())
    q = np.clip(Wx, 0, None).sum(2); SY = np.clip(Wy, 0, None).sum(2).sum(1)
    b = np.digitize(ysh, ED) - 1; nb = len(ED) - 1; o = {k: np.full((nb, 2), np.nan) for k in FLOOR}
    def put(k, j, v):
        if len(v) >= 25: o[k][j] = (v.mean(), v.std() / np.sqrt(len(v)))
    for j in range(nb):
        mb = b == j
        put('A_Y', j, a1[mb & yp]); put('A_X', j, a1[mb & ~yp])
        if mb.sum() >= 25: o['fY'][j] = (yp[mb].mean(), np.sqrt(yp[mb].mean() * (1 - yp[mb].mean()) / mb.sum()))
        m = mb & yp & (SY > 0)
        put('m1', j, q[m, 3] / SY[m]); put('p1', j, q[m, 5] / SY[m]); put('c', j, q[m, 4] / SY[m])
    A = Wx.max(2)
    o['marg'] = dict(n=int(len(a1)), a1=float(a1.mean()), coin=float((A[:, 5] > A[:, 3]).mean()), ysh_med=float(np.median(ysh)))
    return o
def fp_data(sel):
    E = R.ev(sel)
    return fp_arrays(E['Wx'], E['Wy'], E['px'], E['py'])
def score(m, d, keys=tuple(FLOOR)):
    c, n = 0., 0
    for k in keys:
        ok = np.isfinite(m[k][:, 0]) & np.isfinite(d[k][:, 0])
        c += float((((m[k][ok, 0] - d[k][ok, 0]) ** 2) / (m[k][ok, 1] ** 2 + d[k][ok, 1] ** 2 + FLOOR[k] ** 2)).sum()); n += int(ok.sum())
        # a bin where the data has a class and the model has none at all: penalise as a full miss
        miss = np.isnan(m[k][:, 0]) & np.isfinite(d[k][:, 0]) & (k in ('A_Y', 'A_X'))
        c += 4.0 * int(miss.sum()); n += int(miss.sum())
    return c, n
def show(name, o):
    f = lambda a: ' '.join('  -  ' if not np.isfinite(v) else f'{v:+.2f}' for v in a[:, 0])
    g = lambda a: ' '.join('  -  ' if not np.isfinite(v) else f'{v:.3f}' for v in a[:, 0])
    M = o['marg']
    print(f'{name:24s} n={M["n"]:6d} a1 {M["a1"]:+.3f} coin {M["coin"]:.3f} Ysh {M["ysh_med"]:.2f} | A_Y [{f(o["A_Y"])}] A_X [{f(o["A_X"])}] fY [{f(o["fY"])}]')
    print(f'{"":24s} Y-pad bookkeeping m1 [{g(o["m1"])}] p1 [{g(o["p1"])}] c [{g(o["c"])}]')
def score_delta(m, m0, d, d0, abs_keys=('A_Y', 'A_X'), delta_keys=('fY', 'm1', 'p1')):
    """chi2 of the model EFFECT against the data EFFECT: absolute for the class asymmetries (both ~0 on controls),
    differences (model - sim baseline m0) vs (data - control board d0) for fY and the bookkeeping, whose baselines
    differ between sim and data (sim X central fraction is flat in Y share, data's falls)."""
    c, n = score(m, d, abs_keys)
    for k in delta_keys:
        dm = m[k][:, 0] - m0[k][:, 0]; em = np.hypot(m[k][:, 1], m0[k][:, 1])
        dd = d[k][:, 0] - d0[k][:, 0]; ed = np.hypot(d[k][:, 1], d0[k][:, 1])
        ok = np.isfinite(dm) & np.isfinite(dd)
        c += float((((dm[ok] - dd[ok]) ** 2) / (em[ok] ** 2 + ed[ok] ** 2 + FLOOR[k] ** 2)).sum()); n += int(ok.sum())
    return c, n
