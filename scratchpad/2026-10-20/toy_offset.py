"""SUPERSEDED: collection-only toy, its Y share does not separate pad types like the data; see fwd_variants.py.
Toy forward model (2026-10-20): copy vs X-view offset, through det3's own parity / Y-share analysis.
Checkerboard 0.68 mm pads on 0.78 mm, X pad at (i+j) even. Cloud = 2D Gaussian (sigma drawn per event to span the data's
Y-share range), charge fraction per pad (collection only). Variants:
  copy  a  : X[i-1] += a X[i]   (DAQ n+1 = lower column)
  shift d  : the X view sees the cloud displaced by -d in x (toward the lower column), the Y view does not
  ycopy a  : X[i-1] += a * (charge on the Y pad at (i, j)), only for charge on Y pads (pad-local coupling)
Prints, per Y-share bin: Y-pad / X-pad class asymmetry (X[lower]-X[upper])/X[c] and the Y-class fraction.
Usage: toy_offset.py"""
import numpy as np
from scipy.special import erf
P, H = 780., 340.; rng = np.random.default_rng(3); N = 60000; K = np.arange(-5, 6)
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9])
DATA = dict(ypad=[.26, .27, .25, .24, .19, .09, .00], xpad=[-.40, -.29, -.12, -.00, .05, .16, .20], fy=[.87, .74, .67, .59, .48, .38, .13])
def f1(c, s, x0):  # fraction of N(x0, s) in [c-H, c+H]
    return .5 * (erf((c + H - x0) / (np.sqrt(2) * s)) - erf((c - H - x0) / (np.sqrt(2) * s)))
def sim(kind=None, v=0.):
    x0 = rng.uniform(-P / 2, P / 2, N); y0 = rng.uniform(-P / 2, P / 2, N); s = rng.uniform(150, 650, N)
    xd = x0 - (v if kind == 'shift' else 0.)
    FXx = np.stack([f1(P * i, s, xd) for i in K], 1); FYx = np.stack([f1(P * i, s, x0) for i in K], 1)
    Fy = np.stack([f1(P * j, s, y0) for j in K], 1)
    pad = FYx[:, :, None] * Fy[:, None, :]; padX = FXx[:, :, None] * Fy[:, None, :]
    ev = ((K[:, None] + K[None, :]) % 2 == 0)
    X = (padX * ev).sum(2); Y = (pad * ~ev).sum(1)          # X channel i, Y channel j
    X *= rng.normal(1, .03, X.shape); Y *= rng.normal(1, .03, Y.shape)
    i0 = 5 + np.rint(x0 / P).astype(int); j0 = 5 + np.rint(y0 / P).astype(int); ii = np.arange(N)
    if kind == 'copy':
        X2 = X.copy(); X2[:, :-1] += v * X[:, 1:]; X = X2
    if kind == 'ycopy':
        onY = ((K[i0 - 0] + K[j0]) % 2 == 1) & (np.abs(x0 - P * K[i0]) < H) & (np.abs(y0 - P * K[j0]) < H)
        q = pad[ii, i0, j0] * onY; X = X.copy(); X[ii, i0 - 1] += v * q
    px = 1 + X[:, 1:-1].argmax(1); py = 1 + Y[:, 1:-1].argmax(1)
    c = X[ii, px]; d = (X[ii, px - 1] - X[ii, px + 1]) / c
    ysh = (Y[ii, py - 1] + Y[ii, py + 1]) / Y[ii, py]; yp = (K[px] + K[py]) % 2 == 1
    b = np.digitize(ysh, ED) - 1; r = dict(ypad=[], xpad=[], fy=[])
    for k in range(len(ED) - 1):
        m = b == k
        r['ypad'].append(d[m & yp].mean() if (m & yp).sum() > 30 else np.nan)
        r['xpad'].append(d[m & ~yp].mean() if (m & ~yp).sum() > 30 else np.nan); r['fy'].append(yp[m].mean() if m.sum() else np.nan)
    r['all_y'] = d[yp].mean(); r['all_x'] = d[~yp].mean()
    return r
fmt = lambda v: ' '.join('  -  ' if not np.isfinite(x) else f'{x:+.2f}' for x in v)
print(f'{"DATA det3 June-15":22s} Ypad [{fmt(DATA["ypad"])}]  Xpad [{fmt(DATA["xpad"])}]  fY [{fmt(DATA["fy"])}]')
for kind, vals in ((None, [0]), ('copy', [.2, .35]), ('shift', [50, 100, 150]), ('ycopy', [.3, .6, 1.0])):
    for v in vals:
        r = sim(kind, v)
        print(f'{str(kind)+" "+str(v):22s} Ypad [{fmt(r["ypad"])}]  Xpad [{fmt(r["xpad"])}]  fY [{fmt(r["fy"])}]  (class means {r["all_y"]:+.3f} / {r["all_x"]:+.3f})')
