"""Forward model with the solved weighting fields (../2026-10-19/fwd_parity.py machinery), now through det3's full
Y-share-binned parity analysis, for candidate mechanisms (2026-10-20):
  none          symmetric geometry
  copy a        S_X[k-1] += a S_X[k]   (a plain channel copy, into the 'dn' side = data's n+1)
  shift d       X lines see the charge displaced by d um across the lines (toward k-1), Y lines don't
  ycopy a       S_X[k-1] += a * (charge-weight of the cloud over Y pads) * S_Y[peak]: pad-local Y->X coupling
Cloud sigma drawn uniformly 250..650 um per event (spans the data's Y-share range). Prints per Y-share bin the Y-pad /
X-pad class asymmetry (dn - up)/c and the Y-class fraction, next to det3 June-15. Usage: fwd_variants.py [--n 30000]"""
import os, sys, json, argparse, numpy as np
from scipy.interpolate import RegularGridInterpolator
ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=30000); a = ap.parse_args()
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'bus_par')
P = 780.; d = 20.; EX = np.load(f'{B}/d20_nobus_esl.npy').astype(float); EY0 = EX
nx, ny = EX.shape; nrow = int(round(ny * d / P / 2))
xc = (np.arange(nx) + .5) * d - P / 2; yc = (np.arange(ny) + .5) * d - nrow * P
def mk(E):
    Ex = np.concatenate([E[-10:], E, E[:10]]); xx = np.concatenate([xc[-10:] - 2 * P, xc, xc[:10] + 2 * P])
    f = RegularGridInterpolator((xx, yc), Ex, bounds_error=False, fill_value=0.)
    return lambda x, y: f(np.stack([((x + P / 2) % (2 * P)) - P / 2, y], -1))
fX = fY0 = mk(EX)
g = np.linspace(-3, 3, 13); GX, GY = np.meshgrid(g, g, indexing='ij'); W = np.exp(-.5 * (GX ** 2 + GY ** 2)); W /= W.sum()
K = np.arange(-4, 5)
def ypad(px, py):
    i = np.rint((px - P / 2) / P).astype(int); j = np.rint(py / P).astype(int)
    return (np.abs(px - (P / 2 + P * i)) <= 340) & (np.abs(py - P * j) <= 340) & ((i + j) % 2 == 1)
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9])
DATA = dict(ypad=[.26, .27, .25, .24, .19, .09, .00], xpad=[-.40, -.29, -.12, -.00, .05, .16, .20], fy=[.87, .74, .67, .59, .48, .38, .13])
rng = np.random.default_rng(5)
cx = rng.uniform(-P, P, a.n); cy = rng.uniform(-P / 2, P / 2, a.n); s = rng.uniform(250, 650, a.n)
px = cx[:, None, None] + s[:, None, None] * GX[None]; py = cy[:, None, None] + s[:, None, None] * GY[None]
SY = np.stack([(fY0(py - P * c + P / 2, px + P / 2 - P * c) * W).sum((1, 2)) for c in K], 1)
wY = (ypad(px, py) * W).sum((1, 2))                                     # cloud weight over Y pads
SXcache = {}
def SXf(dsh):
    if dsh not in SXcache:
        SXcache[dsh] = np.stack([(fX(px - P * k, py - P * k + dsh) * W).sum((1, 2)) for k in K], 1)
    return SXcache[dsh]
def analyse(SX):
    ii = np.arange(a.n); SXn = SX * rng.normal(1, .03, SX.shape); SYn = SY * rng.normal(1, .03, SY.shape)
    kx = SXn[:, 1:-1].argmax(1) + 1; ky = SYn[:, 1:-1].argmax(1) + 1
    dd = (SXn[ii, kx - 1] - SXn[ii, kx + 1]) / SXn[ii, kx]
    ysh = (SYn[ii, ky - 1] + SYn[ii, ky + 1]) / SYn[ii, ky]; par = (K[kx] + K[ky]) % 2
    yp = par == (0 if ysh[par == 0].mean() < ysh[par == 1].mean() else 1)
    b = np.digitize(ysh, ED) - 1; r = dict(ypad=[], xpad=[], fy=[], n=[])
    for k in range(len(ED) - 1):
        m = b == k; r['n'].append(int(m.sum()))
        r['ypad'].append(float(dd[m & yp].mean()) if (m & yp).sum() > 30 else np.nan)
        r['xpad'].append(float(dd[m & ~yp].mean()) if (m & ~yp).sum() > 30 else np.nan)
        r['fy'].append(float(yp[m].mean()) if m.sum() > 30 else np.nan)
    r['cls'] = [float(dd[yp].mean()), float(dd[~yp].mean())]
    r['chi'] = float(np.nansum([(np.array(r[k]) - np.array(DATA[k])) ** 2 / .02 ** 2 for k in ('ypad', 'xpad', 'fy')]))
    return r
fmt = lambda v: ' '.join('  -  ' if not np.isfinite(x) else f'{x:+.2f}' for x in v)
print(f'{"DATA det3 June-15":14s} Ypad [{fmt(DATA["ypad"])}] Xpad [{fmt(DATA["xpad"])}] fY [{fmt(DATA["fy"])}]')
OUT = []
for kind, vals in (('none', [0]), ('copy', [.3]), ('shift', [-200]), ('combo', [(-150, .1), (-200, .1), (-250, .05), (-250, .1), (-300, .0)])):
    for v in vals:
        if kind == 'shift': SX = SXf(float(v))
        elif kind == 'combo':
            S0 = SXf(float(v[0])); SX = S0.copy(); SX[:, :-1] += v[1] * S0[:, 1:]
        else:
            SX = SXf(0.).copy(); ii = np.arange(a.n)
            if kind == 'copy': SX[:, :-1] += v * SXf(0.)[:, 1:]
            if kind == 'ycopy':
                # the coupling feeds line k-1 from the Y charge sitting next to line k: every X line k gets wY-weighted share
                # of the Y signal of its own pad column, i.e. into k-1 from the cloud over the Y pads adjacent to line k
                kx0 = SXf(0.)[:, 1:-1].argmax(1) + 1; ky0 = SY[:, 1:-1].argmax(1) + 1
                SX[ii, kx0 - 1] += v * wY * SY[ii, ky0]
        r = analyse(SX); r.update(kind=kind, v=list(v) if isinstance(v, tuple) else v); OUT.append(r)
        print(f'{kind+" "+str(v):14s} Ypad [{fmt(r["ypad"])}] Xpad [{fmt(r["xpad"])}] fY [{fmt(r["fy"])}] cls {r["cls"][0]:+.3f}/{r["cls"][1]:+.3f} chi {r["chi"]:.0f}')
print('n per bin (last):', OUT[-1]['n'])
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fwd_variants.json'), 'w'))
