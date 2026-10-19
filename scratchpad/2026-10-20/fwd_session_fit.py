"""Per-session fit of the forward model (X-view offset d + plain copy a) to det3 June-15 (2026-10-20).
The model's Y-share scale and pad-type purity differ from the data's, so bins are Y-share QUANTILES (sextiles) in
each, and only the class asymmetries are fitted (12 numbers per session). Cloud sigma is drawn per event; one global
sigma range is used, then scaled per session by matching the median Y share. Usage: fwd_session_fit.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_more as A
from scipy.interpolate import RegularGridInterpolator
QS = np.linspace(0, 1, 7)
# ---------------- data per session
rows = M.survey_rows()
SESS = {'6-22/23 1000 V': lambda r: r['build'] == 'det3 June-15' and 'det2_det3_overnight_6-22' in r['run'],
        '6-23/24 600 V': lambda r: r['build'] == 'det3 June-15' and '6-23-26' in r['run'],
        '6-27/28 1000 V': lambda r: r['build'] == 'det3 June-15' and r['feu'] in (7, 8) and r['drift'] == '1000',
        '6-27 drift 100-500 V': lambda r: r['build'] == 'det3 June-15' and r['feu'] in (7, 8) and r['drift'] in ('100', '200', '300', '400', '500'),
        'all det3 June-15': lambda r: r['build'] == 'det3 June-15'}
def stats(d, ysh, yp):
    e = np.quantile(ysh, QS); e[-1] += 1e-9; b = np.clip(np.digitize(ysh, e) - 1, 0, 5); r = []
    for k in range(6):
        for cl in (True, False):
            m = (b == k) & (yp == cl); r.append((d[m].mean() if m.sum() > 15 else np.nan, d[m].std() / np.sqrt(max(m.sum(), 1))))
    return np.array(r), float(np.median(ysh))
DATA = {}
for lab, sel in SESS.items():
    acc = {k: [] for k in ('d', 'ysh', 'yp')}
    for rx, ry in A.pairs(rows, sel):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        par = (dx['p'][ix] + dy['p'][iy]) % 2; ysh = (dy['al'][1] + dy['al'][-1])[iy]
        ok = (dx['p'][ix] % 64 >= 4) & (dx['p'][ix] % 64 <= 59)
        acc['d'].append(dx['a1'][ix][ok]); acc['ysh'].append(ysh[ok]); acc['yp'].append((par == int(ysh[par == 1].mean() < ysh[par == 0].mean()))[ok])
    if not acc['d']: continue
    a = {k: np.concatenate(v) for k, v in acc.items()}
    DATA[lab] = (*stats(a['d'], a['ysh'], a['yp']), len(a['d']))
# ---------------- model (as fwd_variants.py)
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'bus_par')
P = 780.; dd = 20.; EX = np.load(f'{B}/d20_nobus_esl.npy').astype(float)
nx, ny = EX.shape; nrow = int(round(ny * dd / P / 2))
xc = (np.arange(nx) + .5) * dd - P / 2; yc = (np.arange(ny) + .5) * dd - nrow * P
Ex = np.concatenate([EX[-10:], EX, EX[:10]]); xx = np.concatenate([xc[-10:] - 2 * P, xc, xc[:10] + 2 * P])
f_ = RegularGridInterpolator((xx, yc), Ex, bounds_error=False, fill_value=0.)
fX = lambda x, y: f_(np.stack([((x + P / 2) % (2 * P)) - P / 2, y], -1))
g = np.linspace(-3, 3, 13); GX, GY = np.meshgrid(g, g, indexing='ij'); W = np.exp(-.5 * (GX ** 2 + GY ** 2)); W /= W.sum()
K = np.arange(-4, 5); N = 20000; rng = np.random.default_rng(7)
def model(slo, shi, shifts, copies):
    cx = rng.uniform(-P, P, N); cy = rng.uniform(-P / 2, P / 2, N); s = rng.uniform(slo, shi, N)
    px = cx[:, None, None] + s[:, None, None] * GX[None]; py = cy[:, None, None] + s[:, None, None] * GY[None]
    SY = np.stack([(fX(py - P * c + P / 2, px + P / 2 - P * c) * W).sum((1, 2)) for c in K], 1)
    nY = rng.normal(1, .03, SY.shape); out = {}
    for dsh in shifts:
        SX0 = np.stack([(fX(px - P * k, py - P * k + dsh) * W).sum((1, 2)) for k in K], 1); nX = rng.normal(1, .03, SX0.shape)
        for a in copies:
            SX = SX0.copy(); SX[:, :-1] += a * SX0[:, 1:]; SX = SX * nX; SYn = SY * nY; ii = np.arange(N)
            kx = SX[:, 1:-1].argmax(1) + 1; ky = SYn[:, 1:-1].argmax(1) + 1
            d = (SX[ii, kx - 1] - SX[ii, kx + 1]) / SX[ii, kx]; ysh = (SYn[ii, ky - 1] + SYn[ii, ky + 1]) / SYn[ii, ky]
            par = (K[kx] + K[ky]) % 2; yp = par == (0 if ysh[par == 0].mean() < ysh[par == 1].mean() else 1)
            out[(dsh, a)] = stats(d, ysh, yp)
    return out
SH = [0, -50, -100, -150, -200, -250, -300]; CP = [0, .05, .1, .15, .2, .25, .3]
RANGES = [(80, 250), (120, 320), (150, 400), (200, 450), (250, 550)]
MODELS = {r: model(*r, SH, CP) for r in RANGES}
OUT = {}
for lab, (dat, med, n) in DATA.items():
    best = None
    for r, mods in MODELS.items():
        for (dsh, a), (t, mm) in mods.items():
            ok = np.isfinite(t[:, 0]) & np.isfinite(dat[:, 0])
            c = float((((t[ok, 0] - dat[ok, 0]) ** 2) / (t[ok, 1] ** 2 + dat[ok, 1] ** 2 + .02 ** 2)).sum())
            if best is None or c < best[0]: best = (c, r, dsh, a, t, mm)
    c, r, dsh, a, t, mm = best
    cs = [(cc, k) for k, (tt, _) in MODELS[r].items() for cc in [float((((tt[:, 0] - dat[:, 0]) ** 2) / (tt[:, 1] ** 2 + dat[:, 1] ** 2 + .02 ** 2))[np.isfinite(tt[:, 0]) & np.isfinite(dat[:, 0])].sum())]]
    near = sorted(k for cc, k in cs if cc < c + 4)
    pc = min((cc, k) for cc, k in cs if k[0] == 0)   # best pure copy, same sigma range
    print(f'{lab:22s} n={n:5d} median Ysh {med:.2f} | best: sigma {r} (model median Ysh {mm:.2f}), offset {dsh} um, copy {a:.2f}, chi2 {c:.0f}/12'
          f' | pure copy best {pc[1][1]:.2f}: chi2 {pc[0]:.0f} | dchi2<4: {near}')
    print('     data Ypad', ' '.join(f'{v:+.2f}' for v in dat[0::2, 0]), ' Xpad', ' '.join(f'{v:+.2f}' for v in dat[1::2, 0]))
    print('     fit  Ypad', ' '.join(f'{v:+.2f}' for v in t[0::2, 0]), ' Xpad', ' '.join(f'{v:+.2f}' for v in t[1::2, 0]))
    OUT[lab] = dict(n=n, median_ysh=med, chi2=c, sigma=r, offset=dsh, copy=a, pure_copy=dict(a=pc[1][1], chi2=pc[0]), near=near,
                    data=dat.tolist(), model=t.tolist())
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fwd_session_fit.json'), 'w'), default=float)
