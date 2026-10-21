"""Y-view neighbour asymmetry split by pad class (2026-10-21). The X-view field shows up on every June build and is
reported flat in Y, but the Y 2D maps averaged both pad classes. Here: a1_Y = alpha_Y(+1) - alpha_Y(-1) (DAQ order),
per build x pad class x share bin, and 4x4 maps per class with a constant-model chi2, next to the same for a1_X.
Usage: NB_RAW=~/x17/nb_asym_raw python3 yview_class.py"""
import json, os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
ED = np.array([0, .4, .55, .7, .85, 1, 1.3])
def ms(v): return (float(v.mean()), float(v.std() / np.sqrt(len(v))), int(len(v))) if len(v) >= 25 else None
def tcol(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
def cmap(a, cx, cy, m):
    G = []; chi = 0; nd = 0; cells = []
    for i in range(4):
        for j in range(4):
            k = m & (cx // 128 == i) & (cy // 128 == j)
            if k.sum() >= 30: cells.append((a[k].mean(), a[k].std() / np.sqrt(k.sum()), i, j))
    if len(cells) < 3: return None
    v = np.array([c[0] for c in cells]); e = np.array([c[1] for c in cells]); w = 1 / e**2
    mu = (v * w).sum() / w.sum(); chi = float((((v - mu) / e)**2).sum())
    rms = float(np.sqrt(max(np.var(v) - np.mean(e**2), 0)))
    return dict(chi2=chi, ndf=len(cells) - 1, mean=float(mu), field_rms=rms, cells=[list(map(float, c)) for c in cells])
OUT = {}
for b, sel in R.BUILDS.items():
    E = R.ev(sel)
    a1x = E['alx'][:, 5] - E['alx'][:, 3]; a1y = E['aly'][:, 5] - E['aly'][:, 3]
    ysh = E['ysh']; xsh = E['alx'][:, 3] + E['alx'][:, 5]; yp = E['ypad']
    cx, cy = tcol(E['px']), tcol(E['py'])
    o = dict(n=int(len(a1x)))
    for cls, m in (('Ypad', yp), ('Xpad', ~yp), ('all', np.ones_like(yp))):
        o[cls] = dict(a1x=ms(a1x[m]), a1y=ms(a1y[m]),
                      a1y_by_ysh=[ms(a1y[m & (np.digitize(ysh, ED) - 1 == j)]) for j in range(len(ED) - 1)],
                      a1y_by_xsh=[ms(a1y[m & (np.digitize(xsh, ED) - 1 == j)]) for j in range(len(ED) - 1)],
                      a1x_by_ysh=[ms(a1x[m & (np.digitize(ysh, ED) - 1 == j)]) for j in range(len(ED) - 1)],
                      map_x=cmap(a1x, cx, cy, m), map_y=cmap(a1y, cx, cy, m))
    OUT[b] = o
    f = lambda t: '   -   ' if t is None else f'{t[0]:+.3f}±{t[1]:.3f}'
    print(f'== {b}  n={o["n"]}')
    for cls in ('Ypad', 'Xpad', 'all'):
        p = o[cls]; mx, my = p['map_x'], p['map_y']
        mm = lambda q: '-' if q is None else f'{q["chi2"]:.0f}/{q["ndf"]} rms {q["field_rms"]:.3f}'
        print(f'  {cls:4s} a1X {f(p["a1x"])}  a1Y {f(p["a1y"])} | map X {mm(mx)} | map Y {mm(my)}')
        print(f'        a1Y by Ysh [{" ".join(f(t) for t in p["a1y_by_ysh"])}]')
        print(f'        a1Y by Xsh [{" ".join(f(t) for t in p["a1y_by_xsh"])}]')
        print(f'        a1X by Ysh [{" ".join(f(t) for t in p["a1x_by_ysh"])}]')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'yview_class.json'), 'w'), indent=1)
