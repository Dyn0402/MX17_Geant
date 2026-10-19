"""Where on the Y pad does the copy come from? (2026-10-20). det3 June-15 (and det4/det7) Y-pad-class events.
1. copy vs y position inside the Y pad (eta_y = alpha_y(+1) - alpha_y(-1), Y view is copy-free)
2. copy vs how centred on the Y pad (A_x/A_y) and vs X +-1 sum (|x| offset proxy)
3. shape: the excess waveform X(n+1) - X(n-1) against the X central and the Y central (fit e = a X_c + b Y_c + c dY_c/dt)
Writes ypad_probe.json next to this file."""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L, asym_more as A
rows = M.survey_rows(); OUT = {}
def bins(x, y, edges):
    r = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        r.append([float((lo + hi) / 2), float(y[m].mean()) if m.sum() > 20 else None, float(y[m].std() / np.sqrt(max(m.sum(), 1))), int(m.sum())])
    return r
for bu in ('det3 June-15', 'det4', 'det7', 'det3 May', 'det2'):
    E = {k: [] for k in ('a1', 'ey', 'axay', 'xsh', 'ysh', 'Wx', 'Wy', 'px', 'py')}
    for rx, ry in A.pairs(rows, lambda r: r['build'] == bu):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        par = (dx['p'][ix] + dy['p'][iy]) % 2; ysh = (dy['al'][1] + dy['al'][-1])[iy]
        yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean())
        ix, iy = ix[yp], iy[yp]
        E['a1'].append(dx['a1'][ix]); E['ey'].append((dy['al'][1] - dy['al'][-1])[iy])
        E['axay'].append(dx['A'][ix][:, 4] / dy['A'][iy][:, 4]); E['xsh'].append((dx['al'][1] + dx['al'][-1])[ix])
        E['ysh'].append((dy['al'][1] + dy['al'][-1])[iy]); E['Wx'].append(dx['W'][ix]); E['Wy'].append(dy['W'][iy])
        E['px'].append(dx['p'][ix]); E['py'].append(dy['p'][iy])
    E = {k: np.concatenate(v) for k, v in E.items()}
    edge = (E['px'] % 64 >= 4) & (E['px'] % 64 <= 59); E = {k: v[edge] for k, v in E.items()}
    n = len(E['a1']); o = dict(n=n)
    o['vs_eta_y'] = bins(E['ey'], E['a1'], np.array([-1.5, -.6, -.3, -.1, .1, .3, .6, 1.5]))
    o['vs_abs_eta_y'] = bins(np.abs(E['ey']), E['a1'], np.array([0, .1, .2, .35, .6, 1.5]))
    o['vs_axay'] = bins(E['axay'], E['a1'], np.quantile(E['axay'], np.linspace(0, 1, 7)))
    o['vs_xsh'] = bins(E['xsh'], E['a1'], np.quantile(E['xsh'], np.linspace(0, 1, 7)))
    # Y-row parity and Y-chip-position, X column parity: any sub-structure?
    o['vs_py_par'] = bins(E['py'] % 2, E['a1'], np.array([-.5, .5, 1.5]))
    # shape: normalise each event by its Y central amplitude, average
    s = E['Wy'][:, 4].max(1)[:, None]
    xc = (E['Wx'][:, 4] / s).mean(0); yc = (E['Wy'][:, 4] / s).mean(0)
    ex = ((E['Wx'][:, 5] - E['Wx'][:, 3]) / s).mean(0); x1 = ((E['Wx'][:, 5] + E['Wx'][:, 3]) / 2 / s).mean(0)
    dy_ = np.gradient(yc)
    for nm, B in (('Xc', [xc]), ('Yc', [yc]), ('Xc+Yc', [xc, yc]), ('Yc+dYc', [yc, dy_]), ('X+-1 mean', [x1])):
        X = np.stack(B, 1); c, *_ = np.linalg.lstsq(X, ex, rcond=None); r = ex - X @ c
        o[f'fit_{nm}'] = dict(coef=[float(v) for v in c], rms=float(r.std() / ex.std()))
    pk = lambda w: float(np.argmax(w)); t50 = lambda w: float(np.flatnonzero(w > .5 * w.max())[0])
    o['shape'] = dict(xc=xc.tolist(), yc=yc.tolist(), ex=ex.tolist(), x1=x1.tolist(),
                      peak=dict(xc=pk(xc), yc=pk(yc), ex=pk(ex), x1=pk(x1)), t50=dict(xc=t50(xc), yc=t50(yc), ex=t50(ex), x1=t50(x1)))
    OUT[bu] = o
    print(f'== {bu}: n={n}  copy {E["a1"].mean():+.3f}')
    for k in ('vs_eta_y', 'vs_abs_eta_y', 'vs_axay', 'vs_xsh', 'vs_py_par'):
        print(f'  {k:13s}', '  '.join(f'{b[0]:+.2f}:{b[1]:+.3f}±{b[2]:.3f}' if b[1] is not None else f'{b[0]:+.2f}:  -   ' for b in o[k]))
    print('  shape peak/t50 (samples, 60 ns):', o['shape']['peak'], o['shape']['t50'])
    for k in [k for k in o if k.startswith('fit_')]: print(f'  {k:14s} coef {np.round(o[k]["coef"], 3)}  resid/rms {o[k]["rms"]:.3f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ypad_probe.json'), 'w'))
