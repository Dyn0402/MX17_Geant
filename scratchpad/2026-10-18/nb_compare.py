"""Neighbour (±1) shape, chain vs data, no fitted scaling: delay, nb/central peak, leading edge (t = 0 at the central's
leading 50 %), rms over -240..+1140 ns. Same alignment as 2026-10-13/xtalk_fit.py. Usage: nb_compare.py NPZ..."""
import sys, numpy as np
DT = 60.; grid = np.arange(-4, 20, .25)
def cross(w, frac, side):
    ipk = int(np.argmax(w)); lvl = frac * w[ipk]
    b = np.flatnonzero(w[:ipk + 1] < lvl)
    if not len(b): return np.nan
    i = b[-1]; j = i + 1
    return i + (lvl - w[i]) / (w[j] - w[i]) if w[j] != w[i] else float(i)
def aligned(W):
    W = W[W[:, 3].max(1) < 3500]; C, N = [], []
    for w in W:
        l = cross(w[3], .5, -1)
        if np.isnan(l): continue
        pk = w[3].max(); x = np.arange(w.shape[1])
        C.append(np.interp(grid + l, x, w[3] / pk)); N.append(.5 * (np.interp(grid + l, x, w[2] / pk) + np.interp(grid + l, x, w[4] / pk)))
    return np.mean(C, 0), np.mean(N, 0)
tpk = lambda y: grid[np.argmax(y)] * DT
for f in sys.argv[1:]:
    z = np.load(f); name = f.split('/')[-2]
    for v in 'xy':
        dc, dn = aligned(z[f'data_{v}_wf']); sc, sn = aligned(z[f'sim_{v}_wf'])
        print(f'{name:10s} {v}: delay data {tpk(dn)-tpk(dc):+4.0f} sim {tpk(sn)-tpk(sc):+4.0f} | nb/c {dn.max():.3f} {sn.max():.3f} | '
              f'nb(0) {dn[16]:.3f} {sn[16]:.3f} | nb rms {np.sqrt(np.mean((sn-dn)**2)):.4f} | central rms {np.sqrt(np.mean((sc-dc)**2)):.4f}')
