"""Sim fingerprint at fixed strip offset (2026-10-20 evening): matched X/Y sim events (xo0, th00 + thy*, seeds 5..16)
grouped by the X peak column's strip offset delta (sim_delta.py, orientation 'col'): POS = delta in [+150, +350],
NEG = [-350, -150], ZERO = [-75, +75]. Same fingerprints as det2_regions.py (classes, bookkeeping, timing).
Usage: sim_delta_fp.py"""
import os, glob, re, json, numpy as np
import r5lib as R, fingerprint as FP, timing as T
SIM = os.path.expanduser('~/x17/nb_asym_raw/simxoff'); P = 780.; ORG = -(512 - 1) / 2 * P
Wx, Wy, px, py = [], [], [], []
for fx in sorted(glob.glob(f'{SIM}/xo0_*_x.npz')):
    fy = fx[:-6] + '_y.npz'
    if not re.search('_th00_|_thy', fx) or not os.path.exists(fy) or np.load(fx)['wf'].ndim != 3 or np.load(fy)['wf'].ndim != 3: continue
    dx, dy = R.L.load(fx), R.L.load(fy)
    kx = dx['fi'].astype(np.int64) * 10**9 + dx['ev']; ky = dy['fi'].astype(np.int64) * 10**9 + dy['ev']
    _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
    Wx.append(dx['W'][ix]); Wy.append(dy['W'][iy]); px.append(dx['p'][ix]); py.append(dy['p'][iy])
Wx, Wy, px, py = [np.concatenate(v) for v in (Wx, Wy, px, py)]
p = px % 512; col = 64 * (p // 64) + 63 - p % 64; delta = ((ORG + col * P + 400.) % 800.) - 400.
print(f'matched sim events {len(px)}')
OUT = {}
for nm, m in (('sim delta +150..+350', (delta >= 150) & (delta <= 350)), ('sim delta -350..-150', (delta >= -350) & (delta <= -150)),
              ('sim delta -75..+75', np.abs(delta) <= 75)):
    o = FP.fp_arrays(Wx[m], Wy[m], px[m], py[m]); FP.show(nm, o)
    cy = Wy[m][:, 4].astype(float); ny = np.einsum('ij,ij->i', cy, cy)
    ysh = (np.einsum('ij,ij->i', Wy[m][:, 3], cy) + np.einsum('ij,ij->i', Wy[m][:, 5], cy)) / ny; par = (px[m] + py[m]) % 2
    yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean()); t = T.prof(Wx[m], Wy[m], yp); print('   ' + T.line('timing', t))
    OUT[nm] = dict(marg=o['marg'], A_Y=o['A_Y'].tolist(), A_X=o['A_X'].tolist(), fY=o['fY'].tolist(), m1=o['m1'].tolist(), p1=o['p1'].tolist(),
                   timing={k: t[k] for k in ('Yc', 'Xc', 'Xm1', 'exc')})
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_delta_fp.json'), 'w'), indent=1, default=float)
