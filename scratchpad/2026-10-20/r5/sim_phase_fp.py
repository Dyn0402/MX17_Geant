"""Stripe-offset hypothesis in the full sim, phase-resolved (2026-10-20). The sim's 800 um stripes over 780 um pads make
each X column phase (col mod 40) a fixed stripe-to-column offset. Rank the 40 phases by their X-only a1 on seeds
5..10 (+ the original 5/6), then measure the fingerprints on the held-out seeds 11..16 for:
  sim POS = the 10 phases with the most positive a1, sim NEG = the 10 most negative, sim MID = the 10 nearest zero.
Fingerprints as for det2's regions (det2_regions.py): fingerprint.py classes + bookkeeping, timing.py, x-position
conditional. Arms: xo0, th00 + thy*. Usage: sim_phase_fp.py"""
import os, glob, re, json, numpy as np
import r5lib as R, fingerprint as FP, timing as T, xpos_probe as XP
SIM = os.path.expanduser('~/x17/nb_asym_raw/simxoff')
def load(seeds):
    Wx, Wy, px, py = [], [], [], []
    for ang in ('th00', 'thy05', 'thy10', 'thy15', 'thy20', 'thy30'):
        for s in seeds:
            fx, fy = f'{SIM}/xo0_{ang}_s{s}_x.npz', f'{SIM}/xo0_{ang}_s{s}_y.npz'
            if not (os.path.exists(fx) and os.path.exists(fy)) or np.load(fx)['wf'].ndim != 3 or np.load(fy)['wf'].ndim != 3: continue
            dx, dy = R.L.load(fx), R.L.load(fy)
            kx = dx['fi'].astype(np.int64) * 10**9 + dx['ev']; ky = dy['fi'].astype(np.int64) * 10**9 + dy['ev']
            _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
            Wx.append(dx['W'][ix]); Wy.append(dy['W'][iy]); px.append(dx['p'][ix]); py.append(dy['p'][iy])
    return [np.concatenate(v) for v in (Wx, Wy, px, py)]
def phase(p): p = p % 512; return (64 * (p // 64) + 63 - p % 64) % 40
def a1of(Wx):
    c = Wx[:, 4].astype(float); n = np.einsum('ij,ij->i', c, c)
    return (np.einsum('ij,ij->i', Wx[:, 5], c) - np.einsum('ij,ij->i', Wx[:, 3], c)) / n
Wr, _, pr, _ = load((5, 6, 7, 8, 9, 10)); ar = a1of(Wr); phr = phase(pr)
rank = np.array([ar[phr == k].mean() if (phr == k).sum() >= 20 else np.nan for k in range(40)])
print('rank a1 per phase (seeds 5-10):', ' '.join('  -  ' if np.isnan(x) else f'{x:+.2f}' for x in rank))
ok = np.flatnonzero(np.isfinite(rank)); o = ok[np.argsort(rank[ok])]
SETS = {'sim POS': o[-10:], 'sim NEG': o[:10], 'sim MID': o[np.argsort(np.abs(rank[o]))[:10]]}
Wx, Wy, px, py = load(range(11, 17)); ph = phase(px)
print(f'held-out matched events: {len(px)}')
OUT = {}
for nm, S in SETS.items():
    m = np.isin(ph, S); fp = FP.fp_arrays(Wx[m], Wy[m], px[m], py[m]); FP.show(nm + f' ({len(S)} ph)', fp)
    cy = Wy[m][:, 4].astype(float); ny = np.einsum('ij,ij->i', cy, cy)
    ysh = (np.einsum('ij,ij->i', Wy[m][:, 3], cy) + np.einsum('ij,ij->i', Wy[m][:, 5], cy)) / ny; par = (px[m] + py[m]) % 2
    yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean())
    t = T.prof(Wx[m], Wy[m], yp); print('   ' + T.line('timing', t))
    OUT[nm] = dict(phases=S.tolist(), marg=fp['marg'], A_Y=fp['A_Y'].tolist(), A_X=fp['A_X'].tolist(), fY=fp['fY'].tolist(),
                   m1=fp['m1'].tolist(), p1=fp['p1'].tolist(), timing={k: t[k] for k in ('Yc', 'Xc', 'Xm1', 'exc', 'fit')})
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_phase_fp.json'), 'w'), indent=1, default=float)
