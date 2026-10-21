"""Strip-offset hypothesis tested in the full sim (2026-10-20). The sim's ESL strips (800 um) beat against the
780 um pads, so each column phase (col mod 40) is a 'board' whose strip centres sit at a fixed offset from the pad
columns. Group the sim's events (xo0 arms, all 24 files; tilted ones included for statistics) by phase and compute
the X-only marginals (a1, coin) per group of 5 phases, plus the matched X/Y class asymmetries in two Y-share halves.
If some phase reproduces det3 June-15 (a1 ~ +0.14, A_Y ~ +0.25, narrow A_X < 0), a pitch-locked print at that
offset is a viable mechanism. Usage: sim_phase.py"""
import os, glob, numpy as np
import r5lib as R
SIM = os.path.expanduser('~/x17/nb_asym_raw/simxoff')
Xs, Ys = [], []
for fx in sorted(glob.glob(f'{SIM}/xo0_*_x.npz')):
    fy = fx[:-6] + '_y.npz'
    if np.load(fx)['wf'].ndim != 3: continue
    dx = R.L.load(fx); Xs.append(dx)
    if np.load(fy)['wf'].ndim == 3:
        dy = R.L.load(fy); kx = dx['fi'].astype(np.int64) * 10**9 + dx['ev']; ky = dy['fi'].astype(np.int64) * 10**9 + dy['ev']
        _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
        Ys.append((dx['a1'][ix], dx['p'][ix] % 512, dy['p'][iy] % 512, (dy['al'][1] + dy['al'][-1])[iy]))
p = np.concatenate([d['p'] % 512 for d in Xs]); a1 = np.concatenate([d['a1'] for d in Xs])
A = np.concatenate([d['A'] for d in Xs]); col = 64 * (p // 64) + 63 - p % 64; ph = col % 40
print(f'X-only: n={len(a1)}; per phase group (col mod 40, 5 phases each): a1 / coin / n')
for g in range(8):
    m = (ph // 5) == g
    print(f'  phases {5*g:2d}-{5*g+4:2d}: a1 {a1[m].mean():+.3f} ± {a1[m].std()/np.sqrt(m.sum()):.3f}  coin {(A[m,5] > A[m,3]).mean():.3f}  n {m.sum()}')
a1m, px, py, ysh = map(np.concatenate, zip(*Ys)); colm = 64 * (px // 64) + 63 - px % 64; phm = colm % 40
par = (px + py) % 2; yp = par == int(ysh[par == 1].mean() < ysh[par == 0].mean()); med = np.median(ysh)
print(f'matched: n={len(a1m)}; per phase group: A_Y / A_X (narrow half | wide half), n')
for g in range(8):
    m = (phm // 5) == g
    s = lambda k: f'{a1m[k].mean():+.2f}' if k.sum() > 15 else '  -  '
    print(f'  phases {5*g:2d}-{5*g+4:2d}: A_Y {s(m & yp & (ysh < med))} | {s(m & yp & (ysh >= med))}   A_X {s(m & ~yp & (ysh < med))} | {s(m & ~yp & (ysh >= med))}   n {m.sum()}')
