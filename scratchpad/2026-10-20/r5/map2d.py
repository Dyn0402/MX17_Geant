"""2-D map of the X asymmetry a1 over the board (2026-10-20): X position from the X peak's true column, Y position from
the matched Y peak's true row, in 4 x 4 cells (128 columns/rows each). A print scale error moves the strip offset along
x only; a rotation moves it along y. Usage: map2d.py"""
import json, os, numpy as np
import r5lib as R
OUT = {}
for b in ('det3 June-15', 'det4', 'det6', 'det7'):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64; a1 = E['alx'][:, 5] - E['alx'][:, 3]
    G = np.full((4, 4, 3), np.nan)
    for i in range(4):
        for j in range(4):
            m = (cx // 128 == i) & (cy // 128 == j)
            if m.sum() >= 30: G[i, j] = (a1[m].mean(), a1[m].std() / np.sqrt(m.sum()), m.sum())
    OUT[b] = G.tolist()
    print(f'== {b} (rows: X quarter 0..3 along x; columns: Y quarter 0..3 along y)')
    for i in range(4):
        print('   ' + '  '.join('       -        ' if np.isnan(G[i, j, 0]) else f'{G[i,j,0]:+.3f}±{G[i,j,1]:.3f} ({int(G[i,j,2]):4d})' for j in range(4)))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'map2d.json'), 'w'))
