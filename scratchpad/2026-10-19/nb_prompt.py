"""Prompt in-phase content of every neighbour strip (2026-10-19): fit n_k(t) = alpha_k c(t) on the early window
-200 <= t <= +300 ns (central's leading 50 % at 0), where the slow resistive spread has not yet arrived, per side,
k = 1..5, with bootstrap errors. Data minus sim alpha is the prompt share the sim lacks.
Usage: nb_prompt.py NPZ...   (B env, default 150)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_asym import load, HW
from nb_side import aligned, grid, DT
B = int(os.environ.get('B', 150)); rng = np.random.default_rng(3); win = (grid * DT >= -200) & (grid * DT <= 300)
def alphas(A):
    M = np.nanmean(A, 0); c = M[HW]
    return np.array([np.nansum(M[HW + k][win] * c[win]) / np.nansum(c[win] ** 2) for k in range(-HW, HW + 1)])
for f in sys.argv[1:]:
    for v in 'xy':
        A = aligned(load(f, v)['W']); a0 = alphas(A); e = np.std([alphas(A[rng.integers(0, len(A), len(A))]) for _ in range(B)], 0)
        print(f'{os.path.basename(f):18s} {v.upper()} n={len(A):5d} ' + '  '.join(f'{k:+d}:{a0[k+HW]:.4f}±{e[k+HW]:.4f}' for k in range(-HW, HW + 1) if k))
