"""Per-side neighbour timing and the one-sided copy model, on the tagged dumps (2026-10-19).
Copy model: the strip at the LOWER position picks up a in-phase copy of its upper neighbour, so
    n_-1 = s_1 + a c,   n_-2 = s_2 + a s_1 (+ a^2 c),   n_+1 = s_1 + a s_2,   n_+2 = s_2 + a s_3
with s_k the true (symmetric) share. Fits on the aligned mean waveforms: d1 = n_-1 - n_+1 = a (c - n_+2),
d2 = n_-2 - n_+2 = a (n_+1 - n_+3) to first order; each gives an independent a. Bootstrap errors (B resamples).
Then the per-side delays (peak and 50 % edge, vs the central) for k = 1..4, and the copy-corrected folded neighbour
s_1 = n_+1 - a n_+2 against the sim's folded ±1. Same selection as nb_asym.load, same alignment as nb_far.
Usage: nb_side.py NPZ...   (B env, default 200; AMAX env as nb_asym)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_asym import load, HW
DT = 60.; grid = np.arange(-8, 24, .25); B = int(os.environ.get('B', 200)); rng = np.random.default_rng(2)

def cross(w, frac):
    ipk = int(np.nanargmax(w)); lvl = frac * w[ipk]; b = np.flatnonzero(~(w[:ipk + 1] >= lvl))
    if not len(b) or b[-1] + 1 >= len(w): return np.nan
    i = b[-1]; return i + (lvl - w[i]) / (w[i + 1] - w[i])

def aligned(W):
    x = np.arange(W.shape[2]); out = []
    for w in W:
        l = cross(w[HW], .5)
        if np.isnan(l): continue
        out.append([np.interp(grid + l, x, s, left=np.nan, right=np.nan) for s in w])
    return np.array(out)

def stats(A):
    cov = np.isfinite(A[:, HW]).mean(0) >= .5
    M = np.nanmean(A, 0); M[:, ~cov] = np.nan; c = M[HW]; ok = np.isfinite(M[HW - 2:HW + 4]).all(0)
    fit = lambda y, x: float(np.dot(x[ok], y[ok]) / np.dot(x[ok], x[ok]))
    a1 = fit(M[HW - 1] - M[HW + 1], c - M[HW + 2]); a2 = fit(M[HW - 2] - M[HW + 2], M[HW + 1] - M[HW + 3])
    tc = grid[np.nanargmax(c)]; ec = cross(np.nan_to_num(c, nan=-1), .5)
    r = [a1, a2]
    for k in (1, 2, 3, 4):
        for y in (M[HW - k], M[HW + k]):
            yy = np.nan_to_num(y, nan=-1); r += [(grid[np.argmax(yy)] - tc) * DT, (np.interp(cross(yy, .5), np.arange(len(grid)), grid) - np.interp(ec, np.arange(len(grid)), grid)) * DT]
    s1 = M[HW + 1] - a1 * M[HW + 2]; s1 = np.nan_to_num(s1, nan=-1)
    r += [(grid[np.argmax(s1)] - tc) * DT, s1.max()]
    return np.array(r), M

if __name__ == '__main__':
    for f in sys.argv[1:]:
        for v in 'xy':
            D = load(f, v); A = aligned(D['W']); v0, M = stats(A)
            bs = np.array([stats(A[rng.integers(0, len(A), len(A))])[0] for _ in range(B)]); e = np.nanstd(bs, 0)
            print(f'\n{os.path.basename(f)} {v.upper()}  n={len(A)}   copy a from d1 = {v0[0]:+.3f}±{e[0]:.3f}, from d2 = {v0[1]:+.3f}±{e[1]:.3f}')
            for i, k in enumerate((1, 2, 3, 4)):
                j = 2 + 4 * i
                print(f'    -{k}: peak {v0[j]:+5.0f}±{e[j]:3.0f}  edge {v0[j+1]:+5.0f}±{e[j+1]:3.0f}   |   +{k}: peak {v0[j+2]:+5.0f}±{e[j+2]:3.0f}  edge {v0[j+3]:+5.0f}±{e[j+3]:3.0f}')
            print(f'    copy-corrected s1 = n_+1 - a n_+2: peak delay {v0[-2]:+.0f}±{e[-2]:.0f} ns, s1/c {v0[-1]:.3f}±{e[-1]:.3f}')
