"""Time shape of the left/right neighbour asymmetry (2026-10-19). Same selection as nb_asym.py and the same alignment as
nb_far.py (t = 0 at the central's leading 50 %, each event / its central peak). Prints the mean waveforms of the
central c, n_-1, n_+1 and their difference d = n_-1 - n_+1, then fits d = a c + b dc/dt (+ const 0) over the record:
a = an in-phase copy of the central (charge-division crosstalk at the preamp input), b = a derivative term
(capacitive pickup after shaping). Also the -2/+2 difference. Usage: nb_asym_shape.py OUT.png NPZ..."""
import os, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_asym import load, HW
DT = 60.; grid = np.arange(-8, 24, .25)

def cross(w, frac):
    ipk = int(np.nanargmax(w)); lvl = frac * w[ipk]; b = np.flatnonzero(~(w[:ipk + 1] >= lvl))
    if not len(b) or b[-1] + 1 >= len(w): return np.nan
    i = b[-1]; return i + (lvl - w[i]) / (w[i + 1] - w[i])

def aligned_means(W):
    x = np.arange(W.shape[2]); out = []
    for w in W:
        l = cross(w[HW], .5)
        if np.isnan(l): continue
        out.append([np.interp(grid + l, x, s, left=np.nan, right=np.nan) for s in w])
    A = np.array(out); cov = np.isfinite(A[:, HW]).mean(0) >= .5
    M = np.nanmean(A, 0); M[:, ~cov] = np.nan
    return M, len(A)

out, files = sys.argv[1], sys.argv[2:]
fig, ax = plt.subplots(2, len(files), figsize=(4.2 * len(files), 7), squeeze=False)
for j, f in enumerate(files):
    for i, v in enumerate('xy'):
        D = load(f, v); M, n = aligned_means(D['W'])
        c, l1, r1, l2, r2 = M[HW], M[HW - 1], M[HW + 1], M[HW - 2], M[HW + 2]
        d1, d2 = l1 - r1, l2 - r2; dc = np.gradient(c, grid * DT)
        ok = np.isfinite(c) & np.isfinite(d1) & np.isfinite(d2) & np.isfinite(dc)
        X = np.stack([c[ok], dc[ok] * 100], 1)                     # b per 100 ns
        (a, b), *_ = np.linalg.lstsq(X, d1[ok], rcond=None); res = d1[ok] - X @ [a, b]
        (a2, b2), *_ = np.linalg.lstsq(X, d2[ok], rcond=None)
        rms = lambda y: np.sqrt(np.nanmean(y ** 2))
        tp = lambda y: grid[np.nanargmax(np.abs(y))] * DT
        print(f'{os.path.basename(f):18s} {v.upper()} n={n:5d}  d1 = {a:+.4f} c {b:+.4f} dc/dt*100ns  (rms d1 {rms(d1[ok]):.4f} -> {rms(res):.4f}),'
              f' |d1| max at {tp(d1):+.0f} ns;  d2 = {a2:+.4f} c {b2:+.4f} dc  | peaks: -1 {np.nanmax(l1):.3f} @{tp(l1):+.0f}  +1 {np.nanmax(r1):.3f} @{tp(r1):+.0f}'
              f'  -2 {np.nanmax(l2):.3f} @{tp(l2):+.0f}  +2 {np.nanmax(r2):.3f} @{tp(r2):+.0f}')
        a_ = ax[i, j]; t = grid * DT
        a_.plot(t, l1, label='-1'); a_.plot(t, r1, label='+1'); a_.plot(t, d1, 'k', lw=2, label='(-1)-(+1)')
        a_.plot(t, a * c + b * dc * 100, 'r--', label=f'{a:+.3f} c {b:+.3f} dc'); a_.plot(t, d2, color='0.5', label='(-2)-(+2)')
        a_.set_title(f'{os.path.basename(f)} {v.upper()} (n {n})', fontsize=9); a_.axhline(0, color='0.7', lw=.6); a_.grid(alpha=.3)
        a_.legend(fontsize=7)
for a_ in ax[-1]: a_.set_xlabel('t − central 50 % edge [ns]')
fig.tight_layout(); fig.savefig(out, dpi=80)
