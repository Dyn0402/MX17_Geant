"""One-sided copy a vs the central's peak amplitude (2026-10-19): a linear induced copy is amplitude independent;
an electronics effect near the rail grows. a from d1 = n_-1 - n_+1 = a (c - n_+2) on aligned means, per amplitude bin,
bootstrap errors. Usage: nb_copy_vs_amp.py NPZ...  (AMAX env as nb_asym; railed (>= 3500) excluded by load)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nb_asym
from nb_asym import HW
from nb_side import aligned, grid
rng = np.random.default_rng(4); EDGES = [200, 700, 1200, 1800, 2500, 3500]
def a_of(A):
    M = np.nanmean(A, 0); ok = np.isfinite(M[HW - 1:HW + 3]).all(0)
    d1 = M[HW - 1] - M[HW + 1]; x = M[HW] - M[HW + 2]; return np.dot(d1[ok], x[ok]) / np.dot(x[ok], x[ok])
for f in sys.argv[1:]:
    z = np.load(f); W = z['x_wf']; m3 = z['x_m3']; amp = np.nanmax(W, 2); c = amp[:, HW]
    sim = np.all(m3[:, 2:] == 0); base = np.isfinite(amp[:, HW - 2:HW + 3]).all(1)
    if not sim: base &= (np.abs(m3[:, 2]) < nb_asym.AMAX) & (np.abs(m3[:, 3]) < nb_asym.AMAX)
    row = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        m = base & (c >= lo) & (c < hi)
        if m.sum() < 30: row.append(f'{lo}-{hi}: –'); continue
        A = aligned(W[m] / c[m, None, None]); a0 = a_of(A); e = np.std([a_of(A[rng.integers(0, len(A), len(A))]) for _ in range(100)])
        row.append(f'{lo}-{hi}: {a0:+.3f}±{e:.3f} ({m.sum()})')
    m = base & (c >= 3500); row.append(f'railed fraction {m.sum() / base.sum():.2f}')
    print(f'{os.path.basename(f):18s} X  ' + '   '.join(row))
