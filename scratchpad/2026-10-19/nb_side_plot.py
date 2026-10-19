"""Mean aligned neighbour waveforms, + side (clean of the one-sided copy) solid and - side dashed, data vs sims,
k = 1..4, X and Y, each / central peak (2026-10-19). Usage: nb_side_plot.py OUT.png DATA.npz SIM.npz..."""
import os, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_asym import load, HW
from nb_side import aligned, grid, DT
out, files = sys.argv[1], sys.argv[2:]
fig, ax = plt.subplots(2, 5, figsize=(22, 8))
cols = ['k', 'C0', 'C1', 'C2']
for fi, f in enumerate(files):
    for vi, v in enumerate('xy'):
        A = aligned(load(f, v)['W']); cov = np.isfinite(A[:, HW]).mean(0) >= .5
        M = np.nanmean(A, 0); M[:, ~cov] = np.nan; t = grid * DT; lab = os.path.basename(f).replace('.npz', '')
        ax[vi, 0].plot(t, M[HW], color=cols[fi], label=lab); ax[vi, 0].set_title(f'{v.upper()} central')
        for k in range(1, 5):
            a = ax[vi, k]; a.plot(t, M[HW + k], color=cols[fi], lw=2 if fi == 0 else 1.2, label=f'{lab} +{k}')
            a.plot(t, M[HW - k], color=cols[fi], ls='--', lw=1, label=f'{lab} -{k}'); a.set_title(f'{v.upper()} ±{k} (+ solid, - dashed)')
for a in ax.flat: a.axvline(0, color='0.7', lw=.6); a.axhline(0, color='0.7', lw=.6); a.grid(alpha=.3)
ax[0, 1].legend(fontsize=7); ax[0, 0].legend(fontsize=7)
for a in ax[-1]: a.set_xlabel('t − central 50 % edge [ns]')
fig.tight_layout(); fig.savefig(out, dpi=75)
