"""Mean aligned waveforms per strip distance (±k folded, each scaled to its own peak, and absolute), data vs sims,
X and Y. Same alignment as nb_far.py. Usage: nb_far_plot.py OUT.png NPZ...  (first NPZ's data leg is the data)"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open(__file__.replace('nb_far_plot.py', 'nb_far.py')).read().split('rng =')[0])
out, files = sys.argv[1], sys.argv[2:]
fig, ax = plt.subplots(4, 4, figsize=(17, 13), sharex=True)
for f_i, f in enumerate(files):
    z = np.load(f); name = f.split('/')[-2]
    for vi, v in enumerate('xy'):
        legs = ([('data', z[f'data_{v}_wf'])] if f_i == 0 else []) + [(name, z[f'sim_{v}_wf'])]
        for lab, W in legs:
            A = aligned(W); hw = A.shape[1] // 2; cov = np.isfinite(A[:, hw]).mean(0) >= .5
            for k in range(1, 5):
                n = np.where(cov, np.nanmean(np.nanmean(A[:, [hw - k, hw + k]], 1), 0), np.nan)
                st = dict(color='k', lw=2) if lab == 'data' else dict(lw=1.3)
                a = ax[2 * vi + (k > 2), (k - 1) % 2 * 2]; a.plot(grid * DT, n, label=lab, **st); a.set_title(f'{v.upper()} ±{k}  (/central peak)')
                b = ax[2 * vi + (k > 2), (k - 1) % 2 * 2 + 1]; b.plot(grid * DT, n / np.nanmax(n), label=lab, **st); b.set_title(f'{v.upper()} ±{k}  (own peak = 1)')
for a in ax.flat: a.axvline(0, color='0.7', lw=.6); a.axhline(0, color='0.7', lw=.6); a.grid(alpha=.3)
for a in ax[-1]: a.set_xlabel('t − central 50 % edge [ns]')
ax[0, 0].legend(); fig.tight_layout(); fig.savefig(out, dpi=80)
