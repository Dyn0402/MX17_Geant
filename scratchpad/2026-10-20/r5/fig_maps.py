"""Figure: X-view and Y-view neighbour asymmetry over each board (2026-10-20). 8 x 8 cells (matched X true column x
Y true row), a1 = alpha(+1) - alpha(-1) in DAQ order (det3's copy side positive). Diverging blue (-) / gray (0) /
red (+), clamped at +-0.25; cells with < 25 events blank. Titles carry the constant-model chi2 (map_y.json).
Writes design/figures/asym_maps_2026-10-20.png. Usage: fig_maps.py"""
import os, json, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import r5lib as R
HERE = os.path.dirname(os.path.abspath(__file__)); OUTP = os.path.join(HERE, '..', '..', '..', 'design', 'figures', 'asym_maps_2026-10-20.png')
cm = LinearSegmentedColormap.from_list('div', ['#2a78d6', '#f0efec', '#e34948']); cm.set_bad('#ffffff')
CH = json.load(open(os.path.join(HERE, 'map_y.json')))
B = ['det3 June-15', 'det2', 'det4', 'det6', 'det7', 'det3 May', 'det1']
fig, ax = plt.subplots(2, len(B), figsize=(2.3 * len(B), 5.4), constrained_layout=True)
for k, b in enumerate(B):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = (64 * (px // 64) + 63 - px % 64) // 64; cy = (64 * (py // 64) + 63 - py % 64) // 64
    for r, (v, nm) in enumerate(((E['alx'][:, 5] - E['alx'][:, 3], 'X'), (E['aly'][:, 5] - E['aly'][:, 3], 'Y'))):
        G = np.full((8, 8), np.nan)
        for i in range(8):
            for j in range(8):
                m = (cx == i) & (cy == j)
                if m.sum() >= 25: G[j, i] = v[m].mean()
        im = ax[r, k].imshow(G, origin='lower', cmap=cm, vmin=-.25, vmax=.25, extent=(0, 512, 0, 512))
        c, n, rms = CH[b][nm]
        ax[r, k].set_title(f'{b} · {nm}\nχ²/ndf {c:.0f}/{n} · mean {CH[b]["mean" + nm]:+.2f}', fontsize=8, color='#333')
        ax[r, k].set_xticks([0, 256, 512]); ax[r, k].set_yticks([0, 256, 512]); ax[r, k].tick_params(labelsize=6, colors='#666')
        for s in ax[r, k].spines.values(): s.set_color('#bbb')
    ax[1, k].set_xlabel('X column', fontsize=7, color='#555')
ax[0, 0].set_ylabel('Y row', fontsize=7, color='#555'); ax[1, 0].set_ylabel('Y row', fontsize=7, color='#555')
cb = fig.colorbar(im, ax=ax, shrink=.6, pad=.01); cb.set_label('neighbour asymmetry α(n+1) − α(n−1)', fontsize=7, color='#333'); cb.ax.tick_params(labelsize=6)
fig.suptitle('X-view neighbour asymmetry varies smoothly over the board, with both signs (det2: bands in y; det3: positive with a tilted ridge); the Y view is flat', fontsize=10, color='#222')
fig.savefig(OUTP, dpi=150, facecolor='#fcfcfb'); print(os.path.abspath(OUTP))
