"""Copy survey figure + summary (2026-10-19) from survey_table.py's CSV.
Quality flag: alpha_+1 outside 0.06..0.6 (window/pickup problems: the 6-23 short sub-run, det6/det7
beam X), run_55 (sparse, 32-sample beam windows), or the June-17 X of the June-15 board (its central waveform rings, asym_shape_jun17.png).
Top: det3 X (filled) and Y (open) copy a vs sub-run start time, colour = drift voltage, square = beam.
Middle: the 6-27 drift scan in time order. Bottom: other boards (weighted means, flagged sets dropped).
Usage: survey_plot.py SURVEY.csv OUT.png"""
import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = pd.read_csv(sys.argv[1]); d['t'] = pd.to_datetime(d.tstart); d['beam'] = d.run.str.startswith('run_')
d['bad'] = (d.alpha_p1 < .06) | (d.alpha_p1 > .6) | (d.run.str.contains('6-17') & (d.phys == 'det3') & (d.view == 'X'))
d['bad'] |= d.run.str.contains('run_55')
good = d[~d.bad]
fig, ax = plt.subplots(3, 1, figsize=(14, 13), gridspec_kw=dict(height_ratios=[3, 1.6, 1.6]))
x3 = good[good.phys == 'det3']
for v, mfc in (('X', None), ('Y', 'none')):
    s = x3[x3.view == v]
    for b, mk in ((False, 'o'), (True, 's')):
        q = s[s.beam == b]
        sc = ax[0].scatter(q.t, q.a_ch, c=q.drift_V, cmap='viridis', vmin=100, vmax=1200, marker=mk, s=28,
                           facecolors=mfc if mfc else None, edgecolors=None if mfc is None else plt.cm.viridis((q.drift_V - 100) / 1100))
        ax[0].errorbar(q.t, q.a_ch, q.ea_ch, fmt='none', ecolor='0.6', lw=.6)
fig.colorbar(sc, ax=ax[0], label='drift voltage [V]')
ax[0].axvline(pd.Timestamp('2026-06-15'), color='r', ls='--', lw=1); ax[0].text(pd.Timestamp('2026-06-15'), .17, ' June-15 bulk', color='r', fontsize=8)
ax[0].axhline(0, color='k', lw=.6); ax[0].set_ylabel('copy a (channel n+1 ← n)'); ax[0].grid(alpha=.3)
ax[0].set_title('det3 X (filled) / Y (open): May build ≈ 0; June-15 build from 6-22, decaying; n_TOF beam (squares) July', fontsize=10)
ds = good[(good.phys == 'det3') & good['sub'].str.startswith('drift_scan') & (good.view == 'X')].sort_values('t')
ax[1].errorbar(ds.drift_V, ds.a_ch, ds.ea_ch, fmt='o-'); [ax[1].annotate(r.t.strftime('%H:%M'), (r.drift_V, r.a_ch), fontsize=7) for _, r in ds.iterrows()]
lr = good[(good.phys == 'det3') & good['sub'].str.startswith('long_run_resist_490V') & (good.view == 'X')]
ax[1].errorbar(lr.drift_V, lr.a_ch, lr.ea_ch, fmt='s', color='C3', label='long run right after (22:49)')
ax[1].set_xlabel('drift voltage [V] (6-27 drift scan, run 1100 → 100 V in time)'); ax[1].set_ylabel('copy a (X)'); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)
o = good[good.phys != 'det3']; labs = []
for j, (k, g) in enumerate(o.groupby(['phys', 'view', 'beam'])):
    w = 1 / g.ea_ch ** 2; ax[2].errorbar(j, (g.a_ch * w).sum() / w.sum(), 1 / np.sqrt(w.sum()), fmt='o', color='C3' if k[2] else 'k')
    ax[2].scatter([j] * len(g), g.a_ch, s=8, color='0.6'); labs.append(f'{k[0]} {k[1]} {"beam" if k[2] else "bench"}')
ax[2].set_xticks(range(len(labs)), labs, rotation=30, fontsize=8); ax[2].axhline(0, color='k', lw=.6); ax[2].grid(alpha=.3)
ax[2].set_ylabel('copy a'); ax[2].set_title('other boards (grey: each dataset; black/red: weighted mean bench/beam)', fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=80)
pd.set_option('display.width', 220); pd.set_option('display.max_rows', 300)
print(d.sort_values(['phys', 'view', 't'])[['phys', 'view', 'tstart', 'run', 'sub', 'feu', 'mesh_V', 'drift_V', 'n', 'a_ch', 'ea_ch', 'alpha_p1', 'bad']]
      .to_string(index=False, float_format=lambda v: f'{v:.3f}'))
