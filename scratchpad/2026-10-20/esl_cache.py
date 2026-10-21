"""Event-level cache of X/Y-matched events for the ESL-strip study (2026-10-20).
Per build: X and Y alpha_k (k=-4..4), DAQ peak channels, parity, sub-run index, Y-pad flag (per sub-run, smaller Y share).
Writes ~/x17/nb_asym_raw/esl_cache.npz (derived data, desktop)."""
import os, sys, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L, asym_more as A

rows = M.survey_rows()
out = {k: [] for k in ('build', 'sr', 'px', 'py', 'par', 'ypad', 'ax', 'ay', 'Ax', 'Ay', 'drift', 'wx')}
srs = []
for b, bu in sorted({(r['board'], r['build']) for r in rows}):
    for rx, ry in A.pairs(rows, lambda r: r['build'] == bu):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        par = (dx['p'][ix] + dy['p'][iy]) % 2; ysh = (dy['al'][1] + dy['al'][-1])[iy]
        ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean())
        n = len(ix); srs.append((bu, rx['run'], rx['sub'], rx['t'], rx['drift'], rx['feu']))
        out['build'] += [bu] * n; out['sr'].append(np.full(n, len(srs) - 1)); out['drift'] += [rx['drift'] or ''] * n
        out['px'].append(dx['p'][ix]); out['py'].append(dy['p'][iy]); out['par'].append(par); out['ypad'].append(par == ypar)
        out['ax'].append(np.stack([dx['al'][k][ix] for k in range(-4, 5)], 1))
        out['ay'].append(np.stack([dy['al'][k][iy] for k in range(-4, 5)], 1))
        out['Ax'].append(dx['A'][ix][:, 4]); out['Ay'].append(dy['A'][iy][:, 4])
        c = dx['W'][ix][:, 4]; out['wx'].append((c > 0.3 * c.max(1, keepdims=True)).sum(1) * 60)
o = {k: (np.array(v) if k in ('build', 'drift') else np.concatenate(v)) for k, v in out.items()}
o['srs'] = np.array(srs, dtype=object)
np.savez(os.path.expanduser('~/x17/nb_asym_raw/esl_cache.npz'), **o)
print({k: len(v) for k, v in o.items()}); print(np.unique(o['build'], return_counts=True))
