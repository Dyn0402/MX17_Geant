"""Is the 2-D asymmetry pattern a track-angle effect? (2026-10-20) Tagged dumps (dump_tagged.py: position-ordered
+-5 windows with M3 x, y, ax, ay at the detector plane, |angles| < 6 deg). Per event, in position order,
a = alpha(pos-1) - alpha(pos+1) (the det3 copy side positive). 4 x 4 map in M3 (x, y) quantiles, for vertical tracks
(|ax|, |ay| < 3 deg) and for 3-6 deg, plus the mean ax, ay per cell (a lean of ~20 deg would be needed for 0.18).
Usage: map_vertical.py"""
import numpy as np, os, json
T = '/media/dylan/data/claude_tmp/tmp/nb/tagged'; OUT = {}
def alpha(W):
    c = W[:, 5]; n = np.einsum('ij,ij->i', c, c)
    return {k: np.einsum('ij,ij->i', np.nan_to_num(W[:, 5 + k]), c) / n for k in (-1, 1)}
for tag in ('g_det2', 'g_det3_wknd', 'o22', 'g_det4', 'g_det6_long', 'g_det7_long'):
    z = np.load(f'{T}/{tag}.npz', allow_pickle=True); W = np.nan_to_num(z['x_wf']); m3 = z['x_m3']
    A = W[:, 5].max(1); ok = (A > 200) & (A < 3500) & np.isfinite(m3).all(1)
    W, m3 = W[ok], m3[ok]; al = alpha(W); a = al[-1] - al[1]
    qx = np.quantile(m3[:, 0], [.25, .5, .75]); qy = np.quantile(m3[:, 1], [.25, .5, .75])
    ix = np.digitize(m3[:, 0], qx); iy = np.digitize(m3[:, 1], qy)
    vert = (np.abs(m3[:, 2]) < 3) & (np.abs(m3[:, 3]) < 3)
    print(f'== {tag}: n={len(a)} (vertical {vert.sum()}); mean a all {a.mean():+.3f}, vertical {a[vert].mean():+.3f}')
    o = {}
    for nm, sel in (('vertical <3deg', vert), ('3-6deg', ~vert)):
        G = np.full((4, 4), np.nan); E_ = np.full((4, 4), np.nan)
        for i in range(4):
            for j in range(4):
                k = sel & (ix == i) & (iy == j)
                if k.sum() >= 25: G[i, j] = a[k].mean(); E_[i, j] = a[k].std() / np.sqrt(k.sum())
        o[nm] = G.tolist()
        chi = np.nansum(((G - np.nanmean(G)) / E_) ** 2); nn = np.isfinite(G).sum()
        print(f'  {nm:15s} spread chi2 {chi:6.1f}/{nn-1}   ' + ' | '.join(' '.join('  -  ' if np.isnan(v) else f'{v:+.2f}' for v in G[i]) for i in range(4)))
    mx = np.array([[m3[(ix == i) & (iy == j), 2].mean() for j in range(4)] for i in range(4)])
    print(f'  mean ax per cell (deg): ' + ' | '.join(' '.join(f'{v:+.1f}' for v in mx[i]) for i in range(4)))
    OUT[tag] = o
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'map_vertical.json'), 'w'))
