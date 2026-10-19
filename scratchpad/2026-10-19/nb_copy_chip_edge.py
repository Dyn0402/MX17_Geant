"""Does the one-sided X copy cross DREAM chip boundaries? (2026-10-19). Per event the in-phase content
alpha_k = <w_k, w_c>/<w_c, w_c> of the position neighbours k = -1, +1; copy estimate = alpha_-1 - alpha_+1.
Grouped by the peak CHANNEL within its 64-channel chip: with connector-inverted wiring, position -1 of channel p is
channel p+1, except at p % 64 == 63 where position -1 is on the next chip; position +1 is channel p-1 except at
p % 64 == 0. Usage: nb_copy_chip_edge.py NPZ...  (data only; AMAX as nb_asym)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nb_asym
from nb_asym import HW
for f in sys.argv[1:]:
    z = np.load(f)
    for v in 'xy':
        W = z[f'{v}_wf']; m3 = z[f'{v}_m3']; ch = z[f'{v}_ch'].astype(int); amp = np.nanmax(W, 2); c = amp[:, HW]
        m = (c > 200) & (c < 3500) & np.isfinite(W[:, HW - 1:HW + 2]).all((1, 2)) & (np.abs(m3[:, 2]) < nb_asym.AMAX) & (np.abs(m3[:, 3]) < nb_asym.AMAX)
        W, ch = np.nan_to_num(W[m]), ch[m]; wc = W[:, HW]
        al = lambda k: np.einsum('ij,ij->i', W[:, HW + k], wc) / np.einsum('ij,ij->i', wc, wc)
        d = al(-1) - al(1); p = ch[:, HW] % 64
        same_m = (ch[:, HW - 1] // 64 == ch[:, HW] // 64); same_p = (ch[:, HW + 1] // 64 == ch[:, HW] // 64)
        row = []
        for nm, s in (('interior', same_m & same_p), ('pos -1 on other chip', ~same_m), ('pos +1 on other chip', ~same_p)):
            row.append(f'{nm}: {d[s].mean():+.3f}±{d[s].std() / np.sqrt(max(s.sum(), 1)):.3f} (n {s.sum()})')
        print(f'{os.path.basename(f):18s} {v.upper()}  copy (alpha_-1 - alpha_+1): ' + '   '.join(row))
        print(f'      -1 neighbour channel - peak channel (interior): {np.unique((ch[:, HW - 1] - ch[:, HW])[same_m & same_p], return_counts=True)}')
