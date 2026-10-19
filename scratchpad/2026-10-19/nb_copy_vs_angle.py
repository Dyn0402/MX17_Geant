"""One-sided X asymmetry vs the M3 track angle (2026-10-19): if det3 were tilted with respect to M3, the asymmetry
would follow the projected angle along the measured coordinate and vanish at the tilt angle. Per event
copy = alpha_-1 - alpha_+1 (in-phase content, as nb_copy_chip_edge.py), binned in the M3 angle that projects onto the
strips' measured coordinate (X strips ~ ref_y -> ay; Y strips ~ -ref_x -> -ax) and in the other angle.
Usage: nb_copy_vs_angle.py NPZ..."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_asym import HW
E_ = np.arange(-6, 6.1, 1.5)
for f in sys.argv[1:]:
    z = np.load(f)
    for v in 'xy':
        W = z[f'{v}_wf']; m3 = z[f'{v}_m3']; amp = np.nanmax(W, 2); c = amp[:, HW]
        m = (c > 200) & (c < 3500) & np.isfinite(W[:, HW - 1:HW + 2]).all((1, 2))
        W = np.nan_to_num(W[m]); m3 = m3[m]; wc = W[:, HW]
        al = lambda k: np.einsum('ij,ij->i', W[:, HW + k], wc) / np.einsum('ij,ij->i', wc, wc)
        d = al(-1) - al(1)
        meas = m3[:, 3] if v == 'x' else -m3[:, 2]; other = m3[:, 2] if v == 'x' else m3[:, 3]
        for nm, ang in (('measured-coord angle', meas), ('other angle', other)):
            b = np.digitize(ang, E_); s = []
            for i in range(1, len(E_)):
                k = b == i
                s.append(f'[{E_[i-1]:+.1f},{E_[i]:+.1f}): {d[k].mean():+.3f}±{d[k].std() / np.sqrt(max(k.sum(), 1)):.3f}' if k.sum() > 20 else f'[{E_[i-1]:+.1f},{E_[i]:+.1f}): –')
            sl = np.polyfit(ang, d, 1)
            print(f'{os.path.basename(f):18s} {v.upper()} {nm:21s}: ' + '  '.join(s) + f'   slope {sl[0]:+.4f}/deg, zero at {-sl[1] / sl[0]:+.1f} deg')
