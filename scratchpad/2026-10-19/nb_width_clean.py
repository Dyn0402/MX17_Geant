"""Vertical cluster width with and without the one-sided copy (2026-10-19). Per event, strip amplitudes A_k (peak of
each waveform, k = -3..3 around the peak strip, negatives clipped at 0, unsaturated peak strip): RMS width from both
sides (as xy_extent_m3 does, here without its thresholds) and from the + side mirrored onto the - side (copy-free in
X; the copy sits on the - side). Also from the - side mirrored. Pitch 0.78 mm. Means over events, standard errors.
Usage: nb_width_clean.py NPZ...   (data: M3 |ax|,|ay| < AMAX, default 3)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nb_asym
HW = 5; K = np.arange(-3, 4) * .78
def rms(A, k):
    w = A / A.sum(1, keepdims=True); m = (w * k).sum(1); return np.sqrt(np.clip((w * k ** 2).sum(1) - m ** 2, 0, None))
for f in sys.argv[1:]:
    z = np.load(f)
    for v in 'xy':
        W = z[f'{v}_wf']; m3 = z[f'{v}_m3']; amp = np.nanmax(W, 2); c = amp[:, HW]
        m = (c > 200) & (c < 3500) & np.isfinite(amp[:, HW - 3:HW + 4]).all(1)
        if not np.all(m3[:, 2:] == 0): m &= (np.abs(m3[:, 2]) < nb_asym.AMAX) & (np.abs(m3[:, 3]) < nb_asym.AMAX)
        A = np.clip(amp[m, HW - 3:HW + 4], 0, None)
        Ap = A.copy(); Ap[:, :3] = A[:, 4:][:, ::-1]; Am = A.copy(); Am[:, 4:] = A[:, :3][:, ::-1]
        r = [rms(x, K) for x in (A, Ap, Am)]
        print(f'{os.path.basename(f):18s} {v.upper()} n={m.sum():5d}  rms both {r[0].mean():.3f}±{r[0].std()/np.sqrt(m.sum()):.3f}  '
              f'+side mirrored {r[1].mean():.3f}  -side mirrored {r[2].mean():.3f} mm')
