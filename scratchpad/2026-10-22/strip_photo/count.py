"""Ruler-free strip count across the whole film from the overview photo (2026-10-22, FOLLOWUP §14). The strips run
left-right in this photo at 4.6-4.9 px period (resolved). Per 20-px column band: vertical profile, band-pass around the
local strip frequency, analytic signal, unwrapped phase; the patterned region is where the envelope exceeds half its median
(the solid border has none). N = phase advance across it / 2 pi. Stretch, tilt and wrinkles change the local period but not
N; a phase slip (low contrast under a reflection) shows up as columns disagreeing. Same-width designs: 512 strips at 780 um
or 499 at 800 um fill the 399.4 mm active width. Usage: count.py RAW_OVERVIEW.jpg   (-> count.json)"""
import sys, json, numpy as np
from PIL import Image, ImageOps
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
g = np.asarray(ImageOps.exif_transpose(Image.open(sys.argv[1])).convert('L'), float); H, W = g.shape
Y0, Y1 = 500, 3700; out = []
for x0 in range(150, 2750, 50):
    p = g[Y0:Y1, x0:x0 + 20].mean(1); p = p - uniform_filter1d(p, 25)
    F = np.fft.rfft(p); fr = np.fft.rfftfreq(len(p)); band = (fr > 1 / 5.6) & (fr < 1 / 4.0)
    k0 = np.argmax(np.abs(F) * band); f0 = fr[k0]
    F2 = F * np.exp(-0.5 * ((fr - f0) / (0.18 * f0)) ** 2); a = hilbert(np.fft.irfft(F2, len(p)))
    env = uniform_filter1d(np.abs(a), 41); ph = np.unwrap(np.angle(a))
    thr = 0.5 * np.median(env[env > np.percentile(env, 30)])
    # pattern ends = inner edges of the solid black border (brightness), not the envelope
    sm = uniform_filter1d(g[Y0:Y1, x0:x0 + 20].mean(1), 15); mid = len(sm) // 2
    def inner(lo, hi, step):
        seg = sm[lo:hi]; bl = seg.min(); pl = np.median(sm[mid - 300:mid + 300]); t = bl + 0.5 * (pl - bl); jb = lo + int(np.argmin(seg))
        j = jb
        while 0 < j < len(sm) - 1 and sm[j] < t: j += step
        return j
    i0 = inner(200, 700, +1); i1 = inner(len(sm) - 800, len(sm) - 200, -1)
    if i1 - i0 < 1000: continue
    n_low = int((env[i0:i1] < thr).sum())
    N = (ph[i1] - ph[i0]) / (2 * np.pi)
    loc = np.diff(ph[i0:i1]) / (2 * np.pi); slips = int((loc < 0).sum())            # phase running backwards = slip
    out.append(dict(x=x0 + 10, y0=int(Y0 + i0), y1=int(Y0 + i1), period_px=float(1 / f0), N=float(N), n_low=n_low, backsteps=slips))
    print(f'x {x0 + 10:5d}  rows {Y0 + i0:5d}-{Y0 + i1:5d}  period {1 / f0:5.3f} px  N {N:7.2f}   rows below threshold inside {n_low:4d}  backward steps {slips}')
Ns = np.array([o['N'] for o in out]); ok = np.array([o['n_low'] == 0 for o in out])
print(f'\nall columns: N median {np.median(Ns):.2f}, p10-p90 {np.percentile(Ns, 10):.2f}-{np.percentile(Ns, 90):.2f};  '
      f'clean columns (no dips, {ok.sum()}): median {np.median(Ns[ok]) if ok.any() else float("nan"):.2f}, spread {np.std(Ns[ok]) if ok.any() else float("nan"):.2f}')
json.dump(out, open('count.json', 'w'), indent=1)
