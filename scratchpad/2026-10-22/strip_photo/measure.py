"""Resistive-strip print film, phone macro photos with a steel ruler (2026-10-22, FOLLOWUP §14).
Input: the levelled photos from level.py ($LV/<tag>.png). Bands are set by hand from the 1/4-scale previews (BANDS, preview
px; x4 = full resolution) and checked with overlays ($OVL/<tag>.jpg):
 1. ruler band A (top of the mm ticks): tick minima -> tick index; cubic fit index -> x gives mm(x) at band A (distortion);
    ruler band B (lower, only the 5 mm / cm ticks): px per mm at B (reported only). The scale's change with image row
    (camera tilt) is taken from the strips' own pitch in px vs row and extrapolated from band A to the strip band;
 2. strip band (just above the rounded gap ends): dark = printed strips, clear = gaps; edges at the local mid-level, sub-pixel;
 3. per strip: centre and width in mm (band-A map x tilt correction); pitch = slope of centre vs index; fill = width / pitch
    (scale-free). Also the pitch in each third of the frame, and strip pitch in px vs row (a second, ruler-free tilt check).
Usage: measure.py   (-> measure.json)"""
import os, json, numpy as np
from PIL import Image
from scipy.ndimage import uniform_filter1d, median_filter, maximum_filter1d, minimum_filter1d
from scipy.signal import find_peaks
LV = os.environ.get('LV', 'leveled'); OVL = os.environ.get('OVL', 'ovl'); os.makedirs(OVL, exist_ok=True)
#        strip rows   ruler A rows  ruler B rows   x range       (preview px)
BANDS = {'p1': ((228, 262), (350, 372), (382, 392), (10, 990)),
         'p2': ((400, 435), (480, 520), (585, 600), (10, 990)),
         'p3': ((455, 510), (560, 600), (650, 665), (15, 985)),
         'p4': ((450, 505), (545, 600), (690, 710), (10, 990)),
         'p5': ((380, 425), (480, 540), (585, 600), (170, 990))}

def ticks(prof, step_hint=None):
    hp = prof - median_filter(prof, 151)
    if step_hint is None:
        f = np.abs(np.fft.rfft(hp * np.hanning(len(hp)))); fr = np.fft.rfftfreq(len(hp)); k = np.argmax(f[5:]) + 5; step_hint = 1 / fr[k]
    pk, _ = find_peaks(-hp, distance=0.6 * step_hint, prominence=0.3 * np.std(hp))
    c = []
    for p in pk:
        if 1 <= p < len(hp) - 1:
            a, b, d = hp[p - 1], hp[p], hp[p + 1]; den = a - 2 * b + d; c.append(p + (0.5 * (a - d) / den if den != 0 else 0))
    c = np.array(c); dif = np.diff(c); step = np.median(dif)
    idx = np.concatenate([[0], np.cumsum(np.rint(dif / step))]).astype(int)
    return c, idx, step

def fit_ticks(c, idx, deg=3):
    keep = np.ones(len(c), bool)
    for _ in range(3):
        pf = np.polyfit(idx[keep], c[keep], deg); res = c - np.polyval(pf, idx); keep = np.abs(res) < max(4 * np.std(res[keep]), 1.0)
    return pf, keep, res[keep]

def strips(prof, per):
    w = int(round(1.5 * per)); hi = uniform_filter1d(maximum_filter1d(prof, w), w); lo = uniform_filter1d(minimum_filter1d(prof, w), w)
    d = prof - 0.5 * (hi + lo); up, dn = [], []
    for i in range(len(d) - 1):
        if d[i] < 0 <= d[i + 1]: up.append(i + (-d[i]) / (d[i + 1] - d[i]))
        if d[i] >= 0 > d[i + 1]: dn.append(i + d[i] / (d[i] - d[i + 1]))
    up, dn = np.array(up), np.array(dn); S, G = [], []
    for a in dn:
        b = up[up > a]
        if len(b) and 0.3 * per < b[0] - a < 0.95 * per: S.append((a, b[0]))
    for a in up:
        b = dn[dn > a]
        if len(b) and 0.1 * per < b[0] - a < 0.6 * per: G.append((a, b[0]))
    return np.array(S), np.array(G)

def index_of(cen):
    st = np.median(np.diff(cen)); return np.concatenate([[0], np.cumsum(np.rint(np.diff(cen) / st))]).astype(int)

def measure(tag):
    g = np.asarray(Image.open(f'{LV}/{tag}.png'), float); (s0, s1), (a0, a1), (b0, b1), (x0, x1) = [tuple(4 * v for v in b) for b in BANDS[tag]]
    cA, iA, stA = ticks(g[a0:a1, x0:x1].mean(0)); cA += x0; pfA, kA, rA = fit_ticks(cA, iA)
    inv = np.polyfit(np.polyval(pfA, iA[kA]), iA[kA].astype(float), 3); mmA = lambda x: np.polyval(inv, x)
    cB, iB, stB = ticks(g[b0:b1, x0:x1].mean(0), 5 * stA); cB += x0; pfB, kB, rB = fit_ticks(cB, iB, 1)
    yA, yB, yS = 0.5 * (a0 + a1), 0.5 * (b0 + b1), 0.5 * (s0 + s1)
    # px per mm at the frame centre in A (derivative of the cubic) and B (5 mm ticks, linear); linear in y -> strip band
    xc = 0.5 * (x0 + x1); kc = np.polyval(inv, xc); sA = np.polyval(np.polyder(pfA), kc); sB = pfB[0] / 5
    # scale vs image row from the strips themselves (pitch in px in 0.5-mm sub-bands over ~4 mm of rows); the ruler's
    # A/B bands disagree with it (barrel distortion, 5-mm tick indexing), so B is reported only
    rows = []
    for y in range(s0 - 3 * int(sA), s1, int(0.5 * sA)):
        Sy, _ = strips(uniform_filter1d(g[y:y + int(0.5 * sA), x0:x1].mean(0), 3), 0.78 * sA)
        if len(Sy) > 10:
            cy = Sy.mean(1); ky = index_of(cy); rows.append((y, float(np.polyfit(ky, cy, 1)[0])))
    ry = np.array(rows); gl = np.polyfit(ry[:, 0], ry[:, 1], 1); gln = gl[0] / np.polyval(gl, yS)      # d ln(scale) / d row
    sS = sA * (1 + gln * (yS - yA)); tilt = sA / sS
    prof = uniform_filter1d(g[s0:s1, x0:x1].mean(0), 3); S, G = strips(prof, 0.78 * sA); S += x0; G += x0
    S = S[(S[:, 0] > x0 + sA) & (S[:, 1] < x1 - sA)]
    cen = mmA(S.mean(1)) * tilt; wid = (mmA(S[:, 1]) - mmA(S[:, 0])) * tilt; gap = (mmA(G[:, 1]) - mmA(G[:, 0])) * tilt
    k = index_of(cen); pp = np.polyfit(k, cen, 1); rp = cen - np.polyval(pp, k); pitch = pp[0] * 1000
    perr = np.std(rp) / np.sqrt(np.sum((k - k.mean()) ** 2)) * 1000
    thirds = []
    for j in range(3):
        m = (k >= k.min() + j * (np.ptp(k) + 1) / 3) & (k < k.min() + (j + 1) * (np.ptp(k) + 1) / 3)
        thirds.append(float(np.polyfit(k[m], cen[m], 1)[0] * 1000))
    r = dict(tag=tag, px_per_mm_A=float(sA), px_per_mm_B=float(sB), px_per_mm_strips=float(sS), tilt_factor=float(tilt), dlnscale_per_1000rows=float(gln * 1000),
             n_ticks=int(kA.sum()), tick_res_um=float(np.std(rA) / sA * 1000), n_strips=int(len(cen)), span_mm=float(cen[-1] - cen[0]),
             pitch_um=float(pitch), pitch_err_um=float(perr), pitch_thirds_um=thirds, centre_res_um=float(np.std(rp) * 1000), jitter_um=float(np.std(np.diff(rp)[np.diff(k) == 1]) / np.sqrt(2) * 1000),
             width_um=float(np.median(wid) * 1000), width_sd_um=float(np.std(wid) * 1000), gap_um=float(np.median(gap) * 1000), gap_sd_um=float(np.std(gap) * 1000),
             fill=float(np.median(wid) * 1000 / pitch), pitch_px_vs_row=rows, bands=BANDS[tag])
    im = Image.fromarray(np.clip(g[s0 - 40: b1 + 20], 0, 255).astype(np.uint8)).convert('RGB'); px = im.load(); oy = s0 - 40
    for a, b in S:
        for x in (a, b):
            for y in range(s0 - oy, s1 - oy): px[int(x), y] = (255, 60, 60)
    for x in cA[kA]:
        for y in range(a0 - oy, a1 - oy): px[int(x), y] = (60, 220, 60)
    for x in cB[kB]:
        for y in range(b0 - oy, b1 - oy): px[int(x), y] = (60, 120, 255)
    im.save(f'{OVL}/{tag}.jpg', quality=85)
    return r

if __name__ == '__main__':
    out = []
    for tag in BANDS:
        r = measure(tag); out.append(r)
        print(f"{tag}: {r['px_per_mm_A']:5.1f} px/mm (B {r['px_per_mm_B']:5.1f}, strips {r['px_per_mm_strips']:5.1f}, tilt x{r['tilt_factor']:.4f})  ticks {r['n_ticks']:3d} res {r['tick_res_um']:4.1f} um  "
              f"strips {r['n_strips']:3d} / {r['span_mm']:5.1f} mm  pitch {r['pitch_um']:6.1f} ± {r['pitch_err_um']:.2f} um (thirds {', '.join(f'{t:.1f}' for t in r['pitch_thirds_um'])}; centre res {r['centre_res_um']:4.1f}, jitter {r['jitter_um']:4.1f})  "
              f"dark {r['width_um']:4.0f} ± {r['width_sd_um']:3.0f}  gap {r['gap_um']:4.0f} ± {r['gap_sd_um']:3.0f}  fill {r['fill']:.3f}")
    json.dump(out, open('measure.json', 'w'), indent=1)
