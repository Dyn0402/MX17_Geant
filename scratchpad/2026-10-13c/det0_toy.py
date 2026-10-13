#!/usr/bin/env python3
"""det0 pulse-shape test done properly (FRESH_LOOK §13.3, revised 2026-10-13 night).

det0 has a 6 mm drift gap (Feb logbook: "new 6mm drift mounted", runs 69+, never changed),
so at 800 V the drift field is ~1.33 kV/cm. Magboltz (scratchpad/2026-10-13c/magboltz,
extrapolated 1.0 -> 1.33 kV/cm): Ar/iso 95/5 v ~27.5 um/ns (column 6 mm = 218 ns),
He/C2H6 96.5/3.5 v ~15 um/ns (400 ns). He/ethane is sparse (~7 clusters/cm, ~11 e/cm,
Sauli-table values), Ar/iso dense (91 e/cm, 3.83 e/cluster), so the He pulse is NOT a point
response and the two columns differ: model both.

Toy per gas (vertical tracks |tan| < 0.05): clusters along 6 mm, drift at v with D_L,
per-electron gamma(1.5) gain, prompt lateral sigma 0.2 mm onto 0.78 mm strips, kernel =
toy_pulse.electron_kernel(shaper(tp), f, T), 60 ns random-phase sampling, noise 8 ADC,
scale set so the toy's all-event median peak matches data (He 179, Ar 197 ADC), then the
data selection (300-3500 ADC) and the same 50 %-edge-aligned mean pulse as he_template_fit.
Step 1: fit the shaper time scale tp on He (its ions are fast: T 30-60 ns irrelevant).
Step 2: with that tp, scan Ar's (f, T)."""
import sys, math, numpy as np
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
import toy_pulse as TP
import he_template_fit as H
from rise_all import cross
from concurrent.futures import ProcessPoolExecutor
GAS = {'He': dict(v=0.015, dL=0.027, ncl=7.1, csize=11.4 / 7.1, med=179.),
       'Ar': dict(v=0.0275, dL=0.0214, ncl=91.2 / 3.83, csize=3.83, med=197.)}
GAP, PITCH, SIGP, NOISE = 6.0, 0.78, 0.2, 8.0
def sim(gas, tp, f, T, n=3000, seed=1):
    g = GAS[gas]; rng = np.random.default_rng(seed)
    h = TP.shaper('manual_1', tp); k, _ = TP.electron_kernel(h, f, T, 0.0, 230.)
    K = np.fft.rfft(k, 4096)
    strips = np.arange(-4, 5); lo, hi = (strips - .5) * PITCH, (strips + .5) * PITCH
    cdf = lambda u: 0.5 * (1 + np.vectorize(math.erf)(u / (SIGP * math.sqrt(2))))
    peaks, wfs = [], []
    for _ in range(n):
        x0 = rng.uniform(-.5, .5) * PITCH; tan = rng.uniform(-.05, .05)
        nc = rng.poisson(g['ncl'] * GAP / 10)
        if nc == 0: peaks.append(0); wfs.append(None); continue
        z = np.repeat(rng.uniform(0, GAP, nc), rng.geometric(1 / g['csize'], nc))
        t = 200 + z / g['v'] + rng.normal(0, 1, z.size) * g['dL'] * np.sqrt(z / 10) * 10 / g['v']
        x = x0 + z * tan
        w = cdf(hi[:, None] - x[None, :]) - cdf(lo[:, None] - x[None, :])
        q = rng.gamma(1.5, 1 / 1.5, z.size)
        it = np.clip(np.round(t).astype(int), 0, 2047)
        dep = np.array([np.bincount(it, weights=q * w[s], minlength=2048)[:2048] for s in range(len(strips))])
        wf = np.fft.irfft(np.fft.rfft(dep, 4096) * K, 4096)[:, :2048]
        ph = rng.uniform(0, 60); ts = (ph + 60 * np.arange(32)).astype(int)
        S = wf[:, np.clip(ts, 0, 2047)]
        pk = S.max(1).argmax(); peaks.append(S[pk].max()); wfs.append(S[pk])
    peaks = np.array(peaks); scale = (g['med'] - 2.5 * NOISE) / np.median(peaks)  # noise lifts the max ~2.5 sigma
    out = []
    for p, w in zip(peaks, wfs):
        if w is None: continue
        ww = w * scale + rng.normal(0, NOISE, 32)
        if 300 <= ww.max() < 3500: out.append(ww)
    return np.array(out)
def job(c):
    gas, tp, f, T = c
    W = sim(gas, tp, f, T)
    g, a = H.fine_avg(W)
    return c, g, a, len(W)
if __name__ == '__main__':
    He = H.collect('det0__*He_HV_Scan_4-2*', 495, 515, 0, .05); Ar = H.collect('det0__*Ar_Iso*', 440, 480, 0, .05)
    g, dHe = H.fine_avg(He); _, dAr = H.fine_avg(Ar); m = (g >= -150) & (g <= 700)
    rms = lambda a, d: float(np.sqrt(np.mean((a - d)[m] ** 2)))
    with ProcessPoolExecutor(12) as ex:
        R = list(ex.map(job, [('He', tp, 0.9, Th) for tp in (140, 155, 166, 180, 195) for Th in (30., 60.)]))
    print('Step 1, He/C2H6 (data n=%d): shaper scale vs He ion transit' % len(He))
    best = None
    for c, gg, a, n in R:
        r = rms(np.interp(g, gg, a), dHe); print(f'   tp {c[1]:3d} T_He {c[3]:3.0f}: rms {r:.4f}  (toy n={n})')
        if best is None or r < best[0]: best = (r, c[1])
    tp = best[1]; print(f'   -> tp = {tp} ns')
    cfg = [('Ar', tp, f, T) for f in (0.9,) for T in (60., 100., 140., 180., 220., 280., 340.)]
    cfg += [('Ar', tp, f, 340.) for f in (0.0, 0.2, 0.5)]
    cfg += [('Ar', tp2, 0.9, T) for tp2 in (tp - 15, tp + 15) for T in (140., 340.)]
    with ProcessPoolExecutor(12) as ex:
        R = list(ex.map(job, cfg))
    print('Step 2, Ar/iso (data n=%d), He-calibrated shaper tp=%d:' % (len(Ar), tp))
    for c, gg, a, n in R:
        print(f'   tp {c[1]:3d} f {c[2]:.1f} T {c[3]:3.0f}: rms {rms(np.interp(g, gg, a), dAr):.4f}  (toy n={n})')
