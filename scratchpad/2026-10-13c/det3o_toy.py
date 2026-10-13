#!/usr/bin/env python3
"""Old det3 (bulked 5 May, 30 mm drift, 900 V = 300 V/cm), same chamber in Ar/iC4H10 95/5
and Ar/CF4 90/10: fit the mean peak-strip pulse of vertical, ~10 deg and ~20 deg tracks with
the det0_toy machinery generalised to a 30 mm gap and inclined tracks (FRESH_LOOK §15).

Free per gas: drift velocity v (the May gas's water content is unknown; Magboltz dry at
300 V/cm gives Ar/iso 40, Ar/CF4 98 um/ns) and ion transit T (f_ion 0.9). Shared: the shaper
time scale tp. Diffusion from Magboltz (dry). Ionisation ~91 e/cm, 3.83 e/cluster for both
(CF4 adds little at 10 %). The vertical pulse length fixes v; the inclined rises then fix T."""
import sys, math, glob, numpy as np, pandas as pd
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
import toy_pulse as TP
import he_template_fit as H
from concurrent.futures import ProcessPoolExecutor
C = '/home/dylan/.claude/jobs/913ff27a/tmp/c/'
GAP, PITCH, SIGP, NOISE = 30.0, 0.78, 0.2, 9.0
DIFF = {'iso': (0.0263, 0.0480), 'cf4': (0.0168, 0.0211)}   # dL, dT [sqrt cm] at 300 V/cm
BINS = {'vert': (0.0, 0.05), '10': (0.15, 0.22), '20': (0.33, 0.40)}
MED = {}   # filled from data: median peak per (gas, bin)
def sim(gas, tp, f, T, v, tlo, thi, med, n=2500, seed=3):
    rng = np.random.default_rng(seed); dL, dT = DIFF[gas]
    h = TP.shaper('manual_1', tp); k, _ = TP.electron_kernel(h, f, T, 0.0, 230.)
    K = np.fft.rfft(k, 8192)
    strips = np.arange(-24, 25); lo, hi = (strips - .5) * PITCH, (strips + .5) * PITCH
    cdf = lambda u: 0.5 * (1 + np.vectorize(math.erf)(u / (SIGP * math.sqrt(2))))
    peaks, wfs = [], []
    for _ in range(n):
        tan = rng.uniform(tlo, thi) * rng.choice([-1, 1]); x0 = rng.uniform(-.5, .5) * PITCH - 0.5 * GAP * tan
        nc = rng.poisson(91.2 / 3.83 * GAP / 10)
        z = np.repeat(rng.uniform(0, GAP, nc), rng.geometric(1 / 3.83, nc))
        t = 200 + z / v + rng.normal(0, 1, z.size) * dL * np.sqrt(z / 10) * 10 / v
        x = x0 + z * tan + rng.normal(0, 1, z.size) * dT * np.sqrt(z / 10) * 10
        w = cdf(hi[:, None] - x[None, :]) - cdf(lo[:, None] - x[None, :])
        q = rng.gamma(1.5, 1 / 1.5, z.size)
        it = np.clip(np.round(t).astype(int), 0, 4095)
        keep = w.sum(1) > 1e-3
        dep = np.array([np.bincount(it, weights=q * w[s], minlength=4096)[:4096] for s in np.flatnonzero(keep)])
        wf = np.fft.irfft(np.fft.rfft(dep, 8192) * K, 8192)[:, :4096]
        ph = rng.uniform(0, 60); ts = (ph + 60 * np.arange(32)).astype(int)
        S = wf[:, np.clip(ts, 0, 4095)]
        pk = S.max(1).argmax(); peaks.append(S[pk].max()); wfs.append(S[pk])
    scale = (med - 2.5 * NOISE) / np.median(peaks)
    out = [w * scale + rng.normal(0, NOISE, 32) for w in wfs]
    out = np.array([w for w in out if 300 <= w.max() < 3500])
    return out
def collect(pat, vlo, vhi, tlo, thi, fid=(-200, 100, -200, 50)):
    P, A = [], []
    for f in sorted(glob.glob(C + pat + '.parquet')):
        d = pd.read_parquet(f); W = np.load(f.replace('.parquet', '_wf.npy'))[:, 2, :].astype(float)
        t = np.abs(np.where(d.view == 'x', d.tan_y, d.tan_x))
        m0 = ((d.volt >= vlo) & (d.volt <= vhi) & (d.ref_x > fid[0]) & (d.ref_x < fid[1]) & (d.ref_y > fid[2]) & (d.ref_y < fid[3])
              & (t >= tlo) & (t < thi) & (d.view == 'x')).values
        A.append(d.peak_amp.values[m0]); m = m0 & (d.peak_amp >= 300).values & (d.peak_amp < 3500).values
        P.append(np.nan_to_num(W[m]))
    return np.concatenate(P), float(np.median(np.concatenate(A)))
DATA = {'iso': ('det3old__*HV_Scan_5-5*', 460, 490), 'cf4': ('det3old__*ArCF4*', 600, 630)}
def job(c):
    gas, tp, T, v, b, med = c
    W = sim(gas, tp, 0.9, T, v / 1000., *BINS[b], med)
    g, a = H.fine_avg(W, -400, 1500); return c, g, a, len(W)
if __name__ == '__main__':
    data = {}
    for gas, (pat, vlo, vhi) in DATA.items():
        for b, (tlo, thi) in BINS.items():
            W, med = collect(pat, vlo, vhi, tlo, thi); g, a = H.fine_avg(W, -400, 1500)
            data[(gas, b)] = (g, a, med, len(W))
            print(f'data {gas} {b}: n={len(W)} median peak (all) {med:.0f}')
    g = data[('iso', 'vert')][0]; m = (g >= -200) & (g <= 1200)
    grid = {'iso': dict(v=(30, 35, 40), T=(140, 200, 260, 352)), 'cf4': dict(v=(60, 80, 100), T=(80, 140, 200, 260, 352))}
    cfg = [(gas, tp, T, v, b, data[(gas, b)][2]) for gas in grid for tp in (150, 166, 180)
           for v in grid[gas]['v'] for T in grid[gas]['T'] for b in BINS]
    with ProcessPoolExecutor(12) as ex:
        R = list(ex.map(job, cfg))
    res = {}
    for (gas, tp, T, v, b, _), gg, a, n in R:
        d = data[(gas, b)][1]
        res.setdefault((gas, tp, T, v), {})[b] = float(np.sqrt(np.mean((np.interp(g, gg, a) - d)[m] ** 2)))
    rows = [dict(gas=k[0], tp=k[1], T=k[2], v=k[3], **{f'rms_{b}': r[b] for b in BINS}, rms_all=float(np.sqrt(np.mean([r[b] ** 2 for b in BINS]))))
            for k, r in res.items()]
    df = pd.DataFrame(rows).sort_values(['gas', 'rms_all']); df.to_csv('det3o_toy.csv', index=False)
    pd.set_option('display.width', 200)
    for gas in grid:
        print(f'== {gas}: best 8'); print(df[df.gas == gas].head(8).round(4).to_string(index=False))
        for tp in (150, 166, 180):
            s = df[(df.gas == gas) & (df.tp == tp)]
            print(f'   tp {tp}: best per T ->', ' '.join(f'T{T}:{s[s["T"] == T].rms_all.min():.4f}' for T in grid[gas]['T']))
