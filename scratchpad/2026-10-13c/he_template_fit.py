#!/usr/bin/env python3
"""Shaper-free ion test on det0 (FRESH_LOOK §13): the same chamber and electronics ran
He/C2H6 96.5/3.5, whose ions are ~10x faster than any Ar ion, so its mean peak-strip pulse
is an empirical point response (shaper (x) electron signal (x) the short column det0 sees).
Model: Ar/iso mean pulse = He mean pulse (x) K(t) (x) gauss(sigma), with
  K = (1 - f) delta + f/T rect(0, T)   (parallel-plate ion current, transit T)
and sigma a free extra arrival smear (Ar drift column vs He's). Both legs averaged the same way:
peak strip, 300-3500 ADC, normalised to its max, aligned on the leading 50 % crossing
(sub-sample, random-phase sampling fills the grid), track bin |tan| < 0.05 (also 0.15-0.25).
The He response is first unfolded to a fine (5 ns) continuous template by averaging the
random-phase samples on a fine grid."""
import glob, numpy as np, pandas as pd
from scipy.optimize import least_squares
from rise_all import cross
C = '/home/dylan/.claude/jobs/913ff27a/tmp/c/'
FINE = 5.0
def collect(pat, vlo, vhi, tlo, thi, fid=(50, 150, -200, 150)):
    P = []
    for f in sorted(glob.glob(C + pat + '.parquet')):
        d = pd.read_parquet(f); W = np.load(f.replace('.parquet', '_wf.npy'))[:, 2, :].astype(float)
        t = np.abs(np.where(d.view == 'x', d.tan_y, d.tan_x))
        m = ((d.volt >= vlo) & (d.volt <= vhi) & (d.ref_x > fid[0]) & (d.ref_x < fid[1]) & (d.ref_y > fid[2]) & (d.ref_y < fid[3])
             & (d.peak_amp >= 300) & (d.peak_amp < 3500) & (t >= tlo) & (t < thi) & ~d.subrun.str.startswith('final')).values
        P.append(np.nan_to_num(W[m]))
    return np.concatenate(P)
def fine_avg(W, lo=-300, hi=900):
    """continuous mean pulse (ns grid), aligned so the leading 50 % crossing is t=0"""
    g = np.arange(lo, hi, FINE); acc = np.zeros_like(g); cnt = np.zeros_like(g)
    for w in W:
        w = w / w.max(); l = cross(w, .5, -1)
        if np.isnan(l): continue
        t = (np.arange(32) - l) * 60.0
        k = np.searchsorted(g, t); ok = (k > 0) & (k < len(g))
        np.add.at(acc, k[ok], w[ok]); np.add.at(cnt, k[ok], 1)
    m = cnt > 0
    return g, np.interp(g, g[m], acc[m] / cnt[m])
def model(g, tmpl, f, T, sig):
    n = int(max(T, FINE) / FINE)
    K = np.zeros(int(1500 / FINE)); K[0] = 1 - f; K[:n] += f / n
    if sig > 1:
        x = np.arange(-4 * sig, 4 * sig + FINE, FINE); gk = np.exp(-0.5 * (x / sig) ** 2); gk /= gk.sum()
        y = np.convolve(np.convolve(tmpl, K)[:len(tmpl)], gk, mode='same')
    else:
        y = np.convolve(tmpl, K)[:len(tmpl)]
    y = y / y.max()
    # re-align on the leading 50 % crossing
    i = np.argmax(y); j = np.flatnonzero(y[:i] < 0.5)[-1]
    t50 = g[j] + (0.5 - y[j]) / (y[j + 1] - y[j]) * FINE
    return np.interp(g, g - t50, y)
if __name__ == '__main__':
    for tlo, thi, lab in ((0, .05, 'vertical'), (.15, .25, '~10 deg')):
        He = collect('det0__*He_HV_Scan_4-2*', 495, 515, tlo, thi)
        Ar = collect('det0__*Ar_Iso*', 440, 480, tlo, thi)
        g, th = fine_avg(He); _, ta = fine_avg(Ar)
        w = (g >= -150) & (g <= 700)
        print(f'== {lab}: He n={len(He)}, Ar n={len(Ar)}  (fit window -150..700 ns from the 50 % edge)')
        print('   t[ns]   ' + ' '.join(f'{x:5.0f}' for x in g[w][::12]))
        print('   He      ' + ' '.join(f'{x:5.2f}' for x in th[w][::12]))
        print('   Ar/iso  ' + ' '.join(f'{x:5.2f}' for x in ta[w][::12]))
        res = lambda p: (model(g, th, p[0], p[1], p[2]) - ta)[w]
        best = min((least_squares(res, [f0, T0, 20], bounds=([0, 5, 0], [1, 600, 150])) for f0 in (.2, .5, .9) for T0 in (60, 150, 300)), key=lambda r: r.cost)
        f, T, s = best.x; rms = np.sqrt(np.mean(best.fun ** 2))
        print(f'   best: f_ion {f:.2f}  T_ion {T:.0f} ns  extra smear {s:.0f} ns   rms {rms:.4f}   (He alone rms {np.sqrt(np.mean((th - ta)[w] ** 2)):.4f})')
        for f0, T0 in ((0.9, 340), (0.9, 200), (0.9, 140), (0.2, 340), (0.5, 340)):
            r = least_squares(lambda p: res([f0, T0, p[0]]), [20], bounds=([0], [150]))
            print(f'   fixed f {f0:.1f} T {T0:3.0f}: smear {r.x[0]:3.0f} ns  rms {np.sqrt(np.mean(r.fun ** 2)):.4f}')
