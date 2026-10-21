"""Predicted strip beat in the sim (v5a recipe, xo0 arms, all 24 files incl. tilted) vs data, same estimator as
beat_col.py but folded on column mod 40 (the sim draws x over one 31.2 mm superperiod, so columns repeat).
Amplitude of the 1/40 harmonic of per-(col mod 40) means of s1 and d1, null from permuting the 40 phase means.
Usage: beat_sim.py"""
import os, glob, numpy as np
import r5lib as R
SIM = os.path.expanduser('~/x17/nb_asym_raw/simxoff'); rng = np.random.default_rng(7)
def fold(col, v):
    ph = col % 40; n = np.bincount(ph, minlength=40); s = np.bincount(ph, v, minlength=40); m = s / np.maximum(n, 1)
    e = np.sqrt(np.bincount(ph, v * v, minlength=40) / np.maximum(n, 1) - m ** 2) / np.sqrt(np.maximum(n, 1))
    return m, e, n
def harm(m):
    k = np.arange(40); c = ((m - m.mean()) * np.exp(2j * np.pi * k / 40)).sum() * 2 / 40; return np.abs(c), np.angle(c)
def report(name, col, s1, d1):
    out = name
    for nm, v in (('s1', s1), ('d1', d1)):
        m, e, n = fold(col, v); A, ph = harm(m)
        # null: Gaussian with the per-phase errors (events independent within a phase bin)
        null = np.array([harm(m.mean() + e * rng.standard_normal(40))[0] for _ in range(2000)])
        out += f' | {nm}: A1/40 {A:.4f} ph {np.degrees(ph):+4.0f} (null 95% {np.quantile(null, .95):.4f}, p {(null >= A).mean():.3f})'
    print(out)
for arm in ('xo0', 'xom200'):
    p, s1, d1 = [], [], []
    for f in sorted(glob.glob(f'{SIM}/{arm}_*_x.npz')):
        if np.load(f)['wf'].ndim != 3: continue
        d = R.L.load(f); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
        p.append(d['p'][ok] % 512); s1.append((d['al'][1] + d['al'][-1])[ok]); d1.append(d['a1'][ok])
    p, s1, d1 = map(np.concatenate, (p, s1, d1)); col = 64 * (p // 64) + 63 - p % 64
    print(f'{arm}: n={len(col)}, distinct columns {len(np.unique(col))}, phases covered {len(np.unique(col % 40))}')
    report(f'sim {arm}', col, s1, d1)
for b in ('det3 June-15', 'det3 May', 'det2', 'det1'):
    p, s1, d1 = [], [], []
    for r in R.rows():
        if r['build'] != b or r['view'] != 'X': continue
        d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
        p.append(d['p'][ok] % 512); s1.append((d['al'][1] + d['al'][-1])[ok]); d1.append(d['a1'][ok])
    p, s1, d1 = map(np.concatenate, (p, s1, d1)); col = 64 * (p // 64) + 63 - p % 64
    m = (np.abs(s1) < 3) & (np.abs(d1) < 3)
    report(f'data {b}', col[m], s1[m], d1[m])
