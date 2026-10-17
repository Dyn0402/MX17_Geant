#!/usr/bin/env python3
"""How far an inclined track spreads across the strips, det3 data by M3 angle vs sim (2026-10-17).
Same strip selection as xy_balance_m3.py (>= 5 sigma and >= 100 ADC). Per event and view: n strips,
span = last - first strip position [mm] over threshold, and the charge-weighted RMS width [mm].
Usage: xy_extent_m3.py SIMROOT   (env PTS = sim point names, default verdict1)"""
import os, sys, glob
import numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
from qa_config import get_config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
setup_paths()
from M3RefTracking import M3RefTracking
from wft.io import FeuReader, strip_position_map

cfg = get_config('sat_det3'); PM = strip_position_map(cfg)
DEC = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'

def extent(fn, feu):
    out = {}
    pos = PM[feu]
    r = FeuReader(fn); noise = np.where(r.noise > 0, r.noise, np.inf)
    for e, _t, w in r.iter_events(set(int(x) for x in r.event_ids)):
        amp = w.max(axis=1); over = (amp / noise >= 5) & (amp >= 100) & np.isfinite(pos)
        if over.sum() < 1: continue
        p, a = pos[over], amp[over]
        mu = (p * a).sum() / a.sum()
        out[int(e)] = (over.sum(), p.max() - p.min(), np.sqrt((a * (p - mu) ** 2).sum() / a.sum()), a.max() >= 3500)
    return out

def collect(d):
    X, Y = {}, {}
    for fx in sorted(glob.glob(f'{d}/*_07.root')):
        X.update({(fx, k): v for k, v in extent(fx, 7).items()})
        Y.update({(fx, k): v for k, v in extent(fx[:-8] + '_08.root', 8).items()})
    keys = [k for k in X if k in Y]
    return [k[1] for k in keys], np.array([X[k] for k in keys], float), np.array([Y[k] for k in keys], float)

def summ(x, y):
    uns = (x[:, 3] == 0) & (y[:, 3] == 0)
    x, y = x[uns], y[uns]
    return (f'{len(x):5d}  X n {np.median(x[:, 0]):4.1f} span {np.median(x[:, 1]):5.2f} rms {np.median(x[:, 2]):5.2f}'
            f'   Y n {np.median(y[:, 0]):4.1f} span {np.median(y[:, 1]):5.2f} rms {np.median(y[:, 2]):5.2f}')

rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
x0, y0, ev = map(np.asarray, rays.get_xy_positions(702.0)); x1, y1, _ = map(np.asarray, rays.get_xy_positions(802.0))
ax = dict(zip(ev.astype(int), np.degrees(np.arctan((x1 - x0) / 100))))
ay = dict(zip(ev.astype(int), np.degrees(np.arctan((y1 - y0) / 100))))
e, dx, dy = collect(DEC)
g = np.array([k in ax and np.isfinite(ax[k]) and np.isfinite(ay[k]) for k in e])
e = np.array(e)[g]; dx, dy = dx[g], dy[g]
thX = np.array([ay[k] for k in e]); thY = np.array([-ax[k] for k in e])
ANG = {'th00': (0, 0), 'thx10': (10, 0), 'thx20': (20, 0), 'thy10': (0, 10), 'thy20': (0, 20)}
print('unsaturated events only; span/rms in mm')
for ang, (a, b) in ANG.items():
    m = (np.abs(np.abs(thX) - a) < 3) & (np.abs(np.abs(thY) - b) < 3)
    print(f'data {ang:8s} {summ(dx[m], dy[m])}')
    for p in os.environ.get('PTS', 'verdict1').split():
        d = f'{sys.argv[1]}/{p}_{ang}'
        if glob.glob(f'{d}/*_07.root'):
            _, sx, sy = collect(d)
            print(f'  {p:11s} {summ(sx, sy)}', flush=True)
