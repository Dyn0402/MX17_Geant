#!/usr/bin/env python3
"""X/Y charge balance in det3 data binned by the TRUE (M3) angle, matched to each sim point (2026-10-17).
Same observable as xy_balance_sim_vs_data.py (FeuReader, strips >= 5 sigma and >= 100 ADC, q = sum of strip
peaks, pk = max, unsat = no strip >= 3500). Angle convention as ledger_m3.py: X-view angle = M3 ay,
Y-view angle = -M3 ax; pk/q per view = peak-strip fraction of that view's sum (the chamber sits rotated by 90 deg), window |angle - a| < 3 deg in each projection.
Usage: xy_balance_m3.py SIMROOT   (SIMROOT/<point>_<ang>/ holds the sim decoded files; points via env PTS)"""
import os, sys, glob
import numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qa_config import get_config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
setup_paths()
from M3RefTracking import M3RefTracking
from xy_balance_sim_vs_data import per_view

cfg = get_config('sat_det3')
DEC = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'
rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
x0, y0, ev = map(np.asarray, rays.get_xy_positions(702.0)); x1, y1, _ = map(np.asarray, rays.get_xy_positions(802.0))
ax = dict(zip(ev.astype(int), np.degrees(np.arctan((x1 - x0) / 100))))
ay = dict(zip(ev.astype(int), np.degrees(np.arctan((y1 - y0) / 100))))

def collect(files07):
    rows = []
    for fx in files07:
        fy = fx[:-8] + '_08.root'
        (qx, px), (qy, py) = per_view(fx), per_view(fy)
        rows += [(k, qx[k], qy[k], px[k], py[k]) for k in set(qx) & set(qy)]
    a = np.array(rows, dtype=float)
    return a[:, 0].astype(int), a[:, 1], a[:, 2], a[:, 3], a[:, 4]

def summ(qx, qy, px, py):
    out = []
    for m in (np.ones(len(qx), bool), (px < 3500) & (py < 3500)):
        out.append(f'{m.sum():5d} {np.median(qx[m] / (qx[m] + qy[m])):6.3f} {np.median(px[m]) / np.median(py[m]):7.3f}'
                   f' {np.median(px[m] / qx[m]):6.3f} {np.median(py[m] / qy[m]):6.3f}')
    return ' | '.join(out)

e, qx, qy, px, py = collect(sorted(glob.glob(f'{DEC}/*_07.root')))
good = np.array([k in ax and np.isfinite(ax[k]) and np.isfinite(ay[k]) for k in e])
e, qx, qy, px, py = e[good], qx[good], qy[good], px[good], py[good]
thX = np.array([ay[k] for k in e]); thY = np.array([-ax[k] for k in e])
ANG = {'th00': (0, 0), 'thx10': (10, 0), 'thx20': (20, 0), 'thy10': (0, 10), 'thy20': (0, 20)}
root = sys.argv[1]; pts = os.environ.get('PTS', 'verdict1').split()
print(f'{len(e)} data events with an M3 slope and both views over threshold')
print(f'{"":14s} {"all: n   f   pkX/pkY pkX/qX pkY/qY":>36s} | {"unsat: same":>36s}')
for ang, (a, b) in ANG.items():
    m = (np.abs(np.abs(thX) - a) < 3) & (np.abs(np.abs(thY) - b) < 3)
    print(f'data {ang:9s} {summ(qx[m], qy[m], px[m], py[m])}')
    for p in pts:
        d = f'{root}/{p}_{ang}'
        if glob.glob(f'{d}/*_07.root'):
            _, sx, sy, spx, spy = collect(sorted(glob.glob(f'{d}/*_07.root')))
            print(f'  {p:12s} {summ(sx, sy, spx, spy)}', flush=True)
