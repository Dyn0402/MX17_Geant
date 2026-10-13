#!/usr/bin/env python3
"""Unbiased amplitude ledger (2026-10-14). T14 compares sim (quality cut, no angle window) against
data (quality cut + reco-angle window); both cuts depend on amplitude in the data (FRESH_LOOK §6,
peak_vs_theta.py), and the reco angle near 0 is not slope-reliable. Here:
  data = target run, M3-good rays in the fiducial box, binned by the TRUE projected angle from M3
         (|view angle - a| < 3 deg, |other projection| < 6 deg), no det3 reco, no quality cut;
  sim  = every simulated event (true angle known), no reco, no quality cut.
Observable: the T14 peak-strip amplitude (FeuReader, max over strips >= 5 sigma; threshold-free on the
data side, identical for any real signal). Quantiles p25/p50 data/sim and the >=3500 fraction.
Sim peaks are cached per decoded dir in <dir>/peak_all_<feu>.npy. Runs on yoga.
Usage: ledger_m3.py [glob of stageB_w2 sim dirs, default all w2_rho2M_*/th*]"""
import os, sys, glob, re, numpy as np, pandas as pd
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17')]
from wft.io import FeuReader
SB = os.path.expanduser('~/x17/response_sim/stageB_w2')
pk = pd.read_parquet(os.path.expanduser('~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet'))
pk = pk[(pk.ref_x > -200) & (pk.ref_x < 120) & (pk.ref_y.abs() < 170)]
pk['th'] = np.where(pk.view == 'x', pk.ay, -pk.ax); pk['oth'] = np.where(pk.view == 'x', pk.ax, -pk.ay)
def data(view, a):
    d = pk[(pk.view == view) & ((pk.th.abs() - a).abs() < 3) & (pk.oth.abs() < 6)] if a else \
        pk[(pk.view == view) & (pk.th.abs() < 3) & (pk.oth.abs() < 6)]
    return d.peak_amp.values
def sim_peaks(d, feu):
    c = f'{d}/peak_all_{feu:02d}.npy'
    if os.path.exists(c): return np.load(c)
    out = []
    for f in sorted(glob.glob(f'{d}/decoded_root/*_{feu:02d}.root')):
        r = FeuReader(f); noise = np.where(r.noise > 0, r.noise, np.inf)
        for _e, _t, w in r.iter_events(set(int(e) for e in r.event_ids)):
            amp = w.max(axis=1); over = amp / noise >= 5
            out.append(amp[over].max() if over.any() else 0.0)
    out = np.array(out); np.save(c, out); return out
ANG = {'th00': (0, 0), 'thx10': (10, 0), 'thx20': (20, 0), 'thy10': (0, 10), 'thy20': (0, 20)}
dirs = sorted(glob.glob(f'{SB}/{sys.argv[1] if len(sys.argv) > 1 else "w2_rho2M_*/th*"}'))
print(f'{"sim point":28s} {"view":4s} {"ang":>3s} {"n d/s":>10s} {"p50 d/s":>8s} {"p25 d/s":>8s} {"sat d / s":>12s}')
rows = []
for d in dirs:
    ang = os.path.basename(d); kind = os.path.basename(os.path.dirname(d)).replace('w2_rho2M_', '')
    if ang not in ANG or not glob.glob(f'{d}/decoded_root/*_08.root'): continue
    for view, feu, a in (('x', 7, ANG[ang][0]), ('y', 8, ANG[ang][1])):
        s = sim_peaks(d, feu); s = s[s > 0]; x = data(view, a)
        r50 = np.median(x) / np.median(s); r25 = np.quantile(x, .25) / np.quantile(s, .25)
        rows.append(dict(kind=kind, ang=ang, view=view, a=a, n_d=len(x), n_s=len(s), r50=r50, r25=r25,
                         sat_d=(x >= 3500).mean(), sat_s=(s >= 3500).mean(), p50_s=np.median(s), p50_d=np.median(x)))
        print(f'{kind + "/" + ang:28s} {view:4s} {a:3d} {len(x):5d}/{len(s):<5d} {r50:8.3f} {r25:8.3f} {(x >= 3500).mean():6.3f}/{(s >= 3500).mean():5.3f}', flush=True)
out = os.path.expanduser('~/ledger_m3.csv')       # merge: a partial glob must not clobber earlier rows
new = pd.DataFrame(rows)
if os.path.exists(out): new = pd.concat([pd.read_csv(out), new]).drop_duplicates(['kind', 'ang', 'view'], keep='last')
new.to_csv(out, index=False)
