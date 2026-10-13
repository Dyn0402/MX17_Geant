#!/usr/bin/env python3
"""Why does the sim fire twice as many strips >= 5 sigma as the data (X 12 vs 6, Y 13 vs 8 at th00)?
Runs on the laptop (yoga) next to the T14 outputs. For the events each T14 point used (its wf_*.parquet),
take every channel in position order and record, vs signed offset from the peak strip:
  amp/peak (all channels, not just those over threshold) and amp/noise;
plus the 5-sigma exceedance rate on channels >= 40 strips from the peak (noise only).
Usage: profile_tails.py <t14 dir> <sim dir> [data dir]"""
import sys, os, glob, json, numpy as np, pandas as pd
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17'), os.path.expanduser('~/PycharmProjects/nTof_x17/mx17_sim_wft')]
from wft.io import FeuReader, strip_position_map
from t14_compare import _leg_cfg
t14, sim_dir = sys.argv[1], sys.argv[2]
dat_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.expanduser('~/x17/response_sim/stageB_w2/data_leg_satdet3')
cfg = _leg_cfg(sim_dir, 'sat_det3'); pm = strip_position_map(cfg)
OFF = np.arange(-30, 31)
def prof(d, feu, ids):
    pos = pm[feu]; valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]
    R, S, far = [], [], []
    for path in sorted(glob.glob(f'{d}/decoded_root/*_{feu:02d}.root')):
        rdr = FeuReader(path); noise = np.where(rdr.noise > 0, rdr.noise, np.inf)
        want = ids & set(int(e) for e in rdr.event_ids)
        if not want: continue
        for eid, _f, w in rdr.iter_events(want):
            amp = w.max(axis=1)[order]; sig = amp / noise[order]
            ip = int(np.argmax(np.where(sig >= 5, amp, -np.inf)))
            r = np.full(len(OFF), np.nan); s = np.full(len(OFF), np.nan)
            j = ip + OFF; ok = (j >= 0) & (j < len(amp)); r[ok] = amp[j[ok]] / amp[ip]; s[ok] = sig[j[ok]]
            R.append(r); S.append(s)
            k = np.abs(np.arange(len(amp)) - ip) >= 40; far.append((sig[k] >= 5).mean())
    return np.array(R), np.array(S), np.array(far)
cfgv = {'x': cfg.MX17_FEU_X, 'y': cfg.MX17_FEU_Y}
out = {}
for v in 'xy':
    for leg, d in (('sim', sim_dir), ('data', dat_dir)):
        ids = set(pd.read_parquet(f'{t14}/wf_{leg}_{v}.parquet').event_id.astype(int))
        R, S, far = prof(d, cfgv[v], ids)
        out[f'{leg}_{v}_R'] = R; out[f'{leg}_{v}_S'] = S; out[f'{leg}_{v}_far'] = far
        nover = np.nansum(S >= 5, axis=1)
        print(f'{v} {leg:4s} n={len(R)} far-channel 5σ rate {far.mean():.4f}  n_over(±30) med {np.median(nover):.0f}')
        print('   off   ' + ' '.join(f'{o:6d}' for o in range(-8, 9)))
        print('   amp/pk' + ' '.join(f'{x:6.3f}' for x in np.nanmedian(R, 0)[22:39]))
        print('   mean  ' + ' '.join(f'{x:6.3f}' for x in np.nanmean(R, 0)[22:39]))
        print('   σ med ' + ' '.join(f'{x:6.1f}' for x in np.nanmedian(S, 0)[22:39]))
np.savez_compressed(os.path.join(t14, 'profile_tails.npz'), off=OFF, **out)
