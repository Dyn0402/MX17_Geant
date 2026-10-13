#!/usr/bin/env python3
"""Electronics compression vs avalanche space charge. Both would suppress the data's vertical
(dense) events relative to the linear sim. A DREAM response compressing below the rail acts per
channel, so the +-1 neighbour / peak ratio would RISE with peak amplitude in the data (peak
squeezed more than its neighbours), and stay flat in the linear sim. Space charge acts on the
avalanche before sharing, so the ratio stays flat in both. Runs on yoga: <t14 dir> <sim dir>"""
import sys, os, glob, numpy as np, pandas as pd
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17'), os.path.expanduser('~/PycharmProjects/nTof_x17/mx17_sim_wft')]
from wft.io import FeuReader, strip_position_map
from t14_compare import _leg_cfg
t14, sim_dir = sys.argv[1], sys.argv[2]
dat_dir = os.path.expanduser('~/x17/response_sim/stageB_w2/data_leg_satdet3')
cfg = _leg_cfg(sim_dir, 'sat_det3'); pm = strip_position_map(cfg)
def rows(d, feu, ids):
    pos = pm[feu]; valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]
    out = []
    for path in sorted(glob.glob(f'{d}/decoded_root/*_{feu:02d}.root')):
        rdr = FeuReader(path); want = ids & set(int(e) for e in rdr.event_ids)
        if not want: continue
        for eid, _f, w in rdr.iter_events(want):
            W = w[order]; amp = W.max(axis=1); ip = int(np.argmax(amp))
            if ip < 2 or ip > len(amp) - 3: continue
            ipk = int(np.argmax(W[ip]))   # neighbours sampled AT the peak strip's peak time, and at their own max
            out.append((amp[ip], (amp[ip - 1] + amp[ip + 1]) / amp[ip], (W[ip - 1, ipk] + W[ip + 1, ipk]) / amp[ip],
                        (amp[ip - 2] + amp[ip + 2]) / amp[ip], amp[ip - 2:ip + 3].sum()))
    return pd.DataFrame(out, columns=['pk', 'nb1', 'nb1_at_pk', 'nb2', 'q5'])
B = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4100]
for v, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
    for leg, d in (('sim', sim_dir), ('data', dat_dir)):
        ids = set(pd.read_parquet(f'{t14}/wf_{leg}_{v}.parquet').event_id.astype(int))
        t = rows(d, feu, ids); t.to_parquet(os.path.expanduser(f'~/nbratio_{os.path.basename(t14)}_{leg}_{v}.parquet'))
        print(f'{v} {leg}: n {len(t)}   peak bin: n | (±1 max)/pk | (±1 at pk time)/pk | (±2)/pk')
        for lo, hi in zip(B[:-1], B[1:]):
            m = (t.pk >= lo) & (t.pk < hi)
            if m.sum() < 15: continue
            x = t[m]; print(f'   {lo:4d}-{hi:<4d} {m.sum():5d} {x.nb1.median():7.3f} {x.nb1_at_pk.median():7.3f} {x.nb2.median():7.3f}')
