#!/usr/bin/env python3
"""Closure tables from the events parquets closure.py wrote (items 1 and 2), in the form data_flip.py measures det3:
per angle point and view, the fraction of direction flips, the core tan ratio (|ratio - 1| <= 0.25), the core
sigma_theta, and the median |reco line - true line| at depths 3/9/15/21/26 mm (core events).
Usage: closure_table.py CLOSURE_DIR CONFIG TAG [STAGEA_DIR]"""
import sys, os, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_calib import truth_lines, ANG

C, cfg, tag = sys.argv[1:4]
SA = sys.argv[4] if len(sys.argv) > 4 else '/media/dylan/data/claude_tmp/stageA_e17'
DEP = (3, 9, 15, 21, 26)
rs = lambda a: float(1.4826 * np.median(np.abs(a - np.median(a)))) if len(a) else np.nan
out = {}
print(f'| point | view | n | flip | core ratio | core σθ [°] | Δθ med [°] | ' + ' | '.join(f'|r| z{z}' for z in DEP) + ' |')
print('|---' * (7 + len(DEP)) + '|')
for ang, (ax, ay) in ANG.items():
    f = f'{C}/events_{cfg}_{ang}_{tag}.parquet'
    if not os.path.exists(f):
        continue
    d = pd.read_parquet(f)
    L = truth_lines(f'{SA}/mx17_muons_3k_ang_{ang}_e17_t0.root', ax, ay)
    d = d[d.event_id.isin(L.keys())]
    for i, v in enumerate('xy'):
        g = d[d[f'{v}_ok'].astype(bool)]
        tt = np.array([L[e][2 + i] for e in g.event_id]); p0t = np.array([L[e][i] for e in g.event_id])
        tr = g[f'{v}_tan_theta'].to_numpy(); p0 = g[f'{v}_p0'].to_numpy()
        incl = abs(np.median(tt)) > 0.05
        if incl:
            r = tr / tt; core = np.abs(r - 1) <= 0.25
            flip = float(((np.sign(tr) != np.sign(tt)) & ~core).mean()); cr = float(np.median(r[core]))
        else:
            core = np.ones(len(tr), bool); flip = np.nan; cr = np.nan
        dth = np.degrees(np.arctan(tr[core]) - np.arctan(tt[core]))
        res = {z: float(np.median(np.abs((p0[core] - p0t[core]) + z * (tr[core] - tt[core])))) for z in DEP}
        out[f'{ang}_{v}'] = dict(n=int(len(g)), flip=flip, core_ratio=cr, core_sig=rs(dth), dth_med=float(np.median(dth)),
                                 **{f'absres_z{z}': x for z, x in res.items()})
        print(f'| {ang} | {v}{" (tilted)" if incl else ""} | {len(g)} | {"" if np.isnan(flip) else f"{flip:.3f}"} | '
              f'{"" if np.isnan(cr) else f"{cr:.3f}"} | {rs(dth):.2f} | {np.median(dth):+.2f} | '
              + ' | '.join(f'{res[z]:.2f}' for z in DEP) + ' |')
json.dump(out, open(f'{C}/closure_table_{cfg}_{tag}.json', 'w'), indent=1)
