"""Audit data prep (2026-10-21): every number the audit deck plots, recomputed from the raw-ADC survey.
-> audit_data.json (small). Usage: NB_RAW=~/x17/nb_asym_raw python audit_data.py"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import r5lib as R
OUT = {}
def mse(v):
    v = np.asarray(v); return [float(v.mean()), float(v.std() / np.sqrt(len(v))), int(len(v))] if len(v) > 20 else [None, None, int(len(v))]
# 1. survey: per sub-run board means (X and Y), from the survey rows
rows = R.rows(); S = []
for r in rows:
    S.append(dict(build=r['build'], view=r['view'], run=r['run'], sub=r['sub'], feu=r['feu'], a=r['a'], ea=r['ea'], ac=r['ac'],
                  eac=r['eac'], n=int(r['n']), coin=float(r['coin']), mesh=r['mesh'], drift=r['drift'], t=r['t']))
OUT['survey'] = S
# 2. matched-event products per build
OUT['build'] = {}; OUT['maps'] = {}; OUT['width'] = {}
NC = 8
for b in ('det3 June-15', 'det2', 'det4', 'det6', 'det7', 'det3 May', 'det1'):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    a1x = E['alx'][:, 5] - E['alx'][:, 3]; a1y = E['aly'][:, 5] - E['aly'][:, 3]; ysh = E['ysh']; yp = E['ypad']
    OUT['build'][b] = dict(n=int(len(a1x)), a1x=mse(a1x), a1y=mse(a1y), a1x_ypad=mse(a1x[yp]), a1x_xpad=mse(a1x[~yp]),
                           a1y_ypad=mse(a1y[yp]), a1y_xpad=mse(a1y[~yp]), coinx=float((E['alx'][:, 5] > E['alx'][:, 3]).mean()),
                           ysh_med=float(np.median(ysh)))
    i = cx * NC // 512; j = cy * NC // 512; cells = []
    for a in range(NC):
        for c in range(NC):
            m = (i == a) & (j == c)
            if m.sum() >= 25:
                cells.append(dict(i=a, j=c, n=int(m.sum()), ax=mse(a1x[m])[:2], ay=mse(a1y[m])[:2]))
    OUT['maps'][b] = cells
    # width curve: a1X vs Y share (all events: the pad labels migrate), per session (det3) or pooled
    ed = [0, .3, .4, .5, .6, .7, .85, 1.0, 1.3, 2.0]; W = {}
    groups = {s: E['sess'] == s for s in np.unique(E['sess'])} if b == 'det3 June-15' else {'all': np.ones(len(a1x), bool)}
    for g, gm in groups.items():
        if gm.sum() < 300: continue
        pts = []
        for lo, hi in zip(ed[:-1], ed[1:]):
            m = gm & (ysh >= lo) & (ysh < hi)
            if m.sum() >= 40: pts.append(dict(lo=lo, hi=hi, ysh=float(ysh[m].mean()), a1=mse(a1x[m])[:2], n=int(m.sum()),
                                              a1y=mse(a1y[m])[:2]))
        W[str(g)] = pts
    OUT['width'][b] = W
    print(b, OUT['build'][b])
json.dump(OUT, open(os.path.join(HERE, 'audit_data.json'), 'w'), indent=0)
