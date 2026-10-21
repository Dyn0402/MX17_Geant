"""Audit (2026-10-21): det3 June-15 X a1 on the board-centre chips 2-5 per session (survey sub-runs pooled,
weighted), plus the 6-18 / 6-19 EOS extractions in extra/. Merged into audit_data.json['sessions'].
Usage: NB_RAW=~/x17/nb_asym_raw python sessions.py"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import r5lib as R
M, L = R.M, R.L
D = json.load(open(os.path.join(HERE, 'audit_data.json'))); G = {}
for r in R.rows():
    if r['build'] != 'det3 June-15' or r['view'] != 'X' or not (r['nc'] > 50) or r['eac'] != r['eac']: continue
    g = G.setdefault(r['run'], dict(w=0., s=0., n=0, feu=r['feu'], drift=set(), t=r['t']))
    w = 1 / r['eac'] ** 2; g['w'] += w; g['s'] += w * r['ac']; g['n'] += r['nc']; g['drift'].add(r['drift']); g['t'] = min(g['t'] or 'z', r['t'] or 'z')
S = [dict(run=k, a=v['s'] / v['w'], e=v['w'] ** -.5, n=v['n'], feu=v['feu'], drift=sorted(x for x in v['drift'] if x), t=v['t']) for k, v in G.items()]
for f, lab, t in (('d617_long_x', '6-18 long run (EOS)', '2026-06-18'), ('d619_longer_x', '6-19/20 longer run (EOS)', '2026-06-19'),
                  ('d619_longer_feu3_x', '6-19/20, other board on FEU 3', '2026-06-19')):
    d = L.load(os.path.join(L.RAW, 'extra', f + '.npz')); a, e, n = M.central(d)
    S.append(dict(run=lab, a=a, e=e, n=n, feu=None, drift=[], t=t, extra=True))
S.sort(key=lambda s: s['t'] or 'z')
for s in S: print(f"{s['t']}  {s['run'][:60]:60s} a(chips2-5) {s['a']:+.3f}±{s['e']:.3f} n {s['n']} FEU {s['feu']} drift {s['drift']}")
D['sessions'] = S; json.dump(D, open(os.path.join(HERE, 'audit_data.json'), 'w'), indent=0)
