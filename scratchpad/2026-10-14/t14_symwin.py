#!/usr/bin/env python3
"""T14 amplitude ledger with the angle window applied to BOTH legs (T14 applies it to data only).
Reco |theta| grows with amplitude in both legs (peak_vs_theta.py), so a data-only window biases
vertical data/sim low. Here the sim keeps only events whose reco theta lies in the same window
as the data. Peak-amplitude median data/sim per view: raw T14 | sim windowed. Runs on yoga."""
import json, os, glob, re, pandas as pd, numpy as np
SB = os.path.expanduser('~/x17/response_sim/stageB_w2')
def evfile(name):
    m = re.match(r't14_(angW|ionspeed\d+|fullphys\d+|bgtp\d+b\d+)_(th\w+)$', name)
    if not m: return None
    k, a = m.groups()
    d = f'{SB}/w2_rho2M_ang/{a}' if k == 'angW' else f'{SB}/w2_rho2M_{k}/{a}'
    f = glob.glob(f'{d}/events_*{a}.parquet'); return f[0] if f else None
print(f'{"point":26s} {"view":4s} {"win":>13s} {"n sim all/win":>14s} {"data/sim raw":>12s} {"symm":>6s} {"shift":>6s} {"rise sim raw/win, data":>24s}')
for t in sorted(glob.glob(f'{SB}/t14_*/t14_summary.json')):
    name = t.split('/')[-2]; ef = evfile(name)
    if not ef: continue
    s = json.load(open(t)); ev = pd.read_parquet(ef)
    for v in 'xy':
        lo, hi = s['theta_windows_deg'][v]
        ws = pd.read_parquet(f'{SB}/{name}/wf_sim_{v}.parquet').merge(ev[['event_id', f'{v}_theta_deg']], on='event_id')
        wd = pd.read_parquet(f'{SB}/{name}/wf_data_{v}.parquet')
        win = ws[(ws[f'{v}_theta_deg'] >= lo) & (ws[f'{v}_theta_deg'] <= hi)]
        r0 = wd.peak_amp.median() / ws.peak_amp.median(); r1 = wd.peak_amp.median() / win.peak_amp.median()
        print(f'{name[4:]:26s} {v:4s} {lo:6.1f}..{hi:5.1f} {len(ws):6d}/{len(win):<6d} {r0:12.3f} {r1:6.3f} {r1 / r0 - 1:+6.1%} '
              f'{ws.rise_ns.median():7.0f}/{win.rise_ns.median():4.0f}, {wd.rise_ns.median():4.0f}')
