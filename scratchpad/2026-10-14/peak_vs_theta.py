#!/usr/bin/env python3
"""Peak-strip amplitude vs RECONSTRUCTED angle, data and sim, with the T14 extractor itself.
Question: does T14's narrow angle window select a low-amplitude subset of the DATA (a reco
selection bias), which would bias every vertical data/sim amplitude ratio low?
Physics says the peak amplitude is maximal at theta = 0 and falls with |theta| (charge spreads).
Runs on yoga. Usage: peak_vs_theta.py <decoded dir> <events parquet> <tag> [max per view]"""
import sys, os, numpy as np, pandas as pd
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17'), os.path.expanduser('~/PycharmProjects/nTof_x17/mx17_sim_wft')]
from t14_compare import extract_view, _leg_cfg
from wft.io import strip_position_map
d, evp, tag = sys.argv[1], sys.argv[2], sys.argv[3]
nmax = int(sys.argv[4]) if len(sys.argv) > 4 else 12000
cfg = _leg_cfg(os.path.expanduser('~/x17/response_sim/stageB_w2/w2_rho2M_fullphys146/th00'), 'sat_det3')
pm = strip_position_map(cfg)
ev = pd.read_parquet(evp)
B = [0, .25, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15]
rows = []
for v, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
    ok = ev[f'{v}_ok'] & ev[f'{v}_quality_ok'] & ~ev['spark'] & (ev[f'{v}_theta_deg'].abs() < 15)
    sel = ev.loc[ok, ['event_id', f'{v}_theta_deg', f'{v}_slope_reliable', f'{v}_n_strips']]
    sel = sel.sort_values('event_id').head(nmax)
    r = extract_view(d, feu, pm[feu], set(sel.event_id.astype(int)), f'{tag}/{v}')
    t = r['table'].merge(sel, on='event_id')
    t['view'] = v; t['theta'] = t[f'{v}_theta_deg']; rows.append(t[['event_id', 'view', 'theta', 'peak_amp', 'q_event', 'saturated', 'rise_ns', f'{v}_slope_reliable']].rename(columns={f'{v}_slope_reliable': 'rel'}))
    print(f'\n{tag} view {v}: n {len(t)}')
    print(' |theta| bin    n   peak med  peak p25  sat frac  rise med(unsat)  reliable')
    for lo, hi in zip(B[:-1], B[1:]):
        m = (t.theta.abs() >= lo) & (t.theta.abs() < hi)
        if m.sum() < 20: continue
        x = t[m]
        print(f' {lo:5.2f}-{hi:<5.2f} {m.sum():5d} {x.peak_amp.median():9.0f} {x.peak_amp.quantile(.25):9.0f} {x.saturated.mean():9.3f} {x.rise_ns.median():12.0f} {x[f"{v}_slope_reliable"].mean():9.2f}')
pd.concat(rows).to_parquet(os.path.expanduser(f'~/peak_vs_theta_{tag}.parquet'))
