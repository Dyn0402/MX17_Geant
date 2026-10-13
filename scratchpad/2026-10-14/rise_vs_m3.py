#!/usr/bin/env python3
"""Data rise / FWHM / amplitude with the TRUE (M3) angle selection, through the T14 extractor
itself (t14_compare.extract_view), against the T14 data legs (reco window + quality cut).
Target run, fiducial box, |other projection| < 6 deg. Runs on yoga."""
import os, sys, json, numpy as np, pandas as pd
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17'), os.path.expanduser('~/PycharmProjects/nTof_x17/mx17_sim_wft')]
from t14_compare import extract_view, _leg_cfg
from wft.io import strip_position_map
SB = os.path.expanduser('~/x17/response_sim/stageB_w2'); D = f'{SB}/data_leg_satdet3'
cfg = _leg_cfg(f'{SB}/w2_rho2M_fullphys146/th00', 'sat_det3'); pm = strip_position_map(cfg)
pk = pd.read_parquet(os.path.expanduser('~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet'))
pk = pk[(pk.ref_x > -200) & (pk.ref_x < 120) & (pk.ref_y.abs() < 170)]
pk['th'] = np.where(pk.view == 'x', pk.ay, -pk.ax); pk['oth'] = np.where(pk.view == 'x', pk.ax, -pk.ay)
out = []
for v, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
    for a, t14 in ((0, 'th00'), (10, 'thx10' if v == 'x' else 'thy10'), (20, 'thx20' if v == 'x' else 'thy20')):
        sel = pk[(pk.view == v) & ((pk.th.abs() - a).abs() < 3) & (pk.oth.abs() < 6)]
        t = extract_view(D, feu, pm[feu], set(sel.event_id.astype(int)), f'm3/{v}{a}')['table']
        ref = pd.read_parquet(f'{SB}/t14_angW_{t14}/wf_data_{v}.parquet')
        f = lambda x, c: x[c].median()
        print(f'{v} {a:2d}deg  TRUE-angle n {len(t):5d}: rise {f(t,"rise_ns"):5.0f} (raw {f(t,"rise_ns_raw"):5.0f})  fwhm {f(t,"fwhm_ns"):5.0f}  peak {f(t,"peak_amp"):5.0f} sat {t.saturated.mean():.3f}'
              f'  |  T14 data leg n {len(ref):5d}: rise {f(ref,"rise_ns"):5.0f}  fwhm {f(ref,"fwhm_ns"):5.0f}  peak {f(ref,"peak_amp"):5.0f} sat {(ref.peak_amp >= 3500).mean():.3f}', flush=True)
        t['view'] = v; t['a'] = a; out.append(t)
pd.concat(out).to_parquet(os.path.expanduser('~/rise_vs_m3.parquet'))
