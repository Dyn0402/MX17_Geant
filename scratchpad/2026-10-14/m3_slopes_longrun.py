#!/usr/bin/env python3
"""Add M3 reference-track slopes to hv_slope/peaks_longrun.parquet (target run, 490 V, threshold-free
peak strip per view, M3-good rays, NO det3 reco or quality cut), and join the det3 reco angle from
data_leg_satdet3/events_data.parquet. Writes ~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet.
Then: peak amplitude vs TRUE (M3) projected angle, and the reco angle against M3. Runs on yoga."""
import os, sys, numpy as np, pandas as pd
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
from qa_config import get_config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
setup_paths()
from M3RefTracking import M3RefTracking
HS = os.path.expanduser('~/x17/response_sim/hv_slope')
cfg = get_config('sat_det3'); cfg.BASE_PATH = '/home/dylan/x17/cosmic_bench/det3/'
cfg.RUN = 'mx17_det3_saturday_scan_6-27-26'; cfg.SUB_RUN = 'long_run_resist_490V_drift_1000V'
z = 702.0
rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
x0, y0, ev = map(np.asarray, rays.get_xy_positions(z)); x1, y1, ev1 = map(np.asarray, rays.get_xy_positions(z + 100))
assert np.array_equal(ev, ev1)
m3 = pd.DataFrame(dict(event_id=ev.astype(int), tan_x=(x1 - x0) / 100, tan_y=(y1 - y0) / 100)).dropna()
pk = pd.read_parquet(f'{HS}/peaks_longrun.parquet').merge(m3, on='event_id', how='left')
reco = pd.read_parquet(os.path.expanduser('~/x17/response_sim/stageB_w2/data_leg_satdet3/events_data.parquet'))
for v in 'xy':
    pk.loc[pk.view == v, 'reco_theta'] = pk.loc[pk.view == v, 'event_id'].map(reco.set_index('event_id')[f'{v}_theta_deg'])
    pk.loc[pk.view == v, 'quality_ok'] = pk.loc[pk.view == v, 'event_id'].map(reco.set_index('event_id')[f'{v}_ok'] & reco.set_index('event_id')[f'{v}_quality_ok'])
pk['ax'] = np.degrees(np.arctan(pk.tan_x)); pk['ay'] = np.degrees(np.arctan(pk.tan_y))
pk.to_parquet(f'{HS}/peaks_longrun_m3.parquet')
print(f'{len(pk)} rows, {pk.tan_x.notna().mean():.3f} with slopes; {pk.reco_theta.notna().mean():.3f} with reco angle')
for v in 'xy':
    d = pk[(pk.view == v) & pk.reco_theta.notna() & pk.tan_x.notna()]
    print(f'view {v}: corr(reco, M3 ax) {np.corrcoef(d.reco_theta, d.ax)[0,1]:+.3f}  corr(reco, M3 ay) {np.corrcoef(d.reco_theta, d.ay)[0,1]:+.3f}')
