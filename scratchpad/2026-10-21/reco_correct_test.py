"""Does x_corr = x_p0 - K * dx(x, y) * scale flatten the M3 residual maps? (2026-10-21)
dx: raw-ADC centroid bias map (MX17_Geant scratchpad/2026-10-20/r5/corr_maps.json), cells of 64 true columns x 64
rows = 49.92 mm cells in the wft frame (pos_mm = 0.78 x true column, both views; chpos_lookup.json). K = the transfer
coefficient measured in reco_resid_map.py (residual = K dx, K ~ -0.77). scale = the run's session scale (det3) or 1.
Residual cell rms with the alignment terms removed (offset + rotation; the alignment has no x scale), before and after; leave-one-out K for each run.
Usage: reco_correct_test.py (nTof_x17 venv)"""
import os, json, pickle, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
A = '/media/dylan/data/x17/cosmic_bench/Analysis'
CM = json.load(open(os.path.join(HERE, '..', '2026-10-20', 'r5', 'corr_maps.json')))
RR = json.load(open(os.path.join(HERE, 'reco_resid_map.json')))
RUNS = {'det3 6-22 long': ('mx17_det2_det3_overnight_6-22-26/long_run/mx17_3', 'det3 June-15', 0.96),
        'det3 6-27 sanity': ('mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check/mx17_3', 'det3 June-15', 0.44),
        'det3 6-27 490V': ('mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V/mx17_3', 'det3 June-15', 0.44),
        'det2 6-22 long': ('mx17_det2_det3_overnight_6-22-26/long_run/mx17_2', 'det2', 1.0),
        'det2 6-22 longer': ('mx17_det2_det3_overnight_6-22-26/longer_run/mx17_2', 'det2', 1.0),
        'det4 6-24': ('mx17_det4_day_6-24-26/long_run/mx17_4', 'det4', 1.0),
        'det6 6-26': ('mx17_det6_det7_overnight_6-26-26/long_run/mx17_6', 'det6', 1.0),
        'det7 6-26': ('mx17_det6_det7_overnight_6-26-26/long_run/mx17_7', 'det7', 1.0)}
NC = 8; CELL = 49.92
k_all = {n: (RR[n]['vs_raw_map']['fx0fy0']['k'], RR[n]['vs_raw_map']['fx0fy0']['k_err']) for n in RUNS if 'vs_raw_map' in RR.get(n, {})}
def dx_at(G, x, y):
    i = np.clip((x // CELL).astype(int), 0, NC - 1); j = np.clip((y // CELL).astype(int), 0, NC - 1)
    v = G[i, j]; return np.where(np.isfinite(v), v, 0.)
def cell_rms(E, r):
    # what the production alignment absorbs: offset and rotation (x residual ~ y); it has no x-scale term
    X = np.stack([np.ones(len(r)), E.y_p0 - 200], 1); c, *_ = np.linalg.lstsq(X, r, rcond=None); r = r - X @ c
    ix = np.clip((E.x_p0 // CELL).astype(int), 0, NC - 1); iy = np.clip((E.y_p0 // CELL).astype(int), 0, NC - 1)
    M, S = [], []
    for a in range(NC):
        for b in range(NC):
            m = (ix == a) & (iy == b)
            if m.sum() >= 20: M.append(r[m].mean()); S.append(r[m].std() / np.sqrt(m.sum()))
    M, S = np.array(M), np.array(S); return float(np.sqrt(np.mean(M ** 2))), float(np.sum((M / S) ** 2)), len(M)
OUT = {}
for name, (rel, build, sc) in RUNS.items():
    d = f'{A}/{rel}/wft'; C = pickle.load(open(f'{d}/bench_cache.pkl', 'rb'))['events']; E = pd.read_parquet(f'{d}/events.parquet')
    E = E[E.x_ok & E.y_ok & E.x_quality_ok & E.y_quality_ok & (E.n_tracks <= 1) & E.event_id.isin(list(C))].copy()
    E['ref_x'] = [C[int(e)]['truth']['ref_x'] for e in E.event_id]; E['tan_x'] = [C[int(e)]['truth']['tan_x'] for e in E.event_id]
    E = E[np.abs(E.tan_x) < 0.15]; r = (E.x_p0 - E.ref_x).values; ok = np.abs(r - np.median(r)) < 2.0; E = E[ok]; r = r[ok] * 1e3
    G = np.array([[np.nan if v is None else v for v in row] for row in CM['builds'][build]['dx_um']])
    others = [(k, e) for n, (k, e) in k_all.items() if n != name]; w = np.array([1 / e ** 2 for _, e in others])
    K = float(np.sum(w * np.array([k for k, _ in others])) / w.sum())
    rc = r - K * dx_at(G, E.x_p0.values, E.y_p0.values) * sc
    b0, c0, n0 = cell_rms(E, r); b1, c1, n1 = cell_rms(E, rc)
    OUT[name] = dict(K_loo=K, scale=sc, rms_before=b0, chi2_before=c0, rms_after=b1, chi2_after=c1, ncell=n0)
    print(f'{name:18s} K(leave-one-out) {K:+.2f}  scale {sc:.2f} | cell rms {b0:5.1f} -> {b1:5.1f} um, chi2 {c0:5.0f} -> {c1:5.0f} / {n0}')
json.dump(OUT, open(os.path.join(HERE, 'reco_correct_test.json'), 'w'), indent=1)
