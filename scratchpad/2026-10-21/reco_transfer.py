"""Transfer of the X neighbour asymmetry into the production wft position, convention-free (2026-10-21).
From the wft bench cache (nTof_x17 mx_june_wft/bench/build_cache.py), per M3-matched event, in the reco frame:
  raw bias  b = (al[+1] - al[-1]) / (1 + al[+1] + al[-1]) x local pitch [mm], +1 = the next strip in +x (position order),
            al = in-phase content of the +-1 strips on the window's peak strip (integral maximum), production
            pre-processed waveforms. This is where the 3-strip centroid would be pulled, in mm along +x.
  reco residual r = x_p0 - ref_x (M3, raw strip frame), production events.parquet.
Per-event and per-cell (8 x 8 over x_p0, y_p0) regressions r = k b + affine(x, y): k = 1 means the reconstruction
carries the full centroid bias, 0 means the forward-model fit is immune. |tan_x| < 0.15, single track, quality.
Usage: reco_transfer.py (nTof_x17 venv)"""
import os, json, pickle, numpy as np, pandas as pd
A = '/media/dylan/data/x17/cosmic_bench/Analysis'
RUNS = {'det3 6-22 long': 'mx17_det2_det3_overnight_6-22-26/long_run/mx17_3',
        'det3 6-27 sanity': 'mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check/mx17_3',
        'det3 6-27 490V': 'mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V/mx17_3',
        'det2 6-22 long': 'mx17_det2_det3_overnight_6-22-26/long_run/mx17_2',
        'det2 6-22 longer': 'mx17_det2_det3_overnight_6-22-26/longer_run/mx17_2',
        'det4 6-24': 'mx17_det4_day_6-24-26/long_run/mx17_4',
        'det6 6-26': 'mx17_det6_det7_overnight_6-26-26/long_run/mx17_6',
        'det7 6-26': 'mx17_det6_det7_overnight_6-26-26/long_run/mx17_7'}
TANCUT = 0.15; NC = 8
def bias(win):
    W = np.asarray(win['W'], float); pos = np.asarray(win['pos'], float); o = np.argsort(pos); W, pos = W[o], pos[o]
    c = int(np.argmax(np.clip(W, 0, None).sum(1)))
    if c < 1 or c > len(W) - 2 or W[c].max() >= 3500: return np.nan
    w0 = W[c]; n2 = w0 @ w0
    if n2 <= 0: return np.nan
    ap, am = (W[c + 1] @ w0) / n2, (W[c - 1] @ w0) / n2
    return (ap - am) / (1 + ap + am) * 0.5 * (pos[c + 1] - pos[c - 1])
OUT = {}
for name, rel in RUNS.items():
    d = f'{A}/{rel}/wft'; C = pickle.load(open(f'{d}/bench_cache.pkl', 'rb'))['events']; E = pd.read_parquet(f'{d}/events.parquet')
    E = E[E.x_ok & E.y_ok & E.x_quality_ok & E.y_quality_ok & (E.n_tracks <= 1) & E.event_id.isin(list(C))].copy()
    rows = []
    for eid in E.event_id:
        v = C[int(eid)]; wx = v['wins'].get('x') or []
        if not wx: rows.append((eid, np.nan, np.nan, np.nan)); continue
        win = max(wx, key=lambda w: np.clip(np.asarray(w['W']), 0, None).sum())
        t = v['truth']; rows.append((eid, t['ref_x'], t['tan_x'], bias(win)))
    T = pd.DataFrame(rows, columns=['event_id', 'ref_x', 'tan_x', 'b']); E = E.merge(T, on='event_id')
    E = E[(np.abs(E.tan_x) < TANCUT) & np.isfinite(E.b)]
    r = E.x_p0 - E.ref_x; ok = np.abs(r - np.median(r)) < 2.0; E = E[ok]; r = r[ok].values; b = E.b.values
    X = np.stack([b, np.ones(len(r)), E.x_p0 - 200, E.y_p0 - 200], 1); c, res, *_ = np.linalg.lstsq(X, r, rcond=None)
    sig = np.std(r - X @ c); cov = np.linalg.inv(X.T @ X) * sig ** 2; k, ek = c[0], np.sqrt(cov[0, 0])
    # cell level (averages out per-event fluctuation/selection correlations)
    xe = np.linspace(0, 399.36, NC + 1); ix = np.clip(np.digitize(E.x_p0, xe[1:-1]), 0, NC - 1); iy = np.clip(np.digitize(E.y_p0, xe[1:-1]), 0, NC - 1)
    cm, cb, cw, cx, cy = [], [], [], [], []
    for a in range(NC):
        for bb in range(NC):
            m = (ix == a) & (iy == bb)
            if m.sum() >= 20: cm.append(r[m].mean()); cb.append(b[m].mean()); cw.append(m.sum() / r[m].var()); cx.append(a); cy.append(bb)
    cm, cb, cw = map(np.array, (cm, cb, cw)); Xc = np.stack([cb, np.ones(len(cm)), np.array(cx) - 3.5, np.array(cy) - 3.5], 1)
    sw = np.sqrt(cw); cc, *_ = np.linalg.lstsq(Xc * sw[:, None], cm * sw, rcond=None); covc = np.linalg.inv((Xc * cw[:, None]).T @ Xc)
    bc = cb - Xc[:, 1:] @ np.linalg.lstsq(Xc[:, 1:] * sw[:, None], cb * sw, rcond=None)[0]
    OUT[name] = dict(n=int(len(r)), k_event=float(k), k_event_err=float(ek), k_cell=float(cc[0]), k_cell_err=float(np.sqrt(covc[0, 0])),
                     b_mean_um=float(b.mean() * 1e3), b_cell_rms_um=float(np.sqrt(np.mean(bc ** 2)) * 1e3), ncell=int(len(cm)))
    o = OUT[name]
    print(f'{name:18s} n {o["n"]:5d}  raw bias mean {o["b_mean_um"]:+6.1f} um, cell rms (affine removed) {o["b_cell_rms_um"]:5.1f} um | '
          f'k per event {k:+.2f} ± {ek:.2f}   k per cell {o["k_cell"]:+.2f} ± {o["k_cell_err"]:.2f} ({o["ncell"]} cells)')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reco_transfer.json'), 'w'), indent=1)
