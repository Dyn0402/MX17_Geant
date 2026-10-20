"""Does the production wft reconstruction carry the X asymmetry as a position bias? (2026-10-21)
Per bench run with a wft bench cache (nTof_x17 mx_june_wft/bench/build_cache.py: M3 truth in the raw strip frame) and
production events.parquet: residual r = x_p0 - ref_x for M3-matched, single-track, quality events with |tan_x| < TANCUT,
in 8 x 8 cells of (x_p0, y_p0) over the active area; the board mean is removed (alignment absorbs it). The cell rms
above the statistical floor is the position-dependent bias the reconstruction leaves in. Compare the raw-ADC centroid
bias rms of corr_maps.json (det3 June-15 77 um, det2 69, controls ~21).
Usage: reco_resid_map.py (run with the nTof_x17 venv)"""
import os, sys, json, pickle, glob, numpy as np, pandas as pd
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
OUT = {}
for name, rel in RUNS.items():
    d = f'{A}/{rel}/wft'
    if not (os.path.exists(f'{d}/bench_cache.pkl') and os.path.exists(f'{d}/events.parquet')): print(name, 'missing'); continue
    C = pickle.load(open(f'{d}/bench_cache.pkl', 'rb'))['events']
    E = pd.read_parquet(f'{d}/events.parquet')
    E = E[E.x_ok & E.y_ok & E.x_quality_ok & E.y_quality_ok & (E.n_tracks <= 1) & E.event_id.isin(list(C))]
    tr = pd.DataFrame([dict(event_id=k, ref_x=v['truth']['ref_x'] if 'truth' in v else v.get('ref_x', np.nan),
                            ref_y=v['truth']['ref_y'] if 'truth' in v else v.get('ref_y', np.nan),
                            tan_x=v['truth']['tan_x'] if 'truth' in v else v.get('tan_x', np.nan)) for k, v in C.items()])
    E = E.merge(tr, on='event_id'); E = E[np.abs(E.tan_x) < TANCUT]
    r = E.x_p0 - E.ref_x; ok = np.abs(r - np.median(r)) < 2.0; E, r = E[ok], r[ok]
    # remove the affine alignment terms (offset, x scale, rotation): r = a + b x + c y
    X = np.stack([np.ones(len(r)), E.x_p0 - 200, E.y_p0 - 200], 1); cf, *_ = np.linalg.lstsq(X, r.values, rcond=None)
    r = pd.Series((r.values - X @ cf) * 1e3, index=E.index)        # um
    xe = np.linspace(0, 399.36, NC + 1); ye = np.linspace(0, 399.36, NC + 1)   # 64-strip blocks of 0.78 mm
    ix = np.clip(np.digitize(E.x_p0, xe[1:-1]), 0, NC - 1); iy = np.clip(np.digitize(E.y_p0, ye[1:-1]), 0, NC - 1)
    M = np.full((NC, NC), np.nan); S = np.full((NC, NC), np.nan)
    for a in range(NC):
        for b in range(NC):
            m = (ix == a) & (iy == b)
            if m.sum() >= 20: M[a, b] = r[m].mean(); S[a, b] = r[m].std() / np.sqrt(m.sum())
    okc = np.isfinite(M); rms = float(np.sqrt(np.nanmean(M ** 2))); floor = float(np.sqrt(np.nanmean(S ** 2)))
    chi = float(np.nansum((M / S) ** 2))
    # 1D along x (column blocks) for the profile
    px = [float(r[ix == a].mean()) for a in range(NC)]
    OUT[name] = dict(affine=dict(scale_um_per_100mm=float(cf[1] * 1e5), rot_um_per_100mm=float(cf[2] * 1e5)), n=int(len(r)), cell_rms_um=rms, stat_floor_um=floor, excess_rms_um=float(np.sqrt(max(rms ** 2 - floor ** 2, 0))),
                     chi2=chi, ncells=int(okc.sum()), map_um=np.where(okc, M, None).tolist(), x_edges=xe.tolist(), y_edges=ye.tolist(),
                     x_profile_um=px)
    print(f'{name:18s} n {len(r):6d}  cell rms {rms:5.1f} um  stat floor {floor:5.1f}  excess {OUT[name]["excess_rms_um"]:5.1f}  '
          f'chi2 {chi:6.0f}/{int(okc.sum())}  affine b {cf[1]*1e5:+.0f} c {cf[2]*1e5:+.0f} um/100mm  x-profile ' + ' '.join(f'{v:+4.0f}' for v in px))
    # against the raw-ADC centroid bias (corr_maps.json, true column / row cells of 64), 4 orientations
    CM = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-20', 'r5', 'corr_maps.json')))
    bname = {'det3': 'det3 June-15'}.get(name.split()[0], name.split()[0])
    if bname in CM['builds']:
        G = np.array([[np.nan if v is None else v for v in row] for row in CM['builds'][bname]['dx_um']])
        res = {}
        cc_ = (np.arange(NC) + .5) - NC / 2; CX, CY = np.meshgrid(cc_, cc_, indexing='ij')
        def deaff(Z, w):
            ok = np.isfinite(Z) & np.isfinite(w); X = np.stack([np.ones(ok.sum()), CX[ok], CY[ok]], 1)
            sw = np.sqrt(w[ok]); c, *_ = np.linalg.lstsq(X * sw[:, None], Z[ok] * sw, rcond=None)
            out = np.full_like(Z, np.nan); out[ok] = Z[ok] - X @ c; return out
        Mw = 1 / S ** 2; Md = deaff(M, Mw)
        for fx in (False, True):            # fx: the true-column axis runs against +x mm (cell order AND sign flip)
            for fy in (False, True):
                Gm = (-G[::-1]) if fx else G.copy(); Gm = Gm[:, ::-1] if fy else Gm
                Gd = deaff(Gm, np.ones_like(Gm))
                ok = np.isfinite(Gd) & np.isfinite(Md) & np.isfinite(S)
                g, mm_, w = Gd[ok], Md[ok], Mw[ok]
                k = float((w * g * mm_).sum() / (w * g * g).sum()); ek = float(1 / np.sqrt((w * g * g).sum()))
                res[f'fx{int(fx)}fy{int(fy)}'] = dict(k=k, k_err=ek, corr=float(np.corrcoef(g, mm_)[0, 1]), n=int(ok.sum()))
        OUT[name]['vs_raw_map'] = res
        print(f'{"":18s} vs raw-ADC map ({bname}, affine removed from both): ' + '  '.join(f'{o}: k {v["k"]:+.2f}±{v["k_err"]:.2f} (r {v["corr"]:+.2f})' for o, v in res.items()))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reco_resid_map.json'), 'w'), indent=1)
