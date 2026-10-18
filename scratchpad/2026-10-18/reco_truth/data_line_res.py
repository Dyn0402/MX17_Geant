#!/usr/bin/env python3
"""Item 2, data side: |wft line - M3 line| at depths 3..26 mm on det3 (production events.parquet, r06 bundle).

The M3 reference geometry is attached exactly as mx_june_wft/03_angles.py does it (alignment.json, rotated tangents,
ref_mesh positions in the raw strip frame). The residual at depth z is (p0 - ref_mesh) + z (tan_fit - tan_ref).
This carries M3's own pointing error (~0.4 mm at the mesh), which the sim's truth does not.
Also returns the same numbers split by |tan_ref| < 0.06 (head-on) and 0.14-0.45 (inclined).
Usage: data_line_res.py OUT.json"""
import os, sys, json
import numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
sys.path[:0] = [NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis']
from qa_config import get_config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
setup_paths()
import cosmic_micro_tpc_analysis as cm
from M3RefTracking import M3RefTracking, get_xy_angles
from wft import compat

DEPTHS = (3.0, 9.0, 15.0, 21.0, 26.0)


def main():
    cfg = get_config('sat_det3')
    table = os.path.join(cfg.OUT_BASE, 'wft', 'events.parquet')
    params = cm.load_alignment(os.path.join(cfg.OUT_BASE, 'wft', 'alignment', 'alignment.json'))
    df = compat.load_table(table)
    results = compat.as_event_results(df)
    rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
    xa, _ya, an = get_xy_angles(rays.ray_data)
    xa = params.ref_x_sign * np.array(xa)
    cm.attach_reference_positions(results, rays, params, xa, an)
    ref = {}
    for r in results:
        if np.isnan(r.ref_tan_theta_x) or np.isnan(r.ref_mesh_x_mm):
            continue
        tx, ty = cm._rotate_ref_tangents(r, params)
        ref[int(r.event_id)] = (tx, ty, r.ref_mesh_x_mm, r.ref_mesh_y_mm)
    df = df[~df['event_id'].duplicated(keep=False)].set_index('event_id')
    out = {}
    for i, plane in enumerate('xy'):
        e = [k for k in df.index if k in ref and df.loc[k, f'{plane}_ok']]
        tr = np.array([ref[k][i] for k in e]); p0r = np.array([ref[k][2 + i] for k in e])
        tf = df.loc[e, f'{plane}_tan_theta'].to_numpy(); p0 = df.loc[e, f'{plane}_p0'].to_numpy()
        g = np.isfinite(tf) & np.isfinite(tr) & np.isfinite(p0) & np.isfinite(p0r)
        tr, p0r, tf, p0 = tr[g], p0r[g], tf[g], p0[g]
        dp = p0 - p0r
        dp = dp - np.median(dp)                    # alignment offset
        for nm, m in (('all', np.ones(len(tr), bool)), ('headon', np.abs(tr) < 0.06),
                      ('inclined', (np.abs(tr) >= 0.14) & (np.abs(tr) < 0.45))):
            o = {'n': int(m.sum()), 'dp0_rsig': float(1.4826 * np.median(np.abs(dp[m] - np.median(dp[m]))))}
            for z in DEPTHS:
                o[f'abs_res_z{z:g}'] = float(np.median(np.abs(dp[m] + z * (tf[m] - tr[m]))))
            out[f'{plane}_{nm}'] = o
            print(plane, nm, json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in o.items()}))
    json.dump(out, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
