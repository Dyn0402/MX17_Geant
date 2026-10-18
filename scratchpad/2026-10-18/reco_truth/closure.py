#!/usr/bin/env python3
"""wft on the sim, scored against the true line (items 1 and 2 of RECO_TRUTH).

For each rewired angle point: run t13_reco.py (seeds from the sim's own waveforms, nothing from truth) with the given
bundle, then join the events table with the Stage A truth line (sim_calib.truth_lines: position at the mesh and
d position / d depth, in wft coordinates) and report, per view:
  * angle: median and robust sigma of theta_reco - theta_true [deg]; median tan_reco / tan_true (inclined views);
  * the angle constants truth would assign, fitted like production's w0/kw:  w*1e3 = w0 + kw * v * tan_true;
  * position: p0_reco - p0_true at the mesh, and |line_reco(z) - line_true(z)| at depths 3, 9, 15, 21, 26 mm
    (median, as in RECONSTRUCTION_BASIS's table);
  * efficiency: fraction of events with the view ok.
Events whose true track crosses a 64-channel connector boundary are kept (the files are rewired, so they are whole).

Usage: closure.py OUT_DIR BUNDLE ROOT CONFIG STAGEA_DIR [--jobs 8] [--tag name]
"""
import argparse, glob, json, os, subprocess, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_calib import truth_lines, ANG

NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
PY = f'{NT}/.venv/bin/python'
DEPTHS = (3.0, 9.0, 15.0, 21.0, 26.0)


def rsig(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a)))) if len(a) else np.nan


def score(ev, lines, view, v):
    t = pd.DataFrame([(e, *lines[e]) for e in ev.event_id if e in lines],
                     columns=['event_id', 'p0x', 'p0y', 'tx', 'ty']).set_index('event_id')
    d = ev.set_index('event_id').join(t, how='inner')
    ok = d[f'{view}_ok'].astype(bool) & np.isfinite(d[f'{view}_tan_theta'])
    n_all = len(d); d = d[ok]
    tt = d['tx' if view == 'x' else 'ty'].to_numpy()
    tr = d[f'{view}_tan_theta'].to_numpy()
    p0t = d['p0x' if view == 'x' else 'p0y'].to_numpy()
    p0r = d[f'{view}_p0'].to_numpy()
    dth = np.degrees(np.arctan(tr) - np.arctan(tt))
    o = {'n': int(len(d)), 'eff': len(d) / max(n_all, 1),
         'tan_true': float(np.median(tt)),
         'dtheta_med_deg': float(np.median(dth)), 'dtheta_rsig_deg': rsig(dth),
         'dtheta_s68_deg': float(np.percentile(np.abs(dth - np.median(dth)), 68)),
         'dp0_med_mm': float(np.median(p0r - p0t)), 'dp0_rsig_mm': rsig(p0r - p0t)}
    if abs(np.median(tt)) > 0.05:
        o['tan_ratio_med'] = float(np.median(tr / tt))
    for z in DEPTHS:
        r = (p0r + z * tr) - (p0t + z * tt)
        o[f'abs_res_z{z:g}'] = float(np.median(np.abs(r)))
        o[f'res_med_z{z:g}'] = float(np.median(r))
    o['_w'] = d[f'{view}_w'].to_numpy().tolist(); o['_tt'] = tt.tolist()
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('bundle'); ap.add_argument('root'); ap.add_argument('config')
    ap.add_argument('stagea')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--tag', default='reco')
    ap.add_argument('--angles', default=','.join(ANG))
    ap.add_argument('--limit', type=int, default=1200, help='events reconstructed per point')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    bj = json.load(open(os.path.join(a.bundle, 'bundle.json')))
    v = float(bj['v_drift'])
    res = {'bundle': os.path.abspath(a.bundle), 'v_drift': v, 'points': {}}
    for ang in a.angles.split(','):
        ax, ay = ANG[ang]
        d = f'{a.root}/{a.config}_{ang}'
        sa = f'{a.stagea}/mx17_muons_3k_ang_{ang}_e17_t0.root'
        if not glob.glob(f'{d}/*_07.root') or not os.path.exists(sa):
            continue
        # t13_reco wants <RUN>/<SUB_RUN>/decoded_root
        lay = f'{a.out}/layout/{a.config}_{ang}/{a.tag}'
        os.makedirs(lay, exist_ok=True)
        if not os.path.exists(f'{lay}/decoded_root'):
            os.symlink(os.path.abspath(d), f'{lay}/decoded_root')
        evp = f'{a.out}/events_{a.config}_{ang}_{a.tag}.parquet'
        if not os.path.exists(evp):
            subprocess.run([PY, f'{NT}/mx17_sim_wft/t13_reco.py', '--decoded-dir', lay, '--bundle', a.bundle,
                            '--out', evp, '--jobs', str(a.jobs), '--limit', str(a.limit)], check=True,
                           stdout=open(f'{a.out}/reco_{a.config}_{ang}_{a.tag}.log', 'w'), stderr=subprocess.STDOUT)
        ev = pd.read_parquet(evp)
        lines = truth_lines(sa, ax, ay)
        res['points'][ang] = {view: score(ev, lines, view, v) for view in 'xy'}
        p = res['points'][ang]
        print(f'{ang:6s} ' + '  '.join(
            f"{vw}: n {p[vw]['n']} dth {p[vw]['dtheta_med_deg']:+.2f}±{p[vw]['dtheta_rsig_deg']:.2f} "
            f"ratio {p[vw].get('tan_ratio_med', float('nan')):.3f} dp0 {p[vw]['dp0_med_mm']:+.2f}±{p[vw]['dp0_rsig_mm']:.2f} "
            f"|r| {p[vw]['abs_res_z3']:.2f}/{p[vw]['abs_res_z15']:.2f}/{p[vw]['abs_res_z26']:.2f}"
            for vw in 'xy'), flush=True)
    # angle constants truth would assign (production form: w*1e3 = w0 + kw * v * tan_true)
    res['angle_constants'] = {}
    for view in 'xy':
        W, T = [], []
        for ang, p in res['points'].items():
            W += p[view]['_w']; T += p[view]['_tt']
        W, T = np.array(W) * 1e3, np.array(T)
        if len(np.unique(np.round(T, 4))) >= 2:
            # robust: median per true-angle group, then a straight line through the group medians
            g = pd.DataFrame({'w': W, 't': np.round(T, 4)}).groupby('t').w.median()
            kw, w0 = np.polyfit(g.index.to_numpy() * v, g.to_numpy(), 1)
            res['angle_constants'][view] = {'kw': float(kw), 'w0_um_ns': float(w0),
                                            'groups': {str(k): float(x) for k, x in g.items()}}
            print(f'angle constants {view}: kw {kw:.4f}  w0 {w0:+.3f} um/ns  (from {len(g)} true angles)')
    for p in res['points'].values():
        for vw in 'xy':
            p[vw].pop('_w'); p[vw].pop('_tt')
    json.dump(res, open(f'{a.out}/closure_{a.config}_{a.tag}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
