#!/usr/bin/env python3
"""Calibrate wft on the SIMULATION the way det3 data is calibrated, with the true line in place of M3 (item 1).

wft.calibrate is the one stage allowed to use a reference track: it caches waveform windows along the reference
corridor, measures the impulse template from bright inclined strips, then fits the kernel hypers with each event's
(p0, w) pinned to the reference. On data the reference is M3 (~0.4 mm pointing floor). Here it is the Geant4 gun:
a pencil beam in angle and the event's own vertex, so the reference is exact.

Mirrors det3's production bundle (calib_bundle_r06): the RC-ladder kernel (share_lp), c2 slaved to c1 by
c2_over_c1 = 0.6, v pinned. On data v is pinned to 36.6 (wft fit = wet Magboltz); here it is pinned to the sim's own
38.02 um/ns, the analogue. The kernel hypers start from r06 and are refitted on the sim.

Frames (all verified 2026-10-18):
  * decoded files must be REWIRED (response/dream/rewire.py), so wft's det3 map reads position = 0.78 * pad column.
  * pad column k is centred at x_local = PAD_ORIGIN + 0.78 k with PAD_ORIGIN = -199.29 mm (solver.kernels), so
    the wft position of a local x is x + 199.29 mm.
  * local z = -(z_geant - 17.644) is the depth from the mesh (ClusterFile); the gun direction (tx, ty, 1) in Geant
    is (-tx, -ty) per unit local depth. wft's tan is d(position)/d(depth), so tan_wft = -tan(gun).

Usage: sim_calib.py OUT_BUNDLE_DIR WORK_DIR REWIRED_ROOT CONFIG STAGEA_DIR [--per-point 35] [--train 180]
       (REWIRED_ROOT/<CONFIG>_<ang>/*.root for every angle found)
"""
import argparse, glob, json, os, pickle, sys
import numpy as np

NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'):
    sys.path.insert(0, p)
from qa_config import get_config, setup_paths
setup_paths()
from wft import io as wio
from wft.calib import CalibrationBundle, HYPER_NAMES
from wft import calibrate as wc

PAD_ORIGIN_MM = -199.29
Z0_GEANT = 17.64405
V_SIM = 38.02
R06 = os.path.expanduser('~/x17/cosmic_bench/Analysis/mx17_det3_saturday_scan_6-27-26/'
                         'long_run_resist_490V_drift_1000V/mx17_3/wft/calib_bundle_r06')
ANG = {'th00': (0, 0), 'thx05': (5, 0), 'thx10': (10, 0), 'thx15': (15, 0), 'thx20': (20, 0), 'thx30': (30, 0),
       'thy05': (0, 5), 'thy10': (0, 10), 'thy15': (0, 15), 'thy20': (0, 20), 'thy30': (0, 30), 'txy10': (10, 10)}


def truth_lines(stagea_file, ax_deg, ay_deg):
    """{eventID: (p0x, p0y, tan_x, tan_y)} in wft coordinates (position at the mesh, d position / d depth)."""
    import uproot
    e = uproot.open(stagea_file)['EventTree'].arrays(['eventID', 'vertexX', 'vertexY', 'vertexZ'], library='np')
    tx, ty = np.tan(np.radians(ax_deg)), np.tan(np.radians(ay_deg))
    zv = -(e['vertexZ'] - Z0_GEANT)                    # vertex depth (local)
    # local x(z) = vx + (z - zv) * (-tx)  ->  at the mesh z = 0
    xm = e['vertexX'] + zv * tx
    ym = e['vertexY'] + zv * ty
    return {int(i): (float(x - PAD_ORIGIN_MM), float(y - PAD_ORIGIN_MM), float(-tx), float(-ty))
            for i, x, y in zip(e['eventID'], xm, ym)}


def build_cache(points, per_point, pad_mm=5.0, z_lo=-3.0, z_hi=33.0, skip=0):
    """The calibrate.build_cache event dict, from sim + truth. `skip` drops the first events of each point
    (to keep the calibration events disjoint from a later evaluation sample if wanted)."""
    cfg = get_config('sat_det3')
    pos_maps = wio.strip_position_map(cfg)
    events = {}
    for ip, (d, sa, ax, ay) in enumerate(points):
        lines = truth_lines(sa, ax, ay)
        got = {}
        for plane, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
            pm = pos_maps[feu]
            for f in sorted(glob.glob(f'{d}/*_{feu:02d}.root')):
                rdr = wio.FeuReader(f)
                want = sorted(int(x) for x in rdr.event_ids)[skip:skip + 3 * per_point]
                for eid, ftst, wfm in rdr.iter_events(set(want)):
                    if eid not in lines:
                        continue
                    p0x, p0y, tnx, tny = lines[eid]
                    p0, tn = (p0x, tnx) if plane == 'x' else (p0y, tny)
                    a, b = p0 + z_lo * tn, p0 + z_hi * tn
                    lo, hi = min(a, b) - pad_mm, max(a, b) + pad_mm
                    ch = np.where((pm >= lo) & (pm <= hi))[0]
                    ch = ch[np.argsort(pm[ch])]
                    if len(ch) < 4 or wfm[ch].max() < 100:
                        continue
                    key = ip * 100000 + eid
                    ev = got.setdefault(key, dict(eid=key, tan_x=tnx, tan_y=tny, ref_mesh_x=p0x, ref_mesh_y=p0y))
                    ev[plane] = dict(ch=ch.astype(np.int16), pos=pm[ch].astype(np.float32),
                                     W=wfm[ch].astype(np.float32),
                                     noise=np.maximum(rdr.noise[ch], 3.0).astype(np.float32))
                    ev[f'ftst_{plane}'] = ftst
        both = [k for k in sorted(got) if 'x' in got[k] and 'y' in got[k]][:per_point]
        events.update({k: got[k] for k in both})
        print(f'[sim-calib] {os.path.basename(d)}: {len(both)} events', flush=True)
    return events


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('work'); ap.add_argument('root'); ap.add_argument('config')
    ap.add_argument('stagea')
    ap.add_argument('--per-point', type=int, default=35)
    ap.add_argument('--train', type=int, default=180)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--maxiter', type=int, default=130)
    ap.add_argument('--v', type=float, default=V_SIM)
    ap.add_argument('--free-v', action='store_true', help='fit v as well (data production pins it)')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)

    points = []
    for ang, (ax, ay) in ANG.items():
        d = f'{a.root}/{a.config}_{ang}'
        sa = f'{a.stagea}/mx17_muons_3k_ang_{ang}_e17_t0.root'
        if glob.glob(f'{d}/*_07.root') and os.path.exists(sa):
            points.append((d, sa, ax, ay))
    print(f'[sim-calib] {len(points)} angle points: {[os.path.basename(p[0]) for p in points]}')

    cache_path = os.path.join(a.work, 'calib_cache.pkl')
    if os.path.exists(cache_path):
        events = pickle.load(open(cache_path, 'rb'))
    else:
        events = build_cache(points, a.per_point)
        pickle.dump(events, open(cache_path, 'wb'), protocol=4)
    print(f'[sim-calib] {len(events)} cached events')

    grid, tmpl = wc.measure_templates(events)
    seed = CalibrationBundle.load(R06)
    cal = CalibrationBundle(hyper=dict(seed.hyper), v_drift=a.v, grid=grid, tmpl=tmpl,
                            gain={'x': np.ones(512), 'y': np.ones(512)}, share_mode=seed.share_mode,
                            detector='mx17_3', run_key='sat_det3_SIM_' + a.config,
                            conditions=dict(run='sim', sub_run=a.config))
    prov = os.path.join(a.work, 'provisional_bundle')
    cal.save(prov, note='sim: templates measured on the sim, hypers = r06 seed')

    # train on an even spread over the angle points (keys are point*100000 + eid)
    keys = sorted(events)
    by_pt = {}
    for k in keys:
        by_pt.setdefault(k // 100000, []).append(k)
    per = max(1, a.train // len(by_pt))
    train = sorted(k for v in by_pt.values() for k in v[:per])
    x0 = np.array([seed.hyper[k] for k in HYPER_NAMES] + [a.v])
    extra = {k: float(v) for k, v in seed.hyper.items() if k not in HYPER_NAMES}
    print(f'[sim-calib] training on {len(train)} events, extra hypers {extra}', flush=True)
    hj = wc.fit_hypers(cache_path, prov, train, jobs=a.jobs, maxiter=a.maxiter, x0=x0,
                       v_fixed=None if a.free_v else a.v, fixed={'c2': 0.0}, extra_hyper=extra)
    if not a.free_v:
        hj['v'] = a.v
    cal.hyper = {k: hj[k] for k in HYPER_NAMES}
    cal.hyper.update(extra)
    cal.v_drift = hj['v']
    cal.provenance.update(n_train=hj['n_train'], chi2=hj['chi2'], chi2_init=hj['chi2_init'],
                          fitted='MX17_Geant scratchpad/2026-10-18/reco_truth/sim_calib.py (wft.calibrate on sim, truth-pinned)',
                          seed_bundle=R06, gain_map='unit (not fitted)')
    cal.save(a.out, note='sim: hypers fitted truth-pinned')
    cal.dt_xy = wc.measure_dt_xy(events, a.out, cal.hyper, cal.v_drift)
    cal.save(a.out, note='sim: hypers fitted truth-pinned; dt_xy measured')
    json.dump(hj, open(os.path.join(a.work, 'hyper_fit.json'), 'w'), indent=1)
    print('[sim-calib]', cal.summary() if hasattr(cal, 'summary') else cal.hyper)


if __name__ == '__main__':
    main()
