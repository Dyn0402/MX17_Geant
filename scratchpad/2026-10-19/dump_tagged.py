"""Per-event tagged ±5 neighbour waveforms (position-ordered) for every M3-good det3 event with |ax|, |ay| < AMAX,
for the X/Y neighbour asymmetry and timing study (2026-10-19). Same peak-strip choice as dump_pm3_m3.py (5 sigma seed,
largest amplitude). Unlike it, keeps per event: run tag, event id, M3 position at det z and slopes, the peak channel,
its rank in position order, the 11 channel ids, their noise. --sim DIR: a Stage BC point (det3 wiring); the
'm3' columns are then (x_true, y_true, 0, 0) from the truth sidecar, in the sim's frame, and CFG_KEY only supplies the map.
Usage: dump_tagged.py CFG_KEY OUT.npz [--amax 6] [--jobs 6] [--sub SUBRUN]"""
import os, sys, glob, argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
from qa_config import get_config, setup_paths, M3_CHI2_CUT, M3_MIN_NCLUS
HW = 5; NS = 32                                         # samples kept; longer records truncated, shorter NaN-padded

def m3_table(cfg):
    setup_paths()
    from M3RefTracking import M3RefTracking
    rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
    z = cfg.DET_PLANE_Z
    x0, y0, ev = map(np.asarray, rays.get_xy_positions(z)); x1, y1, ev1 = map(np.asarray, rays.get_xy_positions(z + 100))
    assert np.array_equal(ev, ev1)
    return {int(e): (a, b, np.degrees(np.arctan((c - a) / 100)), np.degrees(np.arctan((d - b) / 100)))
            for e, a, b, c, d in zip(ev, x0, y0, x1, y1) if np.isfinite([a, b, c, d]).all()}

def one_file(args):
    path, pos, wanted = args
    from wft.io import FeuReader
    valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]
    rank = {c: i for i, c in enumerate(order)}
    r = FeuReader(path); noise = np.where(r.noise > 0, r.noise, np.inf)
    ids = set(int(e) for e in r.event_ids) & wanted
    E, PK, RK, CH, NZ, W = [], [], [], [], [], []
    for e, _t, w in r.iter_events(ids):
        amp = w.max(axis=1); over = valid & (amp / noise >= 5)
        if not over.any(): continue
        c = np.flatnonzero(over); p = int(c[np.argmax(amp[c])])
        if p not in rank: continue
        i = rank[p]; ch = np.full(2 * HW + 1, -1, np.int16); blk = np.full((2 * HW + 1, NS), np.nan, np.float32); n = min(NS, w.shape[1])
        for k in range(-HW, HW + 1):
            if 0 <= i + k < len(order): ch[k + HW] = order[i + k]; blk[k + HW, :n] = w[order[i + k], :n]
        E.append(int(e)); PK.append(p); RK.append(i); CH.append(ch); W.append(blk)
        NZ.append(np.where(ch >= 0, r.noise[np.maximum(ch, 0)], np.nan).astype(np.float32))
    return path, E, PK, RK, CH, NZ, W

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('key'); ap.add_argument('out')
    ap.add_argument('--amax', type=float, default=6.); ap.add_argument('--jobs', type=int, default=6); ap.add_argument('--sub'); ap.add_argument('--sim')
    a = ap.parse_args()
    from wft.io import strip_position_map
    cfg = get_config(a.key)
    if a.sub: cfg.SUB_RUN = a.sub
    PM = strip_position_map(cfg)
    if a.sim:
        import pandas as pd
        t = pd.read_parquet(glob.glob(f'{a.sim}/*_truth.parquet')[0])
        m3 = {int(e): (x, y, 0., 0.) for e, x, y in zip(t.event_id, t.x_true_mm, t.y_true_mm)}
    else:
        m3 = m3_table(cfg)
    keep = {e for e, v in m3.items() if abs(v[2]) < a.amax and abs(v[3]) < a.amax}
    print(f'{a.key}: {len(m3)} M3 rays, {len(keep)} with |a| < {a.amax}', flush=True)
    dec = a.sim or f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'
    out = {}
    for view, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
        files = sorted(glob.glob(f'{dec}/*_{feu:02d}.root'))
        acc = [[] for _ in range(6)]
        with ProcessPoolExecutor(a.jobs) as ex:
            for path, *res in ex.map(one_file, [(f, PM[feu], keep) for f in files]):
                for l, v in zip(acc, res): l.extend(v)
                print(f'  {view} {os.path.basename(path)}: {len(res[0])}', flush=True)
        E = np.array(acc[0], np.int64); M = np.array([m3[e] for e in E], np.float32).reshape(-1, 4)
        out.update({f'{view}_event': E, f'{view}_peak_ch': np.array(acc[1], np.int16), f'{view}_rank': np.array(acc[2], np.int16),
                    f'{view}_ch': np.array(acc[3], np.int16).reshape(-1, 2 * HW + 1), f'{view}_noise': np.array(acc[4], np.float32).reshape(-1, 2 * HW + 1),
                    f'{view}_wf': np.array(acc[5], np.float32), f'{view}_m3': M, f'{view}_pos': PM[feu]})
        print(view, len(E), flush=True)
    np.savez_compressed(a.out, **out)
