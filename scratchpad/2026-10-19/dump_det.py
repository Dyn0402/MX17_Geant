"""Generic M3-free version of dump_tagged.py (2026-10-19): ±5 position-ordered neighbour waveforms around the peak strip
of one detector in any run (cosmic bench or n_TOF beam), FEUs and wiring read from the run's own run_config.json.
Isolated single clusters only: at most NOVER_MAX strips of the FEU over 5 sigma, and the ±5 strips below 5 sigma
(drops the beam flash, multi-track and noisy events). Same output keys as dump_tagged.py; '<v>_m3' is zero.
Usage: dump_det.py RUN_DIR SUBRUN DET_NAME OUT.npz [--files N] [--jobs J] [--nover-max 20]"""
import os, sys, glob, json, argparse, types
import numpy as np
from concurrent.futures import ProcessPoolExecutor
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
HW = 5; NS = 32

def one_file(args):
    path, pos, nover_max = args
    from wft.io import FeuReader
    valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]; rank = {c: i for i, c in enumerate(order)}
    r = FeuReader(path); noise = np.where(r.noise > 0, r.noise, np.inf)
    E, PK, RK, CH, NZ, W = [], [], [], [], [], []
    for e, _t, w in r.iter_events(None):
        amp = w.max(axis=1); over = valid & (amp / noise >= 5)
        if not over.any() or over.sum() > nover_max: continue
        c = np.flatnonzero(over); p = int(c[np.argmax(amp[c])])
        if p not in rank: continue
        i = rank[p]
        if i - HW < 0 or i + HW >= len(order): continue
        if over[order[i - HW]] or over[order[i + HW]]: continue
        n = min(NS, w.shape[1]); ch = order[i - HW:i + HW + 1].astype(np.int16)
        blk = np.full((2 * HW + 1, NS), np.nan, np.float32); blk[:, :n] = w[ch, :n]
        E.append(int(e)); PK.append(p); RK.append(i); CH.append(ch); W.append(blk); NZ.append(r.noise[ch].astype(np.float32))
    return path, E, PK, RK, CH, NZ, W

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('run_dir'); ap.add_argument('sub'); ap.add_argument('det'); ap.add_argument('out')
    ap.add_argument('--files', type=int, default=8); ap.add_argument('--jobs', type=int, default=6); ap.add_argument('--nover-max', type=int, default=20)
    a = ap.parse_args()
    from wft.io import strip_position_map
    rc = json.load(open(f'{a.run_dir}/run_config.json')); d = [x for x in rc['detectors'] if x['name'] == a.det][0]
    fx, fy = d['dream_feus']['x_1'][0], d['dream_feus']['y_1'][0]
    cfg = types.SimpleNamespace(run_config_path=f'{a.run_dir}/run_config.json', MAP_CSV_PATH=f'{NT}/mx17_m1_map.csv', DET_NAME=a.det,
                                MX17_FEU_X=fx, MX17_FEU_Y=fy)
    PM = strip_position_map(cfg); dec = f'{a.run_dir}/{a.sub}/decoded_root'; out = {}
    print(f'{a.run_dir} {a.sub} {a.det}: FEU x {fx} y {fy}, orientation {d.get("dream_feu_orientation", {}).get("x_1")}', flush=True)
    for view, feu in (('x', fx), ('y', fy)):
        files = [f for f in sorted(glob.glob(f'{dec}/*_{feu:02d}.root')) if '_pedestals_' not in f][:a.files]
        acc = [[] for _ in range(6)]
        with ProcessPoolExecutor(a.jobs) as ex:
            for path, *res in ex.map(one_file, [(f, PM[feu], a.nover_max) for f in files]):
                for l, v in zip(acc, res): l.extend(v)
        E = np.array(acc[0], np.int64)
        out.update({f'{view}_event': E, f'{view}_peak_ch': np.array(acc[1], np.int16), f'{view}_rank': np.array(acc[2], np.int16),
                    f'{view}_ch': np.array(acc[3], np.int16).reshape(-1, 2 * HW + 1), f'{view}_noise': np.array(acc[4], np.float32).reshape(-1, 2 * HW + 1),
                    f'{view}_wf': np.array(acc[5], np.float32).reshape(-1, 2 * HW + 1, NS), f'{view}_m3': np.zeros((len(E), 4), np.float32), f'{view}_pos': PM[feu]})
        print(f'  {view} (FEU {feu}): {len(files)} files, {len(E)} isolated clusters', flush=True)
    np.savez_compressed(a.out, **out)
