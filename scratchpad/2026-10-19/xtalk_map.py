"""Crosstalk map over the whole FEU (2026-10-19). For vertical M3 events, alpha_ch = <w_ch, w_p> / <w_p, w_p> of
every channel against the peak channel p (in-phase content), averaged by channel-number offset (ch - p, same chip
only) and by position-rank offset. Big off-cluster entries point at routing / connector coupling.
Usage: xtalk_map.py CFG_KEY [--sub S] [--max-ev N]   writes xtalk_map_<key>.npz in cwd and prints the largest entries."""
import os, sys, glob, argparse, numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qa_config import get_config
from wft.io import FeuReader, strip_position_map
from dump_tagged import m3_table
ap = argparse.ArgumentParser(); ap.add_argument('key'); ap.add_argument('--sub'); ap.add_argument('--max-ev', type=int, default=4000); a = ap.parse_args()
cfg = get_config(a.key)
if a.sub: cfg.SUB_RUN = a.sub
PM = strip_position_map(cfg); m3 = m3_table(cfg); vert = {e for e, v in m3.items() if abs(v[2]) < 3 and abs(v[3]) < 3}
dec = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'; out = {}
for view, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
    pos = PM[feu]; valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]
    rank = np.full(512, -10 ** 6); rank[order] = np.arange(len(order))
    sch, nch = np.zeros(127), np.zeros(127); srk, nrk = np.zeros(1023), np.zeros(1023); n = 0
    for path in sorted(glob.glob(f'{dec}/*_{feu:02d}.root')):
        r = FeuReader(path); noise = np.where(r.noise > 0, r.noise, np.inf)
        for _e, _t, w in r.iter_events(set(int(e) for e in r.event_ids) & vert):
            w = w[:, :32]
            amp = w.max(1); over = valid & (amp / noise >= 5)
            if not over.any(): continue
            c = np.flatnonzero(over); p = int(c[np.argmax(amp[c])])
            if not (200 < amp[p] < 3500): continue
            al = w @ w[p] / np.dot(w[p], w[p])
            chip = np.arange(512) // 64 == p // 64; d = np.arange(512) - p
            np.add.at(sch, d[chip] + 63, al[chip]); np.add.at(nch, d[chip] + 63, 1)
            ok = valid; dr = (rank - rank[p])[ok]
            np.add.at(srk, dr + 511, al[ok]); np.add.at(nrk, dr + 511, 1); n += 1
            if n >= a.max_ev: break
        if n >= a.max_ev: break
    ach = sch / np.maximum(nch, 1); ark = srk / np.maximum(nrk, 1)
    out[f'{view}_ach'], out[f'{view}_nch'], out[f'{view}_ark'], out[f'{view}_nrk'] = ach, nch, ark, nrk
    print(f'\n{a.key} {view.upper()} (FEU {feu}), n={n}')
    print('  by channel offset (same chip), |d| <= 6: ' + ' '.join(f'{d:+d}:{ach[d+63]:+.3f}' for d in range(-6, 7)))
    print('  by position offset, |d| <= 6:            ' + ' '.join(f'{d:+d}:{ark[d+511]:+.3f}' for d in range(-6, 7)))
    far = [(abs(ach[d + 63]), d) for d in range(-63, 64) if abs(d) > 6 and nch[d + 63] > 50]
    print('  largest |alpha| at channel offsets |d| > 6: ' + ', '.join(f'{d:+d}:{ach[d+63]:+.4f}' for _, d in sorted(far)[::-1][:8]))
    farr = [(abs(ark[d + 511]), d) for d in range(-511, 512) if abs(d) > 6 and nrk[d + 511] > 50]
    print('  largest |alpha| at position offsets |d| > 6: ' + ', '.join(f'{d:+d}:{ark[d+511]:+.4f}' for _, d in sorted(farr)[::-1][:8]))
np.savez(f'xtalk_map_{a.key}.npz', **out)
