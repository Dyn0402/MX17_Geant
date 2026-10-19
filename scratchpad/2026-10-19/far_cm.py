"""Common-mode test (2026-10-19): for vertical M3 events, the mean waveform (/ central peak) of strips far from the
peak strip, (a) on the same 64-channel DREAM chip as the peak strip, >= DMIN strips away in position, (b) on the
other chips of the same FEU, against the central. A chip-level ground bounce shows as an in-phase copy in (a) only.
Also reports whether FeuReader subtracts any common mode (it does not: per-channel pedestal only).
Usage: far_cm.py CFG_KEY [--sub S] [--max-ev N] [--dmin 10]   prints and writes far_cm_<key>.npz in cwd."""
import os, sys, glob, argparse, numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qa_config import get_config
from wft.io import FeuReader, strip_position_map
from dump_tagged import m3_table
ap = argparse.ArgumentParser(); ap.add_argument('key'); ap.add_argument('--sub'); ap.add_argument('--max-ev', type=int, default=3000)
ap.add_argument('--dmin', type=int, default=10); a = ap.parse_args()
cfg = get_config(a.key)
if a.sub: cfg.SUB_RUN = a.sub
PM = strip_position_map(cfg); m3 = m3_table(cfg)
vert = {e for e, v in m3.items() if abs(v[2]) < 3 and abs(v[3]) < 3}
dec = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'; out = {}
for view, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
    pos = PM[feu]; valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]; rank = np.full(512, -1); rank[order] = np.arange(len(order))
    C, S, O, n = 0, 0, 0, 0
    for path in sorted(glob.glob(f'{dec}/*_{feu:02d}.root')):
        r = FeuReader(path); noise = np.where(r.noise > 0, r.noise, np.inf)
        for _e, _t, w in r.iter_events(set(int(e) for e in r.event_ids) & vert):
            w = w[:, :32]
            if w.shape[1] < 32: continue
            amp = w.max(1); over = valid & (amp / noise >= 5)
            if not over.any(): continue
            c = np.flatnonzero(over); p = int(c[np.argmax(amp[c])]); pk = amp[p]
            if not (200 < pk < 3500): continue
            far = valid & (np.abs(rank - rank[p]) >= a.dmin) & (rank >= 0)
            same = far & (np.arange(512) // 64 == p // 64); other = far & (np.arange(512) // 64 != p // 64)
            if not same.any(): continue
            C = C + w[p] / pk; S = S + w[same].mean(0) / pk; O = O + w[other].mean(0) / pk; n += 1
            if n >= a.max_ev: break
        if n >= a.max_ev: break
    C, S, O = C / n, S / n, O / n
    k_s = np.dot(S, C) / np.dot(C, C); k_o = np.dot(O, C) / np.dot(C, C)
    print(f'{a.key} {view.upper()} (FEU {feu}) n={n}: same-chip far strips = {k_s:+.4f} x central (rms resid {np.sqrt(np.mean((S-k_s*C)**2)):.4f}); '
          f'other chips = {k_o:+.4f} x central; max |same| {np.abs(S).max():.4f}, max |other| {np.abs(O).max():.4f}')
    print('   t sample   ' + ' '.join(f'{i:6d}' for i in range(0, 32, 3)))
    for nm, y in (('central', C), ('same far', S), ('other', O)): print(f'   {nm:9s}  ' + ' '.join(f'{v:+.4f}' for v in y[::3]))
    out[f'{view}_c'], out[f'{view}_same'], out[f'{view}_other'] = C, S, O
np.savez(f'far_cm_{a.key}.npz', **out)
