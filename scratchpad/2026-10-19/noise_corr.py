"""Electronics vs detector for the one-sided X copy (2026-10-19): adjacent-channel NOISE correlation.
An electronic coupling x_{n+1} += a x_n copies channel n's noise into n+1 -> corr(n, n+1) ~ a; a copy induced by charge
on the resistive layer cannot correlate preamp noise. Uses the first NPRE samples of every event (before any signal),
pedestal-subtracted by FeuReader, with the per-chip, per-sample median across channels removed (common mode).
Also the directional regression (cancels symmetric correlations). Reports the mean correlation of channel pairs (n, n+d) inside a chip, d = 1..4, and by position offset.
Usage: noise_corr.py CFG_KEY [--sub S] [--nev 3000]"""
import os, sys, glob, argparse, numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
from qa_config import get_config
from wft.io import FeuReader, strip_position_map
ap = argparse.ArgumentParser(); ap.add_argument('key'); ap.add_argument('--sub'); ap.add_argument('--nev', type=int, default=3000); a = ap.parse_args()
NPRE = 4
cfg = get_config(a.key)
if a.sub: cfg.SUB_RUN = a.sub
PM = strip_position_map(cfg); dec = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'
for view, feu in (('x', cfg.MX17_FEU_X), ('y', cfg.MX17_FEU_Y)):
    path = sorted(glob.glob(f'{dec}/*_{feu:02d}.root'))[0]; r = FeuReader(path)
    ids = set(int(e) for e in r.event_ids[:a.nev]); X = []
    for _e, _t, w in r.iter_events(ids):
        s = w[:, :NPRE].astype(float)
        s = s - np.median(s.reshape(8, 64, NPRE), axis=1).repeat(64, axis=0)          # per-chip common mode
        X.append(s)
    X = np.concatenate(X, axis=1)                                                     # 512 x (nev*NPRE)
    good = X.std(1) > 0; Xc = (X - X.mean(1, keepdims=True)); sd = np.where(good, Xc.std(1), np.inf)
    corr = lambda i, j: np.mean(Xc[i] * Xc[j]) / (sd[i] * sd[j])
    print(f'{a.key} {view.upper()} (FEU {feu}) {X.shape[1]} samples:')
    for d in (1, 2, 3, 4):
        cs = [corr(n, n + d) for n in range(512 - d) if (n // 64 == (n + d) // 64) and good[n] and good[n + d]]
        print(f'   channel pairs (n, n+{d}) same chip: mean corr {np.mean(cs):+.3f}  (median {np.median(cs):+.3f}, n {len(cs)})')
    # directional test: x_n = b_lo x_{n-1} + b_hi x_{n+1}; a copy x_{n} += a x_{n-1} makes b_lo - b_hi ~ a
    dl = []
    for n in range(1, 511):
        if n // 64 != (n - 1) // 64 or n // 64 != (n + 1) // 64 or not (good[n - 1] and good[n] and good[n + 1]): continue
        M = np.stack([Xc[n - 1], Xc[n + 1]], 1); (bl, bh), *_ = np.linalg.lstsq(M, Xc[n], rcond=None); dl.append(bl - bh)
    print(f'   directional b(n-1) - b(n+1) in channel number: mean {np.mean(dl):+.3f}  median {np.median(dl):+.3f}  (n {len(dl)})')
    pos = PM[feu]; valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]
    cs = [corr(order[i], order[i + 1]) for i in range(len(order) - 1) if order[i] // 64 != order[i + 1] // 64]
    print(f'   position-adjacent pairs ACROSS chips: mean corr {np.mean(cs):+.3f} (n {len(cs)})')
