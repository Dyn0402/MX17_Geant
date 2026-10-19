"""Validation of noise_corr.py's directional test (2026-10-19): read the same pre-signal noise samples, inject a
synthetic electronic copy x_n += A x_{n-1} (within a chip, sequentially as a multiplexer leak would), and recover
b(n-1) - b(n+1) (this fails: a copy correlates both ways), and the excess correlation at channel distance 1
over distances 2-3, which does respond. Usage: noise_corr_inject.py CFG_KEY [A ...]"""
import os, sys, glob, numpy as np
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis'): sys.path.insert(0, p)
from qa_config import get_config
from wft.io import FeuReader
cfg = get_config(sys.argv[1]); dec = f'{cfg.BASE_PATH}{cfg.RUN}/{cfg.SUB_RUN}/decoded_root'
r = FeuReader(sorted(glob.glob(f'{dec}/*_{cfg.MX17_FEU_X:02d}.root'))[0]); X = []
for _e, _t, w in r.iter_events(set(int(e) for e in r.event_ids[:3000])):
    s = w[:, :4].astype(float); X.append(s - np.median(s.reshape(8, 64, 4), axis=1).repeat(64, axis=0))
X0 = np.concatenate(X, axis=1)
for A in [0.] + [float(a) for a in sys.argv[2:]]:
    X = X0.copy()
    for n in range(512):
        if n % 64: X[n] = X[n] + A * X[n - 1]                      # sequential: n sees the already-leaked n-1
    Xc = X - X.mean(1, keepdims=True); good = Xc.std(1) > 0; dl = []
    for n in range(1, 511):
        if n // 64 != (n - 1) // 64 or n // 64 != (n + 1) // 64 or not good[n - 1:n + 2].all(): continue
        (bl, bh), *_ = np.linalg.lstsq(np.stack([Xc[n - 1], Xc[n + 1]], 1), Xc[n], rcond=None); dl.append(bl - bh)
    sd = np.where(good, Xc.std(1), np.inf); cr = lambda i, j: np.mean(Xc[i] * Xc[j]) / (sd[i] * sd[j])
    cd = {d: np.median([cr(n, n + d) for n in range(512 - d) if n // 64 == (n + d) // 64 and good[n] and good[n + d]]) for d in (1, 2, 3)}
    print(f'{sys.argv[1]} X injected A = {A:.2f}: b(n-1) - b(n+1) mean {np.mean(dl):+.3f};  median corr d1 {cd[1]:+.3f} d2 {cd[2]:+.3f} d3 {cd[3]:+.3f}'
          f'  -> d1 excess {cd[1] - (cd[2] + cd[3]) / 2:+.3f}')
