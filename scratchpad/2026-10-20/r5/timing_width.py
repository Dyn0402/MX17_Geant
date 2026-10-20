"""Excess timing vs central pulse width (2026-10-21). The H4/bench cloud_basics study (nTof_x17
sps_beam_test_26/analysis/cloud_basics/FINDINGS.md §42) finds the det3 -1/+1 difference = 0.13 x centre "delayed 114 ns"
for HEAD-ON tracks (first-order fit D = a W0 + b W0', delay = -b/a); delta_timing.py over all angles gives -14 ns with
the same estimator. No track angles in the raw survey, so the central pulse width (samples >= half max) stands in:
head-on tracks give short pulses. Per width quartile: a, delay, and the excess peak vs the central, data and sim.
Usage: timing_width.py [<sim prefix> <dir> ...]  (with arms: sim only, written to timing_width_<prefixes>.json)"""
import os, sys, glob, re, json, numpy as np
import r5lib as R
import timing as T
L = R.L; HW = L.HW; Z = 8; NT = 40
def est(c, e):
    t = int(np.argmax(c)); sl = slice(t - 6, t + 14); W0 = c[sl]; D = e[sl]
    A = np.vstack([W0, np.gradient(W0)]).T; (a, b), *_ = np.linalg.lstsq(A, D, rcond=None)
    return float(a), float(-b / a * 60) if a else float('nan')
def split(W, sgn, nq=4):
    c = W[:, HW]; A = c.max(1); wid = (c >= .5 * A[:, None]).sum(1) + 1e-3 * np.random.default_rng(0).random(len(c))
    q = np.quantile(wid, np.linspace(0, 1, nq + 1)); out = []
    for k in range(nq):
        m = (wid >= q[k]) & (wid <= q[k + 1]); pk = c[m].argmax(1)
        X = T._align(W[m] / A[m, None, None], pk, Z, NT)
        cc = X[:, HW].mean(0); e = (sgn[m, None] * (X[:, HW + 1] - X[:, HW - 1])).mean(0)
        a, dl = est(cc, e)
        out.append(dict(q=k, wid=float(wid[m].mean()), n=int(m.sum()), a=a, delay=dl, epk=T._pk(e) - T._pk(cc)))
    return out
P = 780.; ORG = -(512 - 1) / 2 * P
ARMS = [(sys.argv[i], sys.argv[i + 1]) for i in range(1, len(sys.argv) - 1, 2)] or [('xo0', '~/x17/nb_asym_raw/simxoff')]
OUT = {}
for pre, dd in ARMS:
    Ws, sg = [], []
    for f in sorted(glob.glob(os.path.expanduser(f'{dd}/{pre}_*_x.npz'))):
        if not re.search('_th00_|_thy', f) or np.load(f)['wf'].ndim != 3: continue
        d = L.load(f); p = d['p'] % 512; ok = (p % 64 >= 4) & (p % 64 <= 59)
        col = 64 * (p // 64) + 63 - p % 64; delta = ((ORG + col * P + 400.) % 800.) - 400.; m = ok & (np.abs(delta) >= 150)
        Ws.append(d['W'][m].astype(np.float32)); sg.append(np.sign(delta[m]))
    OUT[f'sim {pre} |delta|>=150'] = split(np.concatenate(Ws), np.concatenate(sg))
if len(sys.argv) == 1:
    for b in ('det3 June-15',):
        W = np.concatenate([L.load(r['f'])['W'].astype(np.float32) for r in R.rows() if r['build'] == b and r['view'] == 'X'])
        OUT[f'data {b} X'] = split(W, np.ones(len(W)))
for k, rows in OUT.items():
    print(k)
    for r in rows: print(f'   width q{r["q"]} ({r["wid"]:.1f} samples) n {r["n"]:6d}  a {r["a"]:+.3f}  delay {r["delay"]:+6.0f} ns  excess peak vs central {r["epk"]:+5.0f} ns')
tag = '' if len(sys.argv) == 1 else '_' + '_'.join(a for a, _ in ARMS)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'timing_width{tag}.json'), 'w'), indent=1)
