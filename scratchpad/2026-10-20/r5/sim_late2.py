"""Late +-2 neighbour content of sim arms with the data's metric (strip_orient2.py: mean waveforms normalised to the
central, aligned on its peak sample; +-2 mean over samples 10..15 after the peak = 600-960 ns), X view (across the
strips). Data: June builds 12-31 x1e-3, May/earlier 1-9, the v5a sim -2. (2026-10-21)
Usage: sim_late2.py <dir> <prefix> [<prefix> ...]   (SIM_VIEW=y: the Y view = along the strips; data June 63-80, May 24-37)"""
import os, sys, glob, re, json, numpy as np
import r5lib as R
L = R.L
def metric(files):
    Ws = []
    for f in files:
        d = L.load(f); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59); Ws.append(d['W'][ok])
    W = np.concatenate(Ws); c = W[:, 4]; A = c.max(1); pk = c.argmax(1); Z = 24; al = np.zeros((len(W), 9, 64))
    for i in range(len(W)):
        s = Z - pk[i]; al[i, :, s:s + 32] = W[i] / A[i]
    late = slice(Z + 10, Z + 16); n2 = lambda m: ((m[2] + m[6]) / 2)[late].mean()
    rng = np.random.default_rng(1); bs = [n2(al[rng.integers(0, len(W), len(W))].mean(0)) for _ in range(50)]
    return len(W), n2(al.mean(0)), float(np.std(bs))
d = os.path.expanduser(sys.argv[1]); JS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_late2.json')
OUT = json.load(open(JS)) if os.path.exists(JS) else {}
for pre in sys.argv[2:]:
    V = os.environ.get('SIM_VIEW', 'x')
    fs = [f for f in sorted(glob.glob(f'{d}/{pre}_*_{V}.npz')) if re.search('_th00_|_thy', f) and np.load(f)['wf'].ndim == 3]
    n, v, e = metric(fs); print(f'{pre:10s} files {len(fs):3d} events {n:6d}  {V.upper()} late +-2 {1e3*v:+6.1f} ± {1e3*e:.1f} x1e-3')
    OUT[pre if V == 'x' else f'{pre}_y'] = dict(files=len(fs), n=int(n), late2=float(v), err=float(e))
    json.dump(OUT, open(JS, 'w'), indent=1)
