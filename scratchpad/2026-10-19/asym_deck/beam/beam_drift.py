"""July n_TOF drift scans (2026-10-19 evening): the X copy of every chamber vs drift voltage, from raw_asym.py --cm
dumps (condor 4502894, beamjob.sh) in ~/x17/nb_asym_raw/beam/<run>/<sub>__<det>_<view>.npz.
A = det3 (June-15), B = det2, C = det6, D = det7. Beam cut: 60 < A_c < 3500 ADC (beam gains are low), peak sample 4..22,
|mean(samples 0-2)| < 0.15 A_c, hot channels vetoed. X/Y are matched by (file index, eventId) for the parity split
(charge over an X pad vs over a Y pad, from (p_x + p_y) % 2 as on the bench).
Usage: beam_drift.py [run ...]  -> table on stdout, beam_drift.json next to this script."""
import os, re, sys, glob, json, collections, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import asym_lib as L
B = os.path.expanduser('~/x17/nb_asym_raw/beam'); HW = 4
BOARD = {'mx17_A': 'det3', 'mx17_B': 'det2', 'mx17_C': 'det6', 'mx17_D': 'det7'}

def load(f):
    z = np.load(f); W = z['wf']
    if len(W) == 0: return None
    c = W[:, HW]; A = c.max(1); pk = c.argmax(1); pre = c[:, :3].mean(1)
    g = (A > 60) & (A < 3500) & (pk >= 4) & (pk <= 22) & (np.abs(pre) < 0.15 * A)
    p = z['p'].astype(int); hv = L.hot(p[g]) if g.sum() else []; g &= ~np.isin(p, hv)
    W = W[g]; c = W[:, HW]; n2 = np.einsum('ij,ij->i', c, c)
    al = {k: np.einsum('ij,ij->i', W[:, HW + k], c) / n2 for k in (-1, 1)}
    fid = np.array([int(re.search(r'_(\d{3})_\d\d\.root$', str(x)).group(1)) for x in z['files']])
    return dict(a1=al[1] - al[-1], sym=al[1] + al[-1], p=p[g], key=fid[z['fi'][g]].astype(np.int64) * 10**9 + z['ev'][g], A=c.max(1))

def drift_of(sub):
    m = re.search(r'dr(\d+)', sub); return int(m.group(1)) if m else None

runs = sys.argv[1:] or sorted(os.listdir(B))
acc = collections.defaultdict(lambda: collections.defaultdict(list))
for run in runs:
    if not os.path.isdir(f'{B}/{run}'): continue
    subs = sorted({os.path.basename(f).split('__')[0] for f in glob.glob(f'{B}/{run}/*.npz')})
    for sub in subs:
        dr = drift_of(sub)
        for det, board in BOARD.items():
            fx, fy = f'{B}/{run}/{sub}__{det}_x.npz', f'{B}/{run}/{sub}__{det}_y.npz'
            if not os.path.exists(fx): continue
            dx = load(fx); dy = load(fy) if os.path.exists(fy) else None
            if dx is None: continue
            k = (run, board, dr); a = acc[k]
            a['a1'].append(dx['a1']); a['sym'].append(dx['sym']); a['A'].append(dx['A'])
            if dy is not None:
                a['ya1'].append(dy['a1']); a['ysym'].append(dy['sym'])
                com, ix, iy = np.intersect1d(dx['key'], dy['key'], return_indices=True)
                par = (dx['p'][ix] + dy['p'][iy]) % 2
                a['a1_Ypad'].append(dx['a1'][ix][par == 1]); a['a1_Xpad'].append(dx['a1'][ix][par == 0])
def ms(v):
    v = np.concatenate(v) if v else np.array([])
    return (float(v.mean()), float(v.std() / np.sqrt(len(v))), int(len(v))) if len(v) > 10 else (np.nan, np.nan, int(len(v)))
res = []
for (run, board, dr), a in sorted(acc.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2] or 0)):
    r = dict(run=run, board=board, drift=dr, **{k: ms(a[k]) for k in ('a1', 'sym', 'ya1', 'ysym', 'a1_Ypad', 'a1_Xpad')},
             medA=float(np.median(np.concatenate(a['A']))) if a['A'] else np.nan)
    res.append(r)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'beam_drift.json'), 'w'), indent=0)
f = lambda t: f'{t[0]:+.3f}±{t[1]:.3f}' if t[2] > 10 else '    --     '
print(f'{"run":7s} {"board":5s} drift   n    X copy         Y-pad copy    X-pad copy    Y copy        Xsym  Ysym  medA')
for r in res:
    print(f'{r["run"]:7s} {r["board"]:5s} {r["drift"]!s:>5} {r["a1"][2]:5d} {f(r["a1"])}  {f(r["a1_Ypad"])}  {f(r["a1_Xpad"])}  {f(r["ya1"])}  '
          f'{r["sym"][0]:.2f}  {r["ysym"][0]:.2f}  {r["medA"]:5.0f}')
