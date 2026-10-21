"""a1 vs time within and across runs (2026-10-21, audit next step 2: charging vs cloud width).
Raw survey sub-runs (X view; matched Y for the width proxy). Time order = (sub-run start, decoded file index).
Per chunk of CH events: X a1 = alpha(+1) - alpha(-1), median X central amplitude, and the matched Y share (width
proxy, X-independent). Usage: NB_RAW=~/x17/nb_asym_raw time_a1.py BOARD [CH]"""
import os, sys, csv, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'))
import asym_lib as L, asym_more as A
RAW = L.RAW; HERE = os.path.dirname(os.path.abspath(__file__))
BOARD = sys.argv[1]; CH = int(sys.argv[2]) if len(sys.argv) > 2 else 800
meta = {(r['run'], r['sub'], r['det'], r['view']): r for r in csv.DictReader(open(f'{HERE}/../../2026-10-19/copy_survey.csv'))}
idx = [l.rstrip('\n').split('\t') for l in open(f'{RAW}/survey/index.tsv')]
rows = {}
for f, rd, sub, det, v, feu in idx:
    m = meta.get((rd.split('/')[-1], sub, det, v.upper()))
    if not m or m['phys'] != BOARD: continue
    if BOARD == 'det3' and not ('6-2' in rd or 'run_79' in rd): continue
    rows.setdefault((m['tstart'], rd, sub), {})[v] = f'{RAW}/survey/{os.path.basename(f)}'
OUT = []
for (t, rd, sub), fv in sorted(rows.items()):
    if 'x' not in fv: continue
    dx = L.load(fv['x']); m = meta[(rd.split('/')[-1], sub, [k for k in meta if k[0] == rd.split('/')[-1] and k[1] == sub][0][2], 'X')]
    ysh = np.full(dx['n'], np.nan)
    if 'y' in fv:
        try:
            _, dy, ix, iy = A.matched(fv['x'], fv['y']); ysh[ix] = (dy['al'][1] + dy['al'][-1])[iy]
        except Exception: pass
    p = dx['p'] % 64; ok = (p >= 4) & (p <= 59)
    order = np.lexsort((dx['ev'], dx['fi'])); order = order[ok[order]]
    a1, A5 = dx['a1'][order], dx['A'][order, L.HW]; ys = ysh[order]
    nch = max(1, len(order) // CH); ch = []
    for k in range(nch):
        s = slice(k * len(order) // nch, (k + 1) * len(order) // nch); y = ys[s]; y = y[np.isfinite(y)]
        ch.append(dict(k=k, n=int(len(a1[s])), a1=float(a1[s].mean()), a1_e=float(a1[s].std() / np.sqrt(len(a1[s]))),
                       amp=float(np.median(A5[s])), ysh=float(np.median(y)) if len(y) > 30 else None, fi0=int(dx['fi'][order][s][0])))
    r = dict(t=t, run=rd.split('/')[-1], sub=sub, mesh=m['mesh_V'], drift=m['drift_V'], nfiles=int(dx['fi'].max() + 1), chunks=ch)
    OUT.append(r)
    print(f'{t} {r["run"][:34]:34s} {sub[:24]:24s} {m["mesh_V"]}/{m["drift_V"]:>4s} files {r["nfiles"]:3d} | ' +
          '  '.join(f'{c["a1"]:+.3f}±{c["a1_e"]:.3f}' + (f' [{c["ysh"]:.2f},{c["amp"]:.0f}]' if c['ysh'] else f' [-,{c["amp"]:.0f}]') for c in ch))
json.dump(OUT, open(f'{HERE}/time_a1_{BOARD}.json', 'w'), indent=1)
