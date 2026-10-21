"""Even term on the same PCB (2026-10-21): det3 May vs det3 June-15 (one PCB, bulked twice). X sharing s1X vs the
X-independent Y share s1Y, in s1Y bins, all events (no pad labels). An attractor that makes June-15's a1 (+0.14) would
raise s1X at fixed s1Y by ~2.2 |a1| (sim) or >= 0.1 (toy). Also det1 / det2 / det4 / det6 / det7 for reference.
Usage: NB_RAW=~/x17/nb_asym_raw same_pcb.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
ED = np.array([0, .3, .45, .6, .75, .9, 1.05, 1.2, 1.4, 2.0]); OUT = {}
for b in ('det3 May', 'det3 June-15', 'det1', 'det2', 'det4', 'det6', 'det7'):
    E = R.ev(R.BUILDS[b]); sx = E['alx'][:, 5] + E['alx'][:, 3]; sy = E['ysh']; a1 = E['alx'][:, 5] - E['alx'][:, 3]
    j = np.digitize(sy, ED) - 1; row = []
    for k in range(len(ED) - 1):
        m = j == k
        row.append((float(sx[m].mean()), float(sx[m].std() / np.sqrt(m.sum())), float(a1[m].mean()), int(m.sum())) if m.sum() >= 40 else None)
    OUT[b] = row
    print(f'{b:13s} ' + ' '.join('     -      ' if r is None else f'{r[0]:.3f}±{r[1]:.3f}' for r in row) + f'   (a1 {a1.mean():+.3f})')
print('s1Y bins:', ' '.join(f'{ED[k]:.2f}-{ED[k+1]:.2f}' for k in range(len(ED) - 1)))
json.dump(dict(edges=ED.tolist(), builds=OUT), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'same_pcb.json'), 'w'), indent=1)
