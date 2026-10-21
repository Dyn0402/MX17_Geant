"""6-19/20 FEU 7 vs det2 / det3 (2026-10-21): hit fraction per 64-ch chip, mean central pulse width (w30, ns),
X +-1 sharing s1 and a1, for the 6-19 extractions and the 6-22 det2 / det3 X sub-runs. Usage: NB_RAW=... board_props.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'))
import asym_lib as L, asym_more as A
RAW = L.RAW
idx = [l.rstrip('\n').split('\t') for l in open(f'{RAW}/survey/index.tsv')]
F = {'6-19 FEU7 (labelled det3)': f'{RAW}/extra/d619_longer_x.npz', '6-19 FEU3': f'{RAW}/extra/d619_longer_feu3_x.npz',
     '6-17 FEU7 det3': f'{RAW}/extra/d617_long_x.npz'}
for f, rd, sub, det, v, feu in idx:
    if v == 'x' and '6-22' in rd and sub in ('long_run', 'longer_run', 'short_run'):
        F[f'6-22 {sub} {"det3" if det == "mx17_3" else "det2"} FEU{feu}'] = f'{RAW}/survey/{os.path.basename(f)}'
    if v == 'x' and '6-17' in rd and det == 'mx17_1': F[f'6-17 {sub} det2 FEU{feu}'] = f'{RAW}/survey/{os.path.basename(f)}'
for k, f in F.items():
    d = L.load(f); ch = (d['p'] % 512) // 64; h = np.bincount(ch, minlength=8) / len(ch)
    c = d['W'][:, L.HW]; w = A.w30(c[d['A'][:, L.HW] > 300]) if (d['A'][:, L.HW] > 300).sum() > 50 else float('nan')
    s1 = d['al'][1] + d['al'][-1]
    print(f'{k:34s} n {d["n"]:5d} | hits/chip ' + ' '.join(f'{x:.2f}' for x in h) + f' | w30 {w:5.0f} ns | s1 {s1.mean():.3f} | a1 {d["a1"].mean():+.3f} | median A {np.median(d["A"][:, L.HW]):.0f}')
