"""Validate the noise-pattern fingerprint (2026-10-21): mean correlation between sub-runs of the same board/view
in different sessions and on different FEU slots, vs between different boards. Usage: NB_RAW=... fp_validate.py"""
import os, sys, csv, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fingerprint619 import ref
lab = [(n.split()[0], n.split()[1], n.split()[2], n.split()[3]) for n, _ in ref]  # board, view, FEU, run
V = np.array([r for _, r in ref]); C = np.corrcoef(V)
same = {}; diff = []
for i, j in itertools.combinations(range(len(ref)), 2):
    if lab[i][3] == lab[j][3] and lab[i][2] == lab[j][2]: continue            # same run and FEU: trivial
    if lab[i][:2] == lab[j][:2]:
        k = (lab[i][0], lab[i][1], 'sameFEU' if lab[i][2] == lab[j][2] else 'otherFEU'); same.setdefault(k, []).append(C[i, j])
    else: diff.append(C[i, j])
for k, v in sorted(same.items()): print(f'{" ".join(k):22s} n {len(v):4d}  median r {np.median(v):.2f}  min {np.min(v):.2f}')
print(f'different board/view   n {len(diff)}  median r {np.median(diff):.2f}  99% {np.quantile(diff, .99):.2f}  max {np.max(diff):.2f}')
