"""Per-chip (true-column order) X-only asymmetry a1 for every board (2026-10-20). A copper/bus mechanism has one sign
everywhere; a stripe-registration offset that drifts across the board (pitch/scale error of the print) can change
sign. Chips in true-column order (chip = col // 64), all X sub-runs of each build pooled, edge channels excluded.
Also per 16-column block for det3 June-15. Usage: chip_profile.py"""
import json, os, numpy as np
import r5lib as R
OUT = {}
for b in ('det3 June-15', 'det3 May', 'det2', 'det1', 'det4', 'det6', 'det7'):
    p, a1 = [], []
    for r in R.rows():
        if r['build'] != b or r['view'] != 'X': continue
        d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59); p.append(d['p'][ok] % 512); a1.append(d['a1'][ok])
    p, a1 = map(np.concatenate, (p, a1)); col = 64 * (p // 64) + 63 - p % 64; ch = col // 64
    row = []
    for c in range(8):
        m = ch == c; row.append((float(a1[m].mean()), float(a1[m].std() / np.sqrt(m.sum())), int(m.sum())) if m.sum() >= 40 else (np.nan, np.nan, int(m.sum())))
    OUT[b] = row
    print(f'{b:13s} ' + '  '.join('   -        ' if np.isnan(x[0]) else f'{x[0]:+.3f}±{x[1]:.3f}' for x in row))
    if b in ('det3 June-15', 'det4', 'det7', 'det2'):
        blk = col // 16; s = []
        for k in range(32):
            m = blk == k
            s.append(f'{a1[m].mean():+.2f}' if m.sum() >= 40 else '  -  ')
        print(f'{"":13s} 16-col blocks: ' + ' '.join(s))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chip_profile.json'), 'w'))
