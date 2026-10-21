"""det2 crossings with a non-parametric width control (2026-10-21). Per event, residual r = s1X - F(s1Y), F = det2's whole-board
mean s1X in fine s1Y bins (the field averages to ~0 over the board, so F is close to the no-pull relation). Per y band:
mean r, and the same restricted to a common s1Y window (0.5-0.9) so bands are compared at equal width.
Also det3 June-15 by x band (one-signed: a1 never crosses zero, so no crossing test there). Usage: NB_RAW=... det2_crossings_np.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
def tc(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
ED = np.quantile(np.r_[0, 2.5], [0, 1]); OUT = {}
for b, axis, NB in (('det2', 'y', 16), ('det3 June-15', 'x', 8), ('det4', 'x', 4), ('det6', 'x', 4), ('det3 May', 'x', 8), ('det1', 'x', 4), ('det2', 'x', 8), ('det3 June-15', 'y', 8)):
    E = R.ev(R.BUILDS[b]); al = E['alx']; a1 = al[:, 5] - al[:, 3]; sx = al[:, 5] + al[:, 3]; sy = E['ysh']
    edges = np.quantile(sy, np.linspace(0, 1, 26)); j = np.clip(np.digitize(sy, edges) - 1, 0, 24)
    F = np.array([sx[j == k].mean() for k in range(25)]); r = sx - F[j]
    pos = tc(E['py'] if axis == 'y' else E['px']); band = pos * NB // 512; rows = []
    print(f'== {b} by {axis}')
    for k in range(NB):
        m = band == k; w = m & (sy > .5) & (sy < .9)
        if m.sum() < 60: continue
        row = dict(band=k, n=int(m.sum()), a1=float(a1[m].mean()), a1_e=float(a1[m].std() / np.sqrt(m.sum())), r=float(r[m].mean()),
                   r_e=float(r[m].std() / np.sqrt(m.sum())), s1Y=float(sy[m].mean()),
                   r_w=float(r[w].mean()) if w.sum() > 40 else None, r_w_e=float(r[w].std() / np.sqrt(w.sum())) if w.sum() > 40 else None, n_w=int(w.sum()))
        rows.append(row)
        rw = '      -       ' if row['r_w'] is None else f'{row["r_w"]:+.3f}±{row["r_w_e"]:.3f}'
        print(f'  band {k:2d} n {row["n"]:4d} a1 {row["a1"]:+.3f}±{row["a1_e"]:.3f}  s1Y {row["s1Y"]:.2f}  resid {row["r"]:+.3f}±{row["r_e"]:.3f}   resid(s1Y .5-.9, n {row["n_w"]:3d}) {rw}')
    OUT[f'{b} by {axis}'] = rows
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'det2_crossings_np.json'), 'w'), indent=1)
