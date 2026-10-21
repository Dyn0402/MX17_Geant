"""Shape of the odd transfer within one board (2026-10-21). det2's field is two-signed in y bands, so per cell
(8 y-bands x 4 x-quarters, true rows/columns, all events) regress the X neighbour content alpha(k), k = -3..3, on
the cell's a1 and |a1|, with the cell's Y share and Y +-2 content as width controls:
  alpha(k) = c + b_a1 * a1 + b_abs * |a1| + c_y * s1Y + c_y2 * s2Y.
Odd transfer: b_a1(k) antisymmetric in k; even response to the mechanism: b_abs(k). Also det3 June-15 (one-signed,
so b_a1 and b_abs are degenerate there: only their sum is reported). Usage: NB_RAW=... det2_profile.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
def tc(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
OUT = {}
for b, (nx, ny) in (('det2', (4, 8)), ('det3 June-15', (8, 4)), ('det4', (4, 4)), ('det6', (4, 4))):
    E = R.ev(R.BUILDS[b]); cx, cy = tc(E['px']), tc(E['py']); al = E['alx']; aly = E['aly']
    a1 = al[:, 5] - al[:, 3]; C = []
    for i in range(nx):
        for j in range(ny):
            m = (cx * nx // 512 == i) & (cy * ny // 512 == j)
            if m.sum() < 60: continue
            C.append(np.concatenate([[a1[m].mean(), (aly[m, 3] + aly[m, 5]).mean(), (aly[m, 2] + aly[m, 6]).mean(), m.sum()],
                                     al[m][:, 1:8].mean(0), al[m][:, 1:8].std(0) / np.sqrt(m.sum())]))
    C = np.array(C); A1, SY, S2, N = C[:, 0], C[:, 1], C[:, 2], C[:, 3]; P, PE = C[:, 4:11], C[:, 11:18]
    X = np.stack([np.ones(len(C)), A1, np.abs(A1), SY, S2], 1); res = []
    for k in range(7):
        if k == 3: res.append(dict(k=0, b_a1=0., b_a1_e=0., b_abs=0., b_abs_e=0.)); continue
        w = 1 / PE[:, k] ** 2; Xw = X * np.sqrt(w)[:, None]
        cf, *_ = np.linalg.lstsq(Xw, P[:, k] * np.sqrt(w), rcond=None); cov = np.linalg.pinv(Xw.T @ Xw)
        r = P[:, k] - X @ cf; s2 = max((w * r ** 2).sum() / max(len(C) - 5, 1), 1); e = np.sqrt(np.diag(cov) * s2)
        res.append(dict(k=k - 3, b_a1=float(cf[1]), b_a1_e=float(e[1]), b_abs=float(cf[2]), b_abs_e=float(e[2])))
    OUT[b] = dict(cells=len(C), a1_range=[float(A1.min()), float(A1.max())], fit=res)
    print(f'== {b}: {len(C)} cells, a1 {A1.min():+.2f} .. {A1.max():+.2f}')
    print('   k      ' + ' '.join(f'{r["k"]:+12d}' for r in res))
    print('   b_a1   ' + ' '.join(f'{r["b_a1"]:+.3f}±{r["b_a1_e"]:.3f}' for r in res))
    print('   b_abs  ' + ' '.join(f'{r["b_abs"]:+.3f}±{r["b_abs_e"]:.3f}' for r in res))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'det2_profile.json'), 'w'), indent=1)
