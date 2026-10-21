"""Audit test (2026-10-21): does the X +-1 sharing follow |a1| across each board, as strip registration predicts?
Sim (sim_delta.json, xo0): s1X = alpha(+1)+alpha(-1) is EVEN in the strip offset delta: 0.63 at delta 0, 0.70-0.77 at
|delta| 75-225, 0.8-1.1 next to a gap, while a1 is odd (+-0.05). So cells with large |a1| (either sign) should share more.
Per cell (NC x NC, true columns/rows, matched X/Y events, Y-pad and X-pad classes pooled and separately):
fit s1X = c0 + c1 s1Y + k_abs |a1| + k_lin a1 (weighted). Sim reference slope from sim_delta.json.
Usage: NB_RAW=~/x17/nb_asym_raw python share_vs_a1.py [NC]"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 8
OUT = {}
S = json.load(open(os.path.join(os.path.dirname(R.__file__), 'sim_delta.json')))['col']
S = [r for r in S if r and abs(r['d']) <= 275]
sa, ss = np.array([r['a1'] for r in S]), np.array([r['s1'] for r in S])
X = np.stack([np.ones_like(sa), np.abs(sa), sa], 1); c, *_ = np.linalg.lstsq(X, ss, rcond=None)
OUT['sim'] = dict(k_abs=float(c[1]), k_lin=float(c[2]), s1_range=[float(ss.min()), float(ss.max())])
print(f'sim (|delta|<=275): s1X = {c[0]:.3f} + {c[1]:+.2f}|a1| + {c[2]:+.2f} a1   (s1 {ss.min():.2f}..{ss.max():.2f})')
for b in ('det3 June-15', 'det2', 'det4', 'det6', 'det7', 'det3 May', 'det1'):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    a1 = E['alx'][:, 5] - E['alx'][:, 3]; sx = E['alx'][:, 5] + E['alx'][:, 3]; sy = E['ysh']
    i = cx * NC // 512; j = cy * NC // 512; res = {}
    for cls, cm in (('all', np.ones(len(a1), bool)), ('ypad', E['ypad']), ('xpad', ~E['ypad'])):
        C = []
        for a in range(NC):
            for cc in range(NC):
                m = (i == a) & (j == cc) & cm
                if m.sum() >= 30:
                    C.append((a1[m].mean(), a1[m].std() / np.sqrt(m.sum()), sx[m].mean(), sx[m].std() / np.sqrt(m.sum()), sy[m].mean(), m.sum()))
        if len(C) < 6: continue
        C = np.array(C); va, ea, vs, es, vy, n = C.T
        # deconvolve noise from |a1|: use sqrt(max(a1^2 - e^2, 0)) for the abs regressor
        aa = np.sqrt(np.clip(va ** 2 - ea ** 2, 0, None))
        w = 1 / np.maximum(es, 0.005) ** 2
        Xm = np.stack([np.ones_like(va), vy, aa, va], 1); Xw = Xm * np.sqrt(w)[:, None]
        cf, *_ = np.linalg.lstsq(Xw, vs * np.sqrt(w), rcond=None); cov = np.linalg.inv(Xw.T @ Xw)
        r = vs - Xm @ cf; chi = float((w * r ** 2).sum()); s2 = max(chi / (len(vs) - 4), 1)
        ce = np.sqrt(np.diag(cov) * s2)
        rho = float(np.corrcoef(aa, vs - cf[1] * vy)[0, 1])
        res[cls] = dict(cells=int(len(va)), k_abs=float(cf[2]), k_abs_e=float(ce[2]), k_lin=float(cf[3]), k_lin_e=float(ce[3]),
                        c_sY=float(cf[1]), a1_rms=float(np.sqrt(np.mean(aa ** 2))), s1X_rms_resid=float(r.std()), chi2=chi,
                        pred_ds1=float(OUT['sim']['k_abs'] * np.sqrt(np.mean(aa ** 2))))
        print(f'{b:13s} {cls:5s} cells {len(va):2d}: k_abs {cf[2]:+.2f}±{ce[2]:.2f}  k_lin {cf[3]:+.2f}±{ce[3]:.2f}'
              f'  (|a1| rms {res[cls]["a1_rms"]:.3f} -> sim-predicted s1X swing {res[cls]["pred_ds1"]:.3f}; resid rms {r.std():.3f})')
    OUT[b] = res
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'share_vs_a1_nc{NC}.json'), 'w'), indent=1)
