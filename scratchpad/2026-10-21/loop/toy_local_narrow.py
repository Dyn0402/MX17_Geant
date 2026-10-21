"""Per-electron pull with a narrow electron cloud and wide induction (2026-10-21). The observed sharing fixes only the
combination of cloud width and induction spread; a per-electron pull acts on the landing positions, before induction.
Scan cloud sigma 100-300 um x induction 150-400 um; keep combos whose no-pull s1 is 0.55-0.9 (data range); report the best
a1 with |ds1| <= 0.03 for per-electron pulls f 0.3-1.0 onto 50/200 um lines. Usage: toy_local_narrow.py"""
import json, os, numpy as np
from toy_loop import resp, run as run_cloud
from toy_local import run_local
OUT = []
for sig_ind in (150., 250., 350.):
    Rm = resp(sig_ind)
    for sigma in (100., 150., 200., 300.):
        _, s0 = run_cloud(sigma, 0., 50., 0., Rm, n=120)
        if not (0.5 <= s0 <= 0.95): print(f'ind {sig_ind:.0f} sigma {sigma:.0f}: no-pull s1 {s0:.2f} outside data range, skipped'); continue
        for f in (0.3, 0.6, 1.0):
            for sl in (50., 200.):
                D = np.arange(0, 391, 30.); L = [run_local(sigma, f, sl, d, Rm, n=120) for d in D]
                a = np.array([l[0] for l in L]); s = np.array([l[1] for l in L]) - s0
                ok = np.abs(s) <= 0.03; best = float(a[ok].max()) if ok.any() else float('nan')
                OUT.append(dict(sig_ind=sig_ind, sigma=sigma, s1_nopull=s0, f=f, sig_line=sl, a1=a.tolist(), ds1=s.tolist(), best=best))
                print(f'ind {sig_ind:3.0f} sigma {sigma:3.0f} (s1 {s0:.2f}) f {f:.1f} line {sl:3.0f}: a1 max {a.max():+.3f} (ds1 {s[np.argmax(a)]:+.3f}); best with |ds1|<=0.03 {best:+.3f}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'toy_local_narrow.json'), 'w'), indent=1)
