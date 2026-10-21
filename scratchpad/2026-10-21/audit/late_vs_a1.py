"""Audit (2026-10-21): does the June boards' late across-strip charge (gap conduction) follow the asymmetry field
inside a board? Per band (det2: 64-row bands; det3 June-15 / det4: 64-column chips), X view, matched events:
a1 and the X +-2 late content (mean of +-2 over 600-960 ns after the central peak, per event normalised to the
central peak, the strip_orient2 metric). Output late_vs_a1.json. Usage: NB_RAW=~/x17/nb_asym_raw python late_vs_a1.py"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import r5lib as R
OUT = {}
def late(W):
    c = W[:, 4]; z = np.argmax(c, 1); A = c[np.arange(len(c)), z]; ok = (z + 16 <= W.shape[2]) & (A > 0)
    W, z, A = W[ok], z[ok], A[ok]; idx = z[:, None] + np.arange(10, 16)[None, :]
    n2 = 0.5 * (np.take_along_axis(W[:, 2], idx, 1) + np.take_along_axis(W[:, 6], idx, 1)).mean(1) / A
    return n2, ok
for b, axis in (('det2', 'y'), ('det3 June-15', 'x'), ('det4', 'x'), ('det6', 'x'), ('det3 May', 'x')):
    E = R.ev(R.BUILDS[b]); px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    n2, ok = late(E['Wx']); a1 = (E['alx'][:, 5] - E['alx'][:, 3])[ok]; k = ((cy if axis == 'y' else cx) // 64)[ok]
    rows = []
    for q in range(8):
        m = k == q
        if m.sum() < 40: continue
        rows.append(dict(band=q, n=int(m.sum()), a1=float(a1[m].mean()), a1e=float(a1[m].std() / np.sqrt(m.sum())),
                         late=float(1e3 * n2[m].mean()), latee=float(1e3 * n2[m].std() / np.sqrt(m.sum()))))
    r = np.array([[x['a1'], x['late']] for x in rows]); rho = float(np.corrcoef(np.abs(r[:, 0]), r[:, 1])[0, 1])
    OUT[b] = dict(axis=axis, rows=rows, corr_abs_a1_late=rho, board_late=float(1e3 * n2.mean()))
    print(f'== {b} by {axis}: board late {1e3*n2.mean():.1f}e-3, corr(|a1|, late) = {rho:+.2f}')
    for x in rows: print(f"  {x['band']} n {x['n']:5d} a1 {x['a1']:+.3f}±{x['a1e']:.3f}  late {x['late']:5.1f}±{x['latee']:.1f}")
json.dump(OUT, open(os.path.join(HERE, 'late_vs_a1.json'), 'w'), indent=1)
