"""Offline emulation (60 ns samples) of the non-conserving KFREE term on top of a chain point: each strip gets
kf * (its two neighbours). Same alignment and observables as nb_compare.py. Usage: emul_kfree.py NPZ"""
import sys, numpy as np
sys.argv, npz = sys.argv[:1], sys.argv[1]
exec(open(__file__.replace('emul_kfree.py', 'nb_compare.py')).read().split('tpk = lambda')[0])
tpk = lambda y: grid[np.argmax(y)] * DT
z = np.load(npz)
for v in 'xy':
    dc, dn = aligned(z[f'data_{v}_wf'])
    for kf in (0, .02, .04, .06, .08):
        W = z[f'sim_{v}_wf'].astype(float); A = W.copy()
        A[:, 1:-1] += kf * (W[:, :-2] + W[:, 2:])
        sc, sn = aligned(A)
        print(f'{v} kf {kf:.2f}: delay {tpk(sn)-tpk(sc):+4.0f} (data {tpk(dn)-tpk(dc):+.0f}) nb/c {sn.max():.3f} ({dn.max():.3f}) '
              f'nb(0) {sn[16]:.3f} ({dn[16]:.3f}) nb rms {np.sqrt(np.mean((sn-dn)**2)):.4f}')
