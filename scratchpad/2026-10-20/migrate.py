"""How often does a Y-pad-like cluster (narrow Y share, copy-free view) get its X peak on the X-pad parity? (2026-10-20)
For a constant copy a this 'migration' fraction rises with a; compare det3 June-15 to copy-free boards with injected a.
Also: for migrated events, the ratio X(old centre) / X(peak) — how far above the true centre the n+1 channel got."""
import os, sys, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uncopy_parity as UP
def mig(bu, gen=None, ycut=(0, .55)):
    W, px0, py, ysh = UP.events(bu)
    if gen: W = W.copy(); W[:, 1:] += gen * W[:, :-1].copy()
    A_ = W.max(2); sh = A_[:, 3:6].argmax(1) - 1; px = px0 + sh; par = (px + py) % 2
    ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean())
    m = (ysh >= ycut[0]) & (ysh < ycut[1]) & (px % 64 >= 4) & (px % 64 <= 59)
    wrong = par[m] != ypar
    # among wrong-parity events: is the parity restored by stepping to p-1 (the copy direction) or p+1?
    ii = np.flatnonzero(m)[wrong]; c0 = 4 + sh[ii]
    lo = A_[ii, c0 - 1]; hi = A_[ii, c0 + 1]
    return m.sum(), wrong.mean(), (lo > hi).mean() if len(ii) else np.nan
for ycut in ((0, .4), (0, .55), (.55, .85), (.85, 9)):
    print('Ysh in', ycut)
    for bu in ('det3 June-15', 'det4', 'det7'):
        n, f, d = mig(bu, None, ycut); print(f'   {bu:13s}        n={n:5d}  X-peak on the other parity {100*f:5.1f} %   (true centre at p-1: {100*d:.0f} %)')
    for bu in ('det2', 'det3 May', 'det1'):
        s = f'   {bu:13s}'
        for g in (0, .2, .3, .5, .8):
            n, f, d = mig(bu, g, ycut); s += f'  a={g}: {100*f:5.1f}%'
        print(s, f' n={n}')
