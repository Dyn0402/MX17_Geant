"""Per-event distribution of the n+1 vs n-1 imbalance for narrow (Y-pad-like, Y share < 0.55) clusters, peak-free:
r = (X[p0+1] - X[p0-1]) / X[p0] with p0 = the X column of Y-pad parity among the raw peak and its two neighbours that
is largest (so a migrated event is put back on the Y-pad column). det3 June-15 vs copy-free boards with an injected
uniform copy. A constant copy shifts the distribution; a variable / pad-local one also widens or skews it.
Usage: dist_narrow.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uncopy_parity as UP
def r_narrow(bu, gen=None, ycut=.55):
    W, px0, py, ysh = UP.events(bu); W = W.astype(float)
    if gen: W2 = W.copy(); W2[:, 1:] += gen * W[:, :-1]; W = W2
    A_ = W.max(2); ii = np.arange(len(W))
    par0 = (px0 + py) % 2; ypar = int(ysh[par0 == 1].mean() < ysh[par0 == 0].mean())
    # candidate centres among 3,4,5 with Y-pad parity; take the largest
    best = np.full(len(W), -1); bv = np.full(len(W), -np.inf)
    for c in (3, 4, 5):
        ok = ((px0 + c - 4 + py) % 2 == ypar) & (A_[:, c] > bv); best[ok] = c; bv[ok] = A_[ok, c]
    c0 = best; m = (ysh < ycut) & ((px0 + c0 - 4) % 64 >= 4) & ((px0 + c0 - 4) % 64 <= 59)
    r = (A_[ii, c0 + 1] - A_[ii, c0 - 1]) / A_[ii, c0]
    return r[m]
OUT = {}
q = [.1, .25, .5, .75, .9]
def show(lab, r):
    qs = np.quantile(r, q); OUT[lab] = dict(n=int(len(r)), mean=float(r.mean()), sd=float(r.std()), q=qs.tolist(),
                                            f_gt1=float((r > .5).mean()))
    print(f'{lab:22s} n={len(r):5d} mean {r.mean():+.3f} sd {r.std():.3f}  q10..90 {np.round(qs, 2)}  frac(r>0.5) {100*(r > .5).mean():4.1f}%')
show('det3 June-15', r_narrow('det3 June-15'))
show('det4', r_narrow('det4')); show('det7', r_narrow('det7'))
for bu in ('det2', 'det3 May', 'det1'):
    for g in (0, .2, .3, .4):
        show(f'{bu} +{g}', r_narrow(bu, g))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dist_narrow.json'), 'w'))
