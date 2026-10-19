"""Two-copy forward fit (2026-10-20): copy a_Y applied to events whose charge is over a Y pad, a_X over an X pad.
Templates: copy-free boards (det2, det3 May, det1), where the parity split is clean (no copy to migrate peaks). Inject
n+1 += a * n per event by its own (clean) class, then run det3's analysis (raw peak, parity, Y-share bins) and compare
with det3 June-15: per Y-share bin and class, mean asymmetry and class population (migration). chi2 over a grid.
Template events are reweighted to det3 June-15's Y-share distribution. Usage: fit_two_copy.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uncopy_parity as UP
ED = np.array([0, .4, .55, .7, .85, 1.0, 1.3, 9])
def analyse(W, px0, py, ysh, wts):
    A_ = W.max(2); sh = A_[:, 3:6].argmax(1) - 1; c0 = 4 + sh; ii = np.arange(len(W))
    c = W[ii, c0]; n2 = np.einsum('ij,ij->i', c, c)
    d = (np.einsum('ij,ij->i', W[ii, c0 + 1], c) - np.einsum('ij,ij->i', W[ii, c0 - 1], c)) / n2
    px = px0 + sh; par = (px + py) % 2; ok = (px % 64 >= 4) & (px % 64 <= 59)
    ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean()); yp = par == ypar
    b = np.digitize(ysh, ED) - 1; r = []
    for k in range(len(ED) - 1):
        for cl in (True, False):
            m = ok & (b == k) & (yp == cl); w = wts[m]; sw = w.sum()
            mu = (w * d[m]).sum() / sw if sw > 0 else np.nan
            se = np.sqrt((w * (d[m] - mu) ** 2).sum() / sw / max(m.sum(), 1)) if sw > 0 else np.nan
            r.append((mu, se))
        mb = ok & (b == k); fy = (wts[mb & yp]).sum() / wts[mb].sum(); r.append((fy, np.sqrt(fy * (1 - fy) / max(mb.sum(), 1))))
    return np.array(r), yp
Wd, pxd, pyd, yshd = UP.events('det3 June-15')
data, _ = analyse(Wd.astype(float), pxd, pyd, yshd, np.ones(len(Wd)))
hd = np.histogram(yshd, ED)[0] / len(yshd)
OUT = dict(data=data.tolist(), fits={})
for tb in ('det2', 'det3 May', 'det1'):
    W, px, py, ysh = UP.events(tb); W = W.astype(float)
    ht = np.histogram(ysh, ED)[0] / len(ysh); wts = (hd / np.maximum(ht, 1e-9))[np.clip(np.digitize(ysh, ED) - 1, 0, len(ED) - 2)]
    _, yp0 = analyse(W, px, py, ysh, wts)          # clean class of each template event (no copy present)
    best = []; grid = {}
    for aY in np.arange(0, 1.01, .05):
        for aX in np.arange(-.1, .41, .05):
            a = np.where(yp0, aY, aX)[:, None, None]; Wi = W.copy(); Wi[:, 1:] += a * W[:, :-1]
            t, _ = analyse(Wi, px, py, ysh, wts)
            ok = np.isfinite(t[:, 0]) & np.isfinite(data[:, 0])
            chi = (((t[ok, 0] - data[ok, 0]) ** 2) / (t[ok, 1] ** 2 + data[ok, 1] ** 2)).sum()
            grid[(round(aY, 2), round(aX, 2))] = (chi, int(ok.sum())); best.append((chi, aY, aX, t))
    best.sort(key=lambda v: v[0]); chi, aY, aX, t = best[0]
    # also the best uniform copy (aY == aX)
    uni = min((v for v in best if abs(v[1] - v[2]) < 1e-9), key=lambda v: v[0])
    print(f'template {tb:9s}: best aY={aY:.2f} aX={aX:.2f} chi2={chi:.0f}/{grid[(round(aY,2),round(aX,2))][1]}'
          f'   | best uniform a={uni[1]:.2f} chi2={uni[0]:.0f}')
    near = sorted(((c, k) for k, (c, n) in grid.items() if c < chi + 4))
    print('   within dchi2<4:', [k for c, k in near][:12])
    lab = [f'{lo:.2f}-{hi:.2f}' for lo, hi in zip(ED[:-1], ED[1:])]
    for k, l in enumerate(lab):
        dY, dX, fy = data[3 * k], data[3 * k + 1], data[3 * k + 2]; mY, mX, my = t[3 * k], t[3 * k + 1], t[3 * k + 2]
        print(f'   Ysh {l}: Ypad data {dY[0]:+.2f}±{dY[1]:.2f} fit {mY[0]:+.2f} | Xpad data {dX[0]:+.2f}±{dX[1]:.2f} fit {mX[0]:+.2f}'
              f' | Y-class frac data {dY[0]*0+fy[0]:.2f} fit {my[0]:.2f} | uniform: {uni[3][3*k][0]:+.2f} {uni[3][3*k+1][0]:+.2f} {uni[3][3*k+2][0]:.2f}')
    OUT['fits'][tb] = dict(aY=aY, aX=aX, chi2=chi, ndf=grid[(round(aY, 2), round(aX, 2))][1] - 2, model=t.tolist(),
                           uniform=dict(a=uni[1], chi2=uni[0], model=uni[3].tolist()),
                           grid=[[k[0], k[1], v[0]] for k, v in grid.items()])
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fit_two_copy.json'), 'w'), default=float)
