"""Wide-grid linear-mixing fit scored on the class asymmetries AND the class +-1 sums per Y-share bin (class fraction
left out: it is board-dependent at zero mixing). Reports the best mixing and its narrow X-pad bins. Usage: mix_reach2.py"""
import itertools, json, os, numpy as np, mix_fit as F
Wd, pxd, pyd, yshd = F.UP.events('det3 June-15'); ED = F.ED
data = F.analyse(Wd.astype(float), pxd, pyd, yshd, np.ones(len(Wd))); hd = np.histogram(yshd, ED)[0] / len(yshd)
idx = np.array([5 * k + j for k in range(7) for j in range(4)])
def chi(t):
    ok = np.isfinite(t[idx, 0]) & np.isfinite(data[idx, 0])
    return float((((t[idx][ok, 0] - data[idx][ok, 0]) ** 2) / (t[idx][ok, 1] ** 2 + data[idx][ok, 1] ** 2)).sum()), int(ok.sum())
OUT = {}
for tb in ('det2', 'det3 May', 'det1'):
    W, px, py, ysh = F.UP.events(tb); W = W.astype(float)
    ht = np.histogram(ysh, ED)[0] / len(ysh); wts = (hd / np.maximum(ht, 1e-9))[np.clip(np.digitize(ysh, ED) - 1, 0, len(ED) - 2)]
    res = []
    for mm1, mp1, mm2, mp2 in itertools.product(np.arange(0, 1.01, .1), np.arange(-.5, .51, .1), (-.2, -.1, 0, .1, .2), (-.2, -.1, 0, .1, .2)):
        t = F.analyse(F.mix(W, {-1: mm1, 1: mp1, -2: mm2, 2: mp2}), px, py, ysh, wts); c, n = chi(t); res.append((c, n, (mm1, mp1, mm2, mp2), t))
    res.sort(key=lambda r: r[0]); c, n, m, t = res[0]
    c0 = min(r[0] for r in res if np.allclose(r[2], 0))
    print(f'{tb:9s} zero mixing chi2 {c0:.0f}/{n} | best mixing m(-1,+1,-2,+2) = {np.round(m, 2)}: chi2 {c:.0f}/{n}'
          f' | narrow X-pad asym {t[2,0]:+.2f} {t[7,0]:+.2f} (data {data[2,0]:+.2f} {data[7,0]:+.2f})'
          f' | narrow Y-pad {t[0,0]:+.2f} {t[5,0]:+.2f} (data {data[0,0]:+.2f} {data[5,0]:+.2f})')
    OUT[tb] = dict(chi2_zero=c0, ndf=n, best=list(map(float, m)), chi2=c, model=t.tolist())
OUT['data'] = data.tolist()
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mix_reach2.json'), 'w'))
