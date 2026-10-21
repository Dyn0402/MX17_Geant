"""Which odd operator? (2026-10-21). Same PCB (det3 May -> June-15), matched Y-share bins, all events. Measured change
of the in-phase neighbour content alpha(k), k = -3..3, against two sum-conserving / sum-raising operators applied to
May's mean profile q_k (q_0 = 1):
  antisymmetric nearest-neighbour coupling  dq_n = (eps/2)(q_{n-1} - q_{n+1})   (odd, sum kept, no even term)
  one-sided copy                            dq_n = a q_{n-1}                     (raises the +1 side and the sum)
eps / a fitted per bin to the measured d a1; then predicted d alpha(k) for k = -3..3 vs measured.
(The peak re-selection is ignored; it matters only where the operator flips the peak, i.e. small eps here.)
Usage: NB_RAW=~/x17/nb_asym_raw antisym.py"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
ED = np.array([0, .45, .6, .75, .9, 1.05, 1.2, 1.4]); K = np.arange(-4, 5)
M, J = R.ev(R.BUILDS['det3 May']), R.ev(R.BUILDS['det3 June-15'])
def prof(E, m): v = E['alx'][m]; return v.mean(0), v.std(0) / np.sqrt(m.sum())
def apply(q, op, e):
    qn = q.copy()
    for n in range(1, 8):
        qn[n] += (e / 2 * (q[n - 1] - q[n + 1])) if op == 'anti' else e * q[n - 1]
    return qn / qn[4]
OUT = []
print('per Y-share bin: measured d alpha(k) (June-15 - May), k = -3..+3; then anti / copy predictions with eps fitted to d a1')
for k in range(len(ED) - 1):
    mm = (M['ysh'] >= ED[k]) & (M['ysh'] < ED[k + 1]); mj = (J['ysh'] >= ED[k]) & (J['ysh'] < ED[k + 1])
    if mm.sum() < 100: continue
    qm, em = prof(M, mm); qj, ej = prof(J, mj); d = qj - qm; de = np.hypot(em, ej)
    da1 = d[5] - d[3]; r = dict(lo=ED[k], hi=ED[k + 1], nM=int(mm.sum()), nJ=int(mj.sum()), d=d[1:8].tolist(), de=de[1:8].tolist())
    for op in ('anti', 'copy'):
        e1 = 0.1; p1 = apply(qm, op, e1) - qm; e = e1 * da1 / (p1[5] - p1[3]); p = apply(qm, op, e) - qm
        chi = float((((d - p)[1:8] / de[1:8]) ** 2).sum()); r[op] = dict(eps=float(e), pred=p[1:8].tolist(), chi2=chi)
    OUT.append(r)
    f = lambda v: ' '.join(f'{x:+.3f}' for x in v)
    print(f'Ysh {ED[k]:.2f}-{ED[k+1]:.2f} (n {mm.sum()}/{mj.sum()})  meas [{f(d[1:8])}] ±{np.mean(de[1:8]):.3f}')
    print(f'      anti eps {r["anti"]["eps"]:.3f}  [{f(r["anti"]["pred"])}]  chi2 {r["anti"]["chi2"]:.0f}/6')
    print(f'      copy a   {r["copy"]["eps"]:.3f}  [{f(r["copy"]["pred"])}]  chi2 {r["copy"]["chi2"]:.0f}/6')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'antisym.json'), 'w'), indent=1, default=float)
