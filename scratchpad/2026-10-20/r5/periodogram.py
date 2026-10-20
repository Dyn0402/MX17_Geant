"""Stripe pitch from the data (2026-10-20). Resistive stripes of pitch P over pads of pitch 780 um modulate the X
response with frequency f = frac(780/P) cycles per pad column (aliased into 0..0.5). The design value P = 800 gives
f = 0.025 (40 columns); a stripe pitch locked to the pads gives f = 0 (a constant). Per build, X view, every survey
sub-run: per-event observables against the TRUE column (col = 64*(p//64) + 63 - p%64):
  s1 = alpha(+1) + alpha(-1) (sharing), d1 = alpha(+1) - alpha(-1) (asymmetry, the copy estimator).
Amplitude spectrum |sum (o - mean) exp(2 pi i f col)| * 2/N for f = 0.002..0.5 (step 0.0005), and its noise level from
column-shuffled data (95 % of the max over f in 20 shuffles). Writes periodogram.json.
Usage: periodogram.py"""
import os, json, numpy as np
import r5lib as R
F = np.arange(0.002, 0.5001, 0.0005); OUT = {}
def spec(col, o):
    o = o - o.mean(); ph = np.exp(2j * np.pi * np.outer(F, col))
    return np.abs(ph @ o) * 2 / len(o)
rng = np.random.default_rng(3)
for b in ('det3 June-15', 'det3 May', 'det2', 'det1', 'det4', 'det6', 'det7'):
    p, s1, d1 = [], [], []
    for r in R.rows():
        if r['build'] != b or r['view'] != 'X': continue
        d = R.L.load(r['f']); ok = (d['p'] % 64 >= 4) & (d['p'] % 64 <= 59)
        p.append(d['p'][ok] % 512); s1.append((d['al'][1] + d['al'][-1])[ok]); d1.append(d['a1'][ok])
    p, s1, d1 = map(np.concatenate, (p, s1, d1)); col = 64 * (p // 64) + 63 - p % 64
    # trim outliers so a few railed/odd events don't dominate
    m = (np.abs(s1) < 3) & (np.abs(d1) < 3); col, s1, d1 = col[m], s1[m], d1[m]
    o = dict(n=int(len(col)))
    for nm, v in (('s1', s1), ('d1', d1)):
        A = spec(col, v)
        nz = np.quantile([spec(rng.permutation(col), v).max() for _ in range(20)], .95)
        top = np.argsort(A)[::-1]; peaks = []
        for i in top:
            if all(abs(F[i] - F[j]) > 0.004 for j in peaks): peaks.append(i)
            if len(peaks) == 5: break
        o[nm] = dict(mean=float(v.mean()), A=A.tolist(), noise95=float(nz), peaks=[(float(F[i]), float(A[i])) for i in peaks],
                     A_0025=float(A[np.argmin(abs(F - 0.025))]))
        print(f'{b:13s} {nm} n={len(v):6d} mean {v.mean():+.3f} | noise95 {nz:.4f} | f=0.025: {o[nm]["A_0025"]:.4f} | top: ' +
              ', '.join(f'{F[i]:.4f}({1/F[i]:.1f} col)={A[i]:.4f}' for i in peaks))
    OUT[b] = o
json.dump(dict(F=F.tolist(), boards=OUT), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'periodogram.json'), 'w'))
