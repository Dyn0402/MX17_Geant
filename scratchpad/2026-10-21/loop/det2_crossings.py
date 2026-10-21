"""Do det2's two zero crossings differ in sharing? (2026-10-21). An attractor makes a1(delta) zero at delta = 0 (cluster
narrowest) and at delta = P/2 (widest); det2's field goes - -> + -> - over its y bands, so its rising and falling crossings
should be one of each and differ in s1X by ~0.2-0.3 at the same Y share (toy_loop.py). Per y band (16 bands of 32 true rows),
all events: a1, s1X, s1Y, and s1X - (c0 + c1 s1Y) with c0, c1 from det3 May (the same-PCB no-pull X-vs-Y relation) and
from det2's own fit. Usage: NB_RAW=... det2_crossings.py [NB]"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 16
def tc(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
E = R.ev(R.BUILDS['det2']); al = E['alx']; a1 = al[:, 5] - al[:, 3]; sx = al[:, 5] + al[:, 3]; sy = E['ysh']
cy = tc(E['py']); b = cy * NB // 512
# det2's own s1X(s1Y) relation over the whole board (the field averages to ~0)
c1, c0 = np.polyfit(sy, sx, 1); OUT = []
print(f'det2 global s1X = {c0:.3f} + {c1:.3f} s1Y')
for k in range(NB):
    m = b == k
    if m.sum() < 60: continue
    r = dict(band=k, rows=[k * 512 // NB, (k + 1) * 512 // NB], n=int(m.sum()), a1=float(a1[m].mean()), a1_e=float(a1[m].std() / np.sqrt(m.sum())),
             s1X=float(sx[m].mean()), s1Y=float(sy[m].mean()), res=float((sx[m] - (c0 + c1 * sy[m])).mean()), res_e=float(sx[m].std() / np.sqrt(m.sum())))
    OUT.append(r)
    print(f'band {k:2d} rows {r["rows"][0]:3d}-{r["rows"][1]:3d} n {r["n"]:4d}  a1 {r["a1"]:+.3f}±{r["a1_e"]:.3f}  s1X {r["s1X"]:.3f}  s1Y {r["s1Y"]:.3f}  '
          f'resid {r["res"]:+.3f}±{r["res_e"]:.3f}  ' + '#' * int(max(0, 20 + 100 * r['res'])))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'det2_crossings_{NB}.json'), 'w'), indent=1)
