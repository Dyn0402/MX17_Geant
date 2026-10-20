"""(C) Charge bookkeeping by pad class (2026-10-20). Is the n+1 excess on det3 June-15 taken from somewhere (induction
redistribution conserves the total induced charge) or added on top (a gain-like copy)?
Per event: integrals Q (sum of positive samples) of the 9-channel X and Y windows. Ratios per pad class:
 SX/SY, X central/SY, X(n+1)/SY, X(n-1)/SY, Y central/SY ... medians and means with errors, per build.
Usage: probe_charge.py"""
import numpy as np, json, os
import r5lib as R
OUT = {}
def med(x): x = x[np.isfinite(x)]; return float(np.median(x)), float(1.25 * x.std() / np.sqrt(len(x)))
print(f'{"build":13s} cls    n   | SX/SY        | Xc/SY   X+1/SY  X-1/SY  X±2/SY | Yc/SY  Y±1/SY | SY median')
for b, sel in R.BUILDS.items():
    E = R.ev(sel); o = {}
    for cls, m in (('Ypad', E['ypad']), ('Xpad', ~E['ypad'])):
        qx, qy = E['qx'][m], E['qy'][m]; SX, SY = qx.sum(1), qy.sum(1); ok = SY > 0
        r = lambda v: float(np.mean(v[ok] / SY[ok]))
        o[cls] = dict(n=int(m.sum()), SXSY=med(SX[ok] / SY[ok]), SXSY_mean=r(SX), Xc=r(qx[:, 4]), Xp1=r(qx[:, 5]), Xm1=r(qx[:, 3]),
                      Xpm2=r(qx[:, 2] + qx[:, 6]), Yc=r(qy[:, 4]), Ypm1=r(qy[:, 3] + qy[:, 5]), SY=float(np.median(SY)), SX=float(np.median(SX)))
        p = o[cls]
        print(f'{b:13s} {cls} {p["n"]:6d} | {p["SXSY"][0]:.3f}±{p["SXSY"][1]:.3f} | {p["Xc"]:.3f}  {p["Xp1"]:.3f}  {p["Xm1"]:.3f}  {p["Xpm2"]:.3f} | {p["Yc"]:.3f}  {p["Ypm1"]:.3f} | {p["SY"]:.0f} (SX {p["SX"]:.0f})')
    OUT[b] = o
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'probe_charge.json'), 'w'), indent=1)
