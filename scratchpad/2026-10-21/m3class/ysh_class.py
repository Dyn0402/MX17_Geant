"""X-independent pad class (2026-10-21). M3 points to 0.9-2 mm at the detector (m3_align.py), too coarse for 0.78 mm
pads, so label events by the Y view alone: Y share s = alpha_Y(-1) + alpha_Y(+1). Charge over a Y pad sits on one Y row
(small s); over an X pad it splits between the two neighbouring Y rows (large s). Per build and Y-share bin:
  fY  = fraction of the (p_x + p_y) parity class with the smaller mean share (the old label; clean on quiet boards),
  a1X = mean alpha_X(+1) - alpha_X(-1) over ALL events in the bin (no X-derived label at all).
If the asymmetry lived only over Y pads, a1X would track fY (fall to ~0 where fY -> 0).
Usage: NB_RAW=~/x17/nb_asym_raw ysh_class.py"""
import json, os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-20', 'r5'))
import r5lib as R
ED = np.array([0, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.25, 1.5, 2.0])
def ms(v): return (float(v.mean()), float(v.std() / np.sqrt(len(v))), int(len(v))) if len(v) >= 30 else None
def tcol(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
OUT = {}
SEL = dict(R.BUILDS)
SEL['det2 band+'] = None
for b in SEL:
    E = R.ev(R.BUILDS['det2'] if b == 'det2 band+' else R.BUILDS[b])
    a1 = E['alx'][:, 5] - E['alx'][:, 3]; s = E['ysh']; yp = E['ypad']
    m0 = np.ones(len(a1), bool)
    if b == 'det2 band+': m0 = (tcol(E['py']) >= 128) & (tcol(E['py']) < 256)
    j = np.digitize(s, ED) - 1; rows = []
    for k in range(len(ED) - 1):
        m = m0 & (j == k)
        rows.append(dict(lo=float(ED[k]), hi=float(ED[k + 1]), n=int(m.sum()), fY=float(yp[m].mean()) if m.sum() >= 30 else None,
                         a1=ms(a1[m]), a1_Ylab=ms(a1[m & yp]), a1_Xlab=ms(a1[m & ~yp])))
    OUT[b] = rows
    print(f'== {b}')
    for r in rows:
        if r['a1'] is None: continue
        f = lambda t: '   -          ' if t is None else f'{t[0]:+.3f}±{t[1]:.3f}'
        print(f'  Ysh {r["lo"]:.2f}-{r["hi"]:.2f} n {r["n"]:5d} fY {r["fY"]:.2f}  a1X(all) {f(r["a1"])}  [Y-lab {f(r["a1_Ylab"])} X-lab {f(r["a1_Xlab"])}]')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ysh_class.json'), 'w'), indent=1)
