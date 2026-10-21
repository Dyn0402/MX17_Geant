"""ysh_class.py on sim arms (2026-10-21): does the strip-registration model give an X asymmetry that is the same over
X and Y pads (flat in the Y-view share), as det3 June-15's is? The sim's asymmetry changes sign with the strip offset
delta, so score O = <sgn(delta) a1> (|delta| >= 50 um; paired_odd.py's definition) per Y-share bin, with the parity
Y-pad fraction fY (sim parity labels) and O split by those labels. Matched X/Y events, th00 + thy* only, edge cut.
Usage: ysh_class_sim.py DIR PREFIX [PREFIX ...]"""
import os, sys, glob, re, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '2026-10-19', 'asym_deck'))
import asym_more as A
P = 780.; ORG = -(512 - 1) / 2 * P
ED = np.array([0, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.25, 1.5, 2.0])
def ms(v): return (float(v.mean()), float(v.std() / np.sqrt(len(v))), int(len(v))) if len(v) >= 30 else None
def tc(p): p = p % 512; return 64 * (p // 64) + 63 - p % 64
D = sys.argv[1]; OUT = {}
for pre in sys.argv[2:]:
    S = dict(o=[], s=[], yp=[], d=[])
    for fx in sorted(glob.glob(f'{D}/{pre}_*_x.npz')):
        if not re.search(rf'{re.escape(pre)}_(th00|thy\d+)_s\d+_x\.npz$', fx): continue
        fy = fx[:-6] + '_y.npz'
        if not os.path.exists(fy) or np.load(fy)['wf'].ndim != 3 or np.load(fx)['wf'].ndim != 3: continue
        dx, dy, ix, iy = A.matched(fx, fy)
        if len(ix) < 20: continue
        px, py = dx['p'][ix], dy['p'][iy]
        ok = (px % 64 >= 4) & (px % 64 <= 59) & (py % 64 >= 4) & (py % 64 <= 59)
        delta = ((ORG + tc(px) * P + 400.) % 800.) - 400.
        S['o'].append((np.sign(delta) * dx['a1'][ix])[ok]); S['d'].append(delta[ok])
        S['s'].append((dy['al'][1] + dy['al'][-1])[iy][ok]); S['yp'].append(((px + py) % 2)[ok])
    S = {k: np.concatenate(v) for k, v in S.items()}
    par = S['yp']; yp = par == int(S['s'][par == 1].mean() < S['s'][par == 0].mean())
    m0 = np.abs(S['d']) >= 50; j = np.digitize(S['s'], ED) - 1; rows = []
    print(f'== {pre}  n {m0.sum()}  O(all) {S["o"][m0].mean():+.4f}')
    for k in range(len(ED) - 1):
        m = m0 & (j == k)
        r = dict(lo=float(ED[k]), hi=float(ED[k + 1]), n=int(m.sum()), fY=float(yp[m].mean()) if m.sum() >= 30 else None,
                 O=ms(S['o'][m]), O_Ylab=ms(S['o'][m & yp]), O_Xlab=ms(S['o'][m & ~yp])); rows.append(r)
        if r['O'] is None: continue
        f = lambda t: '   -          ' if t is None else f'{t[0]:+.3f}±{t[1]:.3f}'
        print(f'  Ysh {r["lo"]:.2f}-{r["hi"]:.2f} n {r["n"]:5d} fY {r["fY"]:.2f}  O(all) {f(r["O"])}  [Y-lab {f(r["O_Ylab"])} X-lab {f(r["O_Xlab"])}]')
    OUT[pre] = rows
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ysh_class_sim.json'), 'w'), indent=1)
