"""SUPERSEDED / BIASED: re-centring after an uncopy with the wrong a fakes ~+0.1 on copy-free boards; see uncopy_parity.py.
Pad-type split with the peak found AFTER undoing the copy (2026-10-20).
The raw peak is biased: when a*c makes channel n+1 the largest, (p_x + p_y) % 2 flips and a Y-pad event is filed as an
X-pad one (at narrow clusters the 'X-pad' class then shows alpha(-1) = 0.59, alpha(+1) = 0.19). Here: uncopy the X
window with a trial a, take the peak among p-1..p+1 of the uncopied window, re-centre, re-classify, and measure the
alphas on the RAW window around the new peak (copy = alpha+1 - alpha-1, raw) and on the uncopied one.
Usage: reclass.py [a_uncopy=0.2]"""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L, asym_more as A
AU = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
rows = M.survey_rows(); OUT = {}
def alph(W, c0):
    ii = np.arange(len(W)); c = W[ii, c0]; n2 = np.einsum('ij,ij->i', c, c)
    return {k: np.einsum('ij,ij->i', W[ii, c0 + k], c) / n2 for k in (-2, -1, 1, 2)}
for bu in ('det3 June-15', 'det4', 'det7', 'det3 May', 'det2', 'det1', 'det6'):
    acc = {k: [] for k in ('yp_raw', 'yp_new', 'r1', 'r_1', 'u1', 'u_1', 'ysh', 'r2', 'r_2')}
    for rx, ry in A.pairs(rows, lambda r: r['build'] == bu):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        Wx = dx['W'][ix]; U = L.uncopy(Wx, AU)
        sh = U[:, 3:6].max(2).argmax(1) - 1                     # -1, 0, +1 shift of the peak after uncopy
        c0 = 4 + sh; px = dx['p'][ix] + sh; py = dy['p'][iy]
        ysh = (dy['al'][1] + dy['al'][-1])[iy]
        par_raw = (dx['p'][ix] + py) % 2; par = (px + py) % 2
        ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean())
        ar, au = alph(Wx, c0), alph(U, c0)
        ok = (px % 64 >= 4) & (px % 64 <= 59)
        for k, v in (('yp_raw', par_raw == ypar), ('yp_new', par == ypar), ('r1', ar[1]), ('r_1', ar[-1]), ('u1', au[1]),
                     ('u_1', au[-1]), ('ysh', ysh), ('r2', ar[2]), ('r_2', ar[-2])):
            acc[k].append(v[ok])
    a = {k: np.concatenate(v) for k, v in acc.items()}
    o = dict(n=int(len(a['ysh'])), moved=float((a['yp_raw'] != a['yp_new']).mean()))
    print(f'== {bu}  n={o["n"]}  reclassified {100*o["moved"]:.1f} %  (uncopy a={AU})')
    ed = [0, .4, .55, .7, .85, 1.0, 1.3, 9]
    for cls, m0 in (('Ypad', a['yp_new']), ('Xpad', ~a['yp_new'])):
        cp = a['r1'][m0] - a['r_1'][m0]
        s = f'  {cls}: copy(raw) {cp.mean():+.3f}±{cp.std()/np.sqrt(m0.sum()):.3f}  n={m0.sum()} |'
        o[cls] = dict(copy=[float(cp.mean()), float(cp.std() / np.sqrt(m0.sum()))], bins=[])
        for lo, hi in zip(ed[:-1], ed[1:]):
            m = m0 & (a['ysh'] >= lo) & (a['ysh'] < hi)
            if m.sum() < 40: s += '     -     |'; continue
            o[cls]['bins'].append([lo, hi, int(m.sum()), float(a['r1'][m].mean()), float(a['r_1'][m].mean()), float(a['r2'][m].mean()), float(a['r_2'][m].mean())])
            s += f' {a["r1"][m].mean():.2f}/{a["r_1"][m].mean():.2f} |'
        print(s)
    OUT[bu] = o
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'reclass_a{AU:.2f}.json'), 'w'))
