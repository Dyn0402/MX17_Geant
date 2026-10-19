"""Is the parity split a selection artefact of a plain channel copy? (2026-10-20)
Inject x'(n+1) = x(n+1) + a x(n) (non-recursive, the copy model) into the raw X windows of copy-free boards, re-find the
peak exactly as raw_asym does (largest channel), re-centre, and redo the X/Y parity split and the Y-share binning.
Compare with det3 June-15 (no injection). Usage: inject_parity.py"""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L, asym_more as A
rows = M.survey_rows(); OUT = []
ED = [0, .4, .55, .7, .85, 1.0, 1.3, 9]
def run(bu, a):
    acc = {k: [] for k in ('yp', 'c1', 'r1', 'r_1', 'ysh')}
    for rx, ry in A.pairs(rows, lambda r: r['build'] == bu):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        W = dx['W'][ix].astype(float)
        if a:
            W2 = W.copy(); W2[:, 1:] += a * W[:, :-1]; W = W2
        A_ = W.max(2); sh = A_[:, 3:6].argmax(1) - 1; c0 = 4 + sh; ii = np.arange(len(W))
        # raw_asym's peak is the largest channel over the board; within +-1 of the old peak is enough here
        c = W[ii, c0]; n2 = np.einsum('ij,ij->i', c, c)
        al = {k: np.einsum('ij,ij->i', W[ii, c0 + k], c) / n2 for k in (-1, 1)}
        px = dx['p'][ix] + sh; py = dy['p'][iy]; par = (px + py) % 2
        ysh = (dy['al'][1] + dy['al'][-1])[iy]; ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean())
        ok = (px % 64 >= 4) & (px % 64 <= 59) & (A_[ii, c0] < 3500)
        for k, v in (('yp', par == ypar), ('c1', al[1] - al[-1]), ('r1', al[1]), ('r_1', al[-1]), ('ysh', ysh)):
            acc[k].append(v[ok])
    x = {k: np.concatenate(v) for k, v in acc.items()}
    o = dict(build=bu, a=a, n=int(len(x['c1'])))
    for cls, m0 in (('ypad', x['yp']), ('xpad', ~x['yp'])):
        cp = x['c1'][m0]; o[cls] = [float(cp.mean()), float(cp.std() / np.sqrt(len(cp))), int(len(cp))]
        o[cls + '_bins'] = []
        for lo, hi in zip(ED[:-1], ED[1:]):
            m = m0 & (x['ysh'] >= lo) & (x['ysh'] < hi)
            o[cls + '_bins'].append([lo, hi, int(m.sum()), float(x['r1'][m].mean()) if m.sum() > 30 else None,
                                     float(x['r_1'][m].mean()) if m.sum() > 30 else None])
    o['all'] = float(x['c1'].mean())
    return o
cases = [('det3 June-15', 0)] + [(b, a) for b in ('det3 May', 'det2', 'det1', 'det6') for a in (0, .1, .2, .3)]
for bu, a in cases:
    o = run(bu, a); OUT.append(o)
    f = lambda bb: ' '.join(f'{b[3]:.2f}/{b[4]:.2f}' if b[3] is not None else '    -    ' for b in bb)
    print(f'{bu:13s} a={a:.1f}  all {o["all"]:+.3f} | Ypad {o["ypad"][0]:+.3f}±{o["ypad"][1]:.3f} | Xpad {o["xpad"][0]:+.3f}±{o["xpad"][1]:.3f}')
    print(f'      Ypad by Ysh: {f(o["ypad_bins"])}\n      Xpad by Ysh: {f(o["xpad_bins"])}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inject_parity.json'), 'w'))
