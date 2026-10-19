"""Inverse test (2026-10-20): undo ONE uniform channel copy a on det3 June-15 / det4 / det7 X (asym_lib.uncopy, exact
inverse of x(n+1) += a x(n)), re-find the peak on the uncopied window, redo the parity split, and measure the
asymmetry on the UNCOPIED window. A plain channel copy => some single a makes both classes (and every Y-share bin)
symmetric. A pad-type-local effect => the Y-pad class needs a much larger a than the X-pad class.
Controls get the same treatment with a = their own fitted (~0) value. Usage: uncopy_parity.py"""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L, asym_more as A
rows = M.survey_rows(); ED = [0, .4, .55, .7, .85, 1.0, 1.3, 9]
CACHE = {}
def events(bu):
    if bu not in CACHE:
        Wl, py, px, ysh = [], [], [], []
        for rx, ry in A.pairs(rows, lambda r: r['build'] == bu):
            dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
            if len(ix) < 20: continue
            Wl.append(dx['W'][ix].astype(float)); px.append(dx['p'][ix]); py.append(dy['p'][iy]); ysh.append((dy['al'][1] + dy['al'][-1])[iy])
        CACHE[bu] = [np.concatenate(v) for v in (Wl, px, py, ysh)]
    return CACHE[bu]
def run(bu, a, gen=None):
    W, px0, py, ysh = events(bu)
    if gen: W = W.copy(); W[:, 1:] += gen * W[:, :-1].copy()     # inject (controls) before undoing
    U = L.uncopy(W, a) if a else W
    A_ = U.max(2); sh = A_[:, 3:6].argmax(1) - 1; c0 = 4 + sh; ii = np.arange(len(U))
    c = U[ii, c0]; n2 = np.einsum('ij,ij->i', c, c)
    r1 = np.einsum('ij,ij->i', U[ii, c0 + 1], c) / n2; r_1 = np.einsum('ij,ij->i', U[ii, c0 - 1], c) / n2
    px = px0 + sh; par = (px + py) % 2; ypar = int(ysh[par == 1].mean() < ysh[par == 0].mean()); yp = par == ypar
    ok = (px % 64 >= 4) & (px % 64 <= 59)
    o = dict(build=bu, a=a, gen=gen)
    for cls, m0 in (('ypad', yp & ok), ('xpad', ~yp & ok)):
        d = (r1 - r_1)[m0]; o[cls] = [float(d.mean()), float(d.std() / np.sqrt(len(d))), int(len(d))]
        o[cls + '_bins'] = [[lo, hi, int((m := m0 & (ysh >= lo) & (ysh < hi)).sum()),
                             float((r1 - r_1)[m].mean()) if m.sum() > 30 else None,
                             float((r1 - r_1)[m].std() / np.sqrt(max(m.sum(), 1)))] for lo, hi in zip(ED[:-1], ED[1:])]
    return o
if __name__ == '__main__':
    OUT = []
    for bu, gen in (('det3 June-15', None), ('det4', None), ('det7', None), ('det2', .3), ('det3 May', .3), ('det2', None)):
        for a in (0, .1, .15, .2, .25, .3, .35):
            if gen is None and bu == 'det2' and a: continue
            o = run(bu, a, gen); OUT.append(o)
            f = lambda bb: ' '.join(f'{b[3]:+.2f}' if b[3] is not None else '  -  ' for b in bb)
            print(f'{bu:13s}{" inj%.1f" % gen if gen else "       "} undo a={a:.2f} | Ypad {o["ypad"][0]:+.3f}±{o["ypad"][1]:.3f} [{f(o["ypad_bins"])}]'
                  f' | Xpad {o["xpad"][0]:+.3f}±{o["xpad"][1]:.3f} [{f(o["xpad_bins"])}]')
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uncopy_parity.json'), 'w'))
