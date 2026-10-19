"""Are the May det3 and the June-15 det3 the same PCB? Dead-line fingerprint (2026-10-20).
From the raw survey (all sub-runs, unvetoed loads): per build, view and FEU slot, how often each TRUE pad column/row is
the peak of an isolated event; a line is 'dead' when it is the peak < 15 % as often as the median of its +-8 neighbours
(and the neighbours are busy). Lines dead on two FEU slots of the same build follow the board. Then compare the sets
between builds. Usage: pcb_fingerprint.py"""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L
COL = np.array([64 * (c // 64) + 63 - c % 64 for c in range(512)])
rows = M.survey_rows(); occ = {}
for r in rows:
    z = np.load(r['f']); p = z['p'].astype(int) % 512
    k = (r['build'], r['view'], r['feu']); occ.setdefault(k, np.zeros(512))
    occ[k] += np.bincount(COL[p], minlength=512)
def dead(h):
    out = []
    for c in range(512):
        if c % 64 < 4 or c % 64 > 59: continue          # extractor never keeps peaks there
        lo, hi = max(0, c - 8), min(512, c + 9)
        nb = np.r_[h[lo:c], h[c + 1:hi]]; nb = nb[[(lo + i) % 64 >= 4 and (lo + i) % 64 <= 59 for i in range(len(nb))]] if False else nb
        med = np.median(nb)
        if med >= 8 and h[c] < 0.15 * med: out.append(c)
    return out
D_ = {}
for (bu, v, feu), h in sorted(occ.items()):
    d = dead(h); D_[(bu, v, feu)] = d
    print(f'{bu:13s} {v} FEU {feu}: {int(h.sum()):6d} peaks, dead lines {d}')
OUT = {f'{k[0]}|{k[1]}|{k[2]}': v for k, v in D_.items()}
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pcb_fingerprint.json'), 'w'))
