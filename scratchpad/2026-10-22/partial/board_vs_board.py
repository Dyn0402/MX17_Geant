"""Direct det2 (6-22 clean) vs det3 June (6-27 clean) comparison over the cells both fill (det2_rescore.json): chi2 per point
type and the weighted mean difference of late/peak and r25/peak ratios (2026-10-22, §13). Usage: board_vs_board.py"""
import json, numpy as np
R = json.load(open('det2_rescore.json')); out = {}
for sel in ('dtn', 'm3v'):
    A = {(tuple(r['fw']), tuple(r['s1'])): r for r in R[f'det2 6-22 (s004+s006) [{sel}]']['rows']}
    B = {(tuple(r['fw']), tuple(r['s1'])): r for r in R[f'det3 June 6-27 1000 V (s051+s025) [{sel}]']['rows']}
    ks = sorted(set(A) & set(B)); res = {}
    for p in ('r05', 'r25', 'peak', 'f50'):
        res[p] = float(sum((A[k][p][0] - B[k][p][0]) ** 2 / (A[k][p][1] ** 2 + B[k][p][1] ** 2) for k in ks))
    def ratio(k, D, p):
        v = D[k][p][0] / D[k]['peak'][0]; e = v * np.hypot(D[k][p][1] / D[k][p][0], D[k]['peak'][1] / D[k]['peak'][0]); return v, e
    for p in ('r25', 'f50'):
        d = np.array([ratio(k, A, p)[0] - ratio(k, B, p)[0] for k in ks]); e = np.array([np.hypot(ratio(k, A, p)[1], ratio(k, B, p)[1]) for k in ks])
        w = 1 / e ** 2; res[f'{p}/pk diff'] = [float((d * w).sum() / w.sum()), float(w.sum() ** -.5)]
    out[sel] = dict(cells=[list(map(list, k)) for k in ks], **res)
    print(f"[{sel}] {len(ks)} common cells: chi2 r05 {res['r05']:.1f}  r25 {res['r25']:.1f}  peak {res['peak']:.1f}  f50 {res['f50']:.1f}  (each / {len(ks)});  "
          f"det2 - det3June: r25/pk {res['r25/pk diff'][0]:+.2f} ± {res['r25/pk diff'][1]:.2f}   f50/pk {res['f50/pk diff'][0]:+.2f} ± {res['f50/pk diff'][1]:.2f}")
json.dump(out, open('board_vs_board.json', 'w'), indent=1)
