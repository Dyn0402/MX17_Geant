"""chi2 of each data set vs the sim over the common (FWHM, s1) cells of rescore.json, points r05, r25, peak, f50
(errors in quadrature). Also the early/peak and late/peak ratios averaged over the s1 0.45-0.7 cells. Usage: rescore_chi2.py"""
import json, numpy as np
R = json.load(open('rescore.json')); out = {}
for sel, sim in (('dtn', 'sim glg0 [dtn]'), ('m3v', 'sim glg0 th00 [m3v]')):
    S = {(tuple(r['fw']), tuple(r['s1'])): r for r in R[sim]['rows']}
    for name, v in R.items():
        if not name.endswith(f'[{sel}]') or name.startswith('sim'): continue
        chi, nd, el = 0., 0, []
        for r in v['rows']:
            k = (tuple(r['fw']), tuple(r['s1']))
            if k not in S: continue
            for p in ('r05', 'r25', 'peak', 'f50'):
                chi += (r[p][0] - S[k][p][0]) ** 2 / (r[p][1] ** 2 + S[k][p][1] ** 2); nd += 1
            if k[1][0] == .45: el.append((r['r05'][0] / r['peak'][0], r['f50'][0] / r['peak'][0], S[k]['r05'][0] / S[k]['peak'][0], S[k]['f50'][0] / S[k]['peak'][0]))
        el = np.array(el); out[name] = dict(chi2=chi, ndof=nd, early_pk=el[:, 0].mean() if len(el) else float("nan"), late_pk=el[:, 1].mean() if len(el) else float("nan"),
                                            sim_early_pk=el[:, 2].mean() if len(el) else float("nan"), sim_late_pk=el[:, 3].mean() if len(el) else float("nan"))
        print(f'{name:45s} chi2 {chi:6.1f} / {nd:2d}   s1 0.45-0.7 cells: early/pk {out[name]["early_pk"]:.2f} (sim {out[name]["sim_early_pk"]:.2f})  late/pk {out[name]["late_pk"]:.2f} (sim {out[name]["sim_late_pk"]:.2f})')
json.dump(out, open('rescore_chi2.json', 'w'), indent=1)
