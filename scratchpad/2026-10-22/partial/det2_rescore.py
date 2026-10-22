"""§13 step 1: det2's §8 timing regime ("narrow-prompt / late-wide") re-scored with M3 angles on its clean sessions
(2026-10-22). det2 = mx17_2 on 6-22 (top slot, z 702, 495 V / 1000 V = sim field): survey s004 (long_run) + s006 (longer_run).
The degraded 6-17 det2 (bottom slot, s001) for reference. Same cells and chi2 as rescore.py / rescore_chi2.py.
Usage: det2_rescore.py   (-> det2_rescore.json)"""
import json, numpy as np
import rescore as RS, slot_table as ST
RS.SESS.update({s['tag']: s for s in ST.sessions()})
R = {}
sets = {'det2 6-22 (s004+s006)': RS.data_set(['s004', 's006']), 'det2 6-17 bottom, degraded (s001)': RS.data_set(['s001']),
        'det3 June 6-27 1000 V (s051+s025)': RS.data_set(['s051', 's025'])}
sims = {'dtn': RS.sim_set(f'{RS.RAW}/simgapleak/glg0_th*_x.npz'), 'm3v': RS.sim_set(f'{RS.RAW}/simgapleak/glg0_th00_*_x.npz')}
for sel in ('dtn', 'm3v'):
    S = RS.table(f'sim glg0 [{sel}]', sims[sel], sel); R[f'sim [{sel}]'] = S; Sk = {(tuple(r['fw']), tuple(r['s1'])): r for r in S['rows']}
    for n, D in sets.items():
        T = RS.table(n, D, sel); chi, nd, el = 0., 0, []
        for r in T['rows']:
            k = (tuple(r['fw']), tuple(r['s1']))
            if k not in Sk: continue
            for p in ('r05', 'r25', 'peak', 'f50'): chi += (r[p][0] - Sk[k][p][0]) ** 2 / (r[p][1] ** 2 + Sk[k][p][1] ** 2); nd += 1
            if k[1][0] == .45: el.append((r['r05'][0] / r['peak'][0], r['r25'][0] / r['peak'][0], r['f50'][0] / r['peak'][0],
                                          Sk[k]['r05'][0] / Sk[k]['peak'][0], Sk[k]['r25'][0] / Sk[k]['peak'][0], Sk[k]['f50'][0] / Sk[k]['peak'][0]))
        el = np.array(el) if el else np.full((1, 6), np.nan); T.update(chi2=chi, ndof=nd, ratios=el.mean(0).tolist()); R[f'{n} [{sel}]'] = T
for k, v in R.items():
    if 'chi2' in v: print(f'{k:48s} chi2 {v["chi2"]:6.1f} / {v["ndof"]:2d}   s1 0.45-0.7: r05/pk {v["ratios"][0]:.2f} (sim {v["ratios"][3]:.2f})  r25/pk {v["ratios"][1]:.2f} (sim {v["ratios"][4]:.2f})  f50/pk {v["ratios"][2]:.2f} (sim {v["ratios"][5]:.2f})')
json.dump(R, open('det2_rescore.json', 'w'), indent=1)
