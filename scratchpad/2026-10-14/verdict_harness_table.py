import json, os
SB = os.path.expanduser('~/x17/response_sim/stageB_w2')
rows = [('fullphys175 (Aug noise)', 't14_fullphys175_{a}'), ('noisecm146', 't14_noisecm146_{a}'), ('verdict1', 't14_verdict1_{a}')]
angs = ('th00', 'thx10', 'thx20', 'thy10', 'thy20')
for view in ('x', 'y'):
    print(f'=== view {view}: rise sim/data | FWHM sim/data | amp data/sim')
    print(f'{"":24s}' + ''.join(f'{a:>30s}' for a in angs))
    for lab, pat in rows:
        cells = []
        for a in angs:
            f = os.path.join(SB, pat.format(a=a), 't14_summary.json')
            if not os.path.exists(f):
                cells.append(f'{"-":>30s}'); continue
            v = json.load(open(f))['views']; v = eval(v) if isinstance(v, str) else v
            x = v[view]
            r, w, p = x['rise_ns_med'], x['fwhm_ns_med'], x['peak_amp_med']
            cells.append(f'{r["sim"]:4.0f}/{r["data"]:3.0f} | {w["sim"]:4.0f}/{w["data"]:3.0f} | {p["ratio"]:4.2f}'.rjust(30))
        print(f'{lab:24s}' + ''.join(cells))
