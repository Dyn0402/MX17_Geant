#!/usr/bin/env python3
"""Tabulate the T14 X/Y summaries for the ion-speed / full-physics scans against the s = 1
angled baseline (t14_angW_*). Run on the laptop: python3 t14_ionspeed_table.py > table.txt"""
import json, os
SB = os.path.expanduser('~/x17/response_sim/stageB_w2')
rows = [('s=1 (Ar+), no fb', 't14_angW_{a}')] + [(f's={s / 100:.2f}, no fb', f't14_ionspeed{s}_{{a}}') for s in (146, 175, 245)] \
     + [(f's={s / 100:.2f}, bG 0.45', f't14_fullphys{s}_{{a}}') for s in (146, 175)]
for view in ('x', 'y'):
    print(f'=== view {view}: sim / data  (rise median ns | FWHM median ns | peak-amp median, data/sim)')
    print(f'{"":20s}' + ''.join(f'{a:>34s}' for a in ('th00', 'thx10', 'thx20')))
    for lab, pat in rows:
        cells = []
        for a in ('th00', 'thx10', 'thx20'):
            f = os.path.join(SB, pat.format(a=a), 't14_summary.json')
            if not os.path.exists(f):
                cells.append(f'{"pending":>34s}'); continue
            v = json.load(open(f))['views']; v = eval(v) if isinstance(v, str) else v
            x = v[view]
            r, w, p = x['rise_ns_med'], x['fwhm_ns_med'], x['peak_amp_med']
            cells.append(f'{r["sim"]:4.0f}/{r["data"]:3.0f} | {w["sim"]:4.0f}/{w["data"]:3.0f} | x{1 / p["ratio"]:4.2f}'.rjust(34))
        print(f'{lab:20s}' + ''.join(cells))
