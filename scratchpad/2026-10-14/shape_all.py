#!/usr/bin/env python3
"""Inclined-track shape of every sim point with no reco/quality selection: the T14 extractor
(t14_compare.extract_view, censored rise / FWHM medians) on ALL simulated events, against the data
by TRUE M3 angle (rise_vs_m3.log). Cached per dir as <dir>/shape_all_<view>.json. Runs on yoga.
Usage: shape_all.py "<glob of stageB_w2 sim dirs>" """
import os, sys, glob, json, numpy as np
sys.path[:0] = [os.path.expanduser('~/PycharmProjects/nTof_x17'), os.path.expanduser('~/PycharmProjects/nTof_x17/mx17_sim_wft')]
from t14_compare import extract_view, _leg_cfg
from wft.io import strip_position_map, FeuReader
SB = os.path.expanduser('~/x17/response_sim/stageB_w2')
cfg = _leg_cfg(f'{SB}/w2_rho2M_fullphys146/th00', 'sat_det3'); pm = strip_position_map(cfg)
DATA = {('x', 'thx10'): (187, 321), ('x', 'thx20'): (169, 295), ('y', 'thy10'): (188, 310), ('y', 'thy20'): (160, 286)}
res = {}
for d in sorted(glob.glob(f'{SB}/{sys.argv[1]}')):
    ang = os.path.basename(d); kind = os.path.basename(os.path.dirname(d)).replace('w2_rho2M_', '')
    for (v, a), (dr, dw) in DATA.items():
        if a != ang: continue
        c = f'{d}/shape_all_{v}.json'
        if os.path.exists(c): r = json.load(open(c))
        else:
            feu = cfg.MX17_FEU_X if v == 'x' else cfg.MX17_FEU_Y
            f = glob.glob(f'{d}/decoded_root/*_{feu:02d}.root')
            if not f: continue
            ids = set(int(e) for e in FeuReader(f[0]).event_ids)
            t = extract_view(d, feu, pm[feu], ids, f'{kind}/{ang}/{v}')['table']
            r = dict(n=len(t), rise=float(t.rise_ns.median()), fwhm=float(t.fwhm_ns.median()),
                     rise_raw=float(t.rise_ns_raw.median()), sat=float(t.saturated.mean()))
            json.dump(r, open(c, 'w'))
        res[(kind, v, ang)] = r
        print(f'{kind:16s} {v}{ang[3:]:>3s}  rise {r["rise"]:5.0f} ({r["rise"] - dr:+4.0f})  fwhm {r["fwhm"]:5.0f} ({r["fwhm"] - dw:+4.0f})  sat {r["sat"]:.2f}', flush=True)
