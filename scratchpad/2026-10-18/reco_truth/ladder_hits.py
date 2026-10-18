#!/usr/bin/env python3
"""The ladder of ladder.py, but from the PRODUCTION hit times (combined_hits `time`, and `time_of_max`) on the det3
data. Checks that ladder.py's waveform CFD50 stands in for what the hits chain actually used.
Same selection: >= 5 sigma, >= 100 ADC, no saturated hit in the view, >= 4 hits; same M3 angle windows.
Usage: ladder_hits.py OUT.json"""
import sys, glob, json
import numpy as np
import uproot
sys.argv, OUT = sys.argv[:1], sys.argv[1]
from ladder import cfg, PM, V_DATA, summarise, M3_CHI2_CUT, M3_MIN_NCLUS

def main():
    from M3RefTracking import M3RefTracking
    rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
    x0, y0, ev = map(np.asarray, rays.get_xy_positions(702.0))
    x1, y1, _ = map(np.asarray, rays.get_xy_positions(802.0))
    ax = dict(zip(ev.astype(int), np.degrees(np.arctan((x1 - x0) / 100))))
    ay = dict(zip(ev.astype(int), np.degrees(np.arctan((y1 - y0) / 100))))
    th = {'x': ay, 'y': {k: -v for k, v in ax.items()}}
    rows = {'x': [], 'y': []}
    for f in sorted(glob.glob(cfg.combined_hits_dir + '*combined_hits.root')):
        d = uproot.open(f)['hits'].arrays(['eventId', 'channel', 'feu', 'time', 'time_of_max', 'amplitude',
                                           'saturated', 'significance'], library='pd')
        for view, feu in (('x', 7), ('y', 8)):
            other = 'y' if view == 'x' else 'x'
            g = d[d.feu == feu]
            for e, h in g.groupby('eventId'):
                e = int(e)
                if e not in th[view] or e not in th[other] or h.saturated.any() or (h.amplitude >= 3500).any():
                    continue
                h = h[(h.significance >= 5) & (h.amplitude >= 100)]
                if len(h) < 4:
                    continue
                p = PM[feu][h.channel.to_numpy()]; a = h.amplitude.to_numpy()
                r = {'th': th[view][e], 'tho': th[other][e]}
                for nm, col in (('cfd', 'time'), ('pk', 'time_of_max')):
                    t = h[col].to_numpy(float)
                    xm = (a * p).sum() / a.sum(); ym = (a * t).sum() / a.sum()
                    r[nm] = (a * (p - xm) * (t - ym)).sum() / (a * (p - xm) ** 2).sum()
                    b = (a * (p - xm) * (t - ym)).sum() / (a * (t - ym) ** 2).sum()
                    r[nm + '_xt'] = 1.0 / b if b != 0 else np.nan
                rows[view].append(r)
    res = {}
    for view in 'xy':
        for a in (5, 10, 15, 20):
            sel = [r for r in rows[view] if abs(abs(r['th']) - a) < (2 if a == 5 else 3) and abs(r['tho']) < 3]
            res[f'th{view}{a:02d}'] = summarise(sel, np.tan(np.radians(a)), V_DATA)
            print(view, a, json.dumps(res[f'th{view}{a:02d}']), flush=True)
    json.dump(res, open(OUT, 'w'), indent=1)

main()
