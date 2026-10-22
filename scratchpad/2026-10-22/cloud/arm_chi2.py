"""Rank the cloud arms against det3 May (X and Y) in the absolute-FWHM x s1 cells of fwhm_bins.py (2026-10-22): chi2 of
the time-profile SHAPE (r05/peak, r25/peak, f50/peak) and of the level (peak), over cells where data and arm both
have >= 25 events. Usage: arm_chi2.py [fwhm_bins_dtn80.json]"""
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'fwhm_bins_dtn80.json')))
def cells(rows): return {(tuple(r['fw']), tuple(r['s1'])): r for r in rows}
def q(r, k):  # ratio to peak and its error (errors treated as independent: conservative)
    v = r[k][0] / r['peak'][0]; return v, v * np.hypot(r[k][1] / r[k][0] if r[k][0] else 0, r['peak'][1] / r['peak'][0])
arms = sorted({k.rsplit(' ', 1)[0] for k in R if k.startswith('sim ')})
print('arm             view  cells   chi2 shape/ndf (r05,r25,f50 over peak)   chi2 peak/ndf')
for arm in arms:
    tot = [0, 0, 0, 0]
    for view in ('X', 'Y'):
        d, s = cells(R[f'det3 May {view}']), cells(R.get(f'{arm} {view}', []))
        cs, ns, cp, n = 0., 0, 0., 0
        for key in d.keys() & s.keys():
            for k in ('r05', 'r25', 'f50'):
                (a, ea), (b, eb) = q(d[key], k), q(s[key], k); cs += (a - b) ** 2 / (ea ** 2 + eb ** 2); ns += 1
            cp += (d[key]['peak'][0] - s[key]['peak'][0]) ** 2 / (d[key]['peak'][1] ** 2 + s[key]['peak'][1] ** 2); n += 1
        print(f'{arm:15s} {view}    {n:3d}     {cs:7.1f} / {ns:2d} = {cs / max(ns, 1):5.2f}                    {cp:6.1f} / {n:2d}')
        tot = [tot[0] + cs, tot[1] + ns, tot[2] + cp, tot[3] + n]
    print(f'{arm:15s} X+Y  {tot[3]:3d}     {tot[0]:7.1f} / {tot[1]:2d} = {tot[0] / max(tot[1], 1):5.2f}                    {tot[2]:6.1f} / {tot[3]:2d}\n')

# the same on the X cells common to every arm (the arms differ in which cells reach 25 events)
common = set.intersection(*[set(cells(R[f'{a} X'])) for a in arms]) & set(cells(R['det3 May X']))
print(f'X cells common to all arms: {sorted(common)}')
for arm in arms:
    d, s = cells(R['det3 May X']), cells(R[f'{arm} X']); cs = 0.; ns = 0
    for key in common:
        for k in ('r05', 'r25', 'f50'):
            (a, ea), (b, eb) = q(d[key], k), q(s[key], k); cs += (a - b) ** 2 / (ea ** 2 + eb ** 2); ns += 1
    print(f'  {arm:15s} shape chi2 {cs:6.1f} / {ns}')
