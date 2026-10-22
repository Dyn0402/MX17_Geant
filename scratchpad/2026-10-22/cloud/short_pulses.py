"""What are the data's short central pulses (FWHM < 350 ns, |dt_nb| < 80 ns) that the sim lacks? (2026-10-22) Per FWHM
class: onset (5 % rise) and 50 %-fall sample relative to the trigger (60 ns samples; fixed trigger -> arrival time =
depth), central peak amplitude, s1, and the fraction of all |dt_nb| < 80 events. A partial drift column shows up as an
onset late (deep only) or an early fall (shallow only) compared with the long pulses.
Usage: short_pulses.py   (-> short_pulses.json)"""
import os, sys, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import prompt_share as PS
from fwhm_bins import load, HW
RAW = os.path.expanduser('~/x17/nb_asym_raw')
SETS = {'det3 May X': [f'{RAW}/d3may_x.npz'], 'det3 May Y': [f'{RAW}/d3may_y.npz'], 'det3 June-15 X': [f'{RAW}/d3jun_x.npz'],
        'det3 June-15 Y': [f'{RAW}/d3jun_y.npz'], 'sim glg0 X': sorted(glob.glob(f'{RAW}/simgapleak/glg0_th*_x.npz')),
        'sim glg0 Y': sorted(glob.glob(f'{RAW}/simgapleak/glg0_th00_*_y.npz'))}
import json; OUT = {}
for name, files in SETS.items():
    W, s1, dtn, fw, pk = load(files); c = W[:, HW]; on = PS.cross(c, pk, .05); off = PS.cross(c, pk, .5, False)
    A = c[np.arange(len(c)), pk]; v = np.abs(dtn) < 80
    print(f'\n{name}  (|dt_nb| < 80: {v.sum()})\n   FWHM ns    frac   onset p25/50/75 (samples)   peak p50   f50 p50   amp p50   s1 p50')
    for lo, hi in ((0, 350), (350, 450), (450, 550), (550, 2000)):
        m = v & (fw >= lo) & (fw < hi)
        if m.sum() < 15: continue
        q = lambda x: np.nanpercentile(x[m], [25, 50, 75])
        OUT.setdefault(name, []).append(dict(fw=[lo, hi], frac=float(m.mean() / v.mean()), onset=[float(x) for x in q(on)], peak=float(np.median(pk[m])),
                                             f50=float(np.nanmedian(off[m])), amp=float(np.median(A[m])), s1=float(np.median(s1[m]))))
        print(f'   {lo:4d}-{hi:<5d} {m.mean() / v.mean():5.2f}   {q(on)[0]:5.1f} {q(on)[1]:5.1f} {q(on)[2]:5.1f}            {np.median(pk[m]):5.1f}    '
              f'{np.nanmedian(off[m]):5.1f}    {np.median(A[m]):6.0f}   {np.median(s1[m]):.2f}')
json.dump(OUT, open(os.path.join(HERE, 'short_pulses.json'), 'w'), indent=1)
