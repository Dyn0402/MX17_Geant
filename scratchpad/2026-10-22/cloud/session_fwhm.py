"""Per session (survey rows), det3 June-15 / det3 May / det2, X and Y: FWHM scale (p90 of the central FWHM, the near-
vertical end) and r(f) for the LONGEST pulses of each session (FWHM >= 0.8 x its own p90, |dt_nb| < 80 ns), all s1
(2026-10-22). Tests whether det3 June X's flat r(f) survives once each session is cut to its own vertical tracks.
Usage: session_fwhm.py"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt')); sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-20', 'r5'))
import prompt_share as PS, r5lib
from fwhm_bins import load, HW
SETS = ('det3 June-15', 'det3 May', 'det2')
print('build | run sub | mesh drift | view   n_all  FWHM p50/p90 | n_long  <s1>   r05   r25   peak  f50   early/pk late/pk')
for r in r5lib.rows():
    if r['build'] not in SETS: continue
    try: W, s1, dtn, fw, pk = load([r['f']])
    except Exception: continue
    if len(s1) < 150: continue
    p50, p90 = np.nanpercentile(fw, [50, 90]); m = (fw >= .8 * p90) & (np.abs(dtn) < 80) & (s1 > .1) & (s1 < 1.)
    if m.sum() < 25: continue
    c = W[:, HW]; S = W[:, HW - 1] + W[:, HW + 1]; v = []
    for t in (PS.cross(c, pk, .05), PS.cross(c, pk, .25), pk.astype(float), PS.cross(c, pk, .5, False)):
        nu, de = PS.at(S, t), PS.at(c, t); g = m & np.isfinite(nu) & np.isfinite(de); v.append(nu[g].sum() / de[g].sum())
    print(f"{r['build']:12s} | {r['run'][-14:]} {r['sub'][:20]:20s} | {r['mesh']} {r['drift']:>5} | {r['view']}  {len(s1):5d}  {p50:4.0f}/{p90:4.0f} | {m.sum():4d}  "
          f"{s1[m].mean():.2f}  " + '  '.join(f'{x:.3f}' for x in v) + f'   {v[0] / v[2]:.2f}  {v[3] / v[2]:.2f}')
