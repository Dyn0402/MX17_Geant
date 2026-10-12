# §7.5 rerun at the noise-measured shaper time scale (FRESH_LOOK §10.2: ~166 ns, s = 0.92),
# against the inclined (10.7 / 20 deg) det3 shapes. Same toy, same fit as
# scratchpad/2026-10-12/wf_ion_T_scan.py; adds 166 ns and an f_ion scan.
# Run with the nTof_x17 .venv, both npz files in the cwd.
import sys, numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_tail20 import avg, grid
import toy_pulse as T
za = np.load('waveforms_angled_pm3.npz')
D = {th: avg(za[f'{tag}_data_wf'][:, 3]) for tag, th in (('thx20', 20.0), ('thx10', 10.7))}
from wf_matched import table
R = {th: table(za[f'{tag}_data_wf'][:, 3]).rise.median() for tag, th in (('thx20', 20.0), ('thx10', 10.7))}
_ORIG = T.shaper
def job(c):
    fam, tp, f, Tion = c
    T.shaper = lambda family, t_peak=180.0: _ORIG(family, tp)
    out = []
    for th in (20.0, 10.7):
        B = T.sampled_peak_strips({**T.BASE, 'shaper': fam, 'f_ion': f, 't_ion': Tion, 'theta': th}, 800, seed=4, target_peak=2000)[:, 3]
        a = avg(B); m = grid >= -2
        out.append((float(np.sqrt(np.mean((a[m] - D[th][m]) ** 2))), float(table(B).rise.median())))
    return c, out
cfgs = [('manual_1', tp, 0.908, Ti) for tp in (166, 180) for Ti in (340., 260., 200., 170., 140., 110.)]
cfgs += [('manual_1', 166, f, 340.) for f in (0.0, 0.2, 0.4, 0.6, 0.8)]
with ProcessPoolExecutor(12) as ex: res = list(ex.map(job, cfgs))
print('data rise median 20deg %.0f  10deg %.0f' % (R[20.0], R[10.7]))
print('shaper   tpk  f_ion T_ion | shape-rms20 rise20 | shape-rms10 rise10')
for c, o in sorted(res, key=lambda r: (r[0][1], -r[0][2], -r[0][3])):
    print(f'{c[0]:8s} {c[1]:4d} {c[2]:5.2f} {c[3]:5.0f} |   {o[0][0]:.3f}    {o[0][1]:4.0f}  |   {o[1][0]:.3f}    {o[1][1]:4.0f}')
