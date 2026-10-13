# Protected-FEU input network as a filter in front of the DREAM shaper (FRESH_LOOK §17).
# Same toy + fit as scratchpad/2026-10-13/wf_ion_T_scan166.py (inclined det3 X shapes,
# 10.7 / 20 deg), with the shaper impulse response h replaced by
#   high-pass (series C into R to ground):  h - (1/tau) * (h (x) e^{-t/tau})
#   low-pass  (series R into strip C):      (1/tau) e^{-t/tau} (x) h
# Full ions (f_ion 0.908) at Ar+ (T 340) and hydrocarbon (T 233 = 340/1.46) transit.
# Also reports the undershoot: min of the normalised average after the peak.
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
t = np.arange(0, 3000, T.DT_FINE)

def net(h, kind, tau):
    if kind == 'none': return h
    k = np.exp(-t / tau) / tau * T.DT_FINE
    f = np.convolve(h, k)[:len(h)]
    g = h - f if kind == 'hp' else f
    return g / g.max()

def under(a):
    i = int(np.argmax(a)); return float(a[i:].min())

def job(c):
    tp, kind, tau, f, Tion = c
    T.shaper = lambda family, t_peak=180.0: net(_ORIG(family, tp), kind, tau)
    out = []
    for th in (20.0, 10.7):
        B = T.sampled_peak_strips({**T.BASE, 'shaper': 'manual_1', 'f_ion': f, 't_ion': Tion, 'theta': th}, 800, seed=4, target_peak=2000)[:, 3]
        a = avg(B); m = grid >= -2
        out.append((float(np.sqrt(np.mean((a[m] - D[th][m]) ** 2))), float(table(B).rise.median()), under(a)))
    return c, out

cfgs = [(tp, 'none', 0, 0.908, Ti) for tp in (166, 180) for Ti in (340., 233.)]
cfgs += [(tp, 'hp', tau, 0.908, Ti) for tp in (166, 180) for Ti in (340., 233.)
         for tau in (3000., 1000., 600., 400., 300., 200., 150., 100.)]
cfgs += [(166, 'lp', tau, 0.908, 233.) for tau in (5., 10., 20., 40.)]
with ProcessPoolExecutor(12) as ex: res = list(ex.map(job, cfgs))
print('data rise median 20deg %.0f  10deg %.0f   undershoot 20deg %.3f  10deg %.3f' % (R[20.0], R[10.7], under(D[20.0]), under(D[10.7])))
print(' tpk net   tau  f_ion T_ion | rms20 rise20 under20 | rms10 rise10 under10')
for c, o in res:
    print(f'{c[0]:4d} {c[1]:4s} {c[2]:5.0f} {c[3]:5.2f} {c[4]:5.0f} | {o[0][0]:.3f}  {o[0][1]:4.0f}  {o[0][2]:6.3f} | {o[1][0]:.3f}  {o[1][1]:4.0f}  {o[1][2]:6.3f}')
