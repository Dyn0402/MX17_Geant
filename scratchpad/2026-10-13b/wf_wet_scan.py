# Does the bench's moisture change the inclined-track rise puzzle? (2026-10-13)
# Same toy, data and fit as scratchpad/2026-10-13/wf_ion_T_scan166.py (166 ns shaper),
# but the drift transport is a parameter: v, dT, dL (toy_pulse hardcodes dry 95/5).
# Gas arms at det3's real field (1000 V / 27.9 mm = 358 V/cm), response/params/wet_ariso_drift.json:
#   dry        v 39.14  dT 0.04766  dL 0.02425  (the toy's/chain's dry table)
#   wet 1 %    v 35.91  dT 0.0306   dL 0.0258   (matches the bench's measured v ~35-37)
#   v only     wet v, dry diffusion   (cloud_basics §27: bench widths look dry)
#   dT only    dry v, wet dT
# Run with the nTof_x17 .venv, waveforms_angled_pm3.npz in the cwd.
import sys, math, numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_tail20 import avg, grid
import toy_pulse as T
from wf_matched import table

za = np.load('waveforms_angled_pm3.npz')
D = {th: avg(za[f'{tag}_data_wf'][:, 3]) for tag, th in (('thx20', 20.0), ('thx10', 10.7))}
R = {th: table(za[f'{tag}_data_wf'][:, 3]).rise.median() for tag, th in (('thx20', 20.0), ('thx10', 10.7))}


def simulate(p, n_events, rng, h):
    """toy_pulse.simulate with dT, dL taken from p (identical otherwise)."""
    k_self, _ = T.electron_kernel(h, p["f_ion"], p["t_ion"], 0.0, p["tau_s"])
    v, dT, dL = p["v"], p["dT"], p["dL"]
    tan = math.tan(math.radians(p["theta"]))
    n_half = 6 + int(abs(T.GAP_MM * tan) / T.PITCH)
    strips = np.arange(-n_half, n_half + 1)
    lo, hi = (strips - 0.5) * T.PITCH, (strips + 0.5) * T.PITCH
    rows = []
    for _ in range(n_events):
        x0 = rng.uniform(-0.5, 0.5) * T.PITCH - 0.5 * T.GAP_MM * tan
        ncl = rng.poisson(91.2 / 3.83 * T.GAP_MM / 10)
        z = np.repeat(rng.uniform(0, T.GAP_MM, ncl), rng.geometric(1 / 3.83, ncl))
        x = x0 + z * tan + rng.normal(0, 1, z.size) * dT * np.sqrt(z / 10) * 10
        t = 300.0 + z / v + rng.normal(0, 1, z.size) * dL * np.sqrt(z / 10) * 10 / v
        g = rng.gamma(1.5, 1 / 1.5, z.size)
        sp = p["sigma_p"]
        cdf = lambda u: 0.5 * (1 + T._erf(u / (sp * math.sqrt(2))))
        w = cdf(hi[:, None] - x[None, :]) - cdf(lo[:, None] - x[None, :])
        it = np.clip(np.round(t / T.DT_FINE).astype(int), 0, T.N_FINE - 1)
        dep = np.array([np.bincount(it, weights=g * w[s], minlength=T.N_FINE)[:T.N_FINE]
                        for s in range(len(strips))])
        F = np.fft.rfft(dep, n=2 * T.N_FINE)
        rows.append(np.fft.irfft(F * np.fft.rfft(k_self, 2 * T.N_FINE), 2 * T.N_FINE)[:, :T.N_FINE])
    return rows


T.simulate = simulate
_ORIG = T.shaper
GAS = {'dry': (0.03914, 0.04766, 0.02425), 'wet1%': (0.03591, 0.0306, 0.0258),
       'v only': (0.03591, 0.04766, 0.02425), 'dT only': (0.03914, 0.0306, 0.02425)}


def job(c):
    gas, f, Tion = c
    v, dT, dL = GAS[gas]
    T.shaper = lambda family, t_peak=180.0: _ORIG(family, 166)
    out = []
    for th in (20.0, 10.7):
        B = T.sampled_peak_strips({**T.BASE, 'shaper': 'manual_1', 'f_ion': f, 't_ion': Tion,
                                   'theta': th, 'v': v, 'dT': dT, 'dL': dL}, 800, seed=4, target_peak=2000)[:, 3]
        a = avg(B); m = grid >= -2
        out.append((float(np.sqrt(np.mean((a[m] - D[th][m]) ** 2))), float(table(B).rise.median())))
    return c, out


cfgs = [(g, f, Ti) for g in GAS for (f, Ti) in ((0.908, 340.), (0.908, 200.), (0.908, 140.), (0.2, 340.))]
if __name__ == '__main__':
    with ProcessPoolExecutor(12) as ex:
        res = list(ex.map(job, cfgs))
    print('shaper 166 ns; data rise median 20deg %.0f  10deg %.0f' % (R[20.0], R[10.7]))
    print('gas      f_ion T_ion | shape-rms20 rise20 | shape-rms10 rise10')
    for c, o in res:
        print(f'{c[0]:8s} {c[1]:5.2f} {c[2]:5.0f} |   {o[0][0]:.3f}    {o[0][1]:4.0f}  |   {o[1][0]:.3f}    {o[1][1]:4.0f}')
