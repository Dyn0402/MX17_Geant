#!/usr/bin/env python3
"""
delay_vs_rho.py -- the +-1 neighbour peak delay against the sheet resistivity, from the
chain's own W2 kernels (validated against COMSOL 3D in FRESH_LOOK §10.1).

Data: +54-61 ns, invariant across six SPS arms, two gases, x2 gain (DEEP_DIVE §5).
Sim (full chain, rho_s = 2 MOhm/sq): +93-142 ns. The kernels' timing is right (COMSOL),
so the question is which input sets the delay. The sheet spreading speed scales as
1/(rho_s C), so rho_s is first.

Per rho_s: avalanche clouds (Gaussian, sigma_T, both axes) centred uniformly over one
pad cell; each cloud's Y-row currents = sum over the cloud of the chain operator
((1-f) G + f ion (x) G, the exact ion-transit shape); DREAM 180 ns shaper; the peak
row's shaped peak time against rows +-1 (each event's mean of the two). Mean over events.

Usage: python3 delay_vs_rho.py <greens_w2_rho*.npz> ...   (from the repo root)
"""
import sys
import numpy as np

sys.path.insert(0, ".")
from response.common import constants as C
from response.solver import kernels as K
from response.dream.shaper import DreamShaper

GAP, LZ, T_ION, F0 = 150.0, 14.8, 340.0, 0.901
DT, TMAX = 1.0, 2500.0


def ion_shape(T_ION=T_ION):
    tg = np.arange(0, TMAX, DT)
    v = GAP / T_ION
    zc = np.clip(GAP - v * (tg + 0.5 * DT), 0, GAP)
    cdf = (1 - np.exp(-zc / LZ)) / (1 - np.exp(-GAP / LZ))
    return cdf / cdf.sum()


def run(npz, sigma_um=400.0, n_ev=200, seed=1, t_ion=T_ION, f_ion=F0):
    d = np.load(npz)
    t = d["t"]; x = d["x"]; yY = d["y_Y"]
    rho = float(__import__("json").loads(str(d["meta"]))["rho_s_ohm_sq"])
    cand = np.arange(0, C.N_PAD_PER_SUPER, 2)
    c0 = int(cand[np.argmin(np.abs((K.esl_phase(K.pad_x(cand)) + 400e-6) % 800e-6 - 400e-6))])
    x0c = float(K.pad_x(c0))
    # x window of +-3 mm around c0, every 4th sample (40 um)
    ix = np.where(np.abs(((x - x0c) + C.SUPERPERIOD_M / 2) % C.SUPERPERIOD_M - C.SUPERPERIOD_M / 2) < 3e-3)[0][::4]
    xr = ((x[ix] - x0c) + C.SUPERPERIOD_M / 2) % C.SUPERPERIOD_M - C.SUPERPERIOD_M / 2
    o = np.argsort(xr); ix, xr = ix[o], xr[o]
    iy = np.where(np.abs(yY) < 4e-3)[0]
    Ge = d["G_Y_even"][:, iy][:, :, ix].astype(float)     # (nt, ny, nx): row at y=0, even parity
    Go = d["G_Y_odd"][:, iy][:, :, ix].astype(float)
    yr = yY[iy]
    tg = np.arange(0, TMAX, DT)

    def to_current(G):
        q = np.empty((len(tg) + 1,) + G.shape[1:])
        e = np.append(tg, TMAX) * 1e-9
        flat = G.reshape(len(t), -1)
        # linear interpolation in time for all columns at once
        k = np.clip(np.searchsorted(t, e) - 1, 0, len(t) - 2)
        w = np.clip((e - t[k]) / (t[k + 1] - t[k]), 0, 1)[:, None]
        qi = flat[k] * (1 - w) + flat[k + 1] * w
        qi[0] = 0.0
        return (np.diff(qi, axis=0) / (DT * 1e-9)).reshape((len(tg),) + G.shape[1:])
    Ie, Io = to_current(Ge), to_current(Go)
    ion = ion_shape(t_ion)
    def op(I):   # chain operator along time, for every (y, x)
        fl = I.reshape(len(tg), -1)
        n = 1 << int(np.ceil(np.log2(2 * len(tg))))
        conv = np.fft.irfft(np.fft.rfft(fl, n, axis=0) * np.fft.rfft(ion, n)[:, None], n, axis=0)[:len(tg)]
        return ((1 - f_ion) * fl + f_ion * conv).reshape(I.shape)
    Ie, Io = op(Ie), op(Io)
    sh = DreamShaper(dt_ns=DT)

    rng = np.random.default_rng(seed)
    P = C.PAD_PITCH_M
    delays = []
    for _ in range(n_ev):
        xc, yc = rng.uniform(-P, P), rng.uniform(-P / 2, P / 2)   # one Y-row pad cell pair
        wx = np.exp(-0.5 * ((xr - xc) / (sigma_um * 1e-6)) ** 2)
        rows = {}
        for j in (-1, 0, 1):
            G = Ie if j % 2 == 0 else Io
            # channel row j sits at y = j P; the kernel is in the row's frame: dy = y_src - j P
            wy = np.exp(-0.5 * ((yr + j * P - yc) / (sigma_um * 1e-6)) ** 2)
            w = np.outer(wy, wx); w /= w.sum()
            cur = np.tensordot(G, w, axes=([1, 2], [0, 1]))
            rows[j] = sh.apply(cur, dt_ns=DT)
        amp = {j: rows[j].max() for j in rows}
        jc = max(amp, key=amp.get)
        if jc != 0:
            continue                       # keep events whose peak row is the centre row
        tp = {j: tg[np.argmax(rows[j])] for j in rows}
        delays.append(0.5 * ((tp[-1] - tp[0]) + (tp[1] - tp[0])))
    delays = np.array(delays)
    return rho, float(np.mean(delays)), float(np.median(delays)), len(delays)


def main_ion():
    """the same delay with faster ions / a weaker ion signal, at one rho_s"""
    npz = sys.argv[2]
    for sig in (200.0, 400.0):
        for t_ion, f in ((340.0, F0), (200.0, F0), (140.0, F0), (340.0, 0.2), (340.0, 0.0)):
            rho, m, med, n = run(npz, sigma_um=sig, t_ion=t_ion, f_ion=f)
            print(f"rho_s {rho/1e6:3.1f}  sigma_T {sig:4.0f} um  T_ion {t_ion:4.0f} ns  f_ion {f:.2f}:  "
                  f"+-1 delay mean {m:6.1f} ns, median {med:6.1f} ns ({n})", flush=True)


def main():
    if sys.argv[1] == "--ion":
        return main_ion()
    for sig in (200.0, 400.0, 800.0):
        for npz in sys.argv[1:]:
            rho, m, med, n = run(npz, sigma_um=sig)
            print(f"rho_s {rho/1e6:4.1f} MOhm/sq  sigma_T {sig:4.0f} um:  +-1 delay mean {m:6.1f} ns, "
                  f"median {med:6.1f} ns  ({n} events)", flush=True)


if __name__ == "__main__":
    main()
