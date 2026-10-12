#!/usr/bin/env python3
"""
Fresh-look step 1 (2026-10-12): an independent toy forward model of the det3
peak-strip waveform, to ask whether ANY ion signal of the modelled size can
give the data's 10-90 % rise through a 180 ns DREAM shaper.

Deliberately shares NO code with response/ (numpy only). Inputs are numbers:

  gas      Ar/iC4H10 95/5 dry Magboltz at 333 V/cm: v 39.14 um/ns,
           dT 0.04766, dL 0.02425 sqrt(cm)   (design/gas/drift_velocity_*.json)
  ionis.   91.2 e/cm, 3.83 e/cluster          (primary_ionisation_2026-08-10)
  gap      30 mm drift, pitch 0.78 mm, prompt induction sigma 222 um
  ion      induced charge Q(t) = (1-f) + f*min(t/T, 1)   (parallel plate)
  shaper   impulse (charge) response, three families, each scaled so the
           5 -> 100 % peaking time is 180 ns (DREAM manual Table 9, code 2)
  DAQ      60 ns samples (run_config sample_period 60; cfg WrClk_Div 6),
           32 samples, white noise 11 ADC, rail 3872, shape censored >= 3500
  metric   10-90 % rise and FWHM on the peak strip, by the SAME linear
           interpolation t14_compare.py uses (it is the metric under test,
           so it is reproduced, not reinvented)

Optional resistive spreading on the strip: a fraction c of the charge under a
strip leaks to its two neighbours with time constant tau_s.
"""
from __future__ import annotations

import argparse
import json
import math
import sys

import numpy as np

DT_FINE = 1.0            # ns
N_FINE = 2600            # ns of fine grid
DT_S = 60.0
N_S = 32
NOISE = 11.0
RAIL = 3872.0
SAT = 3500.0
SIGMA = 5.0
PITCH = 0.78             # mm
GAP_MM = 30.0


# ---------------------------------------------------------------- shapers
def _peak_5_100(h, dt=DT_FINE):
    i = int(np.argmax(h))
    lvl = 0.05 * h[i]
    j = np.flatnonzero(h[:i + 1] < lvl)[-1]
    t5 = j + (lvl - h[j]) / (h[j + 1] - h[j])
    return (i - t5) * dt


def _sk_manual(t, tau_r, w0, zeta=0.75):
    """real pole (PZC) x 2 complex poles (Sallen-Key), charge response."""
    wd = w0 * math.sqrt(1 - zeta ** 2)
    sk = (w0 / math.sqrt(1 - zeta ** 2)) * np.exp(-zeta * w0 * t) * np.sin(wd * t)
    ex = np.exp(-t / tau_r)
    return np.convolve(ex, sk)[:len(t)] * DT_FINE


def shaper(family, t_peak=180.0):
    t = np.arange(0, 3000, DT_FINE)

    def make(scale):
        if family == "crrc2":
            return (t / scale) ** 2 * np.exp(-t / scale)
        if family == "crrc4":
            return (t / scale) ** 4 * np.exp(-t / scale)
        if family.startswith("manual"):
            r = float(family.split("_")[1]) if "_" in family else 1.0
            return _sk_manual(t, tau_r=scale * r, w0=1.0 / scale)
        raise ValueError(family)

    lo, hi = 5.0, 300.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _peak_5_100(make(mid)) < t_peak:
            lo = mid
        else:
            hi = mid
    h = make(0.5 * (lo + hi))
    return h / h.max()


# ------------------------------------------------------------ the toy
def electron_kernel(h, f_ion, t_ion, c_spread, tau_s, psi_p=1.0):
    """Peak-strip output for one electron landing at t=0, plus the response it
    leaks onto EACH neighbour. Convolution of the shaper with dQ/dt."""
    t = np.arange(0, len(h)) * DT_FINE
    # ion climbs at constant speed; induced charge = psi(0) - psi(z) with
    # psi(z) = (1 - z/g)**psi_p  (psi_p = 1: parallel plate, a rectangle current)
    u = np.clip(t / t_ion, 0, 1)
    q = (1 - f_ion) + f_ion * (1 - (1 - u) ** psi_p) if f_ion > 0 else np.ones_like(t)
    if c_spread > 0:
        keep = (1 - c_spread) + c_spread * np.exp(-t / tau_s)
        leak = 0.5 * c_spread * (1 - np.exp(-t / tau_s))
    else:
        keep, leak = np.ones_like(t), np.zeros_like(t)
    out = []
    for frac in (keep, leak):
        Q = q * frac
        i = np.diff(Q, prepend=0.0)          # includes the delta at t=0
        out.append(np.convolve(i, h)[:N_FINE])
    return out


def simulate(p, n_events, rng, h):
    k_self, k_nb = electron_kernel(h, p["f_ion"], p["t_ion"], p["c_spread"], p["tau_s"], p.get("psi_p", 1.0))
    v = p["v"]                                 # mm/ns = um/ns * 1e-3
    tan = math.tan(math.radians(p["theta"]))
    n_half = 6 + int(abs(GAP_MM * tan) / PITCH)
    strips = np.arange(-n_half, n_half + 1)
    edges_lo, edges_hi = (strips - 0.5) * PITCH, (strips + 0.5) * PITCH
    rows = []
    t_off = 300.0                              # first arrival on the fine grid
    for _ in range(n_events):
        x0 = rng.uniform(-0.5, 0.5) * PITCH - 0.5 * GAP_MM * tan
        ncl = rng.poisson(91.2 / 3.83 * GAP_MM / 10)
        z = rng.uniform(0, GAP_MM, ncl)
        size = rng.geometric(1 / 3.83, ncl)
        z = np.repeat(z, size)
        # sigma[mm] = d[sqrt(cm)] * sqrt(z[cm]) * 10
        x = x0 + z * tan + rng.normal(0, 1, z.size) * 0.04766 * np.sqrt(z / 10) * 10
        t = t_off + z / v + rng.normal(0, 1, z.size) * 0.02425 * np.sqrt(z / 10) * 10 / v
        g = rng.gamma(1.5, 1 / 1.5, z.size)
        # prompt-induction lateral smear, then strip fractions
        sp = p["sigma_p"]
        cdf = lambda u: 0.5 * (1 + _erf(u / (sp * math.sqrt(2))))
        w = cdf(edges_hi[:, None] - x[None, :]) - cdf(edges_lo[:, None] - x[None, :])
        it = np.clip(np.round(t / DT_FINE).astype(int), 0, N_FINE - 1)
        dep = np.zeros((len(strips), N_FINE))
        for s in range(len(strips)):
            dep[s] = np.bincount(it, weights=g * w[s], minlength=N_FINE)[:N_FINE]
        F = np.fft.rfft(dep, n=2 * N_FINE)
        Ks, Kn = np.fft.rfft(k_self, 2 * N_FINE), np.fft.rfft(k_nb, 2 * N_FINE)
        wf_self = np.fft.irfft(F * Ks, 2 * N_FINE)[:, :N_FINE]
        wf = wf_self.copy()
        if p["c_spread"] > 0:
            nb = np.fft.irfft(F * Kn, 2 * N_FINE)[:, :N_FINE]
            wf[1:] += nb[:-1]
            wf[:-1] += nb[1:]
        rows.append(wf)
    return rows


def _erf(a):
    a = np.asarray(a, float)
    # Abramowitz-Stegun 7.1.26 is 1.5e-7; fine for strip fractions
    s = np.sign(a)
    a = np.abs(a)
    tt = 1 / (1 + 0.3275911 * a)
    y = 1 - (((((1.061405429 * tt - 1.453152027) * tt) + 1.421413741) * tt
              - 0.284496736) * tt + 0.254829592) * tt * np.exp(-a * a)
    return s * y


def sample_and_measure(fine_wfs, scale, rng):
    out = []
    for wf in fine_wfs:
        phase = rng.uniform(0, DT_S)
        ts = 180.0 + phase + DT_S * np.arange(N_S)       # pulse start ~ sample 2
        idx = np.clip((ts / DT_FINE).astype(int), 0, N_FINE - 1)
        w = wf[:, idx] * scale + rng.normal(0, NOISE, (wf.shape[0], N_S))
        w = np.minimum(w, RAIL)
        amp = w.max(axis=1)
        over = amp / NOISE >= SIGMA
        if not over.any():
            continue
        cand = np.flatnonzero(over)
        pk = int(cand[np.argmax(amp[cand])])
        out.append(measure(w[pk]))
    return out


def measure(w):
    """Identical to t14_compare.extract_view's rise / FWHM."""
    ipk = int(np.argmax(w))
    a = float(w[ipk])

    def _cross(frac, side):
        lvl = frac * a
        if side < 0:
            below = np.flatnonzero(w[:ipk + 1] < lvl)
            if len(below) == 0:
                return None
            i = below[-1]
            j = i + 1
        else:
            below = np.flatnonzero(w[ipk:] < lvl)
            if len(below) == 0:
                return None
            j = ipk + below[0]
            i = j - 1
        if w[j] == w[i]:
            return float(i)
        return i + (lvl - w[i]) / (w[j] - w[i])

    r10, r90 = _cross(0.10, -1), _cross(0.90, -1)
    hl, hr = _cross(0.50, -1), _cross(0.50, +1)
    sat = a >= SAT
    rise = (r90 - r10) * DT_S if (r10 is not None and r90 is not None) else np.nan
    fwhm = (hr - hl) * DT_S if (hl is not None and hr is not None) else np.nan
    return dict(peak_amp=a, saturated=sat,
                rise_ns=np.nan if sat else rise, fwhm_ns=np.nan if sat else fwhm)


Q = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90)


def run(p, n_events=1500, seed=1, target_peak=None):
    rng = np.random.default_rng(seed)
    h = shaper(p["shaper"])
    fine = simulate(p, n_events, rng, h)
    # gain: put the median peak-strip amplitude where the data leg's is
    pk = np.median([wf.max() for wf in fine])
    scale = (target_peak or 1000.0) / pk
    res = sample_and_measure(fine, scale, rng)
    r = np.array([x["rise_ns"] for x in res], float)
    fw = np.array([x["fwhm_ns"] for x in res], float)
    sat = np.mean([x["saturated"] for x in res])
    rr = r[~np.isnan(r)]
    return dict(params=p, n=len(res), sat=float(sat),
                rise_q=[float(np.quantile(rr, q)) for q in Q],
                fwhm_med=float(np.nanmedian(fw)),
                undefined=float(np.mean(np.isnan(r)) - sat))


BASE = dict(f_ion=0.908, t_ion=340.0, shaper="crrc2", v=0.03914, theta=0.0,
            sigma_p=0.222, c_spread=0.0, tau_s=230.0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", nargs="*", default=[], help="key=value overrides")
    ap.add_argument("-n", type=int, default=1500)
    ap.add_argument("--target-peak", type=float, default=1000.0)
    a = ap.parse_args()
    p = dict(BASE)
    for kv in a.set:
        k, v = kv.split("=")
        p[k] = v if k == "shaper" else float(v)
    print(json.dumps(run(p, a.n, target_peak=a.target_peak)))


def sampled_peak_strips(p, n_events=1500, seed=1, target_peak=1000.0):
    """Peak-strip waveform with its +-3 neighbours, sampled like the DAQ, for
    shape studies alongside the real legs (waveforms_pm3.npz layout)."""
    rng = np.random.default_rng(seed)
    h = shaper(p["shaper"])
    fine = simulate(p, n_events, rng, h)
    scale = target_peak / np.median([wf.max() for wf in fine])
    out = []
    for wf in fine:
        phase = rng.uniform(0, DT_S)
        ts = 180.0 + phase + DT_S * np.arange(N_S)
        idx = np.clip((ts / DT_FINE).astype(int), 0, N_FINE - 1)
        w = np.minimum(wf[:, idx] * scale + rng.normal(0, NOISE, (wf.shape[0], N_S)), RAIL)
        amp = w.max(axis=1)
        cand = np.flatnonzero(amp / NOISE >= SIGMA)
        if not len(cand):
            continue
        pk = int(cand[np.argmax(amp[cand])])
        blk = np.full((7, N_S), np.nan)
        for k in range(-3, 4):
            if 0 <= pk + k < w.shape[0]:
                blk[k + 3] = w[pk + k]
        out.append(blk)
    return np.array(out)
