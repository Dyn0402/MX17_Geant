#!/usr/bin/env python
"""
tau_g_reinterpretation.py — the resolution of the "A1 problem" (2026-08-07).

THE QUESTION
============
The overnight V4 test refuted plan assumption A1 (ESL strips grounded only at
the two y-ends) by comparing its predicted drain, L²/(π²D) = 4–41 ms, against
"the measured global drain constant τ_g = 5.3–7.3 µs". The hardware is not in
doubt (user 2026-08-07: copper bus strips at both y-ends, nothing in between),
so something had to give.

What gives is the label on τ_g. It comes from nTof_x17
`mx_june_wft/bench/rc_line_step1/2.py`, which fits

    T_Y = T_X ⊛ [δ + d/dt( Gaussian-line-diffusion(t_p) · exp(−t/τ_g) )]

on a ≤1.4 µs window with the measured X template as the electronics reference.
Two consequences:

  1. Anything common to both views cancels by construction — DREAM shaping,
     mesh HV sag, and crucially a REAL global drain of the sheet. The fit is
     blind to A1's drain no matter its value.
  2. An ms-scale drain is a factor 0.9998 across the window — invisible even
     in absolute terms.

So τ_g is the EXTRA decay of the Y view relative to X beyond the Gaussian toy
model: a kernel-shape observable, not a grounding measurement.

THE TEST THIS SCRIPT RUNS
=========================
Compute the true induced-charge evolution G(t) from the S1 solver with
tau_drain_s=None (strict charge conservation — the honest A1 limit) for the
REAL comb channels (checkerboard, 256 pads on 1.56 mm — channel_map.py), then
fit the bench's own toy model to it. If the drain-free kernel produces an
apparent τ_g of the measured order, the "refutation" was a category error.

RESULT (first run, 2026-08-07, laptop, ~10 min)
===============================================
rho_s = 1 MΩ/sq, d_k = 75 µm (D = 2.12 m²/s), NO drain:
  Y comb, charge on-strip (3 phases):  apparent τ_g = 2.46 / 3.02 / 3.31 µs
  Y comb, charge in a 250 µm gap:      apparent τ_g = ∞ (kept 0.70 at 1.4 µs
                                       deep in the gap; 0.22 near its edge)
  X comb (column over a strip):        drops to ~0.36 by 300 ns, then FLAT
                                       ("X drain-free everywhere", as measured)
rho_s = 2 MΩ/sq: on-strip apparent τ_g = 2.22 / 2.70 / 2.96 µs.

Reading: a drain-free detector manufactures a few-µs apparent "drain" in this
observable purely from non-Gaussian sheet relaxation (C(k) dispersion) plus
comb sampling. On-strip alone sits ~×2 below the measured 5.3–7.3 µs; the
measured template averages strip and gap deposits (~31% of the area lands in
gaps and barely decays), and a single-exponential fit to that mixture lands
slower — in the measured band. The mechanism is y-independent (position-flat
along the strip) and geometric (fleet-consistent), both as observed
(nTof_x17 RECO_BENCH_2026-07-29.md: 8.6/8.5/6.5/8.0 µs along the strip;
7.3/5.8/5.4/5.3 µs across det3/det2/det6/det7).

Final closure at waveform level (through the digitizer and the unmodified
rc_line fit) is plan task T13b. The two-component-tail refit of the measured
templates is task Td.

Run:  ~/PycharmProjects/nTof_x17/.venv/bin/python -m response.validation.tau_g_reinterpretation
"""
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import erf

from ..common import constants as C
from ..solver.wpot import WeightingSolver

PAD_P = C.PAD_PITCH_M          # 0.78 mm
PAD_S = C.PAD_SIZE_UM * 1e-6   # 0.68 mm
COMB_P = 2 * PAD_P             # 1.56 mm — the checkerboard comb pitch

# prompt + log-spaced to 2 µs; the fit window matches the bench's 1.4 µs
TIMES = np.concatenate([[0.0], np.logspace(np.log10(5e-9), np.log10(2e-6), 16)])
FIT_MAX = 1.4e-6


def comb_row_pattern(s, y0_m, x_phase=0.5 * PAD_P):
    """Y channel: pads at row y0, every other column (comb along x)."""
    X, Y = np.meshgrid(s.x, s.y, indexing="xy")
    v = np.zeros((s.ny, s.nx))
    for j in range(int(round(s.lx / COMB_P))):
        xc = x_phase + j * COMB_P
        dx = np.abs(((X - xc + s.lx / 2) % s.lx) - s.lx / 2)
        v[(dx <= PAD_S / 2) & (np.abs(Y - y0_m) <= PAD_S / 2)] = 1.0
    return v


def comb_col_pattern(s, xc_m, y_phase=0.0):
    """X channel: pads at column xc, every other row (comb along y)."""
    X, Y = np.meshgrid(s.x, s.y, indexing="xy")
    v = np.zeros((s.ny, s.nx))
    dx = np.abs(((X - xc_m + s.lx / 2) % s.lx) - s.lx / 2)
    yc = (Y - y_phase) / COMB_P
    dyc = np.abs(yc - np.round(yc)) * COMB_P
    v[(dx <= PAD_S / 2) & (dyc <= PAD_S / 2)] = 1.0
    return v


def toy_q0(t, t_p, tau_g, pitch=PAD_P, n_src=9):
    """rc_line_step1.charge_fractions, D=0 term, in SI units."""
    sig = pitch * np.sqrt(np.maximum(t, 1e-12) / t_p)
    y0 = (np.arange(n_src) + 0.5) / n_src * pitch - pitch / 2
    lo = -0.5 * pitch - y0[:, None]
    hi = 0.5 * pitch - y0[:, None]
    z = 1.0 / (np.sqrt(2.0) * sig)[None, :]
    q = 0.5 * (erf(hi * z) - erf(lo * z))
    return q.mean(axis=0) * np.exp(-t / tau_g)


def fit_toy(t, g):
    """Fit (t_p, tau_g) of the toy model to a normalized decay curve g(t)."""
    m = t <= FIT_MAX
    tm, gm = t[m], g[m]

    def loss(v):
        return float(((toy_q0(tm, np.exp(v[0]), np.exp(v[1])) - gm) ** 2).sum())

    best = None
    for tp0 in (50e-9, 150e-9, 500e-9):
        for tg0 in (1e-6, 5e-6, 50e-6):
            r = minimize(loss, np.log([tp0, tg0]), method="Nelder-Mead",
                         options=dict(xatol=1e-4, fatol=1e-12, maxiter=2000))
            if best is None or r.fun < best.fun:
                best = r
    return tuple(np.exp(best.x))


def run(rho_s, d_k=75e-6):
    t0 = time.time()
    s = WeightingSolver(rho_s_ohm_sq=rho_s, d_kapton_m=d_k,
                        nx=3120, ny=1024, ly_m=0.0512, tau_drain_s=None)
    iy0 = s.ny // 2
    assert abs(s.y[iy0]) < 1e-9
    print(f"\n=== rho_s = {rho_s/1e6:.1f} MΩ/sq, d_k = {d_k*1e6:.0f} µm "
          f"(D = {1.0/(rho_s*C.c_prime(s.gap, d_k, s.eps_r)):.2f} m²/s), "
          f"NO drain ===")

    # Y channel: comb row at y=0. Charge positions walk the strip phase
    # (pad centres at (2j+0.5)·0.78 mm advance 20 µm/pad against the 800 µm
    # ESL period; x mod 0.8 mm < 0.55 mm = on-strip).
    psi = s.solve(comb_row_pattern(s, 0.0), TIMES)
    picks = [(j, (0.5 + 2 * j) * PAD_P) for j in range(20)]
    phases = [(j, xc, np.mod(xc, C.ESL_PITCH_M)) for j, xc in picks]
    sel = ([p for p in phases if p[2] < 0.45e-3][:3]
           + [p for p in phases if p[2] >= C.ESL_WIDTH_M][:2])
    print("  Y comb, G(t)/G(0+) for a charge above a driven pad:")
    for j, xc, ph in sel:
        ix = int(np.argmin(np.abs(s.x - xc)))
        g = psi[:, iy0, ix]
        gn = g / g[0]
        t_p, tau_g = fit_toy(TIMES, gn)
        tag = "strip" if ph < C.ESL_WIDTH_M else "GAP  "
        print(f"    x={xc*1e3:6.2f} mm  phase={ph*1e6:3.0f} µm [{tag}]  "
              f"kept@0.3/1.4µs: {np.interp(0.3e-6, TIMES, gn):.3f}/"
              f"{np.interp(1.4e-6, TIMES, gn):.3f}  -> apparent "
              f"t_p={t_p*1e9:5.0f} ns  tau_g="
              + (f"{tau_g*1e6:5.2f} µs" if tau_g < 1e-2 else "  ~inf"))

    # X channel: comb column over a strip — must be drain-free after the
    # early sharing transient, as the bench observes.
    xc = 0.5 * PAD_P
    psix = s.solve(comb_col_pattern(s, xc), TIMES)
    ix = int(np.argmin(np.abs(s.x - xc)))
    gx = psix[:, iy0, ix]
    gxn = gx / gx[0]
    print("  X comb (column over strip), kept at "
          + ", ".join(f"{tt*1e9:.0f} ns: {np.interp(tt, TIMES, gxn):.3f}"
                      for tt in (0.1e-6, 0.3e-6, 0.7e-6, 1.4e-6)))
    print(f"  [{time.time()-t0:.0f} s]")


if __name__ == "__main__":
    for rho in (1e6, 2e6):
        run(rho)
