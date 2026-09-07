#!/usr/bin/env python3
"""ADVERSARIAL RE-DERIVATION of the resistive-sheet screening sizing.

WHY THIS EXISTS
===============
The sheet-screening hypothesis was the sharpest structural candidate for the
f_eff contradiction: the ion never lands on the sheet, so its readout signal is
pure induction from a charge in the gas, yet `apply_longitudinal` couples it
through the SAME surface kernel as charge physically deposited on the sheet. A
resistive sheet is a frequency-dependent screen, so the deposited-charge Green's
function is arguably the wrong operator for an induced source.

It was retired by `nTof_x17 mx17_sim_wft/sheet_screening_sizing.py` at 8 %
suppression where x4.5 was needed, with a rho ladder of 0.81 / 0.86 / 0.92 /
1.02 across rho_s = 0.5 / 1 / 2 / 5 MOhm/sq. That retirement is load-bearing:
it is the last eliminated candidate on the central open question, and it rests
on the assumption that `apply_longitudinal` already carries most of the sheet
dynamics.

This project's record -- Fix 1, Fix 7, A1, the toy that over-predicted amplitude
x6, tonight's 200 ns threshold artifact and tonight's own withdrawn factor-0.74
-- is that a second derivation of a load-bearing sizing is never wasted.

PRE-REGISTERED PREDICTION AND FALSIFIER (written before running)
================================================================
Prediction: I expect to roughly reproduce their ladder, because the underlying
scale argument is sound -- tau_k = L^2 rho_s c' at the CHANNEL scale is
~150-600 ns at the production point against a ~350 ns ion, i.e. the transition
regime, which gives a partial and modest suppression.

  * Central: ratio at rho_s = 2 MOhm/sq in 0.85-0.99.
  * FALSIFIER: if my ratio at 2 MOhm/sq differs from their 0.92 by more than
    +-0.10, one of the two derivations is wrong and I must identify WHICH
    assumption diverges before believing either. I do not get to prefer mine.
  * Second falsifier: my method must reproduce f_ion = 0.9056 in the static
    limit (rho_s -> 0, where the sheet is a perfect conductor and the weighting
    potential is the parallel-plate ramp). If it does not, my derivation is
    broken and its disagreement means nothing.

A DELIBERATELY DIFFERENT METHOD
===============================
Theirs: track a sheet counter-charge c(k,t) relaxing toward the ion's image
exp(-k z(t)), and take the net as img - c.

Mine: Riegler superposition on the time-dependent weighting potential, which is
the same object the S1 kernel already is, so it plugs into the existing
machinery instead of introducing a parallel picture.

  * Psi_sheet(k, tau) = (S(k)/C(k)) exp(-k^2 tau / (rho_s C(k)))
    -- the uniform-sheet closed form the repo already uses as its V1 reference
    (`wpot.solve_uniform_analytic`), i.e. the response to unit charge sitting ON
    the sheet since tau ago. This IS the model's operator.

  * Upward continuation into the gas. The weighting problem has the mesh
    grounded at z = GAP and the sheet plane at Psi_sheet, with no sources in
    between, so the Laplace solution is exact:

        Psi(k, z, tau) = Psi_sheet(k, tau) * sinh(k(GAP - z)) / sinh(k GAP)

    At k -> 0 this is the parallel-plate ramp (GAP - z)/GAP; at large k it tends
    to exp(-k z), which is precisely their image factor. So their exp(-kz) is
    the large-k limit of my exact continuation -- the two methods should agree
    where that limit is good and may diverge at small k, which is exactly the
    regime that decides screening.

  * A charge that appears at height z at time t and stays contributes
    q Psi(z, T - t). A MOVING charge is a sequence of dipoles: at each step,
    +q appears at z_{i+1} and -q at z_i, both switched on at t_{i+1}. So

        Q_true(T) = q Psi(z0, T)
                    + q sum_i [ Psi(z_{i+1}, T - t_{i+1})
                              - Psi(z_i,     T - t_{i+1}) ]

    which is manifestly correct by superposition and needs no appeal to a
    signal theorem.

  * The MODEL's null, in the same language: the ion's charge is injected onto
    the SHEET at a uniform rate and then spreads, i.e. Psi evaluated at z = 0:

        Q_model(T) = q sum_i (1/N) Psi(0, T - t_i)

Both branches use one operator, differing only in the height at which it is
evaluated. That makes the induced-vs-injected distinction the ONLY thing being
compared, which is the cleanest possible statement of the question.
"""
from __future__ import annotations

import json

import numpy as np

from response.common import constants as C
from response.solver.wpot import stack_coeffs

GAP = C.AMP_GAP_M                    # 150 um
Z0 = 15.0e-6                         # ion birth height above the ESL (S3)
T_ION = 306e-9                       # analytic rectangle transit
T_EVAL = 200e-9                      # rising-edge window, as in their sizing
SIG0 = 200e-6                        # transverse charge-cloud sigma
N_STEP = 600

RHOS = (0.5e6, 1.0e6, 2.0e6, 5.0e6)


def psi_sheet(k, tau, rho_s):
    """Response to unit charge sitting ON the sheet since tau ago, per k."""
    Cv, Sv = stack_coeffs(k, GAP, C.KAPTON_THICK_UM * 1e-6, C.KAPTON_EPS_R,
                          C.GLUE_THICK_UM * 1e-6, C.GLUE_EPS_R)
    return (Sv / Cv) * np.exp(-k ** 2 * tau / (rho_s * Cv))


def cont(k, z):
    """Exact upward continuation sinh(k(GAP-z))/sinh(k GAP), k->0 safe."""
    k = np.asarray(k, dtype=float)
    out = np.empty_like(k)
    small = k * GAP < 1e-6
    out[small] = (GAP - z) / GAP
    kk = k[~small]
    # written as a ratio of exponentials to stay finite at large k*GAP
    a, b = kk * (GAP - z), kk * GAP
    out[~small] = np.exp(a - b) * (1 - np.exp(-2 * a)) / (1 - np.exp(-2 * b))
    return out


def make_grid(n=6000, kmax=8e4):
    k = np.linspace(1.0, kmax, n)
    return k


def channel_share(k, F, pitch):
    """Central-channel share: project the k-space profile onto a pad aperture."""
    W = np.sinc(k * pitch / 2.0 / np.pi)
    return float(np.trapz(W * F, k) * pitch / np.pi)


def run(rho_s, pitch, k=None):
    if k is None:
        k = make_grid()
    F0 = np.exp(-k ** 2 * SIG0 ** 2 / 2.0)
    ts = np.linspace(0.0, T_EVAL, N_STEP + 1)
    zs = np.minimum(Z0 + (GAP - Z0) * ts / T_ION, GAP)

    # --- TRUE: induced by a charge moving in the gas -------------------------
    # NOTE (bug found by falsifier 1 on the first run, 2026-08-10): the ion's
    # "appearance" term q Psi(z0, T) must NOT be included here. The ion is not
    # created from nothing -- it appears with its electron, and that prompt term
    # is the ELECTRON's, already accounted for as f_e. Including it made the
    # true branch the total electron+ion signal while the model branch was the
    # ion's arriving charge alone, so the two were not the same observable and
    # their ratio was meaningless (it came out ~0.07, and negative at 0.5M).
    # The ion's own contribution is the dipole sum alone. Psi falls with z, so
    # it is written start-minus-end to come out positive.
    acc = np.zeros_like(k)
    for i in range(N_STEP):
        tau = T_EVAL - ts[i + 1]
        ps = psi_sheet(k, tau, rho_s) * F0
        acc = acc + ps * (cont(k, zs[i]) - cont(k, zs[i + 1]))
    q_true = channel_share(k, acc, pitch)

    # --- MODEL: same operator, charge injected ON the sheet ------------------
    accm = np.zeros_like(k)
    for i in range(N_STEP):
        tau = T_EVAL - ts[i + 1]
        accm = accm + (1.0 / N_STEP) * psi_sheet(k, tau, rho_s) * F0
    q_model = channel_share(k, accm, pitch)

    # arrived-charge fraction, common to both, cancels in the ratio
    q_frac = T_EVAL / T_ION
    return q_true / q_frac, q_model / q_frac, q_true / q_model


def static_check(pitch):
    """Falsifier 2: the static limit must reproduce f_ion = (GAP-Z0)/GAP."""
    k = make_grid()
    F0 = np.exp(-k ** 2 * SIG0 ** 2 / 2.0)
    # perfect conductor: sheet fully relaxed at every tau -> psi is tau-free
    Cv, Sv = stack_coeffs(k, GAP, C.KAPTON_THICK_UM * 1e-6, C.KAPTON_EPS_R,
                          C.GLUE_THICK_UM * 1e-6, C.GLUE_EPS_R)
    ps = (Sv / Cv) * F0
    q_start = channel_share(k, ps * cont(k, Z0), pitch)
    q_end = channel_share(k, ps * cont(k, GAP), pitch)
    q_ref = channel_share(k, ps * cont(k, 0.0), pitch)
    return (q_start - q_end) / q_ref


if __name__ == "__main__":
    THEIRS = {0.5e6: 0.81, 1.0e6: 0.86, 2.0e6: 0.92, 5.0e6: 1.02}
    out = {}
    for pitch, lab in ((C.PAD_PITCH_M, "pad pitch 780 um"),
                       (C.ESL_PITCH_M, "ESL pitch 800 um (theirs)")):
        f_ion_static = static_check(pitch)
        print(f"\n===== aperture: {lab} =====")
        print(f"  falsifier 2 -- static-limit f_ion = {f_ion_static:.4f} "
              f"(expect {(GAP - Z0) / GAP:.4f})")
        print(f"\n  {'rho_s':>7} {'true/Q':>9} {'model/Q':>9} "
              f"{'ratio':>8} {'theirs':>8} {'diff':>8}")
        row = {}
        for rho in RHOS:
            t, m, r = run(rho, pitch)
            d = r - THEIRS[rho]
            print(f"  {rho/1e6:6.1f}M {t:9.4f} {m:9.4f} {r:8.3f} "
                  f"{THEIRS[rho]:8.2f} {d:+8.3f}")
            row[f"{rho/1e6:g}M"] = dict(true=t, model=m, ratio=r,
                                        theirs=THEIRS[rho], diff=d)
        out[lab] = dict(f_ion_static=f_ion_static, ladder=row)

    print("\n=== VERDICT AGAINST THE PRE-REGISTERED FALSIFIERS ===")
    key = "pad pitch 780 um"
    fs = out[key]["f_ion_static"]
    ok2 = abs(fs - (GAP - Z0) / GAP) < 0.01
    print(f"  falsifier 2 (static f_ion): {fs:.4f} vs "
          f"{(GAP - Z0) / GAP:.4f} -> {'PASS' if ok2 else 'FAIL'}")
    r2 = out[key]["ladder"]["2M"]["ratio"]
    ok1 = abs(r2 - 0.92) <= 0.10
    print(f"  falsifier 1 (2 MOhm/sq within +-0.10 of 0.92): {r2:.3f} "
          f"-> {'PASS, their sizing CONFIRMED' if ok1 else 'FAIL, DIVERGENCE'}")
    need = 4.5
    print(f"\n  suppression needed to explain f_eff: x{need:.1f} "
          f"(ratio ~{1/need:.2f}); measured here {r2:.3f}")

    with open("sheet_screening_rederive.json", "w") as f:
        json.dump(out, f, indent=1)
    print("\nwrote sheet_screening_rederive.json")
