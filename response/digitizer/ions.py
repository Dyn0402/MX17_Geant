#!/usr/bin/env python3
"""
ions.py — the analytic ion component of the induced signal (plan §7 step 5).

ANALYTIC *AND* MEASURED. The analytic model below was written when the S3
campaign's `i_elec` / `i_ion` arrays came back identically zero in every slice
(verified over all 56) — a missing `SetWeightingPotential` on the
ComponentConstant, since Garfield's `m_useWeightingPotential` defaults true.
That is fixed and the v2 campaign (all 56 slices, 2026-08-07) carries real
currents, so `measured_longitudinal()` below supersedes the analytic form.

The analytic model is KEPT, not deleted, because it is the cross-check: it
predicts the shape from first principles and the measurement then says how
well. The verdict at 490 V is that the *form* was right and two numbers were
wrong:

                            analytic     measured (S3 v2)
    f_ion                     0.967          0.908
    transit                306 ns rect     ~340 ns, rectangle confirmed
    mean ion birth height     5 um          13.84 um

The transit sits inside the +-30 % mobility band on 306 ns, so the mobility
assumption survives. The f_ion error does not come from mobility at all: it
comes from `z_aval_um`. This file argued the avalanche sits "within microns of
the anode", but the ionization profile is exponential with 1/alpha = 14 um, so
the MEAN ion is born at 14 um, not 5. Two independent routes through the v2
data agree on that to four digits:

    from the current integrals   Qi/(Qe+Qi)                    = 0.9079
    from the alpha_z histogram   1 - <z>/g, <z> = 13.84 um     = 0.9077

and their agreement is also the proof that no ion charge was lost off the end
of the 400 ns S3 signal window — a truncated tail would have pulled the first
number below the second. The histogram's own decay length, 14.8 um, matches
1/alpha = g/ln(G) = 14.1 um from the measured gain, so the profile is internally
consistent too.

WHY THE ION TERM IS NOT A SMALL CORRECTION. In a Micromegas the avalanche
develops toward the anode, so the ions are created within a few microns of the
ESL and then drift the WHOLE gap back to the mesh. In the parallel-plate
weighting field psi = 1 - z/g, an electron created at height z contributes only
z/g of the charge while the ion contributes the remaining (g-z)/g ~ 1. So the
ions carry essentially the entire induced charge, spread over their transit
time, while the electrons deliver a small prompt spike.

    t_ion = g^2 / (mu_i * V)

At the bench point (g = 150 um, V = 490 V, mu ~ 1.5 cm^2/V/s reduced) that is
~300 ns, comparable to the DREAM 180 ns peaking time — so it shapes the pulse
rather than merely adding a tail.

WHAT THIS MODEL DOES AND DOES NOT CAPTURE.

  * It uses the parallel-plate longitudinal weighting psi(z) = 1 - z/g, which
    is exact for the SUM over all readout channels and is where the g^2/(mu V)
    transit time comes from.
  * It does NOT re-solve the lateral kernel at each height. The true statement
    (plan §7 step 6) is that an ion at height z induces through Psi_n(x, y, z),
    which broadens as the ion climbs toward the mesh, so the ion signal is
    laterally WIDER than the surface kernel. S1 solves only the z=0 plane, so
    that broadening is not available yet; here the ion is given the surface
    kernel's lateral shape, frozen at its creation point. That makes this an
    UNDER-estimate of ion-induced lateral sharing, and it is exactly what the
    T10 slow path (Garfield ComponentGrid over the S1 time slices) exists to
    replace.

Ion mobility is flagged in plan §5 as "the single softest parameter", carried
at +-30 %. It is an explicit argument here and is recorded per run.
"""

from __future__ import annotations

import numpy as np

# Reduced mobility of Ar+ in argon at STP [cm^2 / (V s)]. Ar/iC4H10 95/5 is
# argon-dominated and charge transfer moves the charge onto the lowest-IP
# species, so the drifting ion is mostly isobutane+ in argon; the mobilities
# are close at this level and the plan carries +-30 % regardless.
MU_ION_CM2_VS = 1.5
MU_ION_SYS = 0.30


def ion_transit_ns(gap_m, v_mesh, mu_cm2_vs=MU_ION_CM2_VS):
    """
    Time for an ion to cross the amplification gap [ns].

    v_ion = mu E = mu V / g, so t = g / v_ion = g^2 / (mu V).
    """
    g_cm = gap_m * 1e2
    return g_cm ** 2 / (mu_cm2_vs * v_mesh) * 1e9


def ion_current_shape(t_s, gap_m, v_mesh, mu_cm2_vs=MU_ION_CM2_VS,
                      z0_m=0.0):
    """
    Normalised induced-current shape of the ion component, on the time grid t_s.

    An ion released at height z0 drifts at constant speed to the mesh at z = g.
    With psi = 1 - z/g the induced current is CONSTANT while it moves,

        i(t) = q * v_ion / g   for   0 < t < (g - z0)/v_ion,

    i.e. a rectangle, not an exponential. That matters: a decaying tail would
    put most of the ion charge early, and the real thing delivers it flat over
    the full transit. Returned normalised so that sum(i)*dt = 1.
    """
    t_tr = ion_transit_ns(gap_m, v_mesh, mu_cm2_vs) * 1e-9
    t_tr *= (1.0 - z0_m / gap_m)
    if t_tr <= 0:
        out = np.zeros_like(t_s)
        out[0] = 1.0
        return out
    dt = t_s[1] - t_s[0]
    shape = np.where(t_s < t_tr, 1.0, 0.0)
    # Partial last bin, so the total is exact rather than quantised to dt.
    k = int(np.floor(t_tr / dt))
    if 0 <= k < len(shape):
        shape[k] = (t_tr - k * dt) / dt
    s = shape.sum() * dt
    return shape / s if s > 0 else shape


def charge_split(gap_m, z_aval_m):
    """
    Fraction of the induced charge carried by the electron and by the ion.

    With psi = 1 - z/g, a pair created at height z gives the electron (which
    falls to z=0) a weight z/g and the ion (which climbs to z=g) a weight
    1 - z/g. In a Micromegas the avalanche is concentrated within microns of
    the anode, so z/g << 1 and the ions dominate — the opposite of the naive
    "the fast electron signal is the signal" intuition.
    """
    f_e = np.clip(z_aval_m / gap_m, 0.0, 1.0)
    return f_e, 1.0 - f_e


def describe(gap_m, v_mesh, mu_cm2_vs=MU_ION_CM2_VS, z_aval_m=5e-6):
    f_e, f_i = charge_split(gap_m, z_aval_m)
    t = ion_transit_ns(gap_m, v_mesh, mu_cm2_vs)
    return {
        "mu_ion_cm2_Vs": mu_cm2_vs,
        "mu_ion_systematic": MU_ION_SYS,
        "gap_um": gap_m * 1e6, "V_mesh": v_mesh,
        "t_ion_transit_ns": t,
        "t_ion_transit_ns_lo": ion_transit_ns(gap_m, v_mesh,
                                              mu_cm2_vs * (1 + MU_ION_SYS)),
        "t_ion_transit_ns_hi": ion_transit_ns(gap_m, v_mesh,
                                              mu_cm2_vs * (1 - MU_ION_SYS)),
        "z_avalanche_um": z_aval_m * 1e6,
        "f_electron": float(f_e), "f_ion": float(f_i),
        "lateral_model": "surface kernel frozen at creation point "
                         "(UNDER-estimates ion lateral sharing; T10 replaces)",
    }


def measured_longitudinal(pt):
    """
    The MEASURED longitudinal current profile h(t) of one unit of avalanche
    charge, from an S3 v2 calibration point. Returns (h, dt_s, info).

    h is i_elec + i_ion normalised to unit area, so it carries the electron
    prompt spike, the ion rectangle, the charge split between them, and the
    rounding of the rectangle's trailing edge — all measured rather than
    assumed. The trailing edge is rounded because ions are NOT all born at one
    height: those born high finish early, so the current decays over the last
    ~50 ns instead of cutting off. The single-z analytic rectangle cannot
    represent that.

    Sign convention: Garfield returns induced current negative for this
    geometry. h is returned positive-normalised, since the LUT convolution
    only cares about the shape and the digitizer carries its own sign.
    """
    dt = pt["signal_dt_ns"] * 1e-9
    ie = np.asarray(pt["i_elec"], dtype=np.float64)
    ii = np.asarray(pt["i_ion"], dtype=np.float64)
    tot = ie + ii
    area = tot.sum() * dt
    h = tot / area

    qe, qi = ie.sum() * dt, ii.sum() * dt
    c = np.cumsum(ii) / ii.sum()
    t = (np.arange(len(ii)) + 1) * pt["signal_dt_ns"]

    zh = pt.get("alpha_z_hist")
    z_mean_um = None
    if zh:
        cnt = np.asarray(zh["counts"], dtype=np.float64)
        ed = np.asarray(zh["edges"], dtype=np.float64)
        z_mean_um = float((cnt * 0.5 * (ed[:-1] + ed[1:])).sum() / cnt.sum())

    return h, dt, {
        "source": "S3 v2 measured (i_elec + i_ion)",
        "f_ion": float(qi / (qe + qi)),
        "f_electron": float(qe / (qe + qi)),
        "t_ion_50pct_ns": float(np.interp(0.50, c, t)),
        "t_ion_90pct_ns": float(np.interp(0.90, c, t)),
        "t_ion_99pct_ns": float(np.interp(0.99, c, t)),
        "z_birth_mean_um": z_mean_um,
        # Independent cross-check of f_ion from the position histogram; if this
        # disagrees with f_ion the ion tail was clipped by the S3 window.
        "f_ion_from_z_hist": (None if z_mean_um is None
                              else float(1.0 - z_mean_um / 150.0)),
        "n_samples": int(len(h)),
        "span_ns": float(len(h) * pt["signal_dt_ns"]),
        "lateral_model": "surface kernel frozen at creation point "
                         "(UNDER-estimates ion lateral sharing; T10 replaces)",
    }


if __name__ == "__main__":
    import pprint
    pprint.pprint(describe(150e-6, 490.0))
