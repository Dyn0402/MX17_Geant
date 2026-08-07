#!/usr/bin/env python3
"""
ions.py — the analytic ion component of the induced signal (plan §7 step 5).

WHY ANALYTIC RATHER THAN FROM S3. The plan's T9 line asks for "templates from
S1 surface G + analytic ion term, first pass", and that is now also forced: the
S3 campaign's `i_elec` and `i_ion` arrays came back identically zero in every
slice (verified over all 56), so `aval_calib.json` carries no measured ion
shape. The gain, Polya theta and sigma0 from that campaign are unaffected and
still used.

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


if __name__ == "__main__":
    import pprint
    pprint.pprint(describe(150e-6, 490.0))
