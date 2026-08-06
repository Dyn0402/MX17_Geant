#!/usr/bin/env python3
"""
validate.py — the S1 acceptance tests (design/RESPONSE_SIM_PLAN.md §3).

The plan is explicit that these are the definition of done for T3-T5, not a
nice-to-have: "implement with the validation tests as the definition of done".
Nothing downstream may consume S1 output until these pass.

    V1  uniform sheet -> closed-form Riegler/Gaussian limit          (<1 %)
    V1b small-k limit -> sigma^2(t) = 2t/(rho_s c') exactly
    V2  single isolated strip -> 1D telegraph equation
    V3  charge conservation and the two asymptotic limits

    python3 -m response.solver.validate [--quick]
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from ..common import constants as C
from . import wpot as W
from .wpot import WeightingSolver, stack_coeffs, pad_pattern, strip_pattern_y


def _report(name, value, tol, unit="", extra=""):
    ok = value <= tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<52s} "
          f"{value:9.3e}{unit} (tol {tol:.1e}) {extra}")
    return ok


# ── V1 ───────────────────────────────────────────────────────────────────────

def v1_uniform_sheet(nx=780, ny=256, quick=False):
    """
    Hand the full Bloch solver a UNIFORM conductivity and require that it
    reproduces the closed form it should collapse to.

    This is the strongest single test in the suite: the Bloch block assembly,
    the C^-1/2 similarity transform, the eigendecomposition and the propagator
    are all exercised, while the answer is known analytically. Any sign error,
    any mis-indexed Fourier coefficient, any C/S mix-up shows up here.
    """
    print("\nV1  uniform sheet vs closed form")
    s = WeightingSolver(rho_s_ohm_sq=1.0e6,
                        d_kapton_m=75e-6, nx=nx, ny=ny,
                        uniform_sheet=True)
    vpad = pad_pattern(s, x0_m=s.lx / 2, y0_m=0.0)
    times = np.array([0.0, 1e-9, 10e-9, 100e-9, 1e-6])

    t0 = time.time()
    num = s.solve(vpad, times)
    ana = s.solve_uniform_analytic(vpad, times)
    dt = time.time() - t0

    scale = np.abs(ana).max()
    worst = 0.0
    for it, t in enumerate(times):
        err = np.abs(num[it] - ana[it]).max() / scale
        worst = max(worst, err)
        print(f"       t = {t*1e9:8.1f} ns   max|num-ana|/max|ana| = {err:.3e}")
    print(f"       ({dt:.1f} s for {len(times)} times on a "
          f"{ny}x{nx} grid)")
    return _report("V1 uniform sheet vs closed form", worst, 1e-2)


def v1b_gaussian_limit():
    """
    The k -> 0 limit must be an exact Gaussian diffusion kernel with
    sigma^2(t) = 2t/(rho_s c'), which is the form the plan quotes and the
    number (D ~ 2 m^2/s at 1 MΩ/sq) it sanity-checks against.

    Checked in Fourier space, where the statement is exact and grid-free:
    the decay rate of mode k must approach k^2/(rho_s c').
    """
    print("\nV1b small-k limit -> Gaussian diffusion")
    rho_s, d_k = 1.0e6, 75e-6
    cp = C.c_prime(kapton_m=d_k)
    k = np.logspace(0, 4, 40)                       # 1 .. 1e4 rad/m
    Cv, _ = stack_coeffs(k, C.AMP_GAP_M, d_k)
    rate_exact = k ** 2 / (rho_s * Cv)
    rate_limit = k ** 2 / (rho_s * cp)
    # The limit is only claimed for k d << 1 and k g << 1.
    small = k < 0.05 / C.AMP_GAP_M
    err = np.abs(rate_exact[small] / rate_limit[small] - 1).max()
    D = C.sheet_diffusivity(rho_s, kapton_m=d_k)
    print(f"       c' = {cp:.4e} F/m^2   D = 1/(rho_s c') = {D:.3f} m^2/s")
    print(f"       plan §3 quotes c' ~ 5e-7 F/m^2 and D ~ 2.0 m^2/s")
    # A 130 µm feature should relax in ~8 ns and one 800 µm pitch in ~0.3 µs
    for feat_um, want_ns in ((130, 8), (800, 300)):
        tau = (feat_um * 1e-6) ** 2 / (2 * D)
        print(f"       relaxation of a {feat_um:4d} µm feature: "
              f"{tau*1e9:7.1f} ns   (plan says ~{want_ns} ns)")
    return _report("V1b decay rate -> k^2/(rho_s c')", err, 1e-2)


# ── V2 ───────────────────────────────────────────────────────────────────────

def v2_single_strip(nx=3120, ny=64):
    """
    Strip pattern driven uniformly in y, so only k_y = 0 is excited. This is
    the 1D telegraph case (Galan et al., arXiv:1110.6640) and it isolates the
    in-plane transport along x.

    What must be true WITHOUT any drain term:

      * inside a strip, structure relaxes diffusively at D = 1/(rho_s c');
      * the strip's total charge is CONSERVED — with k_y = 0 there is no
        gradient along y, and the 250 µm gaps isolate the strip in x, so the
        charge has nowhere to go. Any drain must come from the y-ends, i.e.
        from assumption A1 (tested in V4), not from this equation.

    The first version of this test asserted that the strips drain to ground
    here. They do not, and should not: that assertion was wrong, not the
    solver. Charge conservation is the correct statement, and it is a
    sharper test.
    """
    print("\nV2  isolated strip, k_y = 0, telegraph behaviour")
    s = WeightingSolver(rho_s_ohm_sq=1.0e6, d_kapton_m=75e-6,
                        nx=nx, ny=ny, ly_m=0.02)
    vpad = strip_pattern_y(s, x0_m=s.lx / 2)
    times = np.array([0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
    psi = s.solve(vpad, times)

    prof = psi[:, ny // 2, :]                       # (nt, nx) along x
    on_strip = s.sigma_x > 0

    # Sheet charge per mode is C(k) V(k); at k_y = 0 the total sheet charge is
    # the k=0 component, which the in-plane operator cannot change.
    q_tot = [np.fft.fft(prof[it]).real[0] for it in range(len(times))]
    # Structure WITHIN each strip, measured strip by strip.
    #
    # The null space of grad.(sigma_s grad) at k_y = 0 is "V constant on each
    # strip, anything in the gaps" — so different strips keep different
    # asymptotic constants and a std taken over all strip samples at once does
    # NOT go to zero. It must be taken per strip and then combined, which is
    # the statement the telegraph equation actually makes: charge spreads to
    # uniformity along the conductor it is on.
    edges = np.flatnonzero(np.diff(on_strip.astype(int)) != 0) + 1
    runs = np.split(np.arange(s.nx), edges)
    strip_runs = [r for r in runs if on_strip[r[0]] and len(r) > 2]
    struct = []
    for it in range(len(times)):
        struct.append(float(np.mean([np.std(prof[it][r]) for r in strip_runs])))

    print(f"       {len(strip_runs)} strips resolved on the grid")
    print("        t [ns]      total (k=0)      structure (mean std WITHIN a strip)")
    for it, t in enumerate(times):
        print(f"       {t*1e9:9.1f}   {q_tot[it]:14.6e}   {struct[it]:14.6e}")

    cons = abs(q_tot[-1] / q_tot[0] - 1.0)
    relax = struct[-1] / struct[0]
    print(f"       charge conservation |q(10 µs)/q(0) - 1| = {cons:.3e}")
    print(f"       structure relaxation std(10 µs)/std(0) = {relax:.3e}")
    print("       Charge in a 250 µm insulating gap, and charge on a strip "
          "with no y-gradient,\n       both have NO DC path to ground. That "
          "asymmetry between x (blocked) and\n       y (conducting) is the "
          "geometric origin of the measured X/Y sharing\n       difference — "
          "plan §1.")
    # Convergence. sigma_s(x) is a step function, so a uniform grid resolves it
    # to first order and the residual in-strip structure falls ~linearly with
    # the sample spacing. Reporting the rate is worth more than a bare
    # pass/fail: it says what resolution a production run has to use.
    print("\n       grid convergence of the residual (must fall ~1/nx):")
    print("         nx     µm/sample   struct(10 µs)/struct(0)")
    resid = []
    for nxc in (390, 780, 1560, 3120):
        sc = WeightingSolver(rho_s_ohm_sq=1.0e6, d_kapton_m=75e-6,
                             nx=nxc, ny=8, ly_m=0.02)
        pc = sc.solve(strip_pattern_y(sc, x0_m=sc.lx / 2), [0.0, 1e-5])[:, 4, :]
        onc = sc.sigma_x > 0
        ec = np.flatnonzero(np.diff(onc.astype(int)) != 0) + 1
        rc = [r for r in np.split(np.arange(nxc), ec) if onc[r[0]] and len(r) > 2]
        st = [np.mean([np.std(pc[i][r]) for r in rc]) for i in (0, 1)]
        resid.append(st[1] / st[0])
        print(f"        {nxc:5d}   {sc.lx/nxc*1e6:8.2f}      {resid[-1]:.4e}")
    rates = [resid[i] / resid[i + 1] for i in range(len(resid) - 1)]
    print(f"       ratio per doubling: "
          f"{', '.join(f'{r:.2f}' for r in rates)}  (first order -> ~2)")
    print("       => a production run needs <=10 µm sampling in x, which is "
          "what plan §2\n          already specifies for the export grid.")

    ok1 = _report("V2 strip charge conserved (no drain term)", cons, 1e-6)
    # One-sided: convergence must be AT LEAST first order. It is slightly
    # better than that on coarse grids (ratio ~3) and approaches 2 asymptotically,
    # so an "equals 2" assertion would fail for the wrong reason.
    ok2 = _report("V2 residual converges at >= first order",
                  max(0.0, 1.8 - min(rates)), 1e-9,
                  extra=f"(ratios {', '.join(f'{r:.2f}' for r in rates)})")
    return ok1 and ok2


# ── V4 ───────────────────────────────────────────────────────────────────────

def v4_drain_vs_measured():
    """
    Plan §3 V4: the late-time drain tail must land within a factor ~2 of the
    measured tau_g = 5.3-7.3 µs, "else revisit A1".

    A1 is that the ESL strips are tied to the bus at both y-ends of the active
    area, so the drain is set by diffusion along a 412 mm strip. This test
    computes that number and compares.
    """
    print("\nV4  global drain: assumption A1 vs the measured tau_g")
    TAU_MEAS = (5.3e-6, 7.3e-6)
    print("        rho_s      tau(A1, L=412 mm)    ratio to measured 5.3-7.3 µs")
    ratios = []
    for rho in C.RHO_S_SCAN_OHM_SQ:
        tau = W.drain_time_from_strip_ends(rho)
        r = tau / np.mean(TAU_MEAS)
        ratios.append(r)
        print(f"       {rho/1e6:4.1f} MΩ/sq   {tau*1e3:10.3f} ms      "
              f"x{r:,.0f}")
    best = min(ratios)
    print()
    print("       Even at the most favourable corner of the resistivity scan, "
          "A1 predicts a\n       drain three orders of magnitude SLOWER than "
          "what is measured.")
    for rho in C.RHO_S_SCAN_OHM_SQ:
        L = W.effective_ground_spacing(np.mean(TAU_MEAS), rho)
        print(f"       to get 6.3 µs at {rho/1e6:4.1f} MΩ/sq the strips would "
              f"have to be grounded\n         every {L*1e3:6.1f} mm — not "
              f"every {412:.0f} mm.")
    print()
    print("       CONCLUSION: A1 as literally stated is refuted by the data.")
    print("       Candidate explanations, in the order worth testing:")
    print("         (a) the strips are not the drain path at all — the "
          "measured tau_g may be\n             the pad-plane/front-end RC "
          "rather than sheet transport;")
    print("         (b) there is a conductive path at a much finer pitch than "
          "the active width\n             (a grounding grid, or the paste "
          "contacting the bus through vias);")
    print("         (c) rho_s is far below the scanned 0.5-5 MΩ/sq.")
    print("       Until this is settled, the solver's tau_drain_s is an "
          "EXPLICIT INPUT and the\n       production scan should carry it as "
          "a nuisance parameter pinned to the\n       measurement, not "
          "derived from A1.")
    # This test is expected to FAIL against A1 — that failure is the result.
    return _report("V4 A1 reproduces measured tau_g (EXPECTED FAIL)",
                   abs(np.log10(best)), 0.3,
                   extra="<- refutes A1, see above")


def v4b_drain_term_works(nx=780, ny=128):
    """
    With tau_drain_s supplied, every mode — including k = 0 and the
    gap-supported null space — must decay at the requested rate. This is what
    makes V3's late-time limit reachable at all.
    """
    print("\nV4b explicit drain term")
    tau = 6.3e-6
    s = WeightingSolver(rho_s_ohm_sq=1.0e6, d_kapton_m=75e-6,
                        nx=nx, ny=ny, tau_drain_s=tau)
    vpad = pad_pattern(s, x0_m=s.lx / 2, y0_m=0.0)
    ts = np.array([0.0, tau, 3 * tau, 10 * tau, 30 * tau])
    psi = s.solve(vpad, ts)
    n0 = np.abs(psi[0]).max()
    print("        t/tau      sup|Psi| / prompt      exp(-t/tau)")
    for t, p in zip(ts, psi):
        print(f"       {t/tau:7.1f}   {np.abs(p).max()/n0:14.6e}   "
              f"{np.exp(-t/tau):12.4e}")
    late = np.abs(psi[-1]).max() / n0

    # The tail is SLOWER than exp(-t/tau), and that is a physical prediction,
    # not a numerical defect. A uniform leak conductance g_leak gives mode k a
    # drain rate g_leak/C(k). The charge stranded in the 250 µm inter-strip
    # gaps lives at k ~ 2*pi/250 µm, where C(k) is roughly twice its k -> 0
    # value c', so those modes drain at roughly HALF the rate of the strips.
    k_gap = 2 * np.pi / (C.ESL_GAP_UM * 1e-6)
    C_gap, _ = stack_coeffs(np.array([k_gap]), C.AMP_GAP_M, 75e-6)
    ratio = float(C_gap[0] / C.c_prime(kapton_m=75e-6))
    print()
    print(f"       C(k) at the 250 µm gap scale is {ratio:.2f}x its k->0 value,")
    print(f"       so gap-stranded charge drains ~{ratio:.1f}x slower than "
          "charge on the strips.")
    print("       PREDICTION: the late-time tail is two-component, not a single")
    print("       exponential — a fast strip component and a slower gap "
          "component. Worth\n       looking for in the measured tau_g fits, "
          "which currently assume one.")
    return _report("V4b drain empties the sheet by t = 30 tau", late, 1e-3)


# ── V3 ───────────────────────────────────────────────────────────────────────

def v3_limits_and_conservation(nx=780, ny=256):
    """
    The two limits the extended Ramo-Shockley construction must satisfy:

      t -> 0+ : the sheet is a plain dielectric. The prompt solution from the
                propagator must equal the independently-computed
                S(k)/C(k) prompt solution to machine precision.
      t -> inf: a UNIFORM sheet is a grounded conductor, so Psi on the sheet
                must go to zero everywhere. (With gaps it does not, and must
                not — see V2.)

    Also checks the propagator is a contraction: no mode may grow.
    """
    print("\nV3  prompt limit, late-time limit, no growing modes")
    # A drain is REQUIRED for the late-time limit to exist at all: on a
    # periodic box with no leak, the k = 0 mode has decay rate k^2 D = 0 and
    # the sheet's mean potential is conserved forever. That is not a solver
    # defect, it is the absence of a ground connection — see V4.
    s = WeightingSolver(rho_s_ohm_sq=1.0e6, d_kapton_m=75e-6,
                        nx=nx, ny=ny, uniform_sheet=True, tau_drain_s=6.3e-6)
    vpad = pad_pattern(s, x0_m=s.lx / 2, y0_m=0.0)

    prompt_direct = s.prompt_potential(vpad)
    prompt_prop = s.solve(vpad, [0.0])[0]
    e_prompt = (np.abs(prompt_prop - prompt_direct).max() /
                np.abs(prompt_direct).max())

    late = s.solve(vpad, [1e-3])[0]
    e_late = np.abs(late).max() / np.abs(prompt_direct).max()

    # Monotonicity: the sup-norm must never increase with time.
    ts = np.array([0.0, 1e-9, 1e-8, 1e-7, 1e-6])
    norms = [np.abs(s.solve(vpad, [t])[0]).max() for t in ts]
    growth = max(0.0, max(norms[i + 1] / norms[i] - 1 for i in range(len(ts) - 1)))
    print("        t [ns]      sup|Psi|")
    for t, n in zip(ts, norms):
        print(f"       {t*1e9:9.1f}   {n:.6e}")

    ok1 = _report("V3 prompt via propagator == S/C direct", e_prompt, 1e-10)
    ok2 = _report("V3 uniform sheet -> grounded conductor", e_late, 1e-6)
    ok3 = _report("V3 no growing modes", growth, 1e-12)
    return ok1 and ok2 and ok3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="coarser grids; use for a smoke test only")
    a = ap.parse_args()

    print("=" * 74)
    print("S1 weighting-potential solver — acceptance tests")
    print(f"gap {C.AMP_GAP_UM:.0f} µm (confirmed) | ESL "
          f"{C.ESL_WIDTH_UM:.0f}/{C.ESL_GAP_UM:.0f} µm on "
          f"{C.ESL_PITCH_UM:.0f} µm | pads {C.PAD_PITCH_UM:.0f} µm | "
          f"superperiod {C.SUPERPERIOD_UM/1000:.1f} mm")
    print("=" * 74)

    nx, ny = (390, 128) if a.quick else (780, 256)
    # Solver correctness: these MUST all pass. A failure here is a code bug.
    correctness = [
        v1_uniform_sheet(nx=nx, ny=ny),
        v1b_gaussian_limit(),
        v2_single_strip(nx=3120 if not a.quick else 780,
                        ny=64 if not a.quick else 32),
        v3_limits_and_conservation(nx=nx, ny=ny),
        v4b_drain_term_works(nx=nx, ny=min(ny, 128)),
    ]
    # Physics confrontation: this compares an ASSUMPTION against measured data.
    # It is reported separately because a failure here is a finding about the
    # detector, not a defect in this code, and it must not be "fixed" by
    # loosening a tolerance.
    a1_ok = v4_drain_vs_measured()

    print("\n" + "=" * 74)
    print(f"solver correctness : {sum(correctness)}/{len(correctness)} passed")
    print(f"physics assumption : A1 {'consistent' if a1_ok else 'REFUTED'} "
          "by the measured drain time")
    if not a1_ok:
        print()
        print("  A1 being refuted is a RESULT, not a bug. It means the ESL")
        print("  strips are not drained through their y-ends alone. The solver")
        print("  is correct either way: it takes tau_drain_s as an explicit")
        print("  input, so the production scan can be pinned to the measured")
        print("  5.3-7.3 µs while the mechanism is worked out.")
    print("=" * 74)
    return 0 if all(correctness) else 1


if __name__ == "__main__":
    sys.exit(main())
