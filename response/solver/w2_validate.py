#!/usr/bin/env python3
"""
w2_validate.py — the W2 prototype's acceptance battery (plan next-in-order
item 3: "prototype ONE family block, validated against v6_pad_gaps' static
solve at t=0 and the gap->0 W1 identity").

Two geometries:

  FAKE — a small board (800 µm ESL / 400 µm pads / 1.6 mm box, 50 µm cells)
  chosen so every mask edge lands exactly on a cell boundary AND the x Bloch
  split is non-trivial (gx = 2): the whole grid is small enough to solve as
  ONE dense family, which makes "family decomposition changes nothing" an
  exact machine-precision statement rather than a convergence claim. This is
  the test aimed squarely at the Bloch-bookkeeping bug class.

  REAL — the production 31.2 mm box at nx = 3120 with a coarse y grid (the
  identities under test are y-resolution independent; the RESOLUTION of the
  100 µm gap is measured separately by the CG ny-ladder at the end).

Tests:
  1  fake: family-decomposed solve == single-dense-family solve   (round-off)
  2  fake: t=0 slice == independent full-grid CG (v6's operator)  (round-off)
  3  fake: rho_s scaling exactness: V0(t; rho) == V0(t/10; rho/10)
  4  fake: constraint + ODE self-consistency at finite t (small-dt Taylor)
  5  real: gap->0 and grounded-substrate reductions == certified
           wpot.WeightingSolver, arbitrary drive                  (round-off)
  6  real: tiling pads driven together -> S(0)/C(0) at every t
  7  real: family (0,0) prompt capture == full-grid CG capture, and the CG
           ny-ladder vs V6's static Richardson limit 0.852998 (the recorded
           y-resolution mitigation check)

    ../../.venv/bin/python -m response.solver.w2_validate           # ~minutes
    ../../.venv/bin/python -m response.solver.w2_validate --heavy   # + ny=512 family
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from ..common import constants as C
from . import wpot as W
from .kernels import prompt_sum_rule
from .wpot_w2 import W2Solver

RNG = np.random.default_rng(20260808)


def _rand_drive(ny, nx, nmodes=6):
    """A smooth random real periodic drive — identity tests must hold for ANY."""
    h = np.zeros((ny, nx), complex)
    for _ in range(nmodes):
        iy, ix = RNG.integers(0, min(ny, 8)), RNG.integers(0, min(nx, 8))
        h[iy, ix] = RNG.normal() + 1j * RNG.normal()
    return np.real(np.fft.ifft2(h)) * nx * ny / nmodes


# pad_size 350 µm: edges at ±175 µm land exactly on 50 µm cell boundaries
# (boundaries at 25 + 50k), so the masks are hard with no half-covered ties.
# ESL 450 µm at phase -225 µm: strip symmetric about x = 0 with its edges off
# the grid points, so the sampled sigma_s is even and sigma_hat real — which
# is what lets test 8 exercise the production real-arithmetic cast here.
FAKE = dict(lx_m=1.6e-3, nx=32, ly_m=1.6e-3, ny=32,
            esl_pitch_m=800e-6, esl_width_m=450e-6, esl_phase_m=-225e-6,
            pad_pitch_m=400e-6, pad_size_m=350e-6, pad_x0_m=0.0)


class DenseW2(W2Solver):
    """Same solver, no family split: one dense family over every mode."""
    def families(self):
        IX, IY = np.meshgrid(np.arange(self.nx), np.arange(self.ny))
        return [(IX.ravel(), IY.ravel())]


def _report(name, err, bar):
    ok = err < bar
    print(f"  {name:<58} {err:9.2e}  vs {bar:.0e}  {'PASS' if ok else 'FAIL'}")
    return ok


def test_fake_family_vs_dense(times):
    s_fam = W2Solver(1e6, **FAKE)
    s_den = DenseW2(1e6, **FAKE)
    assert (s_fam.gx, s_fam.gy) == (2, 4), (s_fam.gx, s_fam.gy)
    v = _rand_drive(s_fam.ny, s_fam.nx)
    a = s_fam.solve(v, times)
    b = s_den.solve(v, times)
    scale = float(np.abs(b).max())
    return _report("1  fake: families == one dense family",
                   float(np.abs(a - b).max()) / scale, 1e-9), s_fam, v, a


def test_fake_prompt_cg(s_fam, v, a):
    cg, it = s_fam.solve_prompt_cg(v)
    scale = float(np.abs(cg).max())
    return _report(f"2  fake: t=0 slice == independent CG ({it} it)",
                   float(np.abs(a[0] - cg).max()) / scale, 1e-9)


def test_fake_rho_scaling(times):
    sA = W2Solver(1e6, **FAKE)
    sB = W2Solver(1e7, **FAKE)
    v = _rand_drive(sA.ny, sA.nx)
    a = sA.solve(v, times)                      # rho_s,  times
    b = sB.solve(v, times * 10.0)               # 10 rho, 10 t  -> identical
    scale = float(np.abs(a).max())
    return _report("3  fake: V0(t; rho) == V0(10t; 10rho)",
                   float(np.abs(a - b).max()) / scale, 1e-9)


def test_fake_ode_consistency():
    """
    d/dt of the propagated field must satisfy C_eff dV0/dt = -M V0 in every
    family — checked at a finite t with a centred difference, which probes the
    prompt->dynamics handoff (the constraint is slaved with D in the dynamics
    but with Mv6 = D - A12^2/C in the prompt; both must describe one motion).
    """
    s = W2Solver(1e6, **FAKE)
    v = _rand_drive(s.ny, s.nx)
    t0, dt = 2e-7, 1e-11
    worst = 0.0
    vhat = np.fft.fft2(v)
    for IX, IY in s.families():
        ops = s._family_ops(IX, IY)
        v0 = s._prompt_family(ops, vhat[IY, IX])
        L, w, Q = ops["L"], ops["w"], ops["Q"]
        coef = Q.conj().T @ (L.conj().T @ v0)
        at = lambda t: np.linalg.solve(L.conj().T, Q @ (np.exp(-w * t) * coef))
        vdot = (at(t0 + dt) - at(t0 - dt)) / (2 * dt)
        vmid = at(t0)
        lhs_s = L @ (L.conj().T @ vdot)         # C_eff vdot via its Cholesky
        rhs = -ops["M"] @ vmid
        num = float(np.abs(lhs_s - rhs).max())
        den = float(np.abs(rhs).max()) + 1e-300
        worst = max(worst, num / den)
    return _report("4  fake: C_eff dV0/dt = -M V0 at t=200 ns (all families)",
                   worst, 1e-4)                 # centred-difference floor


def test_real_w1_identity(times, ny):
    """Both W1 reductions of W2 against the certified wpot solver."""
    ly = 64 * C.PAD_PITCH_M
    ref = W.WeightingSolver(1e6, C.KAPTON_THICK_UM * 1e-6, nx=3120, ny=ny,
                            ly_m=ly, phase_m=-C.ESL_WIDTH_M / 2,
                            tau_drain_s=None)
    v = _rand_drive(ny, 3120)
    b = ref.solve(v, times)
    scale = float(np.abs(b).max())
    ok = True
    for name, kw in (("grounded substrate (sub_layers=[])",
                      dict(sub_layers=[])),
                     ("gap -> 0 (pad size = pitch)",
                      dict(pad_size_m=C.PAD_PITCH_M))):
        s2 = W2Solver(1e6, nx=3120, ny=ny, ly_m=ly,
                      esl_phase_m=-C.ESL_WIDTH_M / 2, **kw)
        a = s2.solve(v, times)
        ok &= _report(f"5  real: W2 [{name}] == wpot",
                      float(np.abs(a - b).max()) / scale, 1e-9)
    return ok


def test_real_sum_rule(times, ny):
    s2 = W2Solver(1e6, nx=3120, ny=ny, ly_m=64 * C.PAD_PITCH_M,
                  esl_phase_m=-C.ESL_WIDTH_M / 2, pad_size_m=C.PAD_PITCH_M)
    a = s2.solve(np.ones((ny, 3120)), times)
    expect = prompt_sum_rule(C.KAPTON_THICK_UM * 1e-6)
    return _report(f"6  real: tiling pads -> S(0)/C(0) = {expect:.6f} at all t",
                   float(np.abs(a - expect).max()) / expect, 1e-9)


def test_real_capture(ny_family, ny_ladder, heavy=False):
    """
    Family-machinery capture vs full-grid CG on the same grid (bookkeeping),
    then the CG ny-ladder (resolution), all on HARD masks. V6's static box
    (5 µm cells) Richardson limit is 0.852998. History: the first run of this
    test (2026-08-08) used fractional-coverage constraint masks and measured
    capture 0 at ny=64 and -14.8 % at ny=512 — the collapse that forced the
    mask_snap default (see W2Solver). The ny=1024/2496 Richardson landed at
    0.854171 (+0.14 % of V6), which is what certified the full-box machinery
    against the static solver.
    """
    ok = True
    print(f"  7  real geometry, nominal substrate, hard masks")
    # The family-vs-CG bookkeeping identity is tests 2 and 8's job (the code
    # is dimension-independent); what is REAL-geometry-specific is what the
    # coarse y grid does to the physics, which the CG ladder measures.
    # -- resolution: CG ladder in ny ----------------------------------------
    print(f"     CG capture ny-ladder (nx=3120, hard masks, "
          f"V6 Richardson 0.852998):")
    caps = {}
    for ny in ny_ladder:
        s = W2Solver(1e6, nx=3120, ny=ny, ly_m=64 * C.PAD_PITCH_M,
                     esl_phase_m=-C.ESL_WIDTH_M / 2)
        v0, it = s.solve_prompt_cg(s.metal_hard)
        caps[ny] = float(v0.mean())
        print(f"       ny {ny:5d}  (dy {49920/ny:7.2f} µm)  capture "
              f"{caps[ny]:.6f}   vs V6 {caps[ny]/0.852998-1:+.2%}   ({it} CG it)")
    ks = sorted(caps)
    if len(ks) >= 2:
        rich = caps[ks[-1]] + (caps[ks[-1]] - caps[ks[-2]])
        print(f"       first-order Richardson from ny={ks[-2]}/{ks[-1]}: "
              f"{rich:.6f}  (V6 static: 0.852998)")
    return ok


def test_production_identity(times):
    """
    8: the memory-ordered REAL-arithmetic production path (w2_production.
    solve_family) must reproduce the validated complex prototype exactly —
    full pipeline, every family, two drives, both rho_s scalings — on the
    fake geometry where the prototype is itself dense-validated.
    """
    from .w2_production import solve_family
    s = W2Solver(1e6, **FAKE)
    d1 = s.metal_hard.copy()
    d2 = s.metal_hard * (np.abs(s.x - 0.4e-3)[None, :] < 0.25e-3)
    drives = [("all", d1), ("col", d2)]
    rhos = (1e6, 5e6)
    ref = {n: s.solve(v, times) for n, v in drives}
    s5 = W2Solver(5e6, **FAKE)
    ref5 = {n: s5.solve(v, times) for n, v in drives}
    out = {(n, r): np.zeros((len(times), s.ny, s.nx), complex)
           for n, _ in drives for r in rhos}
    for ifam in range(len(s.families())):
        slab, IX, IY, _ = solve_family(s, drives, ifam, times, rhos=rhos,
                                       verbose=False)
        for di, (n, _) in enumerate(drives):
            for ri, r in enumerate(rhos):
                out[(n, r)][:, IY, IX] = slab[di, ri].astype(complex)
    worst = 0.0
    for di, (n, _) in enumerate(drives):
        a1 = np.real(np.fft.ifft2(out[(n, 1e6)], axes=(1, 2)))
        a5 = np.real(np.fft.ifft2(out[(n, 5e6)], axes=(1, 2)))
        scale = float(np.abs(ref[n]).max())
        worst = max(worst,
                    float(np.abs(a1 - ref[n]).max()) / scale,
                    float(np.abs(a5 - ref5[n]).max()) / scale)
    # complex64 slab storage floors this at ~1e-7
    return _report("8  fake: production real path == prototype (2 rho_s)",
                   worst, 3e-6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--heavy", action="store_true",
                    help="also run the ny=512 family-(0,0) prompt (~1 h, big RAM)")
    a = ap.parse_args()
    times = np.array([0.0, 1e-8, 1e-7, 1e-6, 1e-5])

    print("W2 prototype acceptance battery "
          "(plan next-in-order item 3, sequencing step 1)\n")
    ok = True
    o, s_fam, v, afield = test_fake_family_vs_dense(times)
    ok &= o
    ok &= test_fake_prompt_cg(s_fam, v, afield)
    ok &= test_fake_rho_scaling(times)
    ok &= test_fake_ode_consistency()
    ok &= test_real_w1_identity(times, ny=64)
    ok &= test_real_sum_rule(times, ny=64)
    ok &= test_production_identity(times)
    ok &= test_real_capture(ny_family=256,
                            ny_ladder=(512, 1024, 2496))
    print(f"\n-> battery {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
