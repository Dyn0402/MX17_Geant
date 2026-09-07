#!/usr/bin/env python3
"""
wpot_w2.py — S1 on the W2 boundary: the dynamic solver with the 100 µm
inter-pad channels FLOATING instead of grounded.

Why this file exists: V6 (`v6_pad_gaps.py`, design/report/V6_PAD_GAPS_2026-08-08.md)
showed that W1's grounded inter-pad channels cost 27 % of the prompt capture and
paint a 4.5× sub-pad amplitude modulation onto the kernels that the real board
does not have. V6 only solved t = 0+. This module carries the same boundary
through the DYNAMICS, so the full time-dependent comb kernels can be re-solved
on W2. Plan: design/RESPONSE_SIM_PLAN.md, next-in-order item 3 (sizing and
sequencing recorded there 2026-08-08).


THE PHYSICS, EXTENDING wpot.py
==============================

Same stack as `wpot.py` (gas / ESL sheet at z = 0 / kapton+glue / pad plane at
z = -d), with ONE change at z = -d: the plane is a Dirichlet pattern only on the
METAL (680 µm pads on a 780 µm pitch); in the 100 µm channels between pads the
potential Vd floats, loaded from below by the board (glue trench + FR4 to L5,
`substrate_admittance`). Per lateral mode k, with the insulator as the
symmetric two-port (A11, A12, A22) of `v6_pad_gaps.stack_2port`:

    sheet free charge   q  = C(k) V0 - A12(k) Vd          C = eps0 k coth(kg) + A11
    gap free charge     0  = (A22 + Y_sub)(k) Vd - A12(k) V0     ON THE GAP ONLY
    metal               Vd = V_pad (the drive)
    sheet conduction    dq/dt = div( sigma_s(x) grad V0 )  =  -M_sheet V0

The gap carries NO free charge at any time — it is buried insulator with no
conduction path — so the second line is a constraint that holds for all t, not
just at t = 0+. Writing Vd = V_pad + u with u supported on the gap, the
constraint slaves u to V0 linearly:

    G u = R [ A12 V0 - (A22 + Y_sub) V_pad ],   G = R (A22 + Y_sub) R

with R = multiplication by the gap mask. Substituting into the sheet charge and
using that V_pad is constant for t > 0:

    C_eff dV0/dt = -M_sheet V0,     C_eff = C - (A12 R) G^+ (R A12)

C_eff is the Schur complement of the floating-gap block in the full capacitance
form: symmetric, positive definite, and — unlike W1's C — NOT diagonal in k,
because R couples modes across the pad lattice. That coupling is the whole
cost of W2: with it, the exact Bloch families are

    x: ESL couples Δix ≡ 0 (mod 39), pads couple Δix ≡ 0 (mod 40); gcd = 1,
       so ALL nx x-modes join one family.
    y: pads couple Δiy ≡ 0 (mod Ly/pitch = 64); the ESL is uniform in y.
       Families are the 64 residue classes of iy.

(Family membership is computed generically from the pitches, not hardcoded —
the fake small geometries the validation battery uses have non-trivial x
classes, which is exactly what exercises the index bookkeeping.)

Initial condition — V6 verbatim. At t = 0+ the sheet also carries no charge,
q = 0 pointwise, so V0 = (A12/C) Vd and the constraint becomes

    R [ (A22 + Y_sub - A12^2/C) (V_pad + u) ] = 0

which is exactly `v6_pad_gaps.PadPlane.solve_prompt`'s operator. The prompt
solved here per family must therefore agree with an independent full-grid CG
solve of that equation (`solve_prompt_cg` below, same structure as V6's, matrix
free) — that is one of the prototype's acceptance tests, and V6's static solve
is the resolution anchor (Richardson: capture 0.852998).

Propagation is exact, no time stepping: Cholesky C_eff = L L^H, then
B = L^-1 M_sheet L^-H is Hermitian PSD, and

    V0(t) = L^-H Q exp(-w t) Q^H L^H V0(0+).

Genuine zero eigenvalues are kept — charge stranded on gap-isolated ESL strips
never decays, same as W1. And because M_sheet ∝ 1/rho_s while C_eff is pure
electrostatics, the eigenvectors are rho_s-INVARIANT and w scales as 1/rho_s:
one eigendecomposition serves the whole rho_s grid (checked in the battery).

LIMITS THAT MUST COME BACK OUT (the battery enforces all of them):
  * gap -> 0 (pad size = pitch): R = 0, C_eff = C, W2 IS W1 identically —
    checked against the certified `wpot.WeightingSolver` to round-off.
  * sub_layers=[] terminated (ground AT the pad plane): W1 again, u = 0.
  * tiling pads driven together: Psi(t) = S(0)/C(0) = 0.881583 at every t.

RESOLUTION WARNING (recorded in the plan): the production y grid (ny = 512 →
97.5 µm cells, ny = 1024 → 48.75 µm) cannot put the 100 µm gap edges on cell
boundaries, so the masks are fractional-coverage there (exact in x at 10 µm).
The metal edge is an r^-1/2 singularity and V6 needed 5 µm cells to sit 0.12 %
below its own Richardson limit — quantify with `solve_prompt_cg` ny-ladders
against V6's static number before trusting any absolute W2 amplitude.

STATUS: prototype (2026-08-08). Correctness first: everything is complex and
nothing exploits the real-symmetric production registration or the y-parity
split yet — those are the recorded 8× / 2× production optimisations, to be
applied only after this passes its battery. Memory/flop budget at production
size (family N = 24,960 at ny = 512) is condor-scale, not laptop-scale.

    python3 -m response.solver.w2_validate        # the acceptance battery
"""

from __future__ import annotations

from math import gcd

import numpy as np
import scipy.linalg as sla

from ..common import constants as C
from . import wpot as W
from .v6_pad_gaps import stack_2port, substrate_admittance, NOMINAL_SUB
from .kernels import _rect_cov, PAD_ORIGIN_M


def _fam_classes(n, steps):
    """Residue classes of mode indices coupled by lattice steps (mod n)."""
    g = n
    for s in steps:
        g = gcd(g, s)
    return g


class W2Solver:
    """
    Dynamic weighting potential with the inter-pad channels floating.

    Grid/attribute conventions are `wpot.WeightingSolver`'s (x from 0, y
    centred, kernels.py drive-pattern helpers work unchanged). The pad METAL
    pattern is part of the geometry here, not just of the drive: it defines
    where the boundary is Dirichlet.

    `sub_layers` follows v6_pad_gaps: None = NOMINAL_SUB (glue trench + FR4 to
    L5, grounded), [] with sub_terminated=True = ground AT the pad plane = W1,
    [] with sub_terminated=False = nothing below at all (open limit).
    """

    def __init__(self, rho_s_ohm_sq, d_kapton_m=None, gap_m=C.AMP_GAP_M,
                 eps_r=C.KAPTON_EPS_R, d_glue_m=None, glue_eps_r=C.GLUE_EPS_R,
                 lx_m=C.SUPERPERIOD_M, nx=3120, ly_m=None, ny=512,
                 esl_pitch_m=C.ESL_PITCH_M, esl_width_m=C.ESL_WIDTH_M,
                 esl_phase_m=0.0,
                 pad_pitch_m=C.PAD_PITCH_M, pad_size_m=C.PAD_SIZE_UM * 1e-6,
                 pad_x0_m=None, pad_y0_m=0.0,
                 sub_layers=None, sub_terminated=True,
                 pinv_rtol=1e-10, mask_snap=True):
        self.rho_s = rho_s_ohm_sq
        self.d_k = C.KAPTON_THICK_UM * 1e-6 if d_kapton_m is None else d_kapton_m
        self.gap = gap_m
        self.eps_r = eps_r
        self.d_g = C.GLUE_THICK_UM * 1e-6 if d_glue_m is None else d_glue_m
        self.glue_eps_r = glue_eps_r
        self.pinv_rtol = pinv_rtol

        ly_m = 64 * pad_pitch_m if ly_m is None else ly_m       # kernels.Y_BOX_M
        self.lx, self.ly, self.nx, self.ny = lx_m, ly_m, nx, ny
        self.x = np.arange(nx) * (lx_m / nx)
        self.y = (np.arange(ny) - ny // 2) * (ly_m / ny)
        self.kx = 2 * np.pi * np.fft.fftfreq(nx, d=lx_m / nx)
        self.ky = 2 * np.pi * np.fft.fftfreq(ny, d=ly_m / ny)

        # ── lattice bookkeeping: everything must be commensurate ────────────
        def _count(box, pitch, name):
            m = box / pitch
            if abs(m - round(m)) > 1e-9:
                raise ValueError(f"{name}: box {box} not a multiple of {pitch}")
            return int(round(m))

        self.n_esl = _count(lx_m, esl_pitch_m, "ESL x")
        self.n_col = _count(lx_m, pad_pitch_m, "pad x")
        self.n_row = _count(ly_m, pad_pitch_m, "pad y")

        # ── ESL conductivity (x only; uniform in y) ─────────────────────────
        sig_x = W.esl_sigma_profile(self.x, rho_s_ohm_sq, width_m=esl_width_m,
                                    pitch_m=esl_pitch_m, phase_m=esl_phase_m)
        self.sigma_x = sig_x
        self.sigma_hat = np.fft.fft(sig_x) / nx

        # ── pad metal mask, fractional coverage (exact where commensurate) ──
        pad_x0_m = (PAD_ORIGIN_M % pad_pitch_m) if pad_x0_m is None else pad_x0_m
        cols = (pad_x0_m + np.arange(self.n_col) * pad_pitch_m) % lx_m
        rows = (np.arange(self.n_row) - self.n_row // 2) * pad_pitch_m + pad_y0_m
        covx = np.zeros(nx)
        covy = np.zeros(ny)
        for x0 in cols:
            covx += _rect_cov(self.x, x0, pad_size_m, lx_m)
        for y0 in rows:
            covy += _rect_cov(self.y, y0, pad_size_m, ly_m)
        self.metal2d = np.outer(covy, covx)
        frac = 1.0 - self.metal2d                          # gap coverage, [0, 1]
        # THE CONSTRAINT SUPPORT MUST BE HARD (measured 2026-08-08, battery
        # test 7 first run). Any cell with r > 0 gets the FULL no-free-charge
        # equation enforced — the r-weighted residual r·[...] = 0 divides
        # through by r — and that equation's only global solution is Vd = 0.
        # With fractional coverage every partially-covered metal-edge cell
        # therefore leaks the Dirichlet clamp: capture read -14.8 % at ny=512,
        # and at ny=64 (every cell fractional) the whole plane floated and
        # capture collapsed to 0. Fractional coverage is the right tool for
        # DRIVES (charge totals); boundary-condition support is a set, not a
        # weight. Snapping at half coverage makes the ny=512 y-mask 7 clamped
        # cells + one pure 97.5 µm gap cell per row period (gap area -2.5 %,
        # a ~0.5 % capture bias by V6's gap scan) — and is a no-op wherever
        # the grid is commensurate (x at 10 µm, the fake battery geometry,
        # ny=2496).
        # Round before comparing: a cell exactly half-covered lands at
        # 0.5 ± 1 ulp differently at different pads, and an inconsistent
        # tie-break breaks the mask's lattice periodicity (caught by
        # _assert_support on the battery's first fake-geometry run).
        self.rgap = ((np.round(frac, 9) >= 0.5).astype(float)
                     if mask_snap else frac)
        self.metal_hard = 1.0 - self.rgap
        self.ghat = np.fft.fft2(self.rgap) / (nx * ny)

        # ── substrate below the gap ─────────────────────────────────────────
        self.sub_layers = NOMINAL_SUB if sub_layers is None else sub_layers
        self.sub_terminated = sub_terminated
        # sub_layers=[] terminated == ground AT the pad plane == W1's boundary
        self.w1_grounded = sub_terminated and not self.sub_layers
        self.has_gap = (not self.w1_grounded) and float(self.rgap.max()) > 1e-12

        # ── Bloch families, from the ACTIVE couplings ───────────────────────
        # The ESL always couples x modes in steps of n_esl. The pad lattice
        # couples (x, y) in steps of (n_col, n_row) — but only when the gap
        # boundary is live; with the gap grounded or absent the mask never
        # enters an operator and including its steps would only bloat the
        # blocks wpot already solves exactly.
        x_steps = [self.n_esl] + ([self.n_col] if self.has_gap else [])
        y_steps = [self.n_row] if self.has_gap else [ny]     # [ny] = no coupling
        self.gx = _fam_classes(nx, x_steps)
        self.gy = _fam_classes(ny, y_steps)
        self._fam_cache = {}

        # Guard: the operators' Fourier support must actually lie on the
        # lattices the family split assumes, else the split silently discards
        # coupling (the Bloch-bookkeeping bug class the plan warns about).
        self._assert_support()

    # ── support checks ──────────────────────────────────────────────────────

    def _assert_support(self):
        # An m-period-per-box pattern has Fourier support at index multiples
        # of m (mod n): m periods -> fundamental harmonic index m.
        s = np.abs(self.sigma_hat)
        keep = np.zeros(self.nx, bool)
        keep[::gcd(self.nx, self.n_esl)] = True
        off = float(s[~keep].max() / max(s.max(), 1e-300)) if (~keep).any() else 0.0
        if off > 1e-10:
            raise RuntimeError(f"sigma_s support off its lattice ({off:.1e}) — "
                               "ESL pattern not commensurate with the grid")
        if float(self.rgap.max()) < 1e-12:
            return                       # no gap: ghat is pure round-off noise
        g = np.abs(self.ghat)
        kx_ok = np.zeros(self.nx, bool)
        kx_ok[::gcd(self.nx, self.n_col)] = True
        ky_ok = np.zeros(self.ny, bool)
        ky_ok[::gcd(self.ny, self.n_row)] = True
        mask = np.outer(ky_ok, kx_ok)
        off = float(g[~mask].max() / max(g.max(), 1e-300)) if (~mask).any() else 0.0
        if off > 1e-10:
            raise RuntimeError(f"gap-mask support off the pad lattice ({off:.1e})")

    # ── per-mode diagonals ──────────────────────────────────────────────────

    def _diags(self, k):
        A11, A12, A22, Cg = stack_2port(k, gap_m=self.gap, d_kapton_m=self.d_k,
                                        d_glue_m=self.d_g, eps_r=self.eps_r,
                                        glue_eps_r=self.glue_eps_r)
        if self.w1_grounded:
            D = None                                  # infinite Y_sub: u == 0
        else:
            D = A22 + substrate_admittance(k, self.sub_layers,
                                           terminated=self.sub_terminated)
        return A12, Cg, D

    # ── families ────────────────────────────────────────────────────────────

    def families(self):
        """[(IX, IY) flat index arrays], one per Bloch family."""
        out = []
        for ry in range(self.gy):
            iys = np.arange(ry, self.ny, self.gy)
            for rx in range(self.gx):
                ixs = np.arange(rx, self.nx, self.gx)
                IX, IY = np.meshgrid(ixs, iys)
                out.append((IX.ravel(), IY.ravel()))
        return out

    def _pinv(self, H):
        """Pseudo-inverse of a Hermitian PSD matrix via eigh, with rank info."""
        w, V = sla.eigh(H)
        w = np.maximum(w, 0.0)
        keep = w > self.pinv_rtol * max(float(w.max()), 1e-300)
        Vk = V[:, keep]
        return (Vk / w[keep]) @ Vk.conj().T, int(keep.sum())

    def _family_ops(self, IX, IY):
        """
        Everything drive-independent for one family: the propagator pieces and
        the prompt gap operator. Cached, so 42 kernels pay for it once.
        """
        key = (int(IX[0]), int(IY[0]))
        if key in self._fam_cache:
            return self._fam_cache[key]
        kx, ky = self.kx[IX], self.ky[IY]
        k = np.hypot(kx, ky)
        A12, Cg, D = self._diags(k)

        # M_sheet: sigma couples same-ky modes across the ESL x lattice
        dIX = (IX[:, None] - IX[None, :]) % self.nx
        same_iy = (IY[:, None] == IY[None, :])
        M = self.sigma_hat[dIX] * same_iy * (np.outer(kx, kx) + np.outer(ky, ky))
        M = 0.5 * (M + M.conj().T)
        del dIX, same_iy

        ops = {"A12": A12, "Cg": Cg, "M": M}
        if self.has_gap:
            dIXg = (IX[:, None] - IX[None, :]) % self.nx
            dIYg = (IY[:, None] - IY[None, :]) % self.ny
            R = self.ghat[dIYg, dIXg]                  # Hermitian: real mask
            del dIXg, dIYg
            Mv6 = D - A12 ** 2 / Cg                    # v6's prompt operator
            Gp, rank_p = self._pinv(R @ (Mv6[:, None] * R))
            Gd, rank_d = self._pinv(R @ (D[:, None] * R))
            Ceff = np.diag(Cg).astype(complex)
            X = A12[:, None] * R
            Ceff -= X @ Gd @ X.conj().T
            del X, Gd
            Ceff = 0.5 * (Ceff + Ceff.conj().T)
            ops.update(R=R, Mv6=Mv6, Gp=Gp, rank=(rank_p, rank_d))
        else:
            Ceff = np.diag(Cg).astype(complex)

        L = sla.cholesky(Ceff, lower=True)
        del Ceff
        B = sla.solve_triangular(L, M.astype(complex), lower=True)
        B = sla.solve_triangular(L, B.conj().T, lower=True).conj().T
        B = 0.5 * (B + B.conj().T)
        w, Q = sla.eigh(B)
        del B
        w = np.maximum(w, 0.0)
        ops.update(L=L, w=w, Q=Q)
        self._fam_cache[key] = ops
        return ops

    # ── prompt (t = 0+), per family ─────────────────────────────────────────

    def _prompt_family(self, ops, vpad_f):
        """V0(0+) family coefficients for the family drive coefficients."""
        if not self.has_gap:
            return (ops["A12"] / ops["Cg"]) * vpad_f
        rhs = -ops["R"] @ (ops["Mv6"] * vpad_f)
        u = ops["R"] @ (ops["Gp"] @ rhs)
        return (ops["A12"] / ops["Cg"]) * (vpad_f + u)

    # ── public API ──────────────────────────────────────────────────────────

    def solve(self, vpad_xy, times_s, progress=False):
        """
        Psi(x, y, z=0, t) for a pad drive held from t = 0 — wpot.solve's
        contract, on the W2 boundary. t = 0 returns the prompt (= V6) field.
        """
        times = np.atleast_1d(np.asarray(times_s, float))
        vhat = np.fft.fft2(np.asarray(vpad_xy, float))
        out_hat = np.zeros((len(times), self.ny, self.nx), complex)
        fams = self.families()
        for i, (IX, IY) in enumerate(fams):
            ops = self._family_ops(IX, IY)
            v0 = self._prompt_family(ops, vhat[IY, IX])
            s0 = ops["L"].conj().T @ v0
            coef = ops["Q"].conj().T @ s0
            u = ops["Q"] @ (np.exp(-np.outer(ops["w"], times)) * coef[:, None])
            vt = sla.solve_triangular(ops["L"].conj().T, u, lower=False)
            out_hat[:, IY, IX] = vt.T
            if progress:
                print(f"  family {i + 1}/{len(fams)} (N={len(IX)})", flush=True)
        return np.real(np.fft.ifft2(out_hat, axes=(1, 2)))

    def solve_prompt_cg(self, vpad_xy, tol=1e-12, maxit=4000):
        """
        The t = 0+ field by full-grid preconditioned CG — v6_pad_gaps'
        solve_prompt, matrix-free on THIS grid and mask. Shares no code with
        the family assembly above (FFT products vs dense Bloch blocks), so
        agreement is a real check on the bookkeeping, and it is cheap at any
        ny — which makes it the resolution-ladder tool the plan asks for.
        """
        KX, KY = np.meshgrid(self.kx, self.ky, indexing="xy")
        k = np.hypot(KX, KY)
        A12, Cg, D = self._diags(k)
        vpad = np.asarray(vpad_xy, float)
        if not self.has_gap:
            return np.real(np.fft.ifft2(A12 / Cg * np.fft.fft2(vpad))), 0
        Mv6 = D - A12 ** 2 / Cg
        r_ = self.rgap
        mul = lambda f: np.real(np.fft.ifft2(Mv6 * np.fft.fft2(f)))
        amul = lambda u: r_ * mul(r_ * u)
        b = -(r_ * mul(vpad))
        n0 = float(np.abs(b).max())
        if n0 == 0.0:
            return np.real(np.fft.ifft2(A12 / Cg * np.fft.fft2(vpad))), 0
        prec = lambda v: r_ * np.real(np.fft.ifft2(np.fft.fft2(r_ * v) / Mv6))
        u = np.zeros_like(b)
        res = b - amul(u)
        z = prec(res)
        p = z.copy()
        rz = float((res * z).sum())
        it = 0
        for it in range(1, maxit + 1):
            Ap = amul(p)
            al = rz / float((p * Ap).sum())
            u += al * p
            res -= al * Ap
            if float(np.abs(res).max()) <= tol * n0:
                break
            z = prec(res)
            rz_new = float((res * z).sum())
            p = z + (rz_new / rz) * p
            rz = rz_new
        else:
            raise RuntimeError(f"W2 prompt CG did not converge in {maxit} it")
        vd = vpad + r_ * u
        return np.real(np.fft.ifft2(A12 / Cg * np.fft.fft2(vd))), it
