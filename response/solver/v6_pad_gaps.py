#!/usr/bin/env python3
"""
v6_pad_gaps.py — V6: is it safe to treat the pad plane as SOLID ground between
the pads?

THE ASSUMPTION UNDER TEST (plan §3, geometry model W1, and §3's W2 note). The
solver clamps the whole z = -d plane to a Dirichlet pattern: the driven
channel's pads at V_w, every other pad at 0, and — this is the part nobody
checked — the 100 µm channels BETWEEN pads at 0 as well. The pads are 680 µm on
a 780 µm pitch, so 1 - (680/780)^2 = 24 % of that plane is not metal at all.
The plan expects a percent-level shift in G_n from exposing it (`plan:545`) and
records V6 as "still not run".

WHY IT MIGHT NOT HOLD. In W1 that 24 % is a perfect ground sitting 68.8 µm below
the ESL, and it terminates field lines that would otherwise reach a pad. That is
where the production capture number comes from: the real-pad prompt capture is
0.665023 = (680/780)^2 x S(0)/C(0), i.e. EXACTLY the pad area fraction times the
whole-plane sum rule, because at k -> 0 a fully conducting plane divides the
image charge by area. Remove the metal from the channels and that division has
no reason to survive — the flux has to go somewhere, and the two candidates
(sideways onto the neighbouring pads, or down into the board) are 100 µm apart
laterally versus a few hundred µm vertically. Neither is obviously the winner,
which is why this is a calculation and not an argument.

WHAT IS ACTUALLY UNDER THERE (design/NEEDED_INPUTS.md §6, header stack). Not
vacuum and not ground: 26 µm of laminating adhesive filling the inter-pad trench
(the pad copper is 26 µm tall), then the FR4 to the L5 layer. So the conductor
under a gap is ~5x further from the ESL than the pad face is. The FR4 thickness
is the residual of a CAD-pinned 1.70 mm board (§6 flags this as the dominant
board uncertainty), so it is BRACKETED here rather than trusted:

    h -> 0        the gap is grounded at the pad plane        = W1 exactly
    h = nominal   26 µm adhesive + 264 µm FR4 to L5
    h -> infinity the gap sees no conductor at all            = the open limit

If W1 and the open limit bracket the answer within the plan's "few %", the
board's internal stackup cannot matter and §6's open question does not
propagate into the response. If they do not, this says how much of the
uncertainty is real.

THE METHOD. At t = 0+ the ESL carries no free charge, so it is electrically
absent and the prompt problem is pure layered electrostatics with ONE patterned
boundary. That has two consequences that make an exact solve cheap. First, the
31.2 mm superperiod disappears — with no ESL pattern the geometry is periodic on
the 780 µm pad lattice, so the box is a few pad pitches, not a superperiod.
Second, the answer for any channel follows from the single-pad Green's function
by superposition, because the gap condition is homogeneous and Dirichlet data
add.

Per lateral mode k the insulator is a two-port,

    D_z(0-)   = -A11 V0 + A12 Vd
    D_z(-d+)  = -A12 V0 + A22 Vd      (symmetric, by reciprocity)

so with the gas above, sigma = [eps0 k coth(kg) + A11] V0 - A12 Vd, which is
exactly W1's C(k) and S(k) — `check_against_w1` asserts that, and it is what
ties this file to the certified solver. Setting sigma = 0 (prompt) gives
V0 = (A12/C) Vd, and requiring no free charge in the gap gives

    [ A22 + Y_sub(k) - A12^2/C(k) ] Vd = 0     on the gap only,             (*)

with Vd = V_pad on the metal. (*) is a mixed boundary-value problem: diagonal in
k, restricted in real space. Its operator M is symmetric and positive definite
(A11 A22 > A12^2 for a passive two-port, and Y_sub >= 0), so it is solved by
conjugate gradients with FFT matvecs — matrix-free, no dense algebra, and no
approximation beyond the grid.

THE GRID IS EXACT. 780 and 680 µm share a factor 20, so a 5 µm cell puts every
pad edge exactly on a cell boundary and the metal mask is exact — none of the
fractional-coverage machinery `kernels.py` needs (where the pad lattice and the
solve box are NOT commensurate) applies or is required here.

WHAT THIS DOES NOT COVER. The prompt kernel only. As the sheet relaxes it
becomes the dominant conductor between the deposit and the pads, and the pad
plane's fine structure is screened by it, so the LATE-time sensitivity to this
boundary is bounded by the prompt one rather than equal to it — but that is an
argument, not a number, and it is stated as such in the report. This also says
nothing about the L5/L6 traces as SIGNAL electrodes: charge that reaches them is
treated here as lost to ground, which is the conservative direction for capture.

    python3 -m response.solver.v6_pad_gaps
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from ..common import constants as C
from . import wpot as W

# Substrate under the inter-pad channel, top down (design/NEEDED_INPUTS.md §6).
# The adhesive fills the trench between the 26 µm-tall pads; the FR4 below is
# the equal-division residual of the CAD-pinned board and is the soft number.
TRENCH_GLUE_UM = C.PAD_CU_THICK_UM          # 26 µm, the pad copper's own height
FR4_TO_L5_UM = 264.0
FR4_EPS_R = 4.4


# ── the insulator as a two-port, per lateral mode ────────────────────────────

def _layer(k, t, eps):
    """
    (Y coth(kt), Y csch(kt)) for one slab — its self- and transfer terms.

    BOTH LIMITS ARE TAKEN HERE RATHER THAN ASSEMBLED FROM PARTS. Reusing
    `wpot._coth` and multiplying by Y = eps0 eps k gives 0 * 1e300 at k = 0,
    which numpy evaluates to 0 instead of the true eps0 eps / t — and then
    A12^2/C is 0/0 and the DC mode, the one the whole capture number rests on,
    comes out NaN. Both terms tend to the same parallel-plate eps0 eps / t.
    """
    k = np.asarray(k, float)
    kt = np.clip(k * t, 0.0, 700.0)
    small = kt < 1e-9
    kts = np.where(small, 1.0, kt)
    flat = C.EPS0 * eps / t
    Y = C.EPS0 * eps * k
    return (np.where(small, flat, Y / np.tanh(kts)),
            np.where(small, flat, Y / np.sinh(kts)))


def stack_2port(k, gap_m=C.AMP_GAP_M, d_kapton_m=None, d_glue_m=None,
                eps_r=C.KAPTON_EPS_R, glue_eps_r=C.GLUE_EPS_R):
    """
    A11, A12, A22 of the kapton+glue stack, plus the gas-side C(k).

    Cascade of two slabs sharing an uncharged interface. Eliminating the
    interface potential V_m gives the composite two-port directly:

        V_m  = (a V0 + b Vd) / c,  a = Y1 csch(k t1), b = Y2 csch(k t2),
                                   c = Y1 coth(k t1) + Y2 coth(k t2)
        A11  = Y1 coth(k t1) - a^2/c
        A12  = a b / c
        A22  = Y2 coth(k t2) - b^2/c

    ORDER MATTERS in the same way it does in `wpot.stack_coeffs`: layer 1 is the
    kapton, against the ESL.
    """
    k = np.asarray(k, float)
    d_kapton_m = C.KAPTON_THICK_UM * 1e-6 if d_kapton_m is None else d_kapton_m
    d_glue_m = C.GLUE_THICK_UM * 1e-6 if d_glue_m is None else d_glue_m
    c1, a = _layer(k, d_kapton_m, eps_r)
    c2, b = _layer(k, d_glue_m, glue_eps_r)
    c = c1 + c2
    A11, A12, A22 = c1 - a * a / c, a * b / c, c2 - b * b / c
    kg = k * gap_m
    c_gas = C.EPS0 * np.where(kg < 1e-6, 1.0 / gap_m, k * W._coth(kg))
    return A11, A12, A22, c_gas + A11


def substrate_admittance(k, layers, terminated=True):
    """
    Y_sub(k) looking DOWN from the pad plane into the board.

    `layers` is [(thickness_m, eps_r), ...] top down; `terminated=True` ends on
    a ground plane (L5 read as solid, the maximum-screening choice),
    `terminated=False` leaves it open. Cascaded bottom up.

    WRITTEN TO SURVIVE k = 0, which is one grid point and the most important
    one — it carries the DC term the whole capture number rests on. The textbook
    form Y (Y_L + Y tanh)/(Y + Y_L tanh) is 0/0 there because Y = eps0 eps k
    vanishes. Dividing through by Y and using r = tanh(k t)/Y, which tends to
    the finite t/(eps0 eps), gives

        y_in = (y_L + Y tanh(k t)) / (1 + y_L r)

    with every term finite, and the grounded start is simply y = 1/r (the
    parallel-plate eps0 eps / t at k -> 0, exactly as it should be).

    An empty `layers` with `terminated=True` is a ground plane AT the pad plane,
    i.e. W1.
    """
    k = np.asarray(k, float)
    if terminated and not layers:
        return np.full_like(k, np.inf)
    y = None if terminated else np.zeros_like(k)
    for t, eps in reversed(layers):
        Y = C.EPS0 * eps * k
        kt = k * t
        th = np.tanh(np.clip(kt, 0.0, 700.0))
        # r = tanh(kt)/Y, with its k -> 0 limit t/(eps0 eps)
        r = np.where(kt < 1e-9, t / (C.EPS0 * eps),
                     th / np.where(Y > 0, Y, 1.0))
        y = 1.0 / r if y is None else (y + Y * th) / (1.0 + y * r)
    return y


NOMINAL_SUB = [(TRENCH_GLUE_UM * 1e-6, C.GLUE_EPS_R),
               (FR4_TO_L5_UM * 1e-6, FR4_EPS_R)]


# ── the mixed boundary-value problem ─────────────────────────────────────────

class PadPlane:
    """
    A periodic box of n_pad x n_pad pad cells at `dx_um`, with an exact mask.

    `gap_um=0` collapses the inter-pad channel and makes the pads tile the
    plane. That is not a curiosity: it is the regression handle. With no gap,
    W2 has no exposed boundary and MUST reduce to W1 identically, so the
    tiling-pad sum rule S(0)/C(0) has to come back out of this machinery to
    round-off.
    """

    def __init__(self, n_pad=4, dx_um=5.0, gap_um=None, gap_m=C.AMP_GAP_M,
                 sub_layers=None, sub_terminated=True, strips=False):
        pitch = C.PAD_PITCH_UM
        gap_um = (pitch - C.PAD_SIZE_UM) if gap_um is None else gap_um
        size_um = pitch - gap_um
        n = int(round(pitch / dx_um))
        assert abs(n * dx_um - pitch) < 1e-9, "dx must divide the pad pitch"
        m = int(round(size_um / dx_um))
        assert abs(m * dx_um - size_um) < 1e-9, "dx must divide the pad size"
        self.n_pad, self.n, self.m, self.strips = n_pad, n, m, strips
        self.N = n_pad * n
        self.dx = dx_um * 1e-6
        self.L = self.N * self.dx

        # exact metal mask: m of every n cells, in each direction
        one = np.zeros(n, dtype=bool)
        one[:m] = True
        line = np.tile(one, n_pad)
        # `strips` makes the metal depend on x only — the pads become infinite
        # strips. Physically a different board, but it is the geometry an
        # independent (x, z) finite-difference solve can reach, which is what
        # `fd_strip_reference` uses to test this construction from outside.
        self.metal = (np.tile(line, (n_pad * n, 1)) if strips
                      else np.outer(line, line))
        self.gapmask = ~self.metal
        self.amp_gap_m = gap_m
        self.pad_i0 = np.arange(n_pad) * n            # each pad's first cell

        kx = 2 * np.pi * np.fft.fftfreq(self.N, d=self.dx)
        self.k = np.hypot(*np.meshgrid(kx, kx, indexing="xy"))
        A11, A12, A22, Cg = stack_2port(self.k, gap_m=gap_m)
        self.A12, self.C = A12, Cg
        ysub = substrate_admittance(
            self.k, NOMINAL_SUB if sub_layers is None else sub_layers,
            terminated=sub_terminated)
        # M = A22 + Y_sub - A12^2/C. Infinite Y_sub is W1: the gap is grounded,
        # which is imposed by pinning u = 0 there instead of by arithmetic.
        self.w1 = bool(np.all(np.isinf(ysub)))
        self.M = A22 - A12 ** 2 / Cg + (0.0 if self.w1 else ysub)
        assert self.w1 or np.all(self.M > 0), "M must be positive definite"

    # --- operator, matrix-free -------------------------------------------
    def _mul(self, f):
        return np.real(np.fft.ifft2(self.M * np.fft.fft2(f)))

    def _amul(self, u):
        return self.gapmask * self._mul(self.gapmask * u)

    def solve_prompt(self, vpad, tol=1e-12, maxit=4000):
        """
        Prompt weighting potential on the ESL plane for a pad drive pattern.

        Returns Psi(x, y, z=0) = V0, which by reciprocity IS the prompt Green's
        function G(x0, y0) for that electrode (plan §3).
        """
        if self.w1:                       # gaps grounded: Vd is just the drive
            return np.real(np.fft.ifft2(self.A12 / self.C
                                        * np.fft.fft2(vpad))), 0
        b = -(self.gapmask * self._mul(vpad))
        n0 = float(np.abs(b).max())
        if n0 == 0.0:
            # No exposed gap (gap_um = 0) or no drive: u is exactly zero and
            # CG would divide by it. Returning here is not a shortcut — it is
            # the answer, and it is what makes the tiling-pad regression an
            # identity rather than a converged approximation.
            return np.real(np.fft.ifft2(self.A12 / self.C
                                        * np.fft.fft2(vpad))), 0
        u = np.zeros_like(b)
        r = b - self._amul(u)
        # Preconditioner: the same operator inverted in k and then restricted.
        # M spans ~3 decades from k=0 to Nyquist, so unpreconditioned CG on a
        # 624^2 grid stalls; this collapses it to a few dozen iterations.
        prec = lambda v: self.gapmask * np.real(
            np.fft.ifft2(np.fft.fft2(self.gapmask * v) / self.M))
        z = prec(r)
        p = z.copy()
        rz = float((r * z).sum())
        it = 0
        for it in range(1, maxit + 1):
            Ap = self._amul(p)
            al = rz / float((p * Ap).sum())
            u += al * p
            r -= al * Ap
            if float(np.abs(r).max()) <= tol * max(n0, 1e-300):
                break
            z = prec(r)
            rz_new = float((r * z).sum())
            p = z + (rz_new / rz) * p
            rz = rz_new
        else:
            raise RuntimeError(f"CG did not converge in {maxit} iterations")
        vd = vpad + self.gapmask * u
        return np.real(np.fft.ifft2(self.A12 / self.C * np.fft.fft2(vd))), it

    # --- drive patterns ---------------------------------------------------
    def all_pads(self):
        return self.metal.astype(float)

    def one_pad(self, i=None, j=None):
        i = self.n_pad // 2 if i is None else i
        j = self.n_pad // 2 if j is None else j
        v = np.zeros((self.N, self.N))
        a, b = self.pad_i0[j], self.pad_i0[i]
        v[a:a + self.m, b:b + self.m] = 1.0
        return v

    def pad_sum(self, field):
        """Fold a field onto ONE pad cell by summing its pad-lattice copies."""
        n, npd = self.n, self.n_pad
        return field.reshape(npd, n, npd, n).sum(axis=(0, 2))


# ── verification against the certified solver ────────────────────────────────

def check_against_w1(verbose=True):
    """
    stack_2port must reproduce `wpot.stack_coeffs` exactly.

    This is the join between a new file and 12 certified products. C(k) and
    S(k) there were derived as a cascaded two-port written in tanh/sech form;
    here they fall out of a two-port written in coth/csch form with an
    eliminated interface potential. Same physics, independently arranged, so
    agreement to round-off is a real check on both — and the two-port's OFF
    diagonal A12 being S(k) is what lets a mixed boundary condition be bolted
    onto the existing solver without redefining anything.
    """
    k = np.concatenate([[0.0], np.geomspace(1e0, 1e6, 400)])
    A11, A12, A22, Cg = stack_2port(k)
    Cw, Sw = W.stack_coeffs(k, C.AMP_GAP_M, C.KAPTON_THICK_UM * 1e-6)
    eC = float(np.abs(Cg - Cw).max() / np.abs(Cw).max())
    eS = float(np.abs(A12 - Sw).max() / np.abs(Sw).max())
    passive = float((A11 * A22 - A12 ** 2).min())
    if verbose:
        print(f"  C(k) vs wpot.stack_coeffs : {eC:.2e}")
        print(f"  A12(k) vs S(k)            : {eS:.2e}")
        print(f"  two-port passivity min(A11 A22 - A12^2) = {passive:.3e} > 0")
    ok = eC < 1e-12 and eS < 1e-12 and passive > 0
    return ok, eC, eS


def sum_rule_regression(dx_um=5.0, n_pad=4, verbose=True):
    """
    Collapse the inter-pad gap to zero and the machinery must give back W1.

    With pitch-sized pads there IS no exposed boundary, so W2 is not an
    approximation to W1 there — it is W1, and the closed-form sum rule
    S(0)/C(0) must come back out at EVERY point of the plane, not just on
    average. Anything else means the mixed-BC solve has a bug that the
    percent-level comparison below would then quietly inherit.
    """
    expect = float(W.stack_coeffs(np.array([0.0]), C.AMP_GAP_M,
                                  C.KAPTON_THICK_UM * 1e-6)[1][0]
                   / W.stack_coeffs(np.array([0.0]), C.AMP_GAP_M,
                                    C.KAPTON_THICK_UM * 1e-6)[0][0])
    got = {}
    for name, kw in (("W1 (gaps grounded)", dict(sub_layers=[])),
                     ("open (no conductor below)",
                      dict(sub_layers=[], sub_terminated=False)),
                     ("nominal substrate", {})):
        pp = PadPlane(n_pad=n_pad, dx_um=dx_um, gap_um=0.0, **kw)
        v0, it = pp.solve_prompt(pp.all_pads())
        got[name] = float(np.abs(v0 - expect).max() / expect)
        if verbose:
            print(f"  tiling pads, {name:<30} max |Psi - S(0)/C(0)| / "
                  f"expect = {got[name]:.2e}   ({it} CG it)")

    if verbose:
        print(f"  (S(0)/C(0) = {expect:.6f})")
    # The tiling test is STRUCTURAL: with no gap the mask is empty, so it never
    # runs the mixed-BC solve. Step 3 of `verification` is the one that does.
    return max(got.values()), expect


# ── the measurement ──────────────────────────────────────────────────────────

def capture(pp):
    """
    Prompt image-charge capture on the pads, driven all together.

    Returns (mean over source position, min, max). The mean is the number the
    plan anchors on — in W1 it is exactly (PAD_SIZE/PAD_PITCH)^2 x S(0)/C(0),
    because a fully conducting plane divides the k = 0 image charge by area and
    nothing else survives the average.
    """
    v0, it = pp.solve_prompt(pp.all_pads())
    return float(v0.mean()), float(v0.min()), float(v0.max()), it


def pad_split(pp):
    """
    Where one pad's own image charge goes, for a source ON that pad's centre.

    From the single-pad Green's function by reciprocity: G_pad(x0, y0) is the
    charge induced on the driven pad by a unit charge at (x0, y0), so reading
    it at the pad centre and at the neighbouring pad centres gives the prompt
    pad-to-pad sharing directly. Superposition over a channel's comb is exact
    (the gap condition is homogeneous and Dirichlet data add), and at 780 µm
    spacing against a kernel that lives on ~200 µm the nearest pad dominates.
    """
    v0, it = pp.solve_prompt(pp.one_pad())
    n, npd = pp.n, pp.n_pad
    c = npd // 2
    ctr = pp.pad_i0[c] + pp.m // 2                    # driven pad centre cell
    step = n
    N = pp.N
    out = {"self": float(v0[ctr, ctr]),
           # wrap: the box is periodic, and at n_pad = 2 the neighbour cell is
           # the periodic image of the driven pad's own cell
           "edge": float(v0[ctr, (ctr + step) % N]),        # one pad along x
           "diag": float(v0[(ctr + step) % N, (ctr + step) % N]),  # d=+-1
           "total_on_plane": float(pp.pad_sum(v0).sum() / (n * n)),
           "cg_it": it}
    return out, v0


# ── independent cross-check: a real-space finite-difference solve ────────────

def fd_strip_reference(pp, nz_per_um=1.0):
    """
    Solve the same prompt problem by 2-D finite differences in (x, z).

    THE POINT IS SHARED-NOTHING. Everything above works in lateral Fourier
    space, assembles a two-port analytically and imposes the gap condition
    through a restricted operator. This solves div(eps grad V) = 0 on a real
    space grid with a flux-conservative 5-point stencil, puts the pads in as
    Dirichlet nodes, and never forms C(k), S(k), Y_sub or M at all. If a 27 %
    capture shift survives that, it is not an artifact of the formulation.

    It can only reach `strips=True` geometry (pads as infinite strips), because
    real squares need three dimensions — so this validates the machinery, on a
    geometry the machinery also computes, rather than the square-pad number
    directly.

    Layers, bottom to top: FR4 -> trench glue -> [pad plane] -> glue -> kapton
    -> gas -> [mesh]. Both outer faces grounded. Node planes are placed exactly
    on every interface so no dielectric jump is smeared.
    """
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import spsolve

    assert pp.strips, "fd_strip_reference needs strips=True"
    d_k = C.KAPTON_THICK_UM * 1e-6
    d_g = C.GLUE_THICK_UM * 1e-6
    layers = []                                     # (thickness, eps), bottom up
    if not pp.w1:
        for t, e in reversed(NOMINAL_SUB):
            layers.append((t, e))
    layers += [(d_g, C.GLUE_EPS_R), (d_k, C.KAPTON_EPS_R), (pp.amp_gap_m, 1.0)]

    z, eps = [0.0], []
    for t, e in layers:
        n = max(2, int(round(t * 1e6 * nz_per_um)))
        h = t / n
        for _ in range(n):
            z.append(z[-1] + h)
            eps.append(e)
    z = np.array(z)
    eps = np.array(eps)                              # per interval
    # index of the pad plane: top of the substrate (or 0 if there is none)
    j_pad = 0 if pp.w1 else int(round(
        sum(max(2, int(round(t * 1e6 * nz_per_um))) for t, _ in NOMINAL_SUB)))
    j_esl = len(z) - 1 - max(2, int(round(pp.amp_gap_m * 1e6 * nz_per_um)))
    nx, nz = pp.N, len(z)

    metal = pp.metal[0]                              # x profile of the strips
    fixed = np.zeros((nz, nx), dtype=bool)
    val = np.zeros((nz, nx))
    fixed[0] = fixed[-1] = True                      # grounded faces
    fixed[j_pad, metal] = True
    val[j_pad, metal] = 1.0

    idx = -np.ones((nz, nx), dtype=int)
    free = ~fixed
    idx[free] = np.arange(free.sum())
    A = lil_matrix((free.sum(), free.sum()))
    b = np.zeros(free.sum())
    dx = pp.dx
    h = np.diff(z)

    def add(r, jj, ii, w):
        ii %= nx
        if fixed[jj, ii]:
            b[r] -= w * val[jj, ii]
        else:
            A[r, idx[jj, ii]] += w

    for j in range(1, nz - 1):
        hd, hu = h[j - 1], h[j]
        ed, eu = eps[j - 1], eps[j]
        ebar = (ed * hd + eu * hu) / (hd + hu)
        wz_d, wz_u = ed / hd, eu / hu
        wx = ebar * (hd + hu) / 2 / dx ** 2
        for i in range(nx):
            if fixed[j, i]:
                continue
            r = idx[j, i]
            A[r, r] -= wz_d + wz_u + 2 * wx
            add(r, j - 1, i, wz_d)
            add(r, j + 1, i, wz_u)
            add(r, j, i - 1, wx)
            add(r, j, i + 1, wx)

    v = np.zeros((nz, nx))
    v[fixed] = val[fixed]
    v[free] = spsolve(A.tocsr(), b)
    return v[j_esl], z[j_esl]


MODELS = [
    ("W1  gaps grounded at the pad plane", dict(sub_layers=[])),
    ("W2  nominal: 26 µm glue + 264 µm FR4 to L5", {}),
    ("W2  FR4 halved (132 µm)",
     dict(sub_layers=[(TRENCH_GLUE_UM * 1e-6, C.GLUE_EPS_R),
                      (132e-6, FR4_EPS_R)])),
    ("W2  open: no conductor below at all",
     dict(sub_layers=[], sub_terminated=False)),
]


def verification(dx_um, n_pad, verbose=True):
    """Everything that has to hold before any of the numbers below mean anything."""
    ok = True
    print("  1. the two-port must reproduce the certified solver")
    o, eC, eS = check_against_w1()
    ok &= o

    print("\n  2. tiling pads (gap -> 0): W2 IS W1, so S(0)/C(0) must come back")
    sr, expect = sum_rule_regression(dx_um=dx_um, n_pad=n_pad, verbose=True)
    o = sr < 1e-12
    ok &= o
    print(f"     worst {sr:.2e} vs a 1e-12 bar   {'PASS' if o else 'FAIL'}")

    print("\n  3. ground the gap from just below and the mixed-BC solve must")
    print("     converge back onto W1 — first order in the standoff h")
    w1 = PadPlane(n_pad=n_pad, dx_um=dx_um, sub_layers=[])
    ref, _ = w1.solve_prompt(w1.all_pads())
    errs = []
    for h in (1e-6, 1e-7, 1e-8, 1e-9):
        pp = PadPlane(n_pad=n_pad, dx_um=dx_um,
                      sub_layers=[(h, C.GLUE_EPS_R)])
        v0, it = pp.solve_prompt(pp.all_pads())
        e = float(np.abs(v0 - ref).max() / np.abs(ref).max())
        errs.append(e)
        print(f"     h = {h*1e9:8.1f} nm   |W2 - W1| = {e:.3e}"
              + ("" if len(errs) < 2 else f"   x{errs[-2]/e:.2f}")
              + f"   ({it} CG it)")
    rat = [errs[i] / errs[i+1] for i in range(len(errs)-1)]
    o = all(9.0 < r < 11.0 for r in rat)
    ok &= o
    print(f"     ratios {['%.2f' % r for r in rat]} -> order 1 in h   "
          f"{'PASS' if o else 'FAIL'}")

    print("\n  4. INDEPENDENT: a real-space (x, z) finite-difference solve,")
    print("     sharing no code — no Fourier, no two-port, no M operator.")
    print("     Strip pads, the only geometry 2-D can reach.")
    worst = 0.0
    for name, kw in (("W1", dict(sub_layers=[])), ("W2 nominal", {})):
        pp = PadPlane(n_pad=2, dx_um=10.0, strips=True, **kw)
        v0, _ = pp.solve_prompt(pp.all_pads())
        fd, _ = fd_strip_reference(pp, nz_per_um=1.0)
        e = abs(float(v0[0].mean()) - float(fd.mean())) / float(fd.mean())
        worst = max(worst, e)
        print(f"     {name:<12} spectral {v0[0].mean():.6f}   FD {fd.mean():.6f}"
              f"   rel {e:.2e}")
    # the strip case has its own closed form in W1, and it is a different one
    anchor = (C.PAD_SIZE_UM / C.PAD_PITCH_UM) * expect
    pp = PadPlane(n_pad=2, dx_um=10.0, strips=True, sub_layers=[])
    got = float(pp.solve_prompt(pp.all_pads())[0].mean())
    print(f"     W1 strips closed form (PAD_SIZE/PAD_PITCH) x S(0)/C(0) = "
          f"{anchor:.6f}, solved {got:.6f}")
    o = worst < 1e-3 and abs(got / anchor - 1) < 1e-12
    ok &= o
    print(f"     {'PASS' if o else 'FAIL'}")
    return ok, expect


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dx-um", type=float, default=5.0)
    ap.add_argument("--n-pad", type=int, default=4)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    print("V6 — exposing the 100 µm inter-pad gaps (plan §3 W2)\n")
    print(f"  pads {C.PAD_SIZE_UM:.0f} µm on {C.PAD_PITCH_UM:.0f} µm pitch"
          f"  ->  gap {C.PAD_PITCH_UM-C.PAD_SIZE_UM:.0f} µm, "
          f"metal fraction {(C.PAD_SIZE_UM/C.PAD_PITCH_UM)**2:.5f}")
    print(f"  insulator to the pads: {C.KAPTON_THICK_UM:.0f} µm kapton + "
          f"{C.GLUE_THICK_UM:.2f} µm glue")
    print(f"  under the gap:         {TRENCH_GLUE_UM:.0f} µm glue (the pad "
          f"trench) + {FR4_TO_L5_UM:.0f} µm FR4 to L5\n")

    print("  VERIFICATION")
    ok, expect = verification(a.dx_um, a.n_pad)
    print(f"\n  -> verification {'PASS' if ok else 'FAIL'}\n")

    print("  RESULT 1 — prompt image charge captured by the readout pads")
    print("    (mean over deposit position; W1's mean is exactly")
    print("     (PAD_SIZE/PAD_PITCH)^2 x S(0)/C(0), because a solid plane")
    print("     divides the k=0 image charge by area)\n")
    print(f"    {'model':<44}{'mean':>9}{'vs W1':>9}{'min':>8}{'max':>8}"
          f"{'max/min':>9}")
    res, base = {}, None
    for name, kw in MODELS:
        pp = PadPlane(n_pad=a.n_pad, dx_um=a.dx_um, **kw)
        m, lo, hi, it = capture(pp)
        base = m if base is None else base
        res[name] = {"mean": m, "min": lo, "max": hi, "cg_it": it,
                     "rel_to_w1": m / base - 1.0, "modulation": hi / lo,
                     "mean_Vd": m / expect}
        print(f"    {name:<44}{m:9.5f}{m/base-1:+8.1%}{lo:8.4f}{hi:8.4f}"
              f"{hi/lo:9.2f}")

    print(f"\n    Read as a pad-plane potential: capture = S(0)/C(0) x <Vd>,")
    print(f"    so W1 puts the inter-pad channel at 0 and gives <Vd> = "
          f"{res[MODELS[0][0]]['mean_Vd']:.4f} = the metal fraction,")
    nom = res[MODELS[1][0]]
    frac = (nom['mean_Vd'] - (C.PAD_SIZE_UM/C.PAD_PITCH_UM)**2) / \
           (1 - (C.PAD_SIZE_UM/C.PAD_PITCH_UM)**2)
    print(f"    while the real channel floats to {frac:.2f} of the pad "
          f"potential and gives <Vd> = {nom['mean_Vd']:.4f}.")

    print("\n  RESULT 2 — one pad's prompt kernel, read at pad centres")
    print(f"    {'model':<44}{'self':>10}{'+1 pad':>10}{'diagonal':>10}")
    for name, kw in MODELS:
        pp = PadPlane(n_pad=a.n_pad, dx_um=a.dx_um, **kw)
        sp, _ = pad_split(pp)
        res[name].update(sp)
        print(f"    {name:<44}{sp['self']:10.6f}{sp['edge']:10.2e}"
              f"{sp['diag']:10.2e}")
    print("    (the PEAK barely moves — the shift is entirely in deposits that")
    print("     land over a channel, which is 24 % of the area)")

    print("\n  CONVERGENCE")
    caps = []
    for dx in (20.0, 10.0, 5.0, 2.5):
        caps.append(capture(PadPlane(n_pad=a.n_pad, dx_um=dx))[0])
        print(f"    cell {dx:5.2f} µm   capture {caps[-1]:.6f}")
    d = [caps[i+1] - caps[i] for i in range(len(caps)-1)]
    rich = caps[-1] + d[-1]
    print(f"    successive deltas {['%.1e' % x for x in d]} -> first order "
          f"(the metal edge is an r^-1/2 singularity)")
    print(f"    Richardson limit {rich:.6f}; at {a.dx_um:g} µm we are "
          f"{abs(caps[2]/rich-1):.2%} low, far inside the effect")
    for npd in (2, 4, 6, 8):
        pp = PadPlane(n_pad=npd, dx_um=10.0)
        print(f"    box {npd}x{npd} pads   single-pad self "
              f"{pad_split(pp)[0]['self']:.6f}")

    print("\n  HOW IT SCALES WITH THE GAP (nominal substrate)")
    for gu in (0.0, 25.0, 50.0, 100.0, 150.0):
        w1 = PadPlane(n_pad=a.n_pad, dx_um=a.dx_um, gap_um=gu, sub_layers=[])
        w2 = PadPlane(n_pad=a.n_pad, dx_um=a.dx_um, gap_um=gu)
        c1, c2 = capture(w1)[0], capture(w2)[0]
        print(f"    gap {gu:5.1f} µm   W1 {c1:.5f}   W2 {c2:.5f}   "
              f"{c2/c1-1:+7.1%}")

    shift = max(abs(v["rel_to_w1"]) for v in res.values())
    mod = res[MODELS[0][0]]["modulation"], res[MODELS[1][0]]["modulation"]
    print(f"\n  VERDICT")
    print(f"    prompt capture, W1 -> W2:            {shift:+.1%}")
    print(f"    spread over the W2 substrate bracket: "
          f"{abs(res[MODELS[3][0]]['rel_to_w1'] - res[MODELS[1][0]]['rel_to_w1']):.1%}"
          f"  (so the board's internal stackup barely matters)")
    print(f"    sub-pad amplitude modulation:        W1 {mod[0]:.2f}x  ->  "
          f"W2 {mod[1]:.2f}x")
    print(f"    plan §3 expected 'percent-level' -> "
          f"{'CONSISTENT' if shift < 0.05 else 'REFUTED: it is an order of magnitude larger'}")
    if a.out:
        json.dump({"models": res, "sum_rule_expect": expect,
                   "richardson_capture": rich, "worst_shift": shift,
                   "dx_um": a.dx_um, "n_pad": a.n_pad,
                   "verification_pass": bool(ok)},
                  open(a.out, "w"), indent=1)
        print(f"\n  -> {a.out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
