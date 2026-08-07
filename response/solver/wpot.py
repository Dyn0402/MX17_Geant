"""
wpot.py — S1, the semi-spectral dynamic weighting-potential solver.

Stage S1 of design/RESPONSE_SIM_PLAN.md §3. Extended Ramo-Shockley for a
resistive electrode (Riegler, JINST 11 (2016) P11002, arXiv:1602.07949).


THE PHYSICS, DERIVED
====================

Stack (model W1), z measured from the ESL surface, +z toward the mesh:

    z = +g      grounded micromesh plane
    0 < z < g   gas, eps = 1                       (g = 150 µm, CONFIRMED)
    z = 0       the ESL sheet: a zero-thickness layer of surface conductivity
                sigma_s(x) = 1/rho_s on the 550 µm strips, 0 in the 250 µm
                gaps; periodic in x with period 800 µm, uniform in y
    -d < z < 0  kapton, eps_r = 3.5                (d unknown -> scan)
    z = -d      the segmented pad plane; a Dirichlet pattern V_pad(x, y)

To get the weighting potential of readout channel n we hold channel n's pads
at V_w for t > 0 and everything else (other pads, mesh) at ground, then follow
the quasi-static relaxation.

In each homogeneous layer the potential obeys Laplace, so a single lateral
Fourier mode of magnitude k behaves as exp(+/- k z). Applying V(z=g) = 0 above
and V(z=-d) = V_pad below, and writing V0 for the potential ON the sheet, the
displacement field on each side of the sheet is

    D_z(0+) = eps0 k coth(k g) V0
    D_z(0-) = -eps0 eps_r k [ V0 coth(k d) - V_pad / sinh(k d) ]

so the free surface charge on the sheet, sigma = D_z(0+) - D_z(0-), is

    sigma = C(k) V0 - S(k) V_pad
    C(k)  = eps0 k coth(k g) + eps0 eps_r k coth(k d)     <- total areal C
    S(k)  = eps0 eps_r k / sinh(k d)                      <- drive coupling

Charge conservation on the sheet, with surface current J = -sigma_s grad_T V,

    d(sigma)/dt = grad_T . [ sigma_s(x) grad_T V0 ]                        (*)

For t > 0 the drive is constant, so (*) becomes a linear ODE for V0 whose
initial condition comes from the fact that the sheet charge cannot jump:
sigma(0-) = 0 gives

    V0(0+) = S(k)/C(k) * V_w        <- the PROMPT (pure-dielectric) solution

and the sheet relaxes from there toward V0 = 0, i.e. toward behaving as a
grounded conductor. That decay IS the signal: a charge sitting on the ESL is
progressively screened by the sheet's own conduction to ground.

Uniform sheet (sigma_s constant) — used for validation V1
---------------------------------------------------------
grad.(sigma_s grad) -> -sigma_s k^2, the modes decouple, and

    V0(k, t) = V_w S(k)/C(k) exp( -t k^2 sigma_s / C(k) )

At small k, C(k) -> eps0 (1/g + eps_r/d) = c', so the decay is exp(-k^2 D t)
with D = 1/(rho_s c') — an exact Gaussian diffusion kernel of width
sigma^2(t) = 2 D t = 2t/(rho_s c'). That is precisely the closed-form limit the
plan quotes in §3, and validate.py checks the solver against it.

Patterned sheet — the Bloch structure
-------------------------------------
sigma_s(x) has period p = 800 µm, so in Fourier it couples lateral modes whose
x-wavenumbers differ by a multiple of 2*pi/p. Expanding V0 in modes and
projecting (*) gives, per k_y and per Bloch family,

    C_a dV_a/dt = -sum_b M_ab V_b ,   M_ab = sigma_hat[j_a - j_b] (kx_a kx_b + ky^2)

with C_a = C(k_a) diagonal. This solver works on a periodic box whose length is
an exact multiple of BOTH the 800 µm ESL pitch and the 780 µm pad pitch — the
31.2 mm superperiod of plan §1 — which makes the coupling *exactly* block
diagonal by residue class of the mode index. No truncation parameter M is
needed at all: the (2M+1)-truncation the plan describes is replaced by the
grid's own Nyquist cut, which is the same approximation stated in a form that
converges by refining one number (N_x).

M is Hermitian and positive semi-definite, so with C diagonal and positive the
system is solved exactly by symmetric eigendecomposition of

    B = C^-1/2 M C^-1/2 ,      V(t) = C^-1/2 Q exp(-w t) Q^H C^1/2 V(0+)

— exact exponentials, no time-stepping error, exactly as §3 requires.

M's null space is physically meaningful and must NOT be regularised away: a
mode supported entirely on the insulating gaps has zero eigenvalue and never
decays. That is the real statement that charge landing in a 250 µm inter-strip
gap has no DC path to ground, and it is a large part of why the X view (across
the strips) shares differently from the Y view (along them).

Reciprocity
-----------
Psi_n(x0, y0, z=0, t) / V_w is exactly G_n(x0, y0, t), the charge induced on
channel n by a unit point charge sitting on the ESL at (x0, y0) since t = 0
(plan §3). So the z = 0 slice this solver returns is already the Green's
function the digitizer wants; no separate reciprocity step is needed.
"""

from __future__ import annotations

import numpy as np

from ..common import constants as C


# ── Layer transfer coefficients, with safe limits ────────────────────────────

def _coth(x):
    """coth(x), stable at both ends. x >= 0."""
    out = np.ones_like(x, dtype=float)
    small = x < 1e-6
    mid = (~small) & (x < 30.0)
    # x -> 0 : coth x = 1/x + x/3 ; x -> inf : coth x -> 1 (the default above)
    out[small] = 1.0 / np.maximum(x[small], 1e-300) + x[small] / 3.0
    out[mid] = 1.0 / np.tanh(x[mid])
    return out


def _csch(x):
    """1/sinh(x), stable at both ends: underflows to 0 for large x, as it should."""
    out = np.zeros_like(x, dtype=float)
    small = x < 1e-6
    mid = (~small) & (x < 700.0)
    out[small] = 1.0 / np.maximum(x[small], 1e-300) - x[small] / 6.0
    out[mid] = 1.0 / np.sinh(x[mid])
    return out


def stack_coeffs(k, gap_m, d_kapton_m, eps_r=C.KAPTON_EPS_R):
    """
    Areal capacitance C(k) and drive coupling S(k) of the W1 stack [F/m^2].

    C(k) = eps0 k coth(k g) + eps0 eps_r k coth(k d)
    S(k) = eps0 eps_r k / sinh(k d)

    At k = 0 these reduce to C = eps0(1/g + eps_r/d) = c' and S = eps0 eps_r/d.
    """
    k = np.asarray(k, dtype=float)
    kg = k * gap_m
    kd = k * d_kapton_m
    c_gas = C.EPS0 * np.where(kg < 1e-6, 1.0 / gap_m, k * _coth(kg))
    c_kap = C.EPS0 * eps_r * np.where(kd < 1e-6, 1.0 / d_kapton_m, k * _coth(kd))
    s_drv = C.EPS0 * eps_r * np.where(kd < 1e-6, 1.0 / d_kapton_m, k * _csch(kd))
    return c_gas + c_kap, s_drv


# ── Conductivity pattern ─────────────────────────────────────────────────────

def esl_sigma_profile(x_m, rho_s_ohm_sq, width_m=C.ESL_WIDTH_M,
                      pitch_m=C.ESL_PITCH_M, phase_m=0.0, uniform=False):
    """
    Surface conductivity sigma_s(x) [S/sq] of the ESL strip pattern.

    `phase_m` shifts the strip pattern relative to x = 0 — this is the
    screen-print registration, which is unknown (plan §12 item 3) and is
    therefore an input, not a constant.

    `uniform=True` replaces the pattern with a uniform sheet of the SAME
    conductivity (not the same average) — that is the V1 reference case, where
    the closed form applies.
    """
    sig = 1.0 / rho_s_ohm_sq
    if uniform:
        return np.full_like(np.asarray(x_m, dtype=float), sig)
    u = np.mod(np.asarray(x_m, dtype=float) - phase_m, pitch_m)
    # Strip occupies [0, width); gap occupies [width, pitch).
    return np.where(u < width_m, sig, 0.0)


# ── The solver ───────────────────────────────────────────────────────────────

class WeightingSolver:
    """
    Dynamic weighting potential on a periodic box commensurate with both
    lateral pitches.

    Parameters
    ----------
    n_super : how many 31.2 mm superperiods the x box spans. 1 is enough for
              the kernel itself (the response decays in a few mm); use more
              only to check that the box is not wrapping.
    nx      : x samples. MUST be divisible by 39*n_super so that the ESL
              period is an exact number of samples and the Bloch blocks are
              exact. 3120 gives a 10 µm step, which is the resolution plan §2
              asks for near the strip edges.
    ly_m,ny : the y box. y is uniform in the model, so its modes decouple
              completely and ny costs almost nothing.
    """

    def __init__(self, rho_s_ohm_sq, d_kapton_m, gap_m=C.AMP_GAP_M,
                 eps_r=C.KAPTON_EPS_R, n_super=1, nx=3120,
                 ly_m=0.0512, ny=1024, phase_m=0.0, uniform_sheet=False,
                 tau_drain_s=None):
        n_esl = C.N_ESL_PER_SUPER * n_super          # ESL periods in the box
        if nx % n_esl:
            raise ValueError(
                f"nx={nx} must be divisible by 39*n_super={n_esl} so the "
                "800 µm ESL period is an exact number of samples")
        self.rho_s = rho_s_ohm_sq
        self.d_k = d_kapton_m
        self.gap = gap_m
        self.eps_r = eps_r
        self.uniform_sheet = uniform_sheet

        # ── The global drain (assumption A1, plan §3 "Boundary/drain") ───────
        # Everything above conserves charge on the sheet: a mode with no
        # in-plane gradient, and a strip that is isolated in x by its gaps, has
        # nowhere for its charge to go. The box is periodic, so there is no
        # ground at infinity either. Something must connect the sheet to the
        # HV/ground bus, and that is what this term is.
        #
        # Modelled as a uniform leak conductance per unit area, g_leak, which
        # adds -g_leak V to the right-hand side of charge conservation. With
        # sigma = C(k) V that gives every mode an extra decay rate g_leak/C(k).
        # g_leak is fixed by requiring the k -> 0 rate to be exactly
        # 1/tau_drain, i.e. g_leak = c'/tau_drain.
        #
        # tau_drain_s=None means NO drain: pure charge conservation on the
        # sheet. That is the honest default, because A1 as literally stated
        # (strips grounded at their two y-ends) does NOT reduce to a uniform
        # leak — see drain_time_from_strip_ends() and validate.py V4.
        self.tau_drain = tau_drain_s
        self.g_leak = (C.c_prime(gap_m, d_kapton_m, eps_r) / tau_drain_s
                       if tau_drain_s else 0.0)

        self.lx = C.SUPERPERIOD_M * n_super
        self.nx, self.ny, self.ly = nx, ny, ly_m
        self.n_esl = n_esl

        self.x = np.arange(nx) * (self.lx / nx)
        self.y = (np.arange(ny) - ny // 2) * (ly_m / ny)

        self.kx = 2 * np.pi * np.fft.fftfreq(nx, d=self.lx / nx)
        self.ky = 2 * np.pi * np.fft.fftfreq(ny, d=ly_m / ny)

        # Fourier coefficients of sigma_s(x). For the patterned sheet these are
        # nonzero only at multiples of n_esl, which is exactly what makes the
        # Bloch blocks exact rather than truncated.
        sig_x = esl_sigma_profile(self.x, rho_s_ohm_sq, phase_m=phase_m,
                                  pitch_m=C.ESL_PITCH_M,
                                  uniform=uniform_sheet)
        self.sigma_x = sig_x
        self.sigma_hat = np.fft.fft(sig_x) / nx

    # ── internals ────────────────────────────────────────────────────────────

    def _blocks(self):
        """Mode-index members of each Bloch residue class (list of arrays)."""
        n_a = self.nx // self.n_esl
        return [np.arange(c, self.nx, self.n_esl) for c in range(self.n_esl)], n_a

    def _solve_block(self, idx, ky):
        """
        Eigendecomposition for one Bloch block at one k_y.

        Returns (Cvec, w, Q) with the propagator
            V(t) = C^-1/2 Q diag(e^-wt) Q^H C^1/2 V(0+)
        """
        kx = self.kx[idx]
        k = np.hypot(kx, ky)
        Cvec, _ = stack_coeffs(k, self.gap, self.d_k, self.eps_r)

        # M_ab = sigma_hat[j_a - j_b] * (kx_a kx_b + ky^2)
        dj = (idx[:, None] - idx[None, :]) % self.nx
        M = self.sigma_hat[dj] * (np.outer(kx, kx) + ky ** 2)

        # Hermitian by construction for a real sigma_s(x); symmetrise to kill
        # round-off asymmetry before eigh.
        M = 0.5 * (M + M.conj().T)
        inv_sqrt_C = 1.0 / np.sqrt(Cvec)
        B = (inv_sqrt_C[:, None] * M) * inv_sqrt_C[None, :]
        # Uniform leak to the bus. The full operator is
        # C dV/dt = -(M + g_leak I) V, so after the C^-1/2 similarity the leak
        # appears as g_leak C^-1, which is diagonal in the MODE basis. It has
        # to go in here, before the eigendecomposition — adding it to the
        # eigenvalues afterwards would be adding a non-constant diagonal in the
        # wrong basis.
        if self.g_leak:
            B = B + np.diag(self.g_leak / Cvec)
        w, Q = np.linalg.eigh(B)
        # Physical eigenvalues are >= 0 (M is PSD); clip round-off negatives so
        # nothing can grow exponentially. Genuine zeros are KEPT — they are the
        # charge stranded in the insulating inter-strip gaps.
        w = np.maximum(w, 0.0)
        return Cvec, w, Q

    # ── public API ───────────────────────────────────────────────────────────

    def prompt_potential(self, vpad_xy):
        """
        The t = 0+ (pure-dielectric) weighting potential on the sheet.

        V0(0+) = S(k)/C(k) * V_pad(k) — the sheet has had no time to conduct,
        so it acts as a plain dielectric interface and every mode decouples.
        """
        vhat = np.fft.fft2(np.asarray(vpad_xy, dtype=float))
        KX, KY = np.meshgrid(self.kx, self.ky, indexing="xy")
        k = np.hypot(KX, KY)
        Cv, Sv = stack_coeffs(k, self.gap, self.d_k, self.eps_r)
        return np.real(np.fft.ifft2(vhat * (Sv / Cv)))

    def solve(self, vpad_xy, times_s, progress=False):
        """
        Psi(x, y, z=0, t) for a drive pattern V_pad(x, y) held from t = 0.

        Parameters
        ----------
        vpad_xy : (ny, nx) array, the Dirichlet pattern on the pad plane
        times_s : times [s] at which to report; t = 0 gives the prompt solution

        Returns
        -------
        (nt, ny, nx) real array. By reciprocity this IS G_n(x0, y0, t) when
        vpad_xy is channel n's pad pattern and V_w = 1.
        """
        times = np.atleast_1d(np.asarray(times_s, dtype=float))
        vhat = np.fft.fft2(np.asarray(vpad_xy, dtype=float))   # (ny, nx)
        out_hat = np.zeros((len(times), self.ny, self.nx), dtype=complex)

        blocks, _ = self._blocks()
        for iy, ky in enumerate(self.ky):
            for idx in blocks:
                Cvec, w, Q = self._solve_block(idx, ky)
                _, Sv = stack_coeffs(np.hypot(self.kx[idx], ky),
                                     self.gap, self.d_k, self.eps_r)
                v0 = (Sv / Cvec) * vhat[iy, idx]        # prompt amplitude
                # u = C^1/2 v0, propagate in the eigenbasis, come back
                u0 = np.sqrt(Cvec) * v0
                coef = Q.conj().T @ u0
                # All requested times at once: exp(-w_a t_k) is a (nb, nt)
                # outer product, so the whole propagation is one matmul. Same
                # arithmetic as the per-time loop, ~50x less python overhead,
                # which matters once nt = 61 and ny = 512.
                u = Q @ (np.exp(-np.outer(w, times)) * coef[:, None])
                out_hat[:, iy, idx] = (u / np.sqrt(Cvec)[:, None]).T
            if progress and iy % 128 == 0:
                print(f"  ky {iy}/{self.ny}", flush=True)

        return np.real(np.fft.ifft2(out_hat, axes=(1, 2)))

    # ── analytic reference (uniform sheet) ───────────────────────────────────

    def solve_uniform_analytic(self, vpad_xy, times_s):
        """
        Closed form for a UNIFORM sheet — the V1 reference.

        V0(k, t) = S(k)/C(k) * V_pad(k) * exp( -k^2 t / (rho_s C(k)) )

        Valid only when the sheet really is uniform; it is the thing the full
        Bloch machinery must reproduce when handed a uniform sigma_s.
        """
        times = np.atleast_1d(np.asarray(times_s, dtype=float))
        vhat = np.fft.fft2(np.asarray(vpad_xy, dtype=float))
        KX, KY = np.meshgrid(self.kx, self.ky, indexing="xy")
        k = np.hypot(KX, KY)
        Cv, Sv = stack_coeffs(k, self.gap, self.d_k, self.eps_r)
        sig = 1.0 / self.rho_s
        out = np.empty((len(times), self.ny, self.nx))
        for it, t in enumerate(times):
            out[it] = np.real(np.fft.ifft2(
                vhat * (Sv / Cv) * np.exp(-k ** 2 * sig * t / Cv)))
        return out


# ── drive patterns ───────────────────────────────────────────────────────────

def drain_time_from_strip_ends(rho_s_ohm_sq, strip_len_m=0.412,
                               gap_m=C.AMP_GAP_M,
                               d_kapton_m=C.KAPTON_THICK_UM_DEFAULT * 1e-6,
                               eps_r=C.KAPTON_EPS_R):
    """
    Drain time constant implied by assumption A1 taken literally.

    A1 (plan §1) is that the ESL strips are tied to the HV/ground bus at BOTH
    y-ends of the active area. A strip of length L grounded at both ends and
    obeying the sheet diffusion equation has slowest mode sin(pi y / L), so

        tau = 1 / (D k^2) = L^2 / (pi^2 D) ,    D = 1/(rho_s c')

    Returned in seconds. Compare against the measured global drain
    tau_g = 5.3-7.3 µs (nTof_x17 rc_line_step2.py) — that comparison is the
    plan's V4 test, and it is the one that decides whether A1 survives.
    """
    D = 1.0 / (rho_s_ohm_sq * C.c_prime(gap_m, d_kapton_m, eps_r))
    return strip_len_m ** 2 / (np.pi ** 2 * D)


def effective_ground_spacing(tau_g_s, rho_s_ohm_sq, **kw):
    """
    Inverse of the above: the grounding pitch L that WOULD give a measured
    drain time tau_g. If this comes out far below the 412 mm active width,
    A1 cannot be the whole story.
    """
    D = 1.0 / (rho_s_ohm_sq * C.c_prime(**kw))
    return np.sqrt(tau_g_s * np.pi ** 2 * D)


def pad_pattern(solver, x0_m, y0_m, size_m=C.PAD_SIZE_UM * 1e-6, v=1.0):
    """A single square pad held at `v`, centred at (x0, y0)."""
    X, Y = np.meshgrid(solver.x, solver.y, indexing="xy")
    dx = np.abs(((X - x0_m + solver.lx / 2) % solver.lx) - solver.lx / 2)
    return np.where((dx <= size_m / 2) & (np.abs(Y - y0_m) <= size_m / 2),
                    v, 0.0)


def strip_pattern_y(solver, x0_m, size_m=C.PAD_SIZE_UM * 1e-6, v=1.0):
    """
    A whole column of pads at x0, running the full length of y.

    This is one X-measuring channel in the limit where the column is long
    compared with the response — the cheap, exact-in-y case, because a
    y-independent drive excites only k_y = 0.
    """
    X = np.broadcast_to(solver.x, (solver.ny, solver.nx))
    dx = np.abs(((X - x0_m + solver.lx / 2) % solver.lx) - solver.lx / 2)
    return np.where(dx <= size_m / 2, v, 0.0)
