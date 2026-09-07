#!/usr/bin/env python3
"""
zextend.py — Psi_n at z > 0, from the S1 surface solve, exactly.

WHY THIS FILE EXISTS. design/RESPONSE_SIM_PLAN.md §2 specifies an S1 product
`wpot_<ch>_<params>.npz` carrying `psi[t,z,y,x]` with "8-16 z-slices to the
mesh", and §7 step 6 (T10, the slow path) needs those slices because the ions
carry ~91 % of the induced charge while climbing the whole 150 µm gap. That
product was never produced: `WeightingSolver.solve()` returns the z = 0 plane
and nothing else, and `greens_comb_*.npz` is the z = 0 comb kernel by
construction. So T10 read as "blocked on re-solving S1 with z retained".

It is not, and no re-solve is needed. In model W1 the gas region 0 < z < g is
charge-free, bounded below by the ESL at a known potential and above by a
GROUNDED plane at z = +g. Quasi-statically the potential there obeys Laplace at
every instant, so each lateral Fourier mode is a pure exponential pair pinned by
those two boundary values:

    Psi(k, z, t) = Psi(k, 0, t) * sinh(k (g - z)) / sinh(k g)                (1)

    dPsi/dz(k, z, t) = -k * Psi(k, 0, t) * cosh(k (g - z)) / sinh(k g)       (2)

(1) is exact — not an interpolation, not a fit — given only that the mesh is a
grounded plane (assumption W1, whose own error bar is V5's exp(-2*pi*g/p) weave
ripple, 7.8e-7) and that the gas holds no space charge on the weighting
problem's own terms. Every t is independent because the relation involves no
time derivative: the sheet's dynamics live entirely in Psi(k, 0, t), which S1
already solved.

WHAT (1) SAYS PHYSICALLY, AND WHY IT IS THE WHOLE OF T10'S PHYSICS.

At k -> 0 the factor becomes (g - z)/g = 1 - z/g, which is exactly the
parallel-plate longitudinal weighting `response/digitizer/ions.py` already uses,
so the existing charge split f_ion = 1 - z/g is the k = 0 limit of this and is
recovered identically. At large k the factor goes as exp(-k z): the kernel is
LOW-PASS FILTERED with a cutoff that tightens as the ion climbs. In real space
that is a broadening of width ~z. So Psi_n does not merely shrink with height,
it SPREADS, and the two effects do not factorise into "one number times the
surface shape" — which is precisely the approximation the fast path makes.

The consequence the fast path cannot represent: the ion's share of the induced
charge is Psi_n(z_birth)/Psi_n(0), and because Psi_n is filtered rather than
scaled, that ratio is DIFFERENT FOR EVERY CHANNEL. A peaked d = 0 kernel loses
its high-k content and so loses more than 1 - z/g; a neighbour kernel that lives
at low k loses less, or gains. The fast path applies one global f_ion to all of
them.

RANGE OF VALIDITY. Eq. (1) describes the gas gap only, 0 <= z <= g. It says
nothing about z < 0 (inside the dielectric stack), which no drifting charge
occupies, and nothing about the drift region above the mesh, where the mesh is
not transparent to the weighting field of a pad.

    python3 -m response.solver.zextend        # self-test, incl. an FD check
"""

from __future__ import annotations

import numpy as np

from ..common import constants as C


# ── The two z factors ────────────────────────────────────────────────────────

def _exp_form(k, z, gap):
    """
    Shared stable rewrite of both factors.

    sinh(k(g-z))/sinh(kg) overflows if evaluated literally: the S1 x grid is
    10 µm, so k_max = pi/10 µm and k_max * g = 47, and cosh of that squared is
    already 1e41. Factor out exp(-kz) instead,

        sinh(k(g-z))/sinh(kg) = e^{-kz} (1 - e^{-2k(g-z)}) / (1 - e^{-2kg})

    which has no growing exponential anywhere and is accurate at both ends.
    Returns (e^{-kz}, u = e^{-2k(g-z)}, w = e^{-2kg}).
    """
    k = np.asarray(k, dtype=float)
    return np.exp(-k * z), np.exp(-2.0 * k * (gap - z)), np.exp(-2.0 * k * gap)


def psi_z_factor(k, z, gap=C.AMP_GAP_M):
    """
    Psi(k, z) / Psi(k, 0) = sinh(k(g-z)) / sinh(k g).

    At k -> 0 this is 1 - z/g (the parallel-plate weighting); the series is used
    below kg = 1e-6, where the ratio of two vanishing sinh loses all precision.
    Exactly 1 at z = 0 and exactly 0 at z = g, both by construction.
    """
    k = np.asarray(k, dtype=float)
    e, u, w = _exp_form(k, z, gap)
    small = k * gap < 1e-6
    out = np.where(small, 0.0, e * (1.0 - u) / np.where(small, 1.0, 1.0 - w))
    return np.where(small, (gap - z) / gap, out)


def dpsi_dz_factor(k, z, gap=C.AMP_GAP_M):
    """
    (dPsi/dz)(k, z) / Psi(k, 0) = -k cosh(k(g-z)) / sinh(k g).

    This is what the moving-charge (extended Ramo) integral actually needs, and
    it is given analytically rather than by differencing `psi_z_factor` on a z
    grid: the factor varies on the scale 1/k, which for the finest modes is
    3 µm, and a z grid fine enough to difference that would be 50x the one the
    plan asks for. At k -> 0 it is -1/g, the flat-weighting value.
    """
    k = np.asarray(k, dtype=float)
    e, u, w = _exp_form(k, z, gap)
    small = k * gap < 1e-6
    out = np.where(small, 0.0,
                   -k * e * (1.0 + u) / np.where(small, 1.0, 1.0 - w))
    return np.where(small, -1.0 / gap, out)


# ── Applying them to a stored surface kernel ─────────────────────────────────

class ZSlicer:
    """
    Lift a stored z = 0 kernel slab G[t, y, x] to any height in the gap.

    The lift is a multiplication in lateral Fourier space, so the x axis must be
    the solver's own periodic 31.2 mm box (it is, exactly) and the y axis must
    either be periodic (the 1.56 mm X box — also exact) or carry enough margin
    that the wrap is negligible. The filter's real-space range is ~z <= 150 µm,
    so a y margin of a few mm puts the wrap error below any float32 the product
    is stored in; `y_margin_m` is checked against `z_max` on construction rather
    than assumed.

    Parameters
    ----------
    G      : (nt, ny, nx) float array, the z = 0 kernel
    x, y   : the axes of G [m]. x must span exactly one superperiod.
    periodic_y : True for the 1.56 mm X box, where the periodic images ARE the
                 rest of the comb and the FFT is exact. False for the tall Y
                 box, where the margin argument above applies instead.
    """

    def __init__(self, G, x, y, gap=C.AMP_GAP_M, periodic_y=False,
                 y_margin_m=2e-3, hat_dtype=np.complex128):
        self.G = np.asarray(G)
        self.hat_dtype = hat_dtype
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
        self.gap = gap
        self.periodic_y = periodic_y
        self.y_margin_m = y_margin_m

        nt, ny, nx = self.G.shape
        assert len(self.x) == nx and len(self.y) == ny, "axes do not match G"
        lx = nx * (self.x[1] - self.x[0])
        ly = ny * (self.y[1] - self.y[0])
        kx = 2 * np.pi * np.fft.fftfreq(nx, d=lx / nx)
        ky = 2 * np.pi * np.fft.fftfreq(ny, d=ly / ny)
        self.k = np.hypot(*np.meshgrid(kx, ky, indexing="xy"))
        self._Ghat = None                       # built lazily, it is large

    def _hat(self):
        if self._Ghat is None:
            # complex128 by default. The product is float32, so this buys no
            # accuracy in the INPUT — but the reductions below sum ~1e5 terms,
            # and in complex64 that accumulation alone costs 3e-5 of the
            # transform's largest amplitude. Far channels carry ~1e-3 of it, so
            # in complex64 their kernels would be quoted to ~3 %. `hat_dtype`
            # exists to trade that back for memory on a slab that will not fit.
            self._Ghat = np.fft.fft2(self.G, axes=(1, 2)).astype(self.hat_dtype)
        return self._Ghat

    def _check_margin(self, z_max):
        """
        Wrap-around bound for the non-periodic y box.

        The filter is a low-pass whose real-space tail falls as exp(-r/z) at
        worst, so contamination of an interior point from the wrapped copy is
        bounded by exp(-margin / z_max). This raises rather than warns: silently
        filtering a windowed slab is exactly the "runs cleanly, plausible but
        wrong" failure the plan's audit keeps finding.
        """
        if self.periodic_y or z_max <= 0:
            return
        bound = np.exp(-self.y_margin_m / z_max)
        if bound > 1e-9:
            raise ValueError(
                f"y margin {self.y_margin_m*1e3:.2f} mm is too small for "
                f"z_max = {z_max*1e6:.0f} µm (wrap bound {bound:.1e}); widen "
                "the y window or lower z_max")

    def free(self):
        """Drop the transform and the source slab; probes already built stay
        valid, but can no longer be `certify`d against this slicer."""
        self._Ghat = None
        self.G = np.empty((0, 0, 0))

    def slab(self, z, deriv=False):
        """G (or dG/dz) on the full (nt, ny, nx) grid at one height z."""
        self._check_margin(z)
        f = (dpsi_dz_factor if deriv else psi_z_factor)(self.k, z, self.gap)
        out = np.fft.ifft2(self._hat() * f.astype(np.complex64), axes=(1, 2))
        return np.real(out)

    def at_exact(self, iy, ix, z_list, deriv=False):
        """
        G (or dG/dz) at ONE (y, x) sample, for a list of heights: (nt, nz).

        The direct definition — reduce over every k with the sample's own
        phase. Correct and slow (one pass over the whole transform per height),
        so it exists to CERTIFY `RadialProbe` rather than to be used.
        """
        z_list = np.atleast_1d(np.asarray(z_list, dtype=float))
        self._check_margin(float(z_list.max()))
        nt, ny, nx = self.G.shape
        a = self._hat() * self._phase(iy, ix)               # (nt, ny, nx)
        fn = dpsi_dz_factor if deriv else psi_z_factor
        out = np.empty((nt, len(z_list)))
        for j, z in enumerate(z_list):
            f = fn(self.k, z, self.gap)
            out[:, j] = np.real(np.einsum("tyx,yx->t", a, f)) / (ny * nx)
        return out

    def _phase(self, iy, ix):
        """e^{i k . r} for the (iy, ix) grid sample."""
        nt, ny, nx = self.G.shape
        return np.outer(np.exp(2j * np.pi * np.fft.fftfreq(ny) * iy),
                        np.exp(2j * np.pi * np.fft.fftfreq(nx) * ix))

    def probe(self, iy, ix, dk_g=1e-3):
        """A `RadialProbe` at one grid sample (see that class)."""
        return RadialProbe(self, iy, ix, dk_g=dk_g)


class RadialProbe:
    """
    One (y, x) sample's kernel, pre-reduced so any height costs a matvec.

    WHY. The trajectory integral wants dPsi/dz at O(700) heights per birth
    height per channel, and `ZSlicer.at_exact` pays a full pass over a 1 M-point
    transform for each one. But BOTH z factors depend on k only through |k|, so
    the phase-weighted transform can be collapsed onto a 1-D |k| axis ONCE, and
    every height after that is a (nt, nbin) x (nbin,) product. Nothing about the
    physics is approximated; the only error is putting neighbouring |k| into a
    common bin, and it is bounded and checked below.

    THE BIN WIDTH IS SET BY THE PHYSICS, NOT BY TASTE. Both factors vary as
    exp(-k z) with z <= g, so across a bin of width dk the factor moves by at
    most dk * g in relative terms. `dk_g` IS that product. The errors are signed
    at random across ~1/dk_g bins, so the realised error comes out near
    dk_g * sqrt(dk_g) rather than dk_g — measured, not assumed: 4e-3 gives
    4.9e-5 and 1e-3 gives ~1e-5 of the largest channel's kernel.

    The default is 1e-3 rather than the cheaper 4e-3 because the SMALLEST
    number this study reports is a 1e-4 per-channel charge residual, and a
    numerical floor of 4.9e-5 underneath it would leave that residual half
    made of binning. `certify_set` measures the floor on every run so it can
    never drift back under a result.

    k = 0 GETS ITS OWN BIN. It carries the DC term, which is the largest single
    amplitude in the transform and the one whose factor (1 - z/g, the flat
    weighting) the whole comparison is anchored on. Averaged into a bin with its
    neighbours it would acquire a spurious k_bar and shift the flat-weighting
    limit by ~dk*g — small, but exactly on the quantity being tested.
    """

    def __init__(self, slicer, iy, ix, dk_g=1e-3):
        # Copy what is needed rather than holding the slicer: a probe is a few
        # MB, the transform behind it is ~1 GB at ny = 1024, and keeping the
        # parent alive through `self.s` pinned every one of them. Fourteen
        # probes then OOM-killed the process on a 16 GB host. `src` is the
        # only path back and `ZSlicer.free()` severs it deliberately.
        self.src = slicer
        self.gap = slicer.gap
        self._margin = (None if slicer.periodic_y else slicer.y_margin_m)
        self.iy, self.ix = iy, ix
        nt, ny, nx = slicer.G.shape
        kf = slicer.k.ravel()
        dk = dk_g / slicer.gap
        # bin 0 is exactly k == 0; bin j >= 1 is [(j-1)dk, j dk)
        idx = np.where(kf > 0, np.floor(kf / dk).astype(np.int64) + 1, 0)
        nb = int(idx.max()) + 1
        cnt = np.bincount(idx, minlength=nb).astype(float)
        self.kbar = np.bincount(idx, weights=kf, minlength=nb) / np.maximum(cnt, 1)
        a = (slicer._hat() * slicer._phase(iy, ix)).reshape(nt, -1)
        self.A = np.empty((nt, nb), dtype=np.complex128)
        for it in range(nt):
            self.A[it] = (np.bincount(idx, weights=a[it].real, minlength=nb)
                          + 1j * np.bincount(idx, weights=a[it].imag,
                                             minlength=nb))
        self.A /= (ny * nx)
        del a

    def values(self, z_list, deriv=False, sigma_m=0.0, k0_only=False):
        """
        (nt, nz) — Psi (or dPsi/dz) at this sample, at every height.

        `sigma_m` smears the SOURCE position by an isotropic Gaussian of that
        width. It is free here and it is not a convenience: the avalanche
        footprint (~34 µm) and transverse diffusion (107-826 µm) are both
        Gaussian smears of the deposit, they are exp(-k^2 sigma^2/2) in exactly
        the k this reduction is already indexed by, and whether a height effect
        survives them is the only form of the question that matters downstream.

        `k0_only` keeps the DC bin alone, which turns both z factors into their
        flat-weighting limits (1 - z/g and -1/g). Used to prove the slow path
        reduces to the fast path when the lateral structure is removed.
        """
        z_list = np.atleast_1d(np.asarray(z_list, dtype=float))
        self._check_margin(float(z_list.max()))
        fn = dpsi_dz_factor if deriv else psi_z_factor
        F = np.stack([fn(self.kbar, z, self.gap) for z in z_list], axis=1)
        if sigma_m:
            F = F * np.exp(-0.5 * (self.kbar * sigma_m) ** 2)[:, None]
        A = self.A
        if k0_only:
            A = np.zeros_like(A)
            A[:, 0] = self.A[:, 0]
        return np.real(A @ F)

    def _check_margin(self, z_max):
        if self._margin is None or z_max <= 0:
            return
        bound = np.exp(-self._margin / z_max)
        if bound > 1e-9:
            raise ValueError(
                f"y margin {self._margin*1e3:.2f} mm is too small for "
                f"z_max = {z_max*1e6:.0f} um (wrap bound {bound:.1e})")

    def _ref_and_got(self, z_list):
        if self.src is None:
            raise RuntimeError("source ZSlicer was freed; cannot certify")
        out = []
        for deriv in (False, True):
            out.append((self.src.at_exact(self.iy, self.ix, z_list,
                                          deriv=deriv),
                        self.values(z_list, deriv=deriv)))
        return out

    def certify(self, z_list=(1e-6, 40e-6, 110e-6), tol=1e-5):
        """Binned vs the direct reduction, on the same heights. Raises if bad."""
        if self.src is None:
            raise RuntimeError("source ZSlicer was freed; cannot certify")
        worst = 0.0
        for deriv in (False, True):
            ref = self.src.at_exact(self.iy, self.ix, z_list, deriv=deriv)
            got = self.values(z_list, deriv=deriv)
            scale = np.abs(ref).max()
            if scale > 0:
                worst = max(worst, float(np.abs(got - ref).max() / scale))
        if worst > tol:
            raise AssertionError(
                f"RadialProbe binning error {worst:.2e} exceeds {tol:.0e}; "
                "lower dk_g")
        return worst


# ── Self-test ────────────────────────────────────────────────────────────────

def certify_set(probes, z_list=(1e-6, 40e-6, 110e-6), tol=1e-8):
    """
    Binned vs direct reduction over a SET of probes, on a COMMON scale.

    Certifying each probe against its own amplitude is the wrong question and
    it fails for the right reason. The binning error is a fixed relative error
    per |k| bin, so on a channel whose kernel is a small residual of large
    cancelling contributions — which is exactly what the checkerboard produces
    for a channel whose view does not own the pad under the deposit — the
    per-probe relative error blows up while the ABSOLUTE error stays put. What
    matters downstream is the error against the largest channel present, since
    every residual in this study is normalised that way.
    """
    pairs = [pr._ref_and_got(z_list) for pr in probes]
    scale = max(float(np.abs(ref).max()) for p in pairs for ref, _ in p)
    err = max(float(np.abs(got - ref).max()) for p in pairs for ref, got in p)
    err /= max(scale, 1e-300)
    if err > tol:
        raise AssertionError(
            f"RadialProbe binning error {err:.2e} exceeds {tol:.0e}; lower dk_g")
    return err


def _fd_reference(v0, lx, gap, nz=1200):
    """
    Independent finite-difference solve of Laplace in (x, z), for the check.

    Deliberately shares NO code with the spectral form above: real-space
    5-point stencil, periodic in x, Dirichlet V = v0(x) at z = 0 and V = 0 at
    z = g, solved by banded elimination mode by mode in x... which would still
    be spectral in x. So it is solved as a genuine 2-D sparse system instead.
    """
    from scipy.sparse import diags, eye, kron
    from scipy.sparse.linalg import spsolve

    nx = len(v0)
    hx, hz = lx / nx, gap / (nz + 1)
    # periodic second difference in x
    off = np.ones(nx - 1)
    Dx = diags([off, -2 * np.ones(nx), off], [-1, 0, 1], format="lil")
    Dx[0, -1] = Dx[-1, 0] = 1.0
    Dx = (Dx.tocsr()) / hx ** 2
    # Dirichlet second difference in z on the nz interior planes
    Dz = diags([np.ones(nz - 1), -2 * np.ones(nz), np.ones(nz - 1)],
               [-1, 0, 1], format="csr") / hz ** 2
    A = kron(eye(nz), Dx) + kron(Dz, eye(nx))
    b = np.zeros(nz * nx)
    b[:nx] = -v0 / hz ** 2                       # the z = 0 boundary row
    v = spsolve(A.tocsr(), b).reshape(nz, nx)
    z = (np.arange(1, nz + 1)) * hz
    return z, v


def main():
    g = C.AMP_GAP_M
    print("zextend — Psi at z > 0 from the S1 surface solve\n")
    print(f"  amplification gap g = {g*1e6:.0f} µm\n")

    # 1. the two limits that must hold identically
    k = np.geomspace(1e1, 1e6, 40)
    assert np.allclose(psi_z_factor(k, 0.0, g), 1.0), "z=0 must be identity"
    assert np.allclose(psi_z_factor(k, g, g), 0.0, atol=1e-12), \
        "z=g must vanish (the mesh is grounded)"
    print("  z = 0 -> identity, z = g -> 0                       OK")

    # 2. k -> 0 must reproduce the parallel-plate weighting ions.py uses
    zs = np.linspace(0, g, 7)
    flat = psi_z_factor(np.full_like(zs, 1e-9), zs, g)
    assert np.allclose(flat, 1 - zs / g, atol=1e-9)
    print("  k -> 0 -> 1 - z/g (ions.py's f_ion)                 OK")

    # 3. the derivative factor against a difference of the value factor
    z0, h = 40e-6, 1e-9
    num = (psi_z_factor(k, z0 + h, g) - psi_z_factor(k, z0 - h, g)) / (2 * h)
    ana = dpsi_dz_factor(k, z0, g)
    rel = np.abs(num - ana) / np.abs(ana)
    print(f"  analytic dPsi/dz vs central difference: max rel "
          f"{rel.max():.2e}                OK")
    assert rel.max() < 1e-5

    # 4. THE INDEPENDENT CHECK: a real-space FD Laplace solve.
    #
    # Stated as a CONVERGENCE test, not as an agreement-to-a-tolerance one.
    # Agreement at one resolution proves nothing on its own — any bar loose
    # enough to pass is loose enough to hide a real discrepancy, and the first
    # version of this check "failed" at 4.4e-4 purely because it evaluated the
    # spectral form at the nominal z while reading the FD at its own nearest
    # plane, up to h_z/2 away. Refining h_x instead shows the FD residual
    # falling as h_x^2 with a flat constant, i.e. the FD converges TO the
    # spectral form and the spectral form is the exact limit. That is the
    # statement worth making, and it cannot be passed by luck.
    #
    # It converges in h_x and not h_z because the z direction is treated
    # exactly in the comparison (both sides read the same plane) while the
    # 3-point x stencil represents -k^2 as -(2-2cos k h_x)/h_x^2.
    lx = 1.6e-3
    print("\n  FD cross-check — real-space 5-point stencil, periodic in x,")
    print("  Dirichlet at both z faces, solved as one sparse 2-D system:")
    print(f"    {'nx':>6} {'h_x [µm]':>10} {'max |spec - FD|':>17} "
          f"{'x h_x^-2':>12}")
    errs = []
    for nx in (128, 256, 512, 1024):
        x = np.arange(nx) * lx / nx
        v0 = np.exp(-((x - lx / 2) ** 2) / (2 * (60e-6) ** 2))
        kx = 2 * np.pi * np.fft.fftfreq(nx, d=lx / nx)
        vh = np.fft.fft(v0)
        zfd, vfd = _fd_reference(v0, lx, g, nz=600)
        j = int(np.argmin(np.abs(zfd - 0.5 * g)))
        # read the spectral form at the FD's OWN plane, not at the nominal z
        spec = np.real(np.fft.ifft(vh * psi_z_factor(np.abs(kx), zfd[j], g)))
        err = float(np.abs(spec - vfd[j]).max() / np.abs(v0).max())
        errs.append(err)
        print(f"    {nx:6d} {lx/nx*1e6:10.3f} {err:17.3e} "
              f"{err*nx**2:12.2f}")
    # second order: each halving of h_x must cut the error by ~4
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    ok = all(3.6 < r < 4.4 for r in ratios)
    print(f"  successive error ratios {['%.2f' % r for r in ratios]} "
          f"-> order 2 in h_x   {'PASS' if ok else 'FAIL'}")
    print("    (a constant last column IS the proof: the FD converges to the")
    print("     spectral form, so the spectral form is the exact answer)")

    # 5. what it means for the kernel: the filter as a lateral low-pass
    print("\n  the filter, as a lateral low-pass (this is T10's physics):")
    scales = [("pad pitch", C.PAD_PITCH_M), ("ESL pitch", C.ESL_PITCH_M),
              ("200 µm", 200e-6), ("50 µm", 50e-6)]
    print(f"    {'z [µm]':>8} {'1 - z/g':>9}" +
          "".join(f"{'@ ' + n:>12}" for n, _ in scales))
    for z in (0, 15, 40, 75, 110, 149):
        zm = z * 1e-6
        row = "".join(
            f"{float(psi_z_factor(2*np.pi/s, zm, g)):12.4f}" for _, s in scales)
        print(f"    {z:8.0f} {1-zm/g:9.3f}{row}")
    # DERIVED, not typed. An earlier version hardcoded "0.047, a factor 11"
    # in this sentence while the table above printed 1e-4 — prose and table
    # from the same run contradicting each other. Anything quantitative here
    # now comes out of the same function that builds the table.
    z_hi = 75e-6
    flat = 1.0 - z_hi / g
    f200 = float(psi_z_factor(2 * np.pi / 200e-6, z_hi, g))
    f50 = float(psi_z_factor(2 * np.pi / 50e-6, z_hi, g))
    print("\n  A single f_ion = 1 - z/g rescales every column together.")
    print(f"  They do not move together — at z = {z_hi*1e6:.0f} µm the flat")
    print(f"  weighting says {flat:.3f}, while 200 µm structure (the scale the")
    print(f"  prompt kernel actually lives on, plan §3) is at {f200:.3f}, a")
    print(f"  factor {flat/f200:.1f}, and 50 µm structure is at {f50:.1e}, a")
    print(f"  factor {flat/f50:.0f}. The fast path applies the first number to")
    print("  a kernel built of the rest.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
