#!/usr/bin/env python3
"""
kernel_lut.py — turn an S1 product into something the digitizer can evaluate
millions of times (plan §7 step 5, the fast path).

The stored kernels are G_n(t, y0, x0): charge induced on channel n by a unit
point charge sitting on the ESL at (x0, y0) since t=0. The digitizer needs the
induced CURRENT, on a uniform time grid, for every channel within reach of an
avalanche.

SIZE IS THE WHOLE DESIGN PROBLEM. A stored product is 2.3 GB, and the obvious
"resample everything onto a 1 ns grid" blows that up rather than shrinking it:
the 40 X kernels alone would be 40 x 3001 x 32 x 3120 float32 = 48 GB. A first
version of this file did exactly that and reached 9.4 GB resident on a 16 GB
laptop before being killed. Three cuts, each justified by what the digitizer
can actually resolve:

  * **x is decimated** (default 4x, 10 -> 40 µm). Every avalanche is smeared by
    ~826 µm of transverse diffusion before it lands, so 10 µm sampling carries
    no information the digitizer can use. It is kept in the *product* because
    the 31.2 mm beat and the ESL strip edges need it; it is not needed here.
  * **y is windowed** to the few pad pitches a channel can reach.
  * **only the current is kept.** G itself is never used downstream — Ramo gives
    the signal as dQ/dt — so it is differentiated and dropped.

The X kernels get one further transformation. Indexed by absolute column they
are unusable (40 columns x every x), but the digitizer only ever wants "the
channel d pads away from this avalanche". Reindexing to that offset collapses
the 40 to 2*n_side+1 and makes the array a thin band instead of a full matrix.

The x axis is used modulo the 31.2 mm superperiod, which is exact rather than
an approximation: the whole stack is periodic in x with that period.
"""

from __future__ import annotations

import json

import numpy as np

from ..common import constants as C


class CombKernelLUT:
    """Windowed, uniform-time lookup of the comb-channel induced currents."""

    # THE STORED WINDOW MUST COVER THE DAQ FRAME (fix 2026-08-07, audit A1).
    # It was 1000 ns. The S1 product is solved to 10 us and the DAQ area window
    # is 32 x 60 = 1.92 us (64 x 60 = 3.84 us in the SPS config), while sheet
    # transport feeds d = +-2 on the ~0.6 us scale and d = +-3 past 1 us
    # (D ~ 1 m^2/s at rho_s = 2 MOhm/sq). Truncating at 1 us therefore
    # GUARANTEED the sim under-delivered the late halo — structurally, in the
    # same direction as the headline §9 "missing broad slow halo" tension.
    #
    # 3000 ns covers the bench 1.92 us window plus the 340 ns ion transit plus
    # margin. An SPS-config run should ask for 4200.
    T_MAX_NS_DEFAULT = 3000.0

    # y is decimated for the same reason x is (2026-08-08). A product solved at
    # ny = 1024 (audit C6) samples y every 48.75 um, and at n_side = 8 that
    # makes the LUT ~10 GB before transients — it was OOM-killed on lxplus. The
    # fine grid exists so the SOLVE resolves the prompt kernel's pad-edge
    # shoulder, which is a property of the stored product; the digitizer cannot
    # use it, because every avalanche arrives smeared by at least ~107 um. So
    # the LUT decimates to a target spacing and the memory becomes independent
    # of the product's ny.
    Y_TARGET_UM = 97.5              # = the ny = 512 spacing, the old behaviour

    def __init__(self, path, dt_ns=1.0, t_max_ns=T_MAX_NS_DEFAULT,
                 y_window_mm=7.02, x_stride=4, n_side=8, y_target_um=None):
        self.path = path
        self.n_side = n_side
        self.t = np.arange(0.0, t_max_ns + dt_ns, dt_ns) * 1e-9
        self.dt = dt_ns * 1e-9
        ds = np.arange(-n_side, n_side + 1)
        self.ds = ds

        with np.load(path) as d:
            self.meta = json.loads(str(d["meta"]))
            t_src = d["t"]
            self.x = d["x"][::x_stride].copy()
            self.dx = self.x[1] - self.x[0]
            y_Y_all = d["y_Y"]
            self.y_X = d["y_X"].copy()          # re-strided below

            # Decimate y to the target spacing, then window. Both boxes share
            # the y grid by construction (kernels.x_ny), so the same stride
            # applies to each and they stay commensurate.
            tgt = (self.Y_TARGET_UM if y_target_um is None else y_target_um)
            dy_um = float(abs(y_Y_all[1] - y_Y_all[0])) * 1e6
            self.y_stride = max(1, int(round(tgt / dy_um)))
            keep = np.zeros(len(y_Y_all), dtype=bool)
            keep[::self.y_stride] = True
            keep &= np.abs(y_Y_all) <= y_window_mm * 1e-3
            self.y_Y = y_Y_all[keep].copy()

            # --- Y: (parity, nt, ny_win, nx) ---------------------------------
            gy = np.stack([d["G_Y_even"][:, keep, ::x_stride],
                           d["G_Y_odd"][:, keep, ::x_stride]])
            self.I_Y = np.stack([self._resample(t_src, gy[0]),
                                 self._resample(t_src, gy[1])])
            del gy
            self.I_Y = self._time_last(self._to_current(self.I_Y, axis=1))

            # --- X: reindex from absolute column to channel OFFSET ------------
            # For an avalanche at x sample ix the nearest pad column is
            # col(ix); the channel d pads away is (col(ix)+d) mod 40. Build
            # band[d, t, jy, ix] once so the digitizer never searches.
            nx = len(self.x)
            col_of_x = np.rint((self.x - np.mod(_pad_origin(), C.SUPERPERIOD_M))
                               / C.PAD_PITCH_M).astype(int)
            # KEPT, because it is the LUT's own definition of "which column is
            # d=0 at this x sample" and the digitizer must book its channels on
            # exactly that and not on a second, independent rounding of the
            # true x. See `col_at` and audit C2.
            self.col_of_x = col_of_x % C.N_PAD_PER_SUPER
            # Same y decimation on the X box.
            xs_keep = np.zeros(len(self.y_X), dtype=bool)
            xs_keep[::self.y_stride] = True
            self.y_X = self.y_X[xs_keep].copy()
            gx = d["G_X"][:, :, xs_keep, ::x_stride]
            band = np.empty((len(ds), len(self.t), gx.shape[2], nx),
                            dtype=np.float32)
            for j, dd in enumerate(ds):
                cols = (col_of_x + dd) % C.N_PAD_PER_SUPER
                # gather the right column's kernel at each x sample
                sel = gx[cols, :, :, np.arange(nx)]      # (nx, nt, jy)
                band[j] = self._resample(t_src, np.transpose(sel, (1, 2, 0)))
            del gx
            self.I_X = self._time_last(self._to_current(band, axis=1))

    @staticmethod
    def _time_last(a):
        """
        Move the time axis from position 1 to LAST, contiguously.

        This is the single biggest lever on digitizer speed and it is pure
        memory layout. With time slowest-varying, reading one channel's time
        series strides ny*nx floats — 491 kB — between consecutive samples, so
        1001 samples touch 1001 separate cache lines spread over 503 MB. The
        induction ran at ~800 us per avalanche, about 22 Mflop/s, which is
        absurd for 18k float operations. Time last makes the same series 3.9 kB
        of sequential memory.
        """
        return np.ascontiguousarray(np.moveaxis(a, 1, -1))

    def _to_current(self, G, axis):
        """
        Induced charge -> induced current, KEEPING THE PROMPT STEP.

        G(0+) is not zero: it is the instantaneous electrostatic image, formed
        before the sheet conducts at all, and on the ns grid it is a step. The
        physical current is

            i(t) = G(0+) delta(t)  +  dG/dt

        and np.gradient represents only the second term. Dropping the first
        silently deletes the single largest component of the signal — with a
        centred gradient the time-integrated charge came out ~0 for the X view
        and NEGATIVE for the Y view, because all that survived was the late
        redistribution.

        A backward difference with a zero prepended is both the delta and the
        derivative at once, and it telescopes: cumsum(i)*dt == G exactly, so
        charge is conserved by construction rather than by luck.
        """
        return (np.diff(G, prepend=0.0, axis=axis) / self.dt).astype(np.float32)

    def _resample(self, t_src, G):
        """
        (nt_src, ny, nx) on the log time axis -> (nt_dst, ny, nx) uniform.

        Vectorised on the shared time axis. np.interp is 1-D, so the obvious
        implementation loops over columns — but there are millions of them and
        that loop runs for minutes per product.
        """
        t_src = np.asarray(t_src, dtype=float)
        j = np.clip(np.searchsorted(t_src, self.t, side="right") - 1,
                    0, len(t_src) - 2)
        span = t_src[j + 1] - t_src[j]
        w = np.clip((self.t - t_src[j]) / span, 0.0, 1.0)[:, None, None]
        lo = G[j].astype(np.float32)
        hi = G[j + 1].astype(np.float32)
        return lo + (hi - lo) * w.astype(np.float32)

    # ── lookup ───────────────────────────────────────────────────────────────

    def ix(self, x0_m):
        """
        x sample index, folded into the 31.2 mm superperiod.

        The metric is CIRCULAR: the axis is periodic, so a source in the last
        half-stride before the seam is nearest to x[0], not to x[-1]. Clipping
        there (the pre-2026-08-07 behaviour, audit C5) mapped 31.19 mm onto the
        last sample, up to a full stride — 40 µm — away from the truth, and did
        it only in that one sliver of x, which is the hardest kind of error to
        notice in an average.
        """
        xs = np.mod(np.asarray(x0_m, dtype=float), C.SUPERPERIOD_M)
        return np.rint(xs / self.dx).astype(int) % len(self.x)

    def col_at(self, x0_m, ix=None):
        """
        ABSOLUTE pad column whose kernel the LUT serves as d = 0 at x0.

        The X band is built per LUT x SAMPLE: `band[d, ..., ix]` holds the
        kernel of column `col_of_x[ix] + d`. The digitizer used to derive its
        channel number independently, by rounding the true x onto the pad
        lattice — two roundings of the same quantity onto grids that do not
        share their boundaries. The LUT samples every 40 µm and pads are
        780 µm apart, so at every second pad boundary the two answers differ by
        a whole pad: ~1.3 % of avalanches were booked one channel off in X,
        with the sign alternating (audit C2). Both tests were blind to it —
        the charge audit sums over channels, and T10 draws its probes ON the
        LUT grid where the two conventions agree by construction.

        Resolved by making the LUT the single authority: take the column
        modulo the superperiod from the LUT sample, and take only the
        superperiod COUNT from the true x, which no rounding can disagree on.
        """
        if ix is None:
            ix = self.ix(x0_m)
        naive = np.rint((np.asarray(x0_m, dtype=float) - _pad_origin())
                        / C.PAD_PITCH_M).astype(int)
        # Snap to the nearest integer congruent to the LUT's local column.
        d = (self.col_of_x[ix] - naive) % C.N_PAD_PER_SUPER
        d = np.where(d > C.N_PAD_PER_SUPER // 2, d - C.N_PAD_PER_SUPER, d)
        return naive + d

    def iy_Y(self, dy_m):
        return np.clip(
            np.rint((np.asarray(dy_m) - self.y_Y[0])
                    / (self.y_Y[1] - self.y_Y[0])).astype(int),
            0, len(self.y_Y) - 1)

    def iy_X(self, y_rel_m):
        return np.clip(
            np.rint((np.asarray(y_rel_m) - self.y_X[0])
                    / (self.y_X[1] - self.y_X[0])).astype(int),
            0, len(self.y_X) - 1)

    def apply_ion_transit(self, f_electron, transit_ns):
        """
        Fold the ion transit into the stored kernels, ONCE.

        Physics (plan §7 step 5, response/digitizer/ions.py): a pair born at
        height z in the amplification gap splits its induced charge z/g to the
        electron and 1 - z/g to the ion, and a Micromegas avalanche sits
        microns off the anode, so the ions carry ~97 %. The ion delivers that
        share not promptly but spread FLAT over its transit — a rectangle, not
        a decaying tail.

        Done here rather than per avalanche because every avalanche shares the
        same shape, so this is one pass over the LUT instead of 150k
        convolutions per file.

        A rectangle convolution is a running mean, so cumsum gives it in O(N)
        with no FFT and no scipy.

        THE LEADING-EDGE DENOMINATOR IS L, NOT n+1. Before the window fills,
        this still divides by the FULL L. Writing it as a mean over the samples
        seen so far — c[n]/(n+1) — is the natural thing to type and it is
        wrong: the ion delivers f_i/T per unit time from the moment it is born,
        so for t < T the answer is (1/T) * integral_0^t x, i.e. c[n]/L. The
        n+1 version was in place 2026-08-07 and inflated the first T = 340 ns
        by up to 340x, put 5.9 units of charge on the readout for every 1
        induced, and — because it is a TIME-dependent distortion while channels
        peak at different times — did not cancel from inter-channel ratios.
        test_longitudinal.py pins this against an impulse, whose answer is
        exact by inspection.

        THE APPROXIMATION, stated again because it is easy to forget once the
        number looks right: the ion keeps the SURFACE kernel's lateral shape.
        The true Psi_n broadens as the ion climbs toward the mesh, so this
        under-estimates ion lateral sharing. T10 (Garfield ComponentGrid over
        the S1 time slices) is what replaces it.

        ON THE "TRAILING-EDGE TRUNCATION" (audit A1, checked 2026-08-07 and NOT
        a bug). The fix order asked for this convolution to be padded to
        nt + L before truncating, on the theory that it was silently losing its
        own tail the way the leading edge once did. It is not. The running mean
        at sample n reads only x[n-L+1..n], so every output sample inside
        [0, nt) is already complete and padding the record is a bit-for-bit
        no-op (verified directly). The same holds for `apply_longitudinal`,
        whose FFT is sized nt + len(h) - 1 and is therefore a genuine LINEAR
        convolution, not a wrapped one.

        What IS lost is the output beyond sample nt — for a 340 ns transit,
        ~16 % of the delivered area on a decaying test record. That charge is
        outside t_max by construction and no amount of padding recovers it;
        only a longer window does. Which is the real content of A1, and is why
        T_MAX_NS_DEFAULT moved 1000 -> 3000 above.
        """
        L = max(1, int(round(transit_ns / (self.dt * 1e9))))
        f_ion = 1.0 - f_electron

        def smear(a):
            # Chunk over the leading axes and stay in float32. Promoting the
            # whole array to float64 for the cumsum tripled a 1.9 GB LUT and
            # got the process OOM-killed on a 16 GB laptop; the running mean
            # needs no extra precision here because the kernel is float32 to
            # begin with.
            flat = a.reshape(-1, a.shape[-1])
            out = np.empty_like(flat)
            step = max(1, int(4e6 // a.shape[-1]))
            for i in range(0, len(flat), step):
                blk = flat[i:i + step]
                c = np.cumsum(blk, axis=-1)
                run = np.empty_like(c)
                run[:, :L] = c[:, :L] / np.float32(L)
                run[:, L:] = (c[:, L:] - c[:, :-L]) / np.float32(L)
                out[i:i + step] = f_electron * blk + f_ion * run
            return out.reshape(a.shape)

        self.I_Y = smear(self.I_Y)
        self.I_X = smear(self.I_X)
        self.ion_transit_ns = transit_ns
        self.f_electron = f_electron
        self.longitudinal_model = "analytic rectangle"

    def apply_longitudinal(self, h, dt_h):
        """
        Fold a MEASURED longitudinal current profile into the LUT (T9 v2).

        Same physics slot as apply_ion_transit, but h is an arbitrary measured
        shape (response/digitizer/ions.measured_longitudinal) instead of the
        analytic delta + rectangle, so it needs a real FIR rather than the
        cumsum running mean.

        The convolution is the right operation because the LUT holds the
        response to charge delivered INSTANTANEOUSLY at the surface, and h is
        how that charge is actually delivered in time. Total signal =
        (impulse response) conv (delivery profile).

        h is resampled onto the LUT's own dt by area-preserving accumulation,
        not by interpolation: h is 0.2 ns Garfield sampling going onto a 1 ns
        grid, and point-sampling a 5x coarser grid would throw away four
        fifths of the prompt electron spike and silently rescale the charge
        split. Renormalised after, so sum(h)*dt == 1 exactly.

        THE APPROXIMATION IS UNCHANGED and is worth restating precisely
        because a measured input makes it easy to believe the whole thing is
        now measured: h is a purely LONGITUDINAL profile, measured with a flat
        (whole-electrode) weighting field, so it carries no lateral
        information at all. The ion still gets the surface kernel's lateral
        shape frozen at z=0. This makes the timing and the charge split
        first-principles; it does not touch T10.
        """
        h = np.asarray(h, dtype=np.float64)
        nt = self.I_Y.shape[-1]

        # --- area-preserving resample onto the LUT grid ----------------------
        ratio = dt_h / self.dt
        edges = (np.arange(len(h) + 1) * ratio)
        idx = np.floor(edges[:-1]).astype(int)
        nb = int(np.ceil(edges[-1]))
        hb = np.bincount(idx, weights=h * dt_h, minlength=nb)[:nb]
        hb = hb / hb.sum()                      # unit AREA on the coarse grid
        if len(hb) > nt:
            hb = hb[:nt]
            hb = hb / hb.sum()

        # --- chunked FFT convolution, truncated to the causal window ---------
        nfft = 1 << int(np.ceil(np.log2(nt + len(hb) - 1)))
        H = np.fft.rfft(hb, nfft)

        def conv(a):
            flat = a.reshape(-1, nt)
            out = np.empty_like(flat)
            step = max(1, int(2e6 // nfft))
            for i in range(0, len(flat), step):
                blk = flat[i:i + step]
                y = np.fft.irfft(np.fft.rfft(blk, nfft, axis=-1) * H,
                                 nfft, axis=-1)
                out[i:i + step] = y[:, :nt].astype(np.float32)
            return out.reshape(a.shape)

        self.I_Y = conv(self.I_Y)
        self.I_X = conv(self.I_X)
        self.longitudinal_model = "S3 v2 measured"
        self.ion_transit_ns = float(np.sum(np.arange(len(hb)) * hb)
                                    * self.dt * 1e9)   # centroid, for the log
        self.f_electron = None

    def nbytes(self):
        return self.I_Y.nbytes + self.I_X.nbytes

    def describe(self):
        return {
            "product": self.path.split("/")[-1],
            "rho_s_MOhm_sq": self.meta["rho_s_ohm_sq"] / 1e6,
            "d_kapton_um": self.meta["d_kapton_m"] * 1e6,
            "dt_ns": self.dt * 1e9,
            "t_max_ns": self.t[-1] * 1e9,
            "dx_um": self.dx * 1e6,
            "y_window_mm": float(np.abs(self.y_Y).max() * 1e3),
            "y_stride": self.y_stride,
            "dy_um": float(abs(self.y_Y[1] - self.y_Y[0]) * 1e6)
                     if len(self.y_Y) > 1 else None,
            "n_side": self.n_side,
            "lut_GB": self.nbytes() / 1e9,
            "solver_git": self.meta.get("git", "")[:12],
        }


def _pad_origin():
    from ..solver.kernels import PAD_ORIGIN_M
    return PAD_ORIGIN_M
