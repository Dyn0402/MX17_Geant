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

    def __init__(self, path, dt_ns=1.0, t_max_ns=1000.0, y_window_mm=3.9,
                 x_stride=4, n_side=4):
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
            self.y_X = d["y_X"].copy()

            keep = np.abs(y_Y_all) <= y_window_mm * 1e-3
            self.y_Y = y_Y_all[keep].copy()

            # --- Y: (parity, nt, ny_win, nx) ---------------------------------
            gy = np.stack([d["G_Y_even"][:, keep, ::x_stride],
                           d["G_Y_odd"][:, keep, ::x_stride]])
            self.I_Y = np.stack([self._resample(t_src, gy[0]),
                                 self._resample(t_src, gy[1])])
            del gy
            self.I_Y = np.gradient(self.I_Y, self.dt, axis=1).astype(np.float32)

            # --- X: reindex from absolute column to channel OFFSET ------------
            # For an avalanche at x sample ix the nearest pad column is
            # col(ix); the channel d pads away is (col(ix)+d) mod 40. Build
            # band[d, t, jy, ix] once so the digitizer never searches.
            nx = len(self.x)
            col_of_x = np.rint((self.x - np.mod(_pad_origin(), C.SUPERPERIOD_M))
                               / C.PAD_PITCH_M).astype(int)
            gx = d["G_X"][:, :, :, ::x_stride]
            band = np.empty((len(ds), len(self.t), gx.shape[2], nx),
                            dtype=np.float32)
            for j, dd in enumerate(ds):
                cols = (col_of_x + dd) % C.N_PAD_PER_SUPER
                # gather the right column's kernel at each x sample
                sel = gx[cols, :, :, np.arange(nx)]      # (nx, nt, jy)
                band[j] = self._resample(t_src, np.transpose(sel, (1, 2, 0)))
            del gx
            self.I_X = np.gradient(band, self.dt, axis=1).astype(np.float32)

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
        """x sample index, folded into the 31.2 mm superperiod."""
        xs = np.mod(np.asarray(x0_m, dtype=float), C.SUPERPERIOD_M)
        return np.clip(np.rint(xs / self.dx).astype(int), 0, len(self.x) - 1)

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
            "n_side": self.n_side,
            "lut_GB": self.nbytes() / 1e9,
            "solver_git": self.meta.get("git", "")[:12],
        }


def _pad_origin():
    from ..solver.kernels import PAD_ORIGIN_M
    return PAD_ORIGIN_M
