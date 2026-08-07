#!/usr/bin/env python3
"""
kernel_lut.py — turn an S1 product into something the digitizer can evaluate
millions of times (plan §7 step 5, the fast path).

The stored kernels are G_n(t, y0, x0): charge induced on channel n by a unit
point charge sitting on the ESL at (x0, y0) since t=0. The digitizer needs the
same thing per avalanche, but evaluated for arbitrary (x0, y0) and on a uniform
1 ns time grid, for every channel within reach.

Three things this class does that a naive np.load does not:

1. **Resamples the log time axis onto a uniform grid, once.** The products
   carry 61 log-spaced times because that is what resolves a relaxation
   spanning 0.1 ns to 10 µs; convolving with electronics needs uniform steps.
   Interpolating per avalanche would dominate the runtime.

2. **Keeps only the channels within reach.** A Y kernel is 1024 y samples over
   ±25 mm, but the digitizer only ever asks about channels within a few pitches
   of the deposit. Slicing to a window turns a 1.5 GB array into a few MB and
   is what makes the fast path fast.

3. **Differentiates to a current.** G is a CHARGE. Ramo gives the induced
   current as dQ/dt, and it is the current the electronics sees. Doing it once
   on the resampled grid is both faster and better conditioned than
   differentiating a log-spaced array per event.

The x axis is the 31.2 mm superperiod and is used modulo that period — which is
exact, not an approximation, because the whole stack is periodic in x with it.
"""

from __future__ import annotations

import json

import numpy as np

from ..common import constants as C


class CombKernelLUT:
    """Uniform-time, windowed lookup of the comb-channel Green's functions."""

    def __init__(self, path, dt_ns=1.0, t_max_ns=3000.0, y_window_mm=6.0):
        self.path = path
        with np.load(path) as d:
            self.meta = json.loads(str(d["meta"]))
            t = d["t"]
            self.x = d["x"].copy()
            y_Y = d["y_Y"].copy()
            self.y_X = d["y_X"].copy()
            # Uniform time grid the digitizer works on.
            self.t = np.arange(0.0, t_max_ns + dt_ns, dt_ns) * 1e-9
            self.dt = dt_ns * 1e-9

            # Y kernels: window in y around the channel row before resampling,
            # otherwise we resample 1024 rows we will never look at.
            keep = np.abs(y_Y) <= y_window_mm * 1e-3
            self.y_Y = y_Y[keep]
            self.G_Y = np.stack([
                self._resample(t, d["G_Y_even"][:, keep, :]),
                self._resample(t, d["G_Y_odd"][:, keep, :])])
            # X kernels: y box is only 1.56 mm, so keep all of it.
            self.G_X = np.stack([self._resample(t, d["G_X"][c])
                                 for c in range(C.N_PAD_PER_SUPER)])

        # dQ/dt on the uniform grid: the induced CURRENT per unit charge.
        self.I_Y = np.gradient(self.G_Y, self.dt, axis=1)
        self.I_X = np.gradient(self.G_X, self.dt, axis=1)

    def _resample(self, t_src, G):
        """(nt_src, ny, nx) on log times -> (nt_dst, ny, nx) on the uniform grid."""
        out = np.empty((len(self.t),) + G.shape[1:], dtype=np.float32)
        flat = G.reshape(G.shape[0], -1).astype(np.float64)
        res = np.empty((len(self.t), flat.shape[1]), dtype=np.float32)
        for j in range(flat.shape[1]):
            res[:, j] = np.interp(self.t, t_src, flat[:, j])
        return res.reshape(out.shape)

    # ── lookup ───────────────────────────────────────────────────────────────

    def ix(self, x0_m):
        """Column index of the x sample nearest x0, folded into the superperiod."""
        xs = np.mod(np.asarray(x0_m, dtype=float), C.SUPERPERIOD_M)
        return np.clip(np.rint(xs / (self.x[1] - self.x[0])).astype(int),
                       0, len(self.x) - 1)

    def y_channel_current(self, x0_m, dy_m, parity):
        """
        Induced current on the Y channel whose row is dy from the deposit.

        `parity` is the row's column parity (row index mod 2). Returns a
        (nt,) array, per unit deposited charge.
        """
        iy = np.clip(np.rint((np.asarray(dy_m) - self.y_Y[0])
                             / (self.y_Y[1] - self.y_Y[0])).astype(int),
                     0, len(self.y_Y) - 1)
        return self.I_Y[parity, :, iy, self.ix(x0_m)]

    def x_channel_current(self, x0_m, y0_m, col):
        """Induced current on the X channel at column `col` (mod 40 phases)."""
        iy = np.clip(np.rint((np.asarray(y0_m) - self.y_X[0])
                             / (self.y_X[1] - self.y_X[0])).astype(int),
                     0, len(self.y_X) - 1)
        return self.I_X[col % C.N_PAD_PER_SUPER, :, iy, self.ix(x0_m)]

    def describe(self):
        return {
            "product": self.path.split("/")[-1],
            "rho_s_MOhm_sq": self.meta["rho_s_ohm_sq"] / 1e6,
            "d_kapton_um": self.meta["d_kapton_m"] * 1e6,
            "dt_ns": self.dt * 1e9,
            "t_max_ns": self.t[-1] * 1e9,
            "y_window_mm": float(np.abs(self.y_Y).max() * 1e3),
            "solver_git": self.meta.get("git", "")[:12],
        }
