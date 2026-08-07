#!/usr/bin/env python3
"""
digitize.py — Stage B, clusters to per-channel induced current (plan §7).

Chain, per ionization cluster from the Geant4 ClusterTree:

    1. electrons        n = nPrimary (Geant4 already applied W and carried the
                        sub-W remainder probabilistically -- do NOT re-derive
                        from edep, see design/NEEDED_INPUTS.md)
    2. drift            t += z/v_d + Gauss(sigma_L/v_d); (x,y) += Gauss(sigma_T)
                        from the Magboltz table (response/digitizer/gas.py)
    3. mesh             binomial thinning by the transparency eps
    4. avalanche        gain ~ Polya(gbar, theta) from aval_calib.json (S3),
                        landing with an extra Gauss(sigma0) of ~33 um
    5. induction        i_n(t) = Q * dG_n/dt from the S1 comb kernels, summed
                        over channels within reach (response/digitizer/kernel_lut.py)

WHAT SETS THE SHARING, and why the packet approximation is not safe here.
The three lateral scales are wildly separated at the bench point:

    avalanche sigma0            33 um     (T7, and nearly voltage independent)
    prompt induction sigma_p   222 um     (T2b)
    drift diffusion, full gap  826 um     (Magboltz, 333 V/cm over 30 mm)

Diffusion dominates and is the ONLY one comparable to the 780 um pad pitch, so
it is what produces the measured c1 = 0.23-0.28. That also means the per-cluster
"packet" approximation the plan offers as a default (§7 step 2) is wrong for
this detector: a packet puts all of a cluster's electrons at one point, which
discards exactly the spread that carries the observable. Electrons are therefore
transported INDIVIDUALLY by default (`packet=False`); the packet path is kept
only to measure how much it costs, which is what plan §12 item 11 asks for.

    python3 -m response.digitizer.digitize clusters.root --kernel <s1.npz>
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from ..common import constants as C
from ..solver import kernels as K
from .gas import DriftGas
from .kernel_lut import CombKernelLUT

# Bench operating point (plan §1, §5). Everything here is recorded per run.
DEFAULT_DRIFT_V = 1000.0
DEFAULT_DRIFT_GAP_MM = 30.0
DEFAULT_MESH_V = 490.0
# Mesh electron transparency. T6 will replace this with a measured curve; until
# then it is an explicit input, NOT a hidden constant.
DEFAULT_TRANSPARENCY = 0.95


def load_calib(path, mesh_v):
    """Pick the S3 point nearest the requested mesh voltage."""
    with open(path) as f:
        cal = json.load(f)["points"]
    best, bestd = None, None
    for key, rec in cal.items():
        try:
            v = float(key.rsplit("@", 1)[1].rstrip("V"))
        except (IndexError, ValueError):
            continue
        d = abs(v - mesh_v)
        if bestd is None or d < bestd:
            best, bestd, best_v = rec, d, v
    if best is None:
        raise ValueError(f"no usable points in {path}")
    if bestd > 1e-6:
        print(f"  [calib] nearest S3 point is {best_v:.0f} V, asked "
              f"{mesh_v:.0f} V — using it and recording the mismatch")
    return best, best_v


def polya_sample(rng, gbar, theta, n):
    """
    Polya-distributed gains. P(g) ~ (g/gbar)^theta exp(-(1+theta) g/gbar) is a
    Gamma with shape (1+theta) and scale gbar/(1+theta), so sample it directly
    rather than by rejection.
    """
    return rng.gamma(shape=1.0 + theta, scale=gbar / (1.0 + theta), size=n)


class Digitizer:
    def __init__(self, kernel_path, calib_path, *, mesh_v=DEFAULT_MESH_V,
                 drift_v=DEFAULT_DRIFT_V, drift_gap_mm=DEFAULT_DRIFT_GAP_MM,
                 transparency=DEFAULT_TRANSPARENCY, v_scale=1.0,
                 n_chan_side=4, seed=12345, packet=False, gas_table=None):
        self.lut = CombKernelLUT(kernel_path)
        self.gas = (DriftGas(gas_table, v_scale=v_scale) if gas_table
                    else DriftGas(v_scale=v_scale))
        self.calib, self.calib_v = load_calib(calib_path, mesh_v)
        self.E_drift = drift_v / (drift_gap_mm * 0.1)      # V/cm
        self.drift_gap_mm = drift_gap_mm
        self.transparency = transparency
        self.n_side = n_chan_side
        self.packet = packet
        self.rng = np.random.default_rng(seed)

        p = self.calib["polya"]
        self.gbar, self.theta = p["gain_mean"], p["theta"]
        self.sigma0_um = self.calib["sigma0_um"]
        self.v_drift = float(np.ravel(
            self.gas.v_drift_um_ns(self.E_drift))[0])

    # ── the physics ──────────────────────────────────────────────────────────

    def transport(self, x_mm, y_mm, z_mm, t_ns, n_e):
        """
        Clusters -> individual avalanche seeds at the ESL.

        Returns (x, y [m], t [ns], gain) arrays, one entry per surviving
        electron (or per cluster in packet mode).
        """
        n_e = np.asarray(n_e, dtype=int)
        if self.packet:
            xs, ys, zs, ts = map(np.asarray, (x_mm, y_mm, z_mm, t_ns))
            w = n_e.astype(float)
        else:
            rep = np.repeat(np.arange(len(n_e)), n_e)
            if rep.size == 0:
                return (np.empty(0),) * 4
            xs, ys = np.asarray(x_mm)[rep], np.asarray(y_mm)[rep]
            zs, ts = np.asarray(z_mm)[rep], np.asarray(t_ns)[rep]
            w = np.ones(len(rep))

        zs = np.clip(zs, 0.0, None)
        sT = self.gas.sigma_T_um(self.E_drift, zs)
        st = self.gas.sigma_t_ns(self.E_drift, zs)

        x = xs * 1e-3 + self.rng.normal(0.0, sT * 1e-6)
        y = ys * 1e-3 + self.rng.normal(0.0, sT * 1e-6)
        t = ts + zs * 1e3 / self.v_drift + self.rng.normal(0.0, st)

        # Mesh transparency, then the avalanche footprint.
        if self.packet:
            surv = self.rng.binomial(n_e, self.transparency).astype(float)
        else:
            surv = (self.rng.random(len(w)) < self.transparency).astype(float)
        keep = surv > 0
        x, y, t, surv = x[keep], y[keep], t[keep], surv[keep]
        x += self.rng.normal(0.0, self.sigma0_um * 1e-6, len(x))
        y += self.rng.normal(0.0, self.sigma0_um * 1e-6, len(y))

        gain = polya_sample(self.rng, self.gbar, self.theta, len(x)) * surv
        return x, y, t, gain

    def induce(self, x, y, t, q, n_samp):
        """
        Sum the induced current of every avalanche onto the channel grid.

        Returns {("X"|"Y", channel_index): current[n_samp]} on the LUT's 1 ns
        grid, in units of elementary charges per second.
        """
        out = {}
        if len(x) == 0:
            return out
        pitch = C.PAD_PITCH_M
        ds = np.arange(-self.n_side, self.n_side + 1)

        # Channel index of the pad nearest each avalanche, on the 0.78 mm
        # lattice: column from x, row from y.
        col = np.rint((x - K.PAD_ORIGIN_M) / pitch).astype(int)
        row = np.rint((y - K.PAD_ORIGIN_M) / pitch).astype(int)
        k0 = np.rint(t / (self.lut.dt * 1e9)).astype(int)

        ok = np.isfinite(q) & (q > 0) & (k0 >= 0) & (k0 < n_samp)
        idx = np.flatnonzero(ok)

        # Pre-resolve the LUT axis indices once per avalanche.
        ix = self.lut.ix(x)
        dy_step = self.lut.y_Y[1] - self.lut.y_Y[0]
        dyx_step = self.lut.y_X[1] - self.lut.y_X[0]

        for i in idx:
            # --- Y view: channels are the ROWS row[i] + d -------------------
            # Kernel argument is the deposit's offset from THAT channel's row.
            rows = row[i] + ds
            dy = y[i] - (K.PAD_ORIGIN_M + rows * pitch)
            iy = np.clip(np.rint((dy - self.lut.y_Y[0]) / dy_step).astype(int),
                         0, len(self.lut.y_Y) - 1)
            curY = self.lut.I_Y[rows % 2, :, iy, ix[i]]        # (nd, nt)

            # --- X view: channels are the COLUMNS col[i] + d ----------------
            # The kernel for column c already knows where column c sits, so
            # the deposit's ABSOLUTE x is what goes in — shifting x by the
            # channel offset (as an earlier version did) double-counts the
            # separation and destroys the 31.2 mm beat the X view carries.
            cols = col[i] + ds
            y_rel = np.mod(y[i] - K.PAD_ORIGIN_M, 2 * pitch)
            jy = int(np.clip(round((y_rel - self.lut.y_X[0]) / dyx_step),
                             0, len(self.lut.y_X) - 1))
            curX = self.lut.I_X[cols % C.N_PAD_PER_SUPER, :, jy, ix[i]]

            for j, d in enumerate(ds):
                self._accumulate(out, ("Y", int(rows[j])), curY[j], q[i],
                                 k0[i], n_samp)
                self._accumulate(out, ("X", int(cols[j])), curX[j], q[i],
                                 k0[i], n_samp)
        return out

    @staticmethod
    def _accumulate(out, key, cur, q, k0, n_samp):
        n = min(len(cur), n_samp - k0)
        if n <= 0:
            return
        buf = out.get(key)
        if buf is None:
            buf = out[key] = np.zeros(n_samp)
        buf[k0:k0 + n] += q * cur[:n]

    # ── bookkeeping ──────────────────────────────────────────────────────────

    def describe(self):
        d = {
            "kernel": self.lut.describe(),
            "gas": self.gas.describe(self.E_drift, self.drift_gap_mm),
            "E_drift_Vcm": self.E_drift,
            "mesh_transparency": self.transparency,
            "transparency_source": "ASSUMED (T6 pending)",
            "gain_mean": self.gbar, "polya_theta": self.theta,
            "sigma0_um": self.sigma0_um,
            "calib_voltage_V": self.calib_v,
            "field_model": self.calib.get("field_model"),
            "packet_mode": self.packet,
            "n_chan_side": self.n_side,
        }
        return d
