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
from . import ions as ION

# Bench operating point (plan §1, §5). Everything here is recorded per run.
DEFAULT_DRIFT_V = 1000.0
DEFAULT_DRIFT_GAP_MM = 30.0
DEFAULT_MESH_V = 490.0
# Mesh electron transparency — MEASURED, T6 (response/meshcell/mesh_transparency.C,
# 2026-08-07). At the bench point (E_drift 333 V/cm over the 30 mm gap against
# 32.7 kV/cm in the 150 um amplification gap, ratio 98) the woven-mesh cell gives
# 0.873. The curve is remarkably flat: 0.854-0.876 across field ratios 33-653,
# because even the lowest ratio here is far above the ~10-20 where transparency
# starts to collapse. Superseded the 0.95 that was assumed before T6 ran.
DEFAULT_TRANSPARENCY = 0.873


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
                 # 8, not 4 (2026-08-07). With the LUT window corrected to
                 # 3000 ns (audit A1) a +-4 channel window LEAKS 8-10 % of the
                 # induced charge at every depth; +-8 holds it to <=0.6 %. The
                 # two fixes are not separable — see test_charge_audit and
                 # design/report/DESKTOP_RUNS_2026-08-07.md.
                 n_chan_side=8, seed=12345, packet=False, gas_table=None,
                 # 13.84 um, not the old hand-waved 5: it is the MEAN ion birth
                 # height measured in the S3 v2 pass (audit C12). Only
                 # --ion-model analytic consumes it; the measured template
                 # remains the production default and is unaffected.
                 with_ions=True, z_aval_um=13.84, mu_ion=ION.MU_ION_CM2_VS,
                 ion_model="measured", kernel_t_max_ns=None,
                 y_window_mm=None):
        # kernel_t_max_ns exists so the LUT window can be varied without
        # editing the default — the A1 before/after (test_window) and any
        # SPS-config run (64 x 60 ns frame -> 4200) both need it.
        # n_side is ONE quantity, not two. The LUT's X band is built with
        # 2*n_side+1 channel offsets and `induce` books that same range, so a
        # Digitizer asking for more channels than the band holds indexes off
        # the end of it. They are wired from the same argument here.
        #
        # y_window_mm must grow WITH n_side or widening does nothing on the Y
        # view: Y channels sit on the 0.78 mm pad pitch, so a 3.9 mm half-window
        # holds only +-5 rows however large n_side is.
        lut_kw = {"n_side": n_chan_side}
        if kernel_t_max_ns is not None:
            lut_kw["t_max_ns"] = float(kernel_t_max_ns)
        if y_window_mm is not None:
            lut_kw["y_window_mm"] = float(y_window_mm)
        elif n_chan_side > 4:
            # Default it to just past the outermost channel asked for.
            lut_kw["y_window_mm"] = (n_chan_side + 1) * 0.78
        self.lut = CombKernelLUT(kernel_path, **lut_kw)
        self.gas = (DriftGas(gas_table, v_scale=v_scale) if gas_table
                    else DriftGas(v_scale=v_scale))
        self.calib, self.calib_v = load_calib(calib_path, mesh_v)
        self.E_drift = drift_v / (drift_gap_mm * 0.1)      # V/cm
        self.drift_gap_mm = drift_gap_mm
        self.transparency = transparency
        self.n_side = n_chan_side
        self.packet = packet
        self.rng = np.random.default_rng(seed)
        # Avalanches landing off the readout, counted for provenance (C14).
        self.n_outside = 0
        self.n_seen = 0

        p = self.calib["polya"]
        self.gbar, self.theta = p["gain_mean"], p["theta"]
        # P(g>0) from the S3 calib; absent in schema <= 2, where it is 1.0.
        self.aval_survival = float(p.get("survival", 1.0))
        self.sigma0_um = self.calib["sigma0_um"]
        self.v_drift = float(np.ravel(
            self.gas.v_drift_um_ns(self.E_drift))[0])

        # The analytic description is built either way: when ion_model is
        # "measured" it is no longer what drives the LUT, but it stays in the
        # record as the first-principles prediction the measurement is
        # checked against (see ions.py).
        self.with_ions = with_ions
        self.ion_model = ion_model
        self.ion = ION.describe(C.AMP_GAP_M, mesh_v, mu_ion, z_aval_um * 1e-6)
        self.ion_measured = None
        if with_ions and ion_model == "measured":
            # PRESENCE IS NOT VALIDITY. The shipped schema-1 calib carries
            # `i_elec` / `i_ion` keys whose arrays are 2000 zeros, so a
            # key-existence check passes, `measured_longitudinal` divides by a
            # zero area, and the entire LUT becomes nan — silently, because nan
            # currents still sum, still write, and still produce a decoded file.
            # Found 2026-08-07 while re-running the charge audit. Check the
            # CONTENT.
            tmpl = [np.asarray(self.calib.get(k, []), dtype=float)
                    for k in ("i_elec", "i_ion")]
            if (any(t.size == 0 for t in tmpl)
                    or not np.isfinite(np.concatenate(tmpl)).all()
                    or abs(float(sum(t.sum() for t in tmpl))) == 0.0):
                raise SystemExit(
                    "--ion-model measured needs an S3 v2 calib whose i_elec / "
                    "i_ion templates are actually populated; this one's are "
                    f"missing, all-zero or non-finite in {calib_path}. Use "
                    "--ion-model analytic, or point --calib at "
                    "aval_calib_v2.json.")
            h, dt_h, self.ion_measured = ION.measured_longitudinal(self.calib)
            self.lut.apply_longitudinal(h, dt_h)
        elif with_ions:
            self.lut.apply_ion_transit(self.ion["f_electron"],
                                       self.ion["t_ion_transit_ns"])

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

        # ── Deposits born INSIDE the amplification gap (audit C7) ────────────
        # `clusters.py` keeps AmpGas as well as DriftGas, and its contract says
        # they "skip the drift and are amplified where they sit". The code did
        # not honour that: it clipped z < 0 to 0 and then handed them the FULL
        # gap gain AND the mesh transparency, i.e. ~9.3x too much charge, all
        # of it prompt and unshared on d = 0. In the production muon file that
        # is 2204 clusters / 4234 electrons (0.90 % / 0.41 %), contributing
        # 0.41 % of the charge where they should contribute 0.04 %.
        #
        # The gap runs z = 0 (ESL/anode) to z = AMP_GAP (mesh), and an electron
        # multiplies over the distance it actually travels to the anode, so its
        # mean gain is exp(alpha z) with alpha = ln(G)/gap. Deposits at z < 0
        # are the ESL groove deposits NEEDED_INPUTS describes; they sit BELOW
        # the anode plane, have no gap to cross, and fall out of the same
        # formula as gain < 1 — i.e. no multiplication — with no special case.
        gap_mm = C.AMP_GAP_M * 1e3
        in_gap = zs < gap_mm
        # Mean gain relative to a full-gap avalanche, 1.0 for drift electrons.
        gain_frac = np.where(
            in_gap, np.exp(np.log(self.gbar) * np.clip(zs, None, gap_mm)
                           / gap_mm) / self.gbar, 1.0)

        # Drift applies only to what actually drifted. An in-gap deposit has no
        # drift length, so no transverse spread and no drift delay: clipping
        # its z to 0 for the gas lookup gives exactly that, and is why the clip
        # is kept rather than removed.
        zs_drift = np.clip(zs, 0.0, None)
        sT = self.gas.sigma_T_um(self.E_drift, zs_drift)
        st = self.gas.sigma_t_ns(self.E_drift, zs_drift)

        x = xs * 1e-3 + self.rng.normal(0.0, sT * 1e-6)
        y = ys * 1e-3 + self.rng.normal(0.0, sT * 1e-6)
        t = ts + zs_drift * 1e3 / self.v_drift + self.rng.normal(0.0, st)

        # Mesh transparency AND drift attachment, as ONE thinning (plan §7
        # step 2, audit A6). Combining them keeps the statistics binomial and
        # the truth accounting sees a single loss channel with two named
        # factors, instead of two Bernoullis whose product is the same thing
        # with more code. p_surv = eps_mesh * exp(-eta z) * P(avalanche > 0);
        # with a table that has no eta column (the dry production table) the
        # second factor is exactly 1 and results are bit-identical at fixed seed.
        #
        # The THIRD factor closes the conditional-Polya decomposition (audit
        # A7): `polya_sample` draws from P(g | g > 0), so seeds that produce no
        # avalanche at all must be removed here or the per-electron charge is
        # high by 1/P(g>0). Folding it into the SAME thinning rather than adding
        # a second Bernoulli keeps the statistics binomial and consumes no extra
        # randoms, so a calib with survival = 1 is bit-identical. Measured over
        # all 56 S3 raw slices, survival IS 1.0 at every voltage (0 of 6400
        # seeds failed to multiply) — design/report/DESKTOP_RUNS_2026-08-07.md —
        # so this is exact bookkeeping today, not a correction. A v2 calib has
        # no `survival` field, hence the 1.0 default.
        #
        # In-gap deposits are already PAST the mesh and never drifted, so
        # neither the transparency nor the attachment applies to them (C7).
        p_surv = np.where(
            in_gap, self.aval_survival,
            self.transparency * self.gas.survival(self.E_drift, zs_drift)
            * self.aval_survival)
        if self.packet:
            surv = self.rng.binomial(n_e, p_surv).astype(float)
        else:
            surv = (self.rng.random(len(w)) < p_surv).astype(float)
        keep = surv > 0
        x, y, t, surv = x[keep], y[keep], t[keep], surv[keep]
        gain_frac = gain_frac[keep]
        x += self.rng.normal(0.0, self.sigma0_um * 1e-6, len(x))
        y += self.rng.normal(0.0, self.sigma0_um * 1e-6, len(y))

        # PACKET MODE DRAWS THE SUM, NOT A SCALED SINGLE DRAW (fix 2026-08-07,
        # audit C3). `surv` electrons sharing one packet each avalanche
        # INDEPENDENTLY, so the packet's gain is a sum of n independent Polyas
        # = Gamma(n(1+theta), gbar/(1+theta)), with variance n*Var(g).
        # Multiplying one single-electron draw by n instead gives n^2*Var(g) —
        # the mean is right, so nothing in the budget noticed, but the
        # avalanche-to-avalanche fluctuation was n times too wide. Production
        # runs packet=False and is unaffected; this matters the moment §12
        # item 11 turns packet mode on for speed.
        if self.packet:
            gain = self.rng.gamma(shape=surv * (1.0 + self.theta),
                                  scale=self.gbar / (1.0 + self.theta))
        else:
            gain = polya_sample(self.rng, self.gbar, self.theta, len(x)) * surv
        # C7: scale to the gap the electron actually crossed. Applied to the
        # DRAW rather than by re-parameterising the Polya, so the fluctuation
        # keeps its shape and the drift electrons (gain_frac == 1) are
        # bit-identical. A shorter avalanche is really somewhat broader in
        # relative terms — fewer generations — which this does not model; at
        # 0.4 % of the electrons that is far below anything observable.
        return x, y, t, gain * gain_frac

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
        nd, nt = len(ds), len(self.lut.t)

        # Channel index of the pad nearest each avalanche, on the 0.78 mm
        # lattice. The COLUMN comes from the LUT (lut.col_at), not from a second
        # rounding of x here: the X band's d = 0 is defined per LUT x sample, and
        # deriving the channel number independently put ~1.3 % of avalanches one
        # pad off (audit C2). The ROW is unaffected — the Y kernels are indexed
        # by a y OFFSET, not by an absolute row, so there is only one rounding.
        ix = self.lut.ix(x)
        col = self.lut.col_at(x, ix=ix)
        row = np.rint((y - K.PAD_ORIGIN_M) / pitch).astype(int)
        k0 = np.rint(t / (self.lut.dt * 1e9)).astype(int)

        # FIDUCIAL: charge landing off the board is not readout charge (audit
        # C14). `induce` returns every channel it computed, but run.py's DAQ
        # writer drops anything outside 0..511 when it lays the dense plane —
        # so the budget normalised over channels the decoded file does not
        # contain, and the two paths disagreed for out-of-area delta rays.
        # Dropping the avalanche here makes them agree, and the count is
        # reported rather than swallowed.
        inside = ((col >= -self.n_side) & (col < C.PAD_N + self.n_side)
                  & (row >= -self.n_side) & (row < C.PAD_N + self.n_side))
        self.n_outside += int((~inside).sum())
        self.n_seen += int(inside.size)

        ok = (np.isfinite(q) & (q > 0) & (k0 >= 0) & (k0 < n_samp)
              & inside)
        idx = np.flatnonzero(ok)
        if idx.size == 0:
            return out

        dy_step = self.lut.y_Y[1] - self.lut.y_Y[0]
        dyx_step = self.lut.y_X[1] - self.lut.y_X[0]

        # Compact channel numbering, so the scatter target is a dense 2-D
        # (channel, time) array that np.bincount can fill in one call.
        rows_all = row[idx][:, None] + ds[None, :]
        cols_all = col[idx][:, None] + ds[None, :]
        keyY = np.unique(rows_all)
        keyX = np.unique(cols_all)
        nY, nX = len(keyY), len(keyX)
        slotY = np.searchsorted(keyY, rows_all)
        slotX = nY + np.searchsorted(keyX, cols_all)
        acc = np.zeros((nY + nX) * n_samp)

        # Chunked so the (chunk, nd, nt) gather stays a few tens of MB. The
        # whole point is that the scatter happens ONCE per chunk via bincount
        # rather than once per avalanche per channel: the arithmetic here is
        # ~3e9 flops for a 500-event file and is memory-bandwidth bound, but a
        # python loop over avalanches spends 3.4M interpreter round-trips doing
        # a few kB of work each, which is what made it minutes instead of
        # seconds.
        step = max(1, int(4e6 // (nd * nt)))
        tgrid = np.arange(nt)
        for s in range(0, len(idx), step):
            sl = idx[s:s + step]
            n = len(sl)
            r = row[sl][:, None] + ds[None, :]
            dy = y[sl][:, None] - (K.PAD_ORIGIN_M + r * pitch)
            iy = np.clip(np.rint((dy - self.lut.y_Y[0]) / dy_step).astype(int),
                         0, len(self.lut.y_Y) - 1)
            curY = self.lut.I_Y[r % 2, iy, ix[sl][:, None], :]      # (n, nd, nt)

            y_rel = np.mod(y[sl] - K.PAD_ORIGIN_M, 2 * pitch)
            jy = np.clip(np.rint((y_rel - self.lut.y_X[0]) / dyx_step
                                 ).astype(int), 0, len(self.lut.y_X) - 1)
            curX = np.moveaxis(self.lut.I_X[:, jy, ix[sl], :], 1, 0)

            # Time index of every sample, clipped INTO the last row and then
            # masked, so the scatter never wraps a late avalanche around to
            # t=0 (which would show up as a phantom prompt signal).
            tt = k0[sl][:, None, None] + tgrid[None, None, :]
            good = tt < n_samp
            tt = np.where(good, tt, 0)
            w = (q[sl][:, None, None] * good)

            for slot, cur in ((slotY[s:s + step], curY),
                              (slotX[s:s + step], curX)):
                flat = (slot[:, :, None] * n_samp + tt).ravel()
                acc += np.bincount(flat, weights=(cur * w).ravel(),
                                   minlength=acc.size)

        acc = acc.reshape(nY + nX, n_samp)
        for j, r in enumerate(keyY):
            out[("Y", int(r))] = acc[j]
        for j, c in enumerate(keyX):
            out[("X", int(c))] = acc[nY + j]
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
            "transparency_source": "T6 measured (2026-08-07), bench ratio 98",
            "gain_mean": self.gbar, "polya_theta": self.theta,
            "aval_survival": self.aval_survival,
            "aval_survival_source": (
                "S3 calib" if "survival" in self.calib.get("polya", {})
                else "absent from calib (schema <= 2), assumed 1.0"),
            "sigma0_um": self.sigma0_um,
            "calib_voltage_V": self.calib_v,
            "field_model": self.calib.get("field_model"),
            "with_ions": self.with_ions,
            "ion_model": self.ion_model,
            "ion": self.ion,
            "ion_measured": self.ion_measured,
            "packet_mode": self.packet,
            "n_chan_side": self.n_side,
            "avalanches_outside_readout": self.n_outside,
            "avalanches_seen": self.n_seen,
            "y_window_mm": float(abs(self.lut.y_Y).max() * 1e3),
        }
        return d
