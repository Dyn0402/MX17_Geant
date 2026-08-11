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
# Mesh electron transparency for a calibration that does NOT contain one.
#
# MEASURED, T6 3-D (response/meshcell/FIELD_MAP_RUNBOOK.md, gate G7, production
# map accepted 2026-08-08, git 5763342): 0.955 -0.045/+0.005 at the bench point
# (E_drift 333 V/cm over the 30 mm gap against ~31.0 kV/cm in the amplification
# gap). Higher than the 1-D model because electrons dodge the wires in both
# lateral dimensions in 3-D.
#
# CORRECTION 2026-08-11: this constant was 0.873 — the 1-D value from
# `mesh_transparency.C` (2026-08-07) — which the runbook's acceptance note said
# outright the 3-D deliverable "supersedes ... on Stage B's next touch". Stage B
# was never touched, so the 1-D number stayed live for three days. See
# design/report/TRANSPARENCY_DOUBLE_COUNT_2026-08-11.md.
DEFAULT_TRANSPARENCY = 0.955


# The pooled meshfield calib is a DIFFERENT container schema, on purpose: all
# 56 of its slices used the same pre-solved 490 V field map regardless of their
# --voltage label, so it is one bench-point measurement rather than a voltage
# scan, and it stores a single `point` instead of a `points` map. Its voltage
# is stated only in its free-text `note`, not in a machine-readable field —
# hence this constant and the hard guard below. Applying a 490 V avalanche
# calibration at some other mesh voltage would be silently wrong, so it raises.
POOLED_SCHEMA_PREFIX = "aval_calib_meshfield_pooled"
POOLED_MESH_V = 490.0


# Height above the mesh at which a MESHFIELD S3 calibration launches its seed
# electrons, so that transparency and funnelling are measured rather than
# assumed: `seed_z0_cm = 180.0e-4` in response/avalanche/mx17_aval_calib.py.
# Its sigma0 therefore already contains the transverse diffusion of that leg,
# and `transport()` must not drift through it again.
# NOT machine-readable in the calib JSON today — same schema gap as the pooled
# file's mesh voltage. If the emitter ever exports `seed_z0_um`, prefer it and
# delete this constant.
MESHFIELD_SEED_Z0_UM = 180.0


def calib_seed_z0_mm(calib_pt):
    """Drift length [mm] already baked into this calibration's sigma0.

    Prefers the calibration's own `seed_z0_um`, which the S3 emitter records
    per point (2026-08-09). Falls back to the field_model heuristic for files
    written before that: uniform_field seeds sit INSIDE the amplification gap
    so there is no drift leg, meshfield seeds sit MESHFIELD_SEED_Z0_UM above
    the mesh. field_model, not the container schema, is the fallback key
    because it is what actually determines the geometry.
    """
    z0 = calib_pt.get("seed_z0_um")
    if z0 is not None:
        return float(z0) * 1e-3
    fm = str(calib_pt.get("field_model", ""))
    if fm.startswith("meshfield"):
        return MESHFIELD_SEED_Z0_UM * 1e-3
    if fm and fm != "uniform_field":
        raise ValueError(
            f"unrecognised calib field_model {fm!r} and no seed_z0_um: "
            "refusing to guess whether its sigma0 includes a drift leg")
    return 0.0


def split_calib_survival(calib_pt):
    """Split the calib's `survival` into (mesh transparency, P(g > 0)).

    THE CONSTRAINT THIS FUNCTION EXISTS TO ENFORCE (fix 2026-08-11): the S3
    calibration writes `survival = (gains > 0).mean()` over its seeds
    (`mx17_aval_calib.py:460`), and what that fraction MEANS depends on where
    the seeds were launched:

      * meshfield  — seeds start MESHFIELD_SEED_Z0_UM above the mesh, so a seed
        absorbed on a wire records gain 0. `survival` is therefore
        eps_mesh * P(g>0), and it IS the transparency: the production file
        `aval_calib_meshfield_pooled.json` reads 0.9559 +- 0.0026 against T6's
        independent 3-D G7 value of 0.955, and P(g>0) is 1.0 to the precision
        of the uniform-field campaign (0 of 6400 seeds failed to multiply, all
        56 slices — design/report/DESKTOP_RUNS_2026-08-07.md). So the external
        DEFAULT_TRANSPARENCY must NOT be applied on top: doing so thins twice.
      * uniform_field — seeds start INSIDE the amplification gap, past a mesh
        that is not in the model at all. `survival` is P(g>0) alone (= 1.0
        measured, and absent entirely from schema <= 2), and the transparency
        has to come from outside, because nothing in the calibration knows
        about a mesh.

    CORRECTION 2026-08-11: before this split the chain applied
    DEFAULT_TRANSPARENCY * survival unconditionally = 0.873 * 0.9559 = 0.8345
    on the production meshfield path, where the physics is a single 0.955.
    Simulated charge was low by x0.874. The `RESPONSE_SIM_PLAN.md` T7 line
    "survival 0.9559 ~= T6 transparency 0.955" was read as an independent
    cross-check of two numbers; it is one number measured twice, and their
    agreement is the duplication. Ledger consequences:
    design/report/TRANSPARENCY_DOUBLE_COUNT_2026-08-11.md.

    Returns (eps_mesh_or_None, p_multiply). A None transparency means "this
    calibration does not contain one, use the external default".
    """
    surv = float(calib_pt.get("polya", {}).get("survival", 1.0))
    fm = str(calib_pt.get("field_model", ""))
    if fm.startswith("meshfield"):
        return surv, 1.0
    return None, surv


def load_calib(path, mesh_v):
    """Pick the S3 point nearest the requested mesh voltage.

    Handles both container schemas: the `points` map keyed `<gas>@<V>V`, and
    the pooled single-`point` meshfield file (see POOLED_SCHEMA_PREFIX).
    """
    with open(path) as f:
        doc = json.load(f)

    if "points" not in doc and "point" in doc:
        schema = str(doc.get("schema", ""))
        if not schema.startswith(POOLED_SCHEMA_PREFIX):
            raise ValueError(
                f"{path} has a single 'point' but an unrecognised schema "
                f"{schema!r}; refusing to guess its mesh voltage")
        pt = doc["point"]
        # Prefer the calibration's own voltage (emitter records `voltage_V`
        # per point since 2026-08-09); fall back to the pinned constant for
        # files written before that. Keyed on the POINT's field, not the
        # container schema: `aval_calib/3` is shared between meshfield and
        # uniform_field files and only meshfield-produced points carry it.
        pooled_v = float(pt.get("voltage_V", POOLED_MESH_V))
        if abs(mesh_v - pooled_v) > 1e-6:
            raise ValueError(
                f"{path} is a pooled {pooled_v:.0f} V bench point with no "
                f"per-voltage data, but {mesh_v:.0f} V was asked for. Use a "
                f"voltage-scan calib for that point, or run at "
                f"{pooled_v:.0f} V.")
        print(f"  [calib] pooled meshfield point, {pooled_v:.0f} V "
              f"({pt.get('n_slices')} slices, {pt.get('nev_total')} events"
              + (f", {pt['gas_file']}" if pt.get("gas_file") else "") + ")")
        return pt, pooled_v

    cal = doc["points"]
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
                 # None = take it from the calibration if the calibration
                 # measured one (meshfield), else DEFAULT_TRANSPARENCY. Pass a
                 # number only to override a calibration deliberately — see
                 # split_calib_survival, and note that on a meshfield calib an
                 # explicit value multiplies the one already inside `survival`.
                 transparency=None, v_scale=1.0,
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
        self.n_side = n_chan_side
        self.packet = packet
        self.rng = np.random.default_rng(seed)
        # Avalanches landing off the readout, counted for provenance (C14).
        self.n_outside = 0
        self.n_seen = 0

        p = self.calib["polya"]
        self.gbar, self.theta = p["gain_mean"], p["theta"]
        # The calib's single `survival` number splits into a mesh term and an
        # avalanche term, and WHICH of the two it is depends on the field model
        # — see split_calib_survival for the constraint and for the 2026-08-11
        # double-count it exists to prevent. Whatever the split, exactly one
        # mesh transparency reaches p_surv below.
        calib_eps, self.aval_survival = split_calib_survival(self.calib)
        self.calib_transparency = calib_eps
        if transparency is not None:
            self.transparency = float(transparency)
            self.transparency_source = "caller override"
        elif calib_eps is not None:
            self.transparency = calib_eps
            self.transparency_source = (
                f"S3 meshfield calib `survival` "
                f"({self.calib.get('field_model')}); T6 3-D G7 measures the "
                f"same quantity at 0.955")
        else:
            self.transparency = DEFAULT_TRANSPARENCY
            self.transparency_source = (
                "T6 3-D G7 measured (production map accepted 2026-08-08); "
                "calib is uniform_field and contains no mesh")
        self.sigma0_um = self.calib["sigma0_um"]
        self.calib_seed_z0_mm = calib_seed_z0_mm(self.calib)
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
        # TRANSVERSE diffusion stops where the calibration's own seeds start,
        # or the overlapping leg is counted twice. A meshfield calib launches
        # its seeds `seed_z0` ABOVE the mesh (mx17_aval_calib.py: 180 µm) so
        # that transparency and funnelling are measured rather than assumed —
        # which means its sigma0 already contains the diffusion of that last
        # leg. A uniform_field calib seeds INSIDE the amplification gap and
        # contains no drift leg at all, hence seed_z0 = 0 there.
        # Only sigma_T is affected: the drift TIME and the attachment survival
        # below are over the full path, and the calib contributes neither
        # (`t_arrival_mean_ns` is not consumed by this class).
        # Since sigma_T^2 is linear in z, dropping the leg from the path length
        # is exactly equivalent to subtracting its variance.
        zs_diff = np.clip(zs_drift - self.calib_seed_z0_mm, 0.0, None)
        sT = self.gas.sigma_T_um(self.E_drift, zs_diff)
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
        # randoms, so a calib with survival = 1 is bit-identical.
        #
        # CORRECTED 2026-08-11. This comment used to say the third factor is
        # 1.0 and therefore exact bookkeeping — true of the UNIFORM-FIELD S3
        # calibs it was written against (0 of 6400 seeds failed to multiply
        # across all 56 slices, design/report/DESKTOP_RUNS_2026-08-07.md), and
        # false of the meshfield calib that has been production since 08-08,
        # whose `survival` = 0.9559 is the MESH TRANSPARENCY, not P(g>0). Read
        # literally, the old comment licensed multiplying the two — which is
        # what the code did, thinning at 0.873 * 0.9559 = 0.8345 for a physics
        # of 0.955. `split_calib_survival` now routes each calibration's number
        # to the factor it actually measures, so `self.transparency` is the ONE
        # mesh term here and `self.aval_survival` is the P(g>0) term, 1.0 on
        # the meshfield path.
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
            "transparency_source": self.transparency_source,
            "gain_mean": self.gbar, "polya_theta": self.theta,
            "aval_survival": self.aval_survival,
            "aval_survival_source": (
                "1.0: meshfield calib `survival` is the mesh term, reported "
                "as mesh_transparency (P(g>0) = 1.0 measured, 0/6400 seeds)"
                if self.calib_transparency is not None
                else "S3 calib" if "survival" in self.calib.get("polya", {})
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
