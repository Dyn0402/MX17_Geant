#!/usr/bin/env python3
"""
gas.py — drift transport parameters for the digitizer (plan §7 step 2).

Reads the Magboltz table computed for THIS gas at THIS pressure by
`nTof_x17/garfield_sim/mm_drift_velocity_table.py`, which is the only source
covering the drift-field range (the S3 gas tables cover only the amplification
grid). Interpolated in E, never extrapolated silently.

Units in the table: v in µm/ns; dL, dT in sqrt(cm), i.e. the diffusion after
drifting z is sigma = d * sqrt(z) with both in cm.

THE DRIFT VELOCITY IS NOT SETTLED, and the digitizer must not pretend it is.
Three numbers are in circulation for the bench:

  * dry Magboltz at the bench field (1000 V / 30 mm = 333 V/cm): 39.1 µm/ns,
    which is this table;
  * "36.6 µm/ns" quoted in plan §1 as the measured value;
  * 28.1 +/- 0.7 µm/ns, the bias-free micro-TPC measurement in nTof_x17
    `mx_june_cosmic_qa/HANDOFF_det3_vdrift_and_kernels.md` (2026-07-01), which
    supersedes an earlier sigma-scan estimate that was inflated by resistive
    charge spreading. That one is for det3's ~19 mm gap, so it is not directly
    the same operating point.

The gap is real and is most likely water contamination (plan §1 flags ~1% H2O).
`DriftGas` therefore takes an explicit `v_scale` and records it, so every
product says which choice it was made under. Do NOT bake a preference in here.
"""

from __future__ import annotations

import json
import os

import numpy as np

TABLE_NAME = "drift_velocity_Ar_iC4H10_95_5_Saclay.json"

# WHERE THE GAS TABLE LIVES IS NOT A CONSTANT (fix 2026-08-07, audit C11).
# The path used to be hardcoded into a checkout under $HOME, which exists on
# this laptop and on no condor worker; a batch job died on it rather than
# falling back. Resolution order: explicit argument, then $MX17_GAS_TABLE, then
# a table shipped beside the repo, then the historical location.
ENV_VAR = "MX17_GAS_TABLE"
_HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = (
    os.path.join(_HERE, "..", "..", "design", "gas", TABLE_NAME),
    os.path.expanduser(f"~/PycharmProjects/nTof_x17/garfield_sim/results/{TABLE_NAME}"),
)


def default_table():
    """First existing candidate, or the historical path so the error names it."""
    env = os.environ.get(ENV_VAR)
    if env:
        return os.path.expanduser(env)
    for p in _CANDIDATES:
        if os.path.exists(os.path.normpath(p)):
            return os.path.normpath(p)
    return os.path.normpath(_CANDIDATES[-1])


# Kept as a module attribute for callers that print it, but resolved lazily by
# DriftGas so that setting the env var after import still works.
DEFAULT_TABLE = default_table()


def _table_hash(path):
    """Content hash, so provenance pins the TABLE and not just its name."""
    import hashlib
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()[:16]
    except OSError:
        return None


class DriftGas:
    """Drift velocity and diffusion vs field, interpolated from Magboltz."""

    def __init__(self, path=None, v_scale=1.0):
        path = os.path.expanduser(path) if path else default_table()
        if not os.path.exists(path):
            raise SystemExit(
                f"gas table not found at {path}. Point ${ENV_VAR} at it, or "
                "pass --gas-table; the response chain will not guess a gas.")
        with open(path) as f:
            d = json.load(f)
        self.table_sha256 = _table_hash(path)
        pts = [p for p in d["points"] if p.get("ok", True)]
        pts.sort(key=lambda p: p["E_Vcm"])
        self.path = path
        self.gas = d.get("gas", "?")
        self.pressure_torr = d.get("pressure_torr")
        self.E = np.array([p["E_Vcm"] for p in pts])
        self._v = np.array([p["v_um_per_ns"] for p in pts])
        self._dL = np.array([p["dL_sqrtcm"] for p in pts])
        self._dT = np.array([p["dT_sqrtcm"] for p in pts])
        # ATTACHMENT (plan §7 step 2, audit A6). Optional: the dry Ar/iC4H10
        # table has no eta key at all, and a table that lacks it means "this
        # gas was characterised without attachment", NOT "attachment is zero" —
        # so the branch taken is recorded in provenance rather than assumed.
        self._eta = (np.array([p["eta_per_cm"] for p in pts])
                     if all("eta_per_cm" in p for p in pts) else None)
        self.has_attachment = self._eta is not None
        # A table CAN carry the key and still be all zeros — the wet-CF4 suite
        # is exactly that (see design/report/WET_CF4_BRACKET_2026-08-07.md), so
        # distinguish "no column" from "column present and identically zero".
        self.attachment_is_zero = bool(self.has_attachment
                                       and not np.any(self._eta > 0))
        # An overall factor on v only. Diffusion is left alone: contamination
        # that slows the drift does not scale sigma the same way, and faking it
        # here would hide the effect we most want to test.
        self.v_scale = float(v_scale)

    def _interp(self, y, E_Vcm):
        E = np.atleast_1d(np.asarray(E_Vcm, dtype=float))
        if (E < self.E[0]).any() or (E > self.E[-1]).any():
            raise ValueError(
                f"drift field {E.min():.1f}-{E.max():.1f} V/cm outside the "
                f"table's {self.E[0]:.1f}-{self.E[-1]:.1f} V/cm. Extend the "
                "Magboltz table rather than extrapolating.")
        return np.interp(E, self.E, y)

    def v_drift_um_ns(self, E_Vcm):
        return self._interp(self._v, E_Vcm) * self.v_scale

    def eta_per_cm(self, E_Vcm):
        """
        Attachment coefficient [1/cm]. Returns 0 when the table has no column
        for it, so the dry chain is bit-for-bit unchanged.
        """
        if self._eta is None:
            return np.zeros_like(np.atleast_1d(np.asarray(E_Vcm, dtype=float)))
        return self._interp(self._eta, E_Vcm)

    def survival(self, E_Vcm, z_mm):
        """
        Fraction of electrons surviving a drift of z, exp(-eta * z).

        z-DEPENDENT BY CONSTRUCTION, which is the whole reason it cannot be
        folded into a flat gain factor: it re-weights the DEPTH MIXTURE, and
        the depth mixture is what produces the c1-vs-z diffusion signature.
        """
        eta = np.atleast_1d(self.eta_per_cm(E_Vcm))
        return np.exp(-eta * np.asarray(z_mm, dtype=float) * 0.1)   # mm -> cm

    def sigma_T_um(self, E_Vcm, z_mm):
        """Transverse spread [µm] after drifting z_mm."""
        return self._interp(self._dT, E_Vcm) * np.sqrt(np.asarray(z_mm) / 10.0) * 1e4

    def sigma_L_um(self, E_Vcm, z_mm):
        """Longitudinal spread [µm] after drifting z_mm."""
        return self._interp(self._dL, E_Vcm) * np.sqrt(np.asarray(z_mm) / 10.0) * 1e4

    def sigma_t_ns(self, E_Vcm, z_mm):
        """Arrival-time spread [ns] — the longitudinal spread over v."""
        return self.sigma_L_um(E_Vcm, z_mm) / self.v_drift_um_ns(E_Vcm)

    def describe(self, E_Vcm, z_mm=30.0):
        v = float(np.ravel(self.v_drift_um_ns(E_Vcm))[0])
        return {
            "gas": self.gas, "table": os.path.basename(self.path),
            # Full path AND content hash: two runs that name the same table are
            # not the same run if the table was regenerated between them.
            "table_path": self.path, "table_sha256": self.table_sha256,
            # Which attachment branch was taken — "no column" and "column of
            # zeros" are different statements and must not look alike later.
            "attachment": ("absent from table" if not self.has_attachment
                           else "present but identically zero"
                           if self.attachment_is_zero else "active"),
            "eta_per_cm": float(np.ravel(self.eta_per_cm(E_Vcm))[0]),
            "survival_full_gap": float(np.ravel(
                self.survival(E_Vcm, z_mm))[0]),
            "E_Vcm": float(E_Vcm), "v_scale": self.v_scale,
            "v_um_per_ns": v,
            "sigma_T_um_full_gap": float(np.ravel(self.sigma_T_um(E_Vcm, z_mm))[0]),
            "sigma_L_um_full_gap": float(np.ravel(self.sigma_L_um(E_Vcm, z_mm))[0]),
            "sigma_t_ns_full_gap": float(np.ravel(self.sigma_t_ns(E_Vcm, z_mm))[0]),
            "t_cross_ns": float(z_mm * 1e3 / v),
        }


if __name__ == "__main__":
    g = DriftGas()
    import pprint
    pprint.pprint(g.describe(1000.0 / 3.0))
