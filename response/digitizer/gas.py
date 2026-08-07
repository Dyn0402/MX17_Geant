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

DEFAULT_TABLE = os.path.expanduser(
    "~/PycharmProjects/nTof_x17/garfield_sim/results/"
    "drift_velocity_Ar_iC4H10_95_5_Saclay.json")


class DriftGas:
    """Drift velocity and diffusion vs field, interpolated from Magboltz."""

    def __init__(self, path=DEFAULT_TABLE, v_scale=1.0):
        with open(path) as f:
            d = json.load(f)
        pts = [p for p in d["points"] if p.get("ok", True)]
        pts.sort(key=lambda p: p["E_Vcm"])
        self.path = path
        self.gas = d.get("gas", "?")
        self.pressure_torr = d.get("pressure_torr")
        self.E = np.array([p["E_Vcm"] for p in pts])
        self._v = np.array([p["v_um_per_ns"] for p in pts])
        self._dL = np.array([p["dL_sqrtcm"] for p in pts])
        self._dT = np.array([p["dT_sqrtcm"] for p in pts])
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
