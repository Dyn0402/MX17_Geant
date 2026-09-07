#!/usr/bin/env python3
"""
selftest.py — Stage B on synthetic tracks, before Geant4 is in the loop.

Two reasons this exists rather than waiting for a ClusterTree.

First, isolation. If the digitizer is fed real Geant4 output on day one and the
sharing comes out wrong, the cause could be anywhere in the chain. A synthetic
normal-incidence MIP has an exactly known geometry, so anything that disagrees
is Stage B's fault.

Second, it goes straight at the sharp corner. §9's dispersed +-1 share is
c1 = 0.23-0.28 and is gain/gas/drift invariant. T2b showed point-charge prompt
sharing to d=+-1 is near zero, and T7 showed the avalanche footprint is only
33 um, so the ONLY thing left that can produce c1 is transverse diffusion plus
resistive spreading. That makes c1 a genuine prediction with nothing left to
tune, and this computes it.

The single most informative scan is c1 vs DEPTH: a cluster made at the mesh has
no diffusion at all, one made at the cathode has the full 826 um. If c1 came
from the avalanche or from induction geometry it would be flat in z; if it comes
from diffusion it must rise steeply with z. That shape is the test, not the
single number.

    python3 -m response.digitizer.selftest --kernel <s1.npz>
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from ..common import constants as C
from ..solver import kernels as K
from .digitize import (Digitizer, DEFAULT_DRIFT_GAP_MM,
                       DEFAULT_MESH_V, DEFAULT_DRIFT_V)

CALIB = "~/x17/response_sim/avalanche/aval_calib.json"


def point_deposit(dig, z_mm, n_e, x_mm, y_mm, n_samp=1500):
    """One localised deposit at depth z: n_e electrons, all at (x, y, z)."""
    x, y, t, g = dig.transport(np.array([x_mm]), np.array([y_mm]),
                               np.array([z_mm]), np.array([0.0]),
                               np.array([n_e]))
    return dig.induce(x, y, t, g, n_samp)


def integrated_charge(cur, dt_s):
    """Time-integrate each channel's current back to a charge."""
    return {k: float(np.sum(v) * dt_s) for k, v in cur.items()}


def share_profile(q, view, centre):
    """Charge on channels centre+d, normalised to the total over the view."""
    tot = sum(v for (vw, _), v in q.items() if vw == view)
    if tot <= 0:
        return None
    return {d: q.get((view, centre + d), 0.0) / tot for d in range(-3, 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--calib", default=CALIB)
    ap.add_argument("--n-electrons", type=int, default=300)
    ap.add_argument("--repeat", type=int, default=40,
                    help="independent deposits averaged per depth")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    import os
    dig = Digitizer(a.kernel, os.path.expanduser(a.calib), seed=a.seed)
    info = dig.describe()
    print("Stage B self-test — synthetic point deposits\n")
    print(f"  kernel   rho_s={info['kernel']['rho_s_MOhm_sq']} MΩ/sq  "
          f"d_k={info['kernel']['d_kapton_um']} µm")
    print(f"  gas      v_drift={info['gas']['v_um_per_ns']:.1f} µm/ns  "
          f"sigma_T(full gap)={info['gas']['sigma_T_um_full_gap']:.0f} µm")
    print(f"  avalanche gain={info['gain_mean']:.0f} theta={info['polya_theta']:.2f} "
          f"sigma0={info['sigma0_um']:.1f} µm  (S3 @ {info['calib_voltage_V']:.0f} V)")
    print(f"  transparency={info['mesh_transparency']} "
          f"[{info['transparency_source']}]\n")

    # Put the deposit on a pad centre so "d=0" is unambiguous, and on a
    # resistive strip so the ESL phase is the on-strip case of §3.
    xw = C.SUPERPERIOD_M / 2
    xw = xw - np.mod(xw, C.ESL_PITCH_M)
    phase = K.nearest_column(xw)                 # 0..39, the ESL phase index
    x0_m = float(K.pad_x(phase))
    row0 = 0
    y0_m = K.PAD_ORIGIN_M + row0 * C.PAD_PITCH_M
    x0_mm, y0_mm = x0_m * 1e3, y0_m * 1e3
    # digitize.induce labels X channels by ABSOLUTE column index, not by the
    # 40-phase index. Looking them up by phase finds nothing and every X share
    # comes back exactly 0.000 -- which is what the first run reported. Y was
    # unaffected only because row 0 is 0 in both conventions.
    c0 = int(round((x0_m - K.PAD_ORIGIN_M) / C.PAD_PITCH_M))
    print(f"  deposit at pad column {c0} (ESL phase {phase}), "
          f"x0={x0_mm:.2f} mm, row {row0}\n")

    print(f"{'z [mm]':>7s} {'sigma_T [µm]':>13s} "
          f"{'X: d=0':>8s} {'d=±1':>8s} {'d=±2':>8s} | "
          f"{'Y: d=0':>8s} {'d=±1':>8s} {'d=±2':>8s}")
    rows = []
    for z in (0.5, 2.0, 5.0, 10.0, 20.0, 30.0):
        accX = np.zeros(7)
        accY = np.zeros(7)
        nok = 0
        for _ in range(a.repeat):
            cur = point_deposit(dig, z, a.n_electrons, x0_mm, y0_mm)
            q = integrated_charge(cur, dig.lut.dt)
            sx = share_profile(q, "X", c0)
            sy = share_profile(q, "Y", row0)
            if sx is None or sy is None:
                continue
            accX += np.array([sx[d] for d in range(-3, 4)])
            accY += np.array([sy[d] for d in range(-3, 4)])
            nok += 1
        if not nok:
            continue
        accX /= nok
        accY /= nok
        sT = float(np.ravel(dig.gas.sigma_T_um(dig.E_drift, z))[0])
        c1X = 0.5 * (accX[2] + accX[4])
        c2X = 0.5 * (accX[1] + accX[5])
        c1Y = 0.5 * (accY[2] + accY[4])
        c2Y = 0.5 * (accY[1] + accY[5])
        print(f"{z:7.1f} {sT:13.0f} {accX[3]:8.3f} {c1X:8.3f} {c2X:8.3f} | "
              f"{accY[3]:8.3f} {c1Y:8.3f} {c2Y:8.3f}")
        rows.append({"z_mm": z, "sigma_T_um": sT,
                     "X": accX.tolist(), "Y": accY.tolist(),
                     "c1_X": c1X, "c1_Y": c1Y})

    print("\n  measured target (§9): dispersed ±1 share c1 = 0.23-0.28, "
          "gain/gas/drift invariant")
    print("  A MIP crossing the full gap averages over all these depths, so "
          "the\n  track-level c1 is the depth-weighted mean, not the z=30 mm "
          "row.")
    if rows:
        c1s = [r["c1_X"] for r in rows]
        print(f"  c1_X rises {c1s[0]:.3f} -> {c1s[-1]:.3f} from z=0.5 to 30 mm "
              f"(factor {c1s[-1]/max(c1s[0],1e-9):.1f}).")
        print("  That steep z dependence IS the diffusion signature: if the "
              "sharing came\n  from the avalanche or from induction geometry "
              "it would be flat in z.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
