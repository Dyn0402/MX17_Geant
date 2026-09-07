#!/usr/bin/env python3
"""
funnel_ion_endpoints.py — S2 deliverables (b) funneling map and (d) ion
endpoints, per FIELD_MAP_RUNBOOK.md's "Cheap follow-ups", on the accepted
meshfield_production.txt.

(b) FUNNELING MAP. Deterministic (no-diffusion) electron field lines from a
grid of (x0, y0) entry points high in the drift bulk (z = 180 um, the same
convention gates_check.C's G7 probe uses — the field there is uniform, so
the exact entry height above the mesh's near-field zone does not matter) down
to the ESL, via DriftLineRKF (macroscopic drift-velocity integration, no
random walk). This isolates the mesh's geometric field-line convergence from
thermal diffusion, which the digitizer's own Gaussian kernel already accounts
for separately (mixing the two would double-count the spread). Records the
entry -> ESL-arrival offset for every grid point, and whether the field line
reaches the ESL at all or terminates on a wire — the noise-free geometric
transparency boundary, as opposed to G7's statistical average over random
entry points.

(d) ION ENDPOINTS. AvalancheMC ion drift lines from birth points sampled
through the amp gap — z from the measured exponential alpha(z) profile
(ions.py: ~14.8 um decay length above the ESL), (x0, y0) uniform over the
unit cell — classified by where they end up: absorbed on the mesh wire (the
normal Micromegas ion-tail endpoint already assumed everywhere else in this
project, e.g. ions.py's psi = 1 - z/g model) vs escaped past the mesh into
the drift bulk (ion backflow, not currently modelled anywhere downstream).

    source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh
    python3 funnel_ion_endpoints.py --field-map meshfield_production.txt
"""
from __future__ import annotations

import argparse
import ctypes
import json
import os
import sys

import numpy as np

WIRE_R_UM = 9.5
Z_OFF_UM = 9.0  # wire-axis |z| offset (1 um kiss overlap); matches gates_check.C
PITCH_CM = 67.0e-4
GAP_UM_DEFAULT = 150.0
ION_Z_DECAY_UM = 14.8  # measured alpha(z) decay length, aval_calib_v2 (ions.py)
# Margin (um, above and below gap_um) treated as "still at the mesh": an ion
# whose drift line stops (StatusLeftDriftMedium) with its height above the
# anode inside [gap_um - MARGIN, gap_um + MARGIN] is classed as absorbed on
# a wire, not merely "escaped". A post-hoc GetMedium() re-query at the
# reported endpoint is NOT reliable for this: RKF sub-steps can trip the
# inactive-region check mid-step while the reported (bisected) endpoint
# itself lands a hair on the valid side, so the *coordinate* reads as valid
# gas even though the drift line genuinely terminated at the wire — caught
# by a debug run where 15/15 ions all reported StatusLeftDriftMedium at
# z = gap_um +0/+23 um, yet GetMedium() at that exact point returned valid
# gas every time.
MESH_BAND_MARGIN_UM = 25.0


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--field-map", default="meshfield_production.txt")
    p.add_argument("--gas-file", default=(
        "/home/dylan/PycharmProjects/nTof_x17/garfield_sim/gas_tables/"
        "Ar_iC4H10_95_5_Saclay_160m.gas"))
    p.add_argument("--gap-um", type=float, default=GAP_UM_DEFAULT)
    p.add_argument("--n-funnel-side", type=int, default=21,
                   help="funnel map is an N x N grid over one unit cell")
    p.add_argument("--n-ions", type=int, default=3000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", default="funnel_ion_endpoints.json")
    return p.parse_args()


def main():
    args = parse_args()
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    if ROOT.gSystem.Load("libGarfield") < 0:
        sys.exit("[funnel] cannot load libGarfield")
    G = ROOT.Garfield

    gas = G.MediumMagboltz()
    if not gas.LoadGasFile(args.gas_file):
        sys.exit(f"[funnel] failed to load gas file {args.gas_file}")
    ion_file = os.path.join(os.environ.get("GARFIELD_INSTALL", ""),
                            "share", "Garfield", "Data",
                            "IonMobility_Ar+_Ar.txt")
    if not gas.LoadIonMobility(ion_file):
        sys.exit(f"[funnel] failed to load ion mobility {ion_file}")

    grid = G.ComponentGrid()
    if not grid.LoadElectricField(args.field_map, "xyz", True, True):
        sys.exit(f"[funnel] failed to load field map {args.field_map}")
    grid.SetMedium(gas)
    gxlo, gylo, gzlo = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gxhi, gyhi, gzhi = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    grid.GetBoundingBox(gxlo, gylo, gzlo, gxhi, gyhi, gzhi)
    grid.EnablePeriodicityX()
    grid.EnablePeriodicityY()

    anode_z_cm = gzlo.value  # ESL, same convention as mx17_aval_calib.py
    gap_cm = args.gap_um * 1e-4

    sensor = G.Sensor()
    sensor.AddComponent(grid)
    sensor.SetArea(-10 * PITCH_CM, -10 * PITCH_CM, gzlo.value,
                   10 * PITCH_CM, 10 * PITCH_CM, gzhi.value)

    # ── (b) funneling map ────────────────────────────────────────────────────
    drift = G.DriftLineRKF()
    drift.SetSensor(sensor)
    drift.EnableSignalCalculation(False)

    n = args.n_funnel_side
    xs = np.linspace(-PITCH_CM / 2, PITCH_CM / 2, n, endpoint=False) \
        + PITCH_CM / (2 * n)
    funnel_rows = []
    n_reached = 0
    for x0 in xs:
        for y0 in xs:
            ok = drift.DriftElectron(x0, y0, 180.0e-4, 0.)
            xe, ye, ze, te = (ctypes.c_double() for _ in range(4))
            status = ctypes.c_int()
            drift.GetEndPoint(xe, ye, ze, te, status)
            reached = bool(ok) and abs(ze.value - anode_z_cm) < 2.e-5
            if reached:
                n_reached += 1
            funnel_rows.append({
                "x0_um": x0 * 1e4, "y0_um": y0 * 1e4,
                "xe_um": xe.value * 1e4, "ye_um": ye.value * 1e4,
                "dx_um": (xe.value - x0) * 1e4, "dy_um": (ye.value - y0) * 1e4,
                "reached_esl": reached, "status": int(status.value),
            })
    geometric_eps = n_reached / len(funnel_rows)
    print(f"[funnel] geometric (no-diffusion) transparency: {geometric_eps:.3f} "
          f"({n_reached}/{len(funnel_rows)} field lines reach the ESL)")

    # ── (d) ion endpoints ────────────────────────────────────────────────────
    rng = np.random.default_rng(args.seed)
    n_mesh, n_drift, n_other = 0, 0, 0
    ion_rows = []
    for _ in range(args.n_ions):
        x0 = rng.uniform(-PITCH_CM / 2, PITCH_CM / 2)
        y0 = rng.uniform(-PITCH_CM / 2, PITCH_CM / 2)
        z_above_anode = min(rng.exponential(ION_Z_DECAY_UM * 1e-4),
                            gap_cm - 1e-7)  # stay inside the gap
        z0 = anode_z_cm + z_above_anode

        drift_ion = G.AvalancheMC()
        drift_ion.SetSensor(sensor)
        drift_ion.SetDistanceSteps(2.e-5)
        drift_ion.EnableSignalCalculation(False)
        drift_ion.DriftIon(x0, y0, z0, 0.)
        if drift_ion.GetNumberOfIonEndpoints() == 0:
            n_other += 1
            continue
        xs_, ys_, zs_, ts_ = (ctypes.c_double() for _ in range(4))
        xe, ye, ze, te = (ctypes.c_double() for _ in range(4))
        status = ctypes.c_int()
        drift_ion.GetIonEndpoint(0, xs_, ys_, zs_, ts_, xe, ye, ze, te, status)
        ze_above_anode_um = (ze.value - anode_z_cm) * 1e4
        if ze_above_anode_um > args.gap_um + MESH_BAND_MARGIN_UM:
            n_drift += 1
            outcome = "escaped_to_drift"
        elif ze_above_anode_um > args.gap_um - MESH_BAND_MARGIN_UM:
            n_mesh += 1
            outcome = "absorbed_on_mesh"
        else:
            n_other += 1
            outcome = "other"
        ion_rows.append({"z0_above_anode_um": z_above_anode * 1e4,
                         "ze_above_anode_um": ze_above_anode_um,
                         "status": int(status.value), "outcome": outcome})

    n_ok = n_mesh + n_drift + n_other
    print(f"[funnel] ion endpoints (n={n_ok}): "
          f"absorbed_on_mesh={n_mesh / n_ok:.3f}, "
          f"escaped_to_drift={n_drift / n_ok:.3f}, other={n_other / n_ok:.3f}")

    out = {
        "schema": "funnel_ion_endpoints/1",
        "field_map": args.field_map,
        "config": vars(args),
        "funnel_map": {
            "geometric_transparency": geometric_eps,
            "n_points": len(funnel_rows),
            "rows": funnel_rows,
        },
        "ion_endpoints": {
            "n_ions": n_ok,
            "frac_absorbed_on_mesh": n_mesh / n_ok,
            "frac_escaped_to_drift": n_drift / n_ok,
            "frac_other": n_other / n_ok,
            "rows": ion_rows,
        },
    }
    with open(args.out, "w") as f:
        json.dump(out, f)
    print(f"[funnel] wrote {args.out}")


if __name__ == "__main__":
    main()
