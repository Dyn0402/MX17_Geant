#!/usr/bin/env python3
"""
wet_cf4_drift.py — Magboltz for the gas the §9 targets were actually taken in.

WHY THIS EXISTS. Every Stage B comparison so far has run Ar/iC4H10 95/5, dry,
at 333 V/cm. The measured area/peak/peak-time table in §9 comes from run_71,
which was Ar/CF4/iC4H10 88/10/2 carrying 1.3-1.7 % water, at 233 V/cm, drifting
at a MEASURED 13-15 um/ns against a dry-Magboltz prediction of 74.7. So the
model has been validated against a detector it is not a model of, and a ~3x
error in drift velocity is a ~3x error in arrival spread -- which is the broad,
late, low-peak halo the simulation is missing.

That is an input mismatch, not missing physics, and fixing it is not tuning:
it is simulating the operating point the data was taken at, which is what the
§9 firewall asks for.

WHAT IS NEEDED AND WHY IT IS NOT ALREADY AVAILABLE. nTof_x17 already has
drift_velocity_beamtest_cf4_wet{1,2}_CERN.json, which fixed the water fraction
by bracketing the measured v_drift between the 1 % and 2 % curves. Those files
carry v_drift ONLY. Stage B needs TRANSVERSE DIFFUSION -- dT sets the lateral
spread at the mesh and therefore the whole sharing budget -- so the table has
to be regenerated with ElectronDiffusion as well.

Structure follows nTof_x17's mm_drift_9010_contam_cern.py, which is the
established pattern for a contaminated-gas suite here; the mixtures and the
pressure are the run_71 ones.

THE WATER FRACTION IS SCANNED, NOT ASSUMED. 1.3-1.7 % is the data-side
bracket, so 1.0/1.5/2.0 % spans it and the dry point is kept as the reference
that shows how large the effect is. Picking a single value would hide the
sensitivity, and the plan carries drift as an open item precisely because it
has been quoted three incompatible ways.

Run on lxplus (Garfield pinned at master 927e5c21, self-built vs LCG_109):
    source setup_garfield.sh
    nohup python3 wet_cf4_drift.py > wet_cf4.log 2>&1 &
"""

import ctypes
import json
import multiprocessing as mp
import os
import time

# 30 mm gap: 233 V/cm is the run_71 operating point (700 V). The grid runs well
# past it so the drift-field lever points of the scan are covered too.
E_MIN, E_MAX, N_GRID = 40.0, 500.0, 14
NCOLL = 5                      # as the 9010 suite; enough for eta as well as v
PRESSURE_TORR = 720.8          # CERN EAR2, 450 m — matches run_71
TEMP_K = 293.15
NWORKERS = 4
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "wet_cf4_drift.json")

# Water displaces argon; CF4 and isobutane stay at the nominal 10/2.
MIXTURES = {
    "dry":     [("ar", 88.0), ("cf4", 10.0), ("ic4h10", 2.0)],
    "H2O_1.0": [("ar", 87.0), ("cf4", 10.0), ("ic4h10", 2.0), ("h2o", 1.0)],
    "H2O_1.5": [("ar", 86.5), ("cf4", 10.0), ("ic4h10", 2.0), ("h2o", 1.5)],
    "H2O_2.0": [("ar", 86.0), ("cf4", 10.0), ("ic4h10", 2.0), ("h2o", 2.0)],
}

# The data-side numbers this run has to reproduce before anything downstream
# is believable. Dry Magboltz already agrees with the first column, so the
# check that matters is the wet one landing on 13-15 um/ns at 233 V/cm.
VALIDATION = {
    "E_Vcm": 233.0,
    "dry_expected_um_per_ns": 74.7,
    "measured_um_per_ns": [13.0, 15.0],
    "source": "nTof_x17 sps_beam_test_26/analysis/RAW_RUN71_REANALYSIS_2026-08-04.md",
}


def worker(args):
    name, comps = args
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    assert ROOT.gSystem.Load("libGarfield") >= 0
    import numpy as np

    gas = ROOT.Garfield.MediumMagboltz()
    gas.SetComposition(*[x for pair in comps for x in pair])
    gas.SetTemperature(TEMP_K)
    gas.SetPressure(PRESSURE_TORR)
    gas.SetFieldGrid(E_MIN, E_MAX, N_GRID, True)
    t0 = time.time()
    gas.GenerateGasTable(NCOLL)

    rows = []
    for e in np.logspace(np.log10(E_MIN), np.log10(E_MAX), 60):
        vx, vy, vz = (ctypes.c_double(0.) for _ in range(3))
        gas.ElectronVelocity(0., 0., -e, 0., 0., 0., vx, vy, vz)
        eta = ctypes.c_double(0.)
        gas.ElectronAttachment(0., 0., -e, 0., 0., 0., eta)
        dl, dt = ctypes.c_double(0.), ctypes.c_double(0.)
        gas.ElectronDiffusion(0., 0., -e, 0., 0., 0., dl, dt)
        rows.append(dict(E_Vcm=float(e),
                         v_um_per_ns=float(vz.value * 1e4),
                         eta_per_cm=float(eta.value),
                         dL_sqrtcm=float(dl.value),
                         dT_sqrtcm=float(dt.value)))
    print(f"{name}: done in {(time.time()-t0)/60:.1f} min", flush=True)
    return name, rows


def main():
    print(f"{len(MIXTURES)} mixtures, Ar/CF4/iC4H10 88/10/2 base @ "
          f"{PRESSURE_TORR} Torr, ncoll={NCOLL}", flush=True)
    with mp.get_context("spawn").Pool(min(NWORKERS, len(MIXTURES))) as pool:
        results = dict(pool.map(worker, list(MIXTURES.items())))

    with open(OUT, "w") as f:
        json.dump(dict(schema="wet_cf4_drift/1",
                       gas_base="Ar/CF4/iC4H10 88/10/2",
                       pressure_torr=PRESSURE_TORR, temp_K=TEMP_K,
                       ncoll=NCOLL, validation=VALIDATION,
                       mixtures=results), f, indent=1)
    print(f"Written {OUT}")

    # Report the validation point immediately rather than making the next
    # reader dig for it: if the wet curves do not bracket 13-15 um/ns the
    # water hypothesis is not reproduced and nothing downstream should be run.
    print(f"\n  v_drift at {VALIDATION['E_Vcm']:.0f} V/cm "
          f"[measured {VALIDATION['measured_um_per_ns'][0]}-"
          f"{VALIDATION['measured_um_per_ns'][1]} µm/ns]")
    for name, rows in results.items():
        r = min(rows, key=lambda q: abs(q["E_Vcm"] - VALIDATION["E_Vcm"]))
        # sigma_T over the 30 mm gap, the quantity Stage B actually consumes.
        sig_um = r["dT_sqrtcm"] * (3.0 ** 0.5) * 1e4
        print(f"    {name:9s} v={r['v_um_per_ns']:7.2f} µm/ns   "
              f"dT={r['dT_sqrtcm']:.4f} /sqrt(cm)   "
              f"sigma_T(30 mm)={sig_um:6.0f} µm")


if __name__ == "__main__":
    main()
