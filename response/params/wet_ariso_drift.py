#!/usr/bin/env python3
"""
wet_ariso_drift.py — Magboltz for det3's own gas, wet, not dry.

WHY THIS EXISTS. The P1 decision (design/RESPONSE_SIM_PLAN.md §0a, 2026-08-07)
retargets T14 from run_71/SPS onto det3 cosmic-bench data, Ar/iC4H10 95/5 at
1000 V drift / 490 V mesh over a 30 mm nominal gap (333 V/cm) -- exactly what
Stage B already simulates, dry. But the same "the gas probably isn't dry"
suspicion that motivated wet_cf4_drift.py for run_71 applies here too:
~/x17/cosmic_bench/Analysis/fleet_gas_survey.csv records h2o_pct=1.46 for the
EXACT long_run_resist_490V_drift_1000V sub-run behind T12's pedestals/noise,
and v_drift for det3 is quoted three incompatible ways (39.1 dry Magboltz /
36.6 measured, plan §1 / 28.1+/-0.7 micro-TPC, a different ~19mm-gap point per
nTof_x17 HANDOFF_det3_vdrift_and_kernels.md). This script runs the same
water-fraction bracket wet_cf4_drift.py ran for run_71, but for Ar/iC4H10
95/5 at det3's own field point, so T_drift/halo-width predictions for det3
rest on the same footing as they now do for the (shelved) run_71 target.

GAP: per user instruction (2026-08-07), use the 30 mm nominal gap for this
run -- the effective-gap dispute (19 / 23 / 27.9 / 30 mm across three
nTof_x17 docs) is a separate P1 sub-item, not resolved here. sigma_T below is
reported at 30 mm; rescale later if the gap question resolves to something
else (sigma_T ~ sqrt(z), so a 27.9 mm gap would be sqrt(27.9/30) = 0.964x
this table's 30 mm number).

THE WATER FRACTION IS SCANNED, NOT ASSUMED, same reasoning as wet_cf4_drift.py:
0.5/1.0/2.0 % brackets a wide range, and 1.46 % is included explicitly because
it is not a guess -- it is the recorded gas-survey value for the very run this
whole comparison targets. Dry is kept as the reference (should reproduce the
existing drift_velocity_Ar_iC4H10_95_5_Saclay.json table).

Run on lxplus (Garfield pinned at master 927e5c21, self-built vs LCG_109):
    source setup_garfield.sh
    nohup python3 wet_ariso_drift.py > wet_ariso.log 2>&1 &
"""

import ctypes
import json
import multiprocessing as mp
import os
import time

# 30 mm gap (per user instruction, 2026-08-07 -- see GAP note above): the det3
# T14 point is 1000 V / 30 mm nominal = 333 V/cm. The grid runs well past it
# so nearby lever points (the 900 V sub-run, etc.) are covered too.
E_MIN, E_MAX, N_GRID = 40.0, 500.0, 14
NCOLL = 5                      # matches drift_velocity_Ar_iC4H10_95_5_Saclay.json
PRESSURE_TORR = 745.83         # Saclay pressure, matching mesh_transparency.C / the S3 campaign
TEMP_K = 293.15
NWORKERS = 5
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "wet_ariso_drift.json")

# Water displaces argon; isobutane stays at the nominal 5 %.
MIXTURES = {
    "dry":     [("ar", 95.0), ("ic4h10", 5.0)],
    "H2O_0.5": [("ar", 94.5), ("ic4h10", 5.0), ("h2o", 0.5)],
    "H2O_1.0": [("ar", 94.0), ("ic4h10", 5.0), ("h2o", 1.0)],
    "H2O_1.46":[("ar", 93.54),("ic4h10", 5.0), ("h2o", 1.46)],
    "H2O_2.0": [("ar", 93.0), ("ic4h10", 5.0), ("h2o", 2.0)],
}

# The numbers this run has to be read against. Unlike run_71 there is no
# single tight measured band -- three numbers are in circulation for det3 and
# they disagree with each other, which is exactly what §0a P1 flags. Report
# against all three rather than picking one.
VALIDATION = {
    "E_Vcm": 333.0,
    "dry_reference_um_per_ns": 39.1,
    "measured_um_per_ns_candidates": {
        "plan_section1": 36.6,
        "micro_tpc_det3_2026-07-01": {
            "value": 28.1, "sigma": 0.7,
            "caveat": "measured over det3's ~19 mm gap per nTof_x17 "
                      "mx_june_cosmic_qa/HANDOFF_det3_vdrift_and_kernels.md, "
                      "not the 30 mm nominal this table assumes -- not "
                      "directly comparable until the gap question resolves",
        },
    },
    "h2o_survey_pct": 1.46,
    "h2o_survey_source": ("~/x17/cosmic_bench/Analysis/fleet_gas_survey.csv, "
                           "run long_run_resist_490V_drift_1000V, mx17_3"),
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
    print(f"{len(MIXTURES)} mixtures, Ar/iC4H10 95/5 base @ "
          f"{PRESSURE_TORR} Torr, ncoll={NCOLL}", flush=True)
    with mp.get_context("spawn").Pool(min(NWORKERS, len(MIXTURES))) as pool:
        results = dict(pool.map(worker, list(MIXTURES.items())))

    with open(OUT, "w") as f:
        json.dump(dict(schema="wet_ariso_drift/1",
                       gas_base="Ar/iC4H10 95/5",
                       pressure_torr=PRESSURE_TORR, temp_K=TEMP_K,
                       ncoll=NCOLL, drift_gap_mm_assumed=30.0,
                       validation=VALIDATION,
                       mixtures=results), f, indent=1)
    print(f"Written {OUT}")

    print(f"\n  v_drift at {VALIDATION['E_Vcm']:.0f} V/cm "
          f"[candidates: dry ref {VALIDATION['dry_reference_um_per_ns']}, "
          f"plan §1 {VALIDATION['measured_um_per_ns_candidates']['plan_section1']}, "
          f"micro-TPC {VALIDATION['measured_um_per_ns_candidates']['micro_tpc_det3_2026-07-01']['value']}"
          f"+/-{VALIDATION['measured_um_per_ns_candidates']['micro_tpc_det3_2026-07-01']['sigma']} "
          f"(different gap, see caveat)]")
    for name, rows in results.items():
        r = min(rows, key=lambda q: abs(q["E_Vcm"] - VALIDATION["E_Vcm"]))
        # sigma_T over the 30 mm assumed gap -- the quantity Stage B consumes.
        sig_um = r["dT_sqrtcm"] * (3.0 ** 0.5) * 1e4
        print(f"    {name:9s} v={r['v_um_per_ns']:7.2f} µm/ns   "
              f"eta={r['eta_per_cm']:.4g} /cm   "
              f"dT={r['dT_sqrtcm']:.4f} /sqrt(cm)   "
              f"sigma_T(30 mm)={sig_um:6.0f} µm")


if __name__ == "__main__":
    main()
