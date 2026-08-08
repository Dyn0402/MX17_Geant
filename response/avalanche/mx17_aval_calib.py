#!/usr/bin/env python3
"""
mx17_aval_calib.py — one point of the S3 avalanche-calibration campaign.

Stage S3 of design/RESPONSE_SIM_PLAN.md. One invocation = one
(gas table, mesh voltage, seed slice) point, so that the whole grid is a
trivially parallel HTCondor submission (mx17_aval.sub).

Physics model, first pass
-------------------------
Uniform field in the amplification gap, as §5 of the plan explicitly permits
("or uniform-field fallback for first pass"). The gap is 150 µm (CONFIRMED,
plan §1 row 2), anode (ESL surface) at z = 0, mesh at z = +gap. The electric
field is E = (0, 0, V/gap), which pushes electrons toward -z, i.e. down onto
the resistive layer.

The realistic woven-mesh field from S2 (T6) replaces ComponentConstant later;
everything downstream of this script reads only the JSON, so that swap does not
propagate. What the uniform-field pass CANNOT give — and does not pretend to —
is (a) electron transparency ε(E_d/E_a), (b) the funnelling map, and (c) the
fraction of ions that escape back through the mesh into the drift gap. Those
are S2 observables; this file records them as null with a `uniform_field`
provenance flag so no downstream consumer can silently mistake a placeholder
for a measurement.

Weighting field
---------------
For the parallel-plate gap the readout weighting potential is ψ(z) = 1 - z/gap,
so the weighting field is uniform, E_w = (0, 0, +1/gap) [cm^-1]. That is exact
here, and it is the right thing for extracting the *shape* of the ion tail:
the resistive-layer and strip-segmentation structure of the true weighting
potential is S1's job (response/solver), and it multiplies this shape rather
than replacing it.

Extracted per point (written to JSON)
-------------------------------------
  gains[]           per-seed-electron avalanche size (raw; the Polya fit is done
                    by the collector over the merged seed slices)
  r_end[]           radial distance of each final electron from the seed axis
                    [µm] -> transverse avalanche spread σ0 at the ESL
  t_end[]           arrival time of each final electron [ns]
  z_ion[]           z at which each ionisation happened [µm] -> α(z) profile
  i_elec[], i_ion[] induced-current time profiles on the readout electrode
                    [fC/ns], sampled on a uniform grid; i_ion is the
                    normalisable ion tail shape the plan calls i_ion(t)
  ion_subsample     ions are drifted with AvalancheMC, which is slow; only this
                    many are drifted and i_ion is scaled up accordingly

Usage
-----
    python3 mx17_aval_calib.py --gas-file Ar_iC4H10_95_5_Saclay_160m.gas \\
        --voltage 490 --nev 500 --seed 1 --out aval_ArIso955_490V_s1.json
"""

from __future__ import annotations

import argparse
import ctypes
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

# ── Fixed geometry (plan §1; keep in sync with shared/MX17ModuleGeometry.hh) ──
GAP_UM_DEFAULT = 150.0      # amplification gap, CONFIRMED
TEMP_K_DEFAULT = 293.15


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gas-file", required=True,
                   help="Magboltz .gas table (amplification field range)")
    p.add_argument("--voltage", type=float, required=True,
                   help="Mesh voltage [V] across the amplification gap")
    p.add_argument("--gap-um", type=float, default=GAP_UM_DEFAULT)
    p.add_argument("--nev", type=int, default=500,
                   help="Seed electrons in this slice")
    p.add_argument("--seed", type=int, default=1,
                   help="RNG seed; also identifies the slice")
    p.add_argument("--penning", default="auto", choices=["auto", "manual", "off"])
    p.add_argument("--penning-rp", type=float, default=0.40)
    p.add_argument("--penning-gas", default="ar")
    p.add_argument("--temp-k", type=float, default=TEMP_K_DEFAULT)
    p.add_argument("--ion-subsample", type=int, default=150,
                   help="Ions drifted per avalanche (AvalancheMC is slow); the "
                        "ion signal is scaled by n_ion/n_drifted")
    p.add_argument("--tmax-ns", type=float, default=400.0,
                   help="Signal time window; must contain the full ion transit")
    p.add_argument("--nbins", type=int, default=2000)
    p.add_argument("--out", required=True)
    p.add_argument("--max-avalanche", type=int, default=0,
                   help="Cap avalanche size (0 = uncapped). Only for smoke tests; "
                        "a cap biases the gain distribution and is recorded.")
    p.add_argument("--field-map", default=None,
                   help="Path to a meshfield_<tag>.txt map from "
                        "solve_fieldmap.py (T6, FIELD_MAP_RUNBOOK.md). When "
                        "given, this replaces the uniform-field "
                        "ComponentConstant with the real woven-mesh field via "
                        "ComponentGrid; seed electrons start above the mesh "
                        "and transparency/funnelling become emergent instead "
                        "of assumed nulls. --voltage/--gap-um must match the "
                        "map's own geometry (see the map's .json).")
    return p.parse_args()


def git_hash(path):
    try:
        return subprocess.check_output(
            ["git", "-C", os.path.dirname(os.path.abspath(path)),
             "rev-parse", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def main():
    args = parse_args()

    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    if ROOT.gSystem.Load("libGarfield") < 0:
        sys.exit("[aval] cannot load libGarfield")
    G = ROOT.Garfield

    gap_cm = args.gap_um * 1e-4
    e_field = args.voltage / gap_cm          # V/cm

    # ── Gas ──────────────────────────────────────────────────────────────────
    gas = G.MediumMagboltz()
    if not gas.LoadGasFile(args.gas_file):
        sys.exit(f"[aval] failed to load gas file {args.gas_file}")
    gas.SetTemperature(args.temp_k)

    penning_applied = None
    if args.penning == "auto":
        penning_applied = bool(gas.EnablePenningTransfer())
        # A False here means Garfield has no built-in table for this mixture and
        # is silently running at rP = 0 — the exact trap documented at length in
        # nTof_x17 garfield_sim/mm_config.py. Refuse rather than produce a
        # quietly-wrong gain.
        if not penning_applied:
            sys.exit("[aval] EnablePenningTransfer() returned False for this "
                     "mixture: no built-in rP. Re-run with --penning manual "
                     "--penning-rp <value>, and record the choice.")
    elif args.penning == "manual":
        penning_applied = bool(
            gas.EnablePenningTransfer(args.penning_rp, 0., args.penning_gas))

    # Ion mobility: needed for the ion tail. Ar+ in Ar is the dominant ion in
    # Ar/iso after charge transfer to isobutane; the iC4H10 file is the other
    # bracket. Using Ar+ is the conservative (faster ion) choice and the plan
    # already flags ion mobility as the softest parameter (±30 % systematic).
    data_dir = os.path.join(os.environ.get("GARFIELD_INSTALL", ""),
                            "share", "Garfield", "Data")
    ion_file = os.path.join(data_dir, "IonMobility_Ar+_Ar.txt")
    if not gas.LoadIonMobility(ion_file):
        sys.exit(f"[aval] failed to load ion mobility {ion_file}")

    # ── Field: either the uniform-field parallel plate, or the T6 woven-mesh
    #    map loaded through ComponentGrid (FIELD_MAP_RUNBOOK.md) ─────────────
    field_map_grid = None  # keep a reference alive for the sensor's lifetime
    if args.field_map:
        field_map_grid = G.ComponentGrid()
        if not field_map_grid.LoadElectricField(args.field_map, "xyz",
                                                 True, True):
            sys.exit(f"[aval] failed to load field map {args.field_map}")
        field_map_grid.SetMedium(gas)
        # Capture the map's own z-extent BEFORE enabling periodicity — the
        # bounding box becomes laterally infinite afterwards (gates_check.C).
        # This is the SOURCE OF TRUTH for the anode/ESL plane: solve_fieldmap
        # exports zlo = Z_ANODE + 0.5 um (a half-step inset from the true
        # ESL), so recomputing Z_ANODE analytically here would put the
        # "reached anode" test 0.5 um past the last real grid point and the
        # electron would leave the map (medium -> nullptr, absorbed at the
        # boundary) just short of it — every avalanche would read gain = 0.
        # (Caught exactly that way on the first field-map smoke run.)
        gxlo, gylo, gzlo = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
        gxhi, gyhi, gzhi = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
        field_map_grid.GetBoundingBox(gxlo, gylo, gzlo, gxhi, gyhi, gzhi)
        field_map_grid.EnablePeriodicityX()
        field_map_grid.EnablePeriodicityY()

        anode_z_cm = gzlo.value
        mesh_under_z_cm = anode_z_cm + gap_cm  # amp gap = pillar height
        pitch_cm_area = 67.0e-4

        # Readout weighting field/potential: same parallel-plate ansatz as
        # the uniform-field pass (S1's real weighting potential is a separate,
        # later job — plan §7 step 6), ψ = 1 at the ESL ramping to ψ = 0 at
        # the mesh underside. Area is confined to the amp gap ONLY (not the
        # drift bulk above): the mesh screens the readout from induction
        # above it, same physical argument as the uniform-field pass, which
        # had no drift region to begin with. Confining the area also gives
        # WeightingField/WeightingPotential a real GetMedium() to fall back
        # on (see the no-area note below) without ever throwing.
        wcmp = G.ComponentConstant()
        wcmp.SetWeightingField(0., 0., 1.0 / gap_cm, "readout")
        wcmp.SetWeightingPotential(0., 0., anode_z_cm, 1.0)
        # wcmp is used ONLY for its weighting field/potential (added via
        # AddElectrode below, never AddComponent), but Garfield still probes
        # every electrode component's GetMedium() along the trajectory; with
        # no area set that falls through to the base Component::GetMedium(),
        # which throws ("geometry is not set") rather than returning nullptr.
        wcmp.SetArea(-10 * pitch_cm_area, -10 * pitch_cm_area, anode_z_cm,
                     10 * pitch_cm_area, 10 * pitch_cm_area, mesh_under_z_cm)

        sensor = G.Sensor()
        sensor.AddComponent(field_map_grid)
        sensor.AddElectrode(wcmp, "readout")
        sensor.SetArea(-10 * pitch_cm_area, -10 * pitch_cm_area, gzlo.value,
                       10 * pitch_cm_area, 10 * pitch_cm_area, gzhi.value)

        # Seed electrons above the mesh, spread over one unit cell, so
        # transparency and funnelling are measured, not assumed (runbook
        # "What S3 does next"). z chosen deep enough in the drift bulk that
        # the mesh's near-field perturbation is negligible (matches
        # gates_check.C's G7 transparency probe).
        seed_z0_cm = 180.0e-4
        pitch_cm = 67.0e-4
        seed_rng = np.random.default_rng(args.seed)

        def sample_xy():
            return (seed_rng.uniform(-pitch_cm / 2, pitch_cm / 2),
                    seed_rng.uniform(-pitch_cm / 2, pitch_cm / 2))
    else:
        # Lateral half-size: generous, so nothing leaves the area sideways.
        half = 0.05   # cm = 500 µm
        cmp = G.ComponentConstant()
        cmp.SetArea(-half, -half, 0.0, half, half, gap_cm)
        cmp.SetMedium(gas)
        cmp.SetElectricField(0., 0., e_field)
        # Readout weighting field, ψ = 1 - z/gap  ->  E_w = -∇ψ = +ẑ/gap.
        cmp.SetWeightingField(0., 0., 1.0 / gap_cm, "readout")
        # ...AND the weighting POTENTIAL, which is not optional. Both
        # AvalancheMicroscopic and AvalancheMC default to m_useWeightingPotential
        # = true, i.e. they compute the induced current from ψ rather than from
        # E_w. ComponentConstant has no ψ until this call, so without it
        # WeightingPotential() returns 0 everywhere and EVERY signal comes out
        # identically zero — silently, with the avalanche and the ion drift both
        # running normally. That is what emptied i_elec/i_ion in the first
        # campaign (56/56 slices). Anchoring ψ = 1 at the readout plane z = 0 with
        # a constant weighting field gives ψ(z) = 1 - z/gap exactly.
        cmp.SetWeightingPotential(0., 0., 0., 1.0)

        sensor = G.Sensor()
        sensor.AddComponent(cmp)
        sensor.AddElectrode(cmp, "readout")

        anode_z_cm = 0.0
        seed_z0_cm = gap_cm - 1.0e-4

        def sample_xy():
            return (0.0, 0.0)

    dt = args.tmax_ns / args.nbins
    sensor.SetTimeWindow(0., dt, args.nbins)

    aval = G.AvalancheMicroscopic()
    aval.SetSensor(sensor)
    aval.EnableSignalCalculation()
    if args.max_avalanche > 0:
        aval.EnableAvalancheSizeLimit(args.max_avalanche)

    drift_ion = G.AvalancheMC()
    drift_ion.SetSensor(sensor)
    drift_ion.SetDistanceSteps(2.e-5)      # 0.2 µm steps
    drift_ion.EnableSignalCalculation()

    ROOT.gRandom.SetSeed(args.seed)
    try:
        G.randomEngine.Seed(args.seed)
    except Exception:
        pass

    # ── Run ──────────────────────────────────────────────────────────────────
    gains, r_end, t_end, z_ion = [], [], [], []
    n_ion_total, n_ion_drifted = 0, 0
    t_start = time.time()

    # The signal buffers are deliberately NOT cleared inside the loop: every
    # avalanche starts at t = 0, so Garfield's running sum over all nev events
    # is exactly nev x (mean single-avalanche signal). Divided by nev below.
    sensor.ClearSignal()

    for iev in range(args.nev):
        x0, y0 = sample_xy()
        aval.AvalancheElectron(x0, y0, seed_z0_cm, 0., 0.1, 0., 0., 0.)

        ne = aval.GetNumberOfElectronEndpoints()
        n_at_anode = 0
        ion_seeds = []
        for i in range(ne):
            xs, ys, zs, ts, es = (ctypes.c_double() for _ in range(5))
            xe, ye, ze, te, ee = (ctypes.c_double() for _ in range(5))
            status = ctypes.c_int()
            aval.GetElectronEndpoint(i, xs, ys, zs, ts, es,
                                     xe, ye, ze, te, ee, status)
            # Where this electron was born = where an ionisation happened
            # (the i = 0 entry is the seed electron itself, excluded). Height
            # is recorded ABOVE THE ANODE, matching the uniform-field
            # convention that ions.py/kernel_lut.py's alpha_z histogram
            # expects (z=0 at anode).
            if i > 0:
                z_ion.append((zs.value - anode_z_cm) * 1e4)
                ion_seeds.append((xs.value, ys.value, zs.value, ts.value))
            # Reached the anode?
            if abs(ze.value - anode_z_cm) < 1.e-5:   # within 0.1 µm
                n_at_anode += 1
                r_end.append(math.hypot(xe.value - x0, ye.value - y0) * 1e4)
                t_end.append(te.value)
        gains.append(n_at_anode)

        # ── Ion tail ─────────────────────────────────────────────────────────
        # One ion is left behind at every ionisation point. Drift a random
        # subsample and scale the resulting signal.
        n_ion_total += len(ion_seeds)
        if ion_seeds:
            k = min(args.ion_subsample, len(ion_seeds))
            idx = np.random.default_rng(args.seed * 100003 + iev).choice(
                len(ion_seeds), size=k, replace=False)
            for j in idx:
                x, y, z, t = ion_seeds[j]
                drift_ion.DriftIon(x, y, z, t)
            n_ion_drifted += k

        if iev == 0 or (iev + 1) % 50 == 0:
            el = time.time() - t_start
            print(f"[aval] {iev+1}/{args.nev} mean_gain="
                  f"{np.mean(gains):.1f} elapsed={el:.0f}s "
                  f"eta={el/(iev+1)*(args.nev-iev-1):.0f}s", flush=True)

    # ── Accumulated signal, per avalanche ────────────────────────────────────
    # Summed over all nev events (never cleared in the loop) -> divide by nev.
    # The ion component additionally carries the subsampling scale, since only
    # ion_scale^-1 of the ions produced were actually drifted.
    i_elec = np.array([sensor.GetElectronSignal("readout", b)
                       for b in range(args.nbins)]) / args.nev
    ion_scale = (n_ion_total / n_ion_drifted) if n_ion_drifted else 0.0
    i_ion = np.array([sensor.GetIonSignal("readout", b)
                      for b in range(args.nbins)]) * ion_scale / args.nev

    gains = np.asarray(gains, dtype=float)
    out = {
        "schema": "mx17_aval_calib/1",
        "provenance": {
            "field_model": ("meshfield:" + os.path.basename(args.field_map))
                           if args.field_map else "uniform_field",
            "garfield_pin": os.environ.get("MX17_GARFIELD_PIN", "unknown"),
            "script_git": git_hash(__file__),
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "host": os.uname().nodename,
            "runtime_s": time.time() - t_start,
        },
        "config": {
            "gas_file": os.path.basename(args.gas_file),
            "voltage_V": args.voltage,
            "gap_um": args.gap_um,
            "e_field_Vcm": e_field,
            "temp_K": args.temp_k,
            "nev": args.nev,
            "seed": args.seed,
            "penning_mode": args.penning,
            "penning_applied": penning_applied,
            "penning_rp": args.penning_rp if args.penning == "manual" else None,
            "ion_mobility_file": os.path.basename(ion_file),
            "max_avalanche": args.max_avalanche,
            "field_map_file": os.path.basename(args.field_map)
                              if args.field_map else None,
        },
        "results": {
            "gains": gains.tolist(),
            "gain_mean": float(gains.mean()),
            "gain_std": float(gains.std(ddof=1)) if len(gains) > 1 else 0.0,
            "survival": float((gains > 0).mean()),
            "r_end_um": r_end,
            "t_end_ns": t_end,
            "z_ion_um": z_ion,
            "signal_dt_ns": dt,
            "i_elec": i_elec.tolist(),
            "i_ion": i_ion.tolist(),
            "ion_scale": ion_scale,
            "n_ion_total": n_ion_total,
            "n_ion_drifted": n_ion_drifted,
        },
        # S2 observables this script does not itself extract. Under
        # uniform_field they are impossible in principle (explicit null).
        # Under a field map they are physically present in this run (seed
        # electrons above the mesh, so "survival" above already reflects
        # transparency losses, and r_end_um already reflects funnelling) but
        # not reduced to a per-point epsilon/funnel-map here — that is the
        # dedicated transparency-curve / funnel-map deliverable, so the
        # explicit-null contract is kept for both modes to avoid a consumer
        # mistaking an un-reduced run for that measurement.
        "s2_pending": {
            "transparency_eps": None,
            "funnel_map": None,
            "ion_fraction_to_mesh": None,
        },
    }
    with open(args.out, "w") as f:
        json.dump(out, f)
    print(f"[aval] wrote {args.out}  gain={gains.mean():.1f}±"
          f"{gains.std(ddof=1) if len(gains)>1 else 0:.1f} "
          f"runtime={time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
