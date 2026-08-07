#!/usr/bin/env python3
"""
run.py — Stage B over a real Geant4 ClusterTree (plan §7, acceptance for T9).

Turns the single-depth self-test numbers into the track-level observable §9
actually compares against: a MIP crossing the full 30 mm gap deposits clusters
at every depth, so what reaches the electronics is the depth-weighted mixture,
not any one row of the depth scan.

The charge budget is centred on the channel with the most charge IN THAT VIEW,
per event, which is what the data analysis does. It is not centred on the true
track position: with a checkerboard the pad under the track belongs to only one
view, so "the channel nearest the track" is ambiguous for the other one, and
using truth would bias the two views differently.

    python3 -m response.digitizer.run clusters.root --kernel <s1.npz>
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from ..common import constants as C
from .clusters import ClusterFile
from .digitize import Digitizer
from ..dream.shaper import DreamShaper

CALIB = "~/x17/response_sim/avalanche/aval_calib.json"
DMAX = 3


def decoded_path(prefix, tag, feu):
    """
    Decoded filename in the shape wft's own file handling requires.

    TWO separate wft conventions ride on this name, and getting either wrong
    produces a file that opens perfectly and yields nothing (audit A8, Fix 8):

      * `wft/io.subrun_files` globs `*_{feu:02d}.root`, so the FEU id IS the
        suffix, zero-padded — and it is also what keys the strip map. The old
        hardcoded 07/08 named a DIFFERENT physical detector (mx17_2) while the
        chain models det3, whose FEUs are 3 and 4.
      * `wft/io.file_tag` pairs the X and Y file of one subrun by regexing
        `_datrun_(.+)_\\d\\d\\.root$`, falling back to the whole basename. A
        name without `_datrun_` therefore gives the two views DIFFERENT tags,
        `reco.py`'s `by_tag` dict never assembles a pair, and the run
        reconstructs zero events with no error anywhere.

    So the name is built as `<prefix>_datrun_<tag>_<feu>.root`, mirroring the
    real `MX17_long_run_datrun_260506_09H28_000_03.root`.
    """
    if "_datrun_" in os.path.basename(prefix):
        return f"{prefix}_{feu:02d}.root"
    return f"{prefix}_datrun_{tag}_{feu:02d}.root"


def event_budget(cur, dt_s, dmax=DMAX):
    """
    Per-view budget, centred on that view's biggest channel. Reports BOTH the
    integrated area and the peak amplitude, because they behave differently
    under the front end and §9 quotes both.

    WHY BOTH, and a trap in reading the area. The ion smear and the DREAM
    shaper are linear time-invariant filters applied IDENTICALLY to every
    channel, and for such a filter integral(h conv f) = integral(h) *
    integral(f). So the inter-channel AREA ratio is mathematically invariant
    under them: switching shaping on cannot change the area budget over an
    infinite window.

    It appeared to — c1_X went 0.201 -> 0.146 — and that was pure truncation:
    the shaped tail runs ~12 tau = 1.2 us past the last avalanche, and the
    NEIGHBOUR channels are fed later than the central one by resistive
    spreading, so a fixed window clips them harder. A control run at the same
    window with shaping off returned exactly 0.201, as the invariance requires.

    Consequence: the area budget is only meaningful against data when
    integrated over the SAME window the DAQ uses (32 samples x 60 ns =
    1.92 us). Peak amplitude has no such invariance, IS genuinely changed by
    shaping, and is what §9's 1.00/0.16-0.19/0.06-0.08/0.03 row compares to.

    ALSO REPORTED: peak TIME, and the width implied by area/peak. Area and
    peak ratios are two numbers about the same waveform, and their quotient is
    a shape statement that survives any per-channel normalisation:

        w_d = (area_d / area_0) / (peak_d / peak_0)

    is how much broader the d-th neighbour's pulse is than the central one.
    The measured pair -- neighbour at ~3/4 the area but ~1/6 the peak -- says
    w_1 ~ 4, i.e. the data's neighbour signal is four times broader than its
    central one. That is the sharpest form of the §9 comparison and it needs
    no absolute calibration at all.

    Peak time gives the same physics independently and is quoted directly in
    §9 (+54-61 ns median shift for d=+-1), so it is a target rather than a
    derived quantity. Both exist because a model can match the area budget and
    the peak budget separately while getting the arrival TIMING wrong, and
    timing is what resistive spreading actually controls.
    """
    out = {}
    for view in ("X", "Y"):
        chans = {k[1]: v for k, v in cur.items() if k[0] == view}
        if not chans:
            continue
        q = {c: float(np.sum(v) * dt_s) for c, v in chans.items()}
        pk = {c: float(np.max(v)) for c, v in chans.items()}
        tp = {c: float(np.argmax(v) * dt_s * 1e9) for c, v in chans.items()}
        tot = sum(q.values())
        if tot <= 0:
            continue
        c0 = max(q, key=q.get)
        pk0 = pk.get(c0, 0.0)
        q0, tp0 = q[c0], tp[c0]
        ds = range(-dmax, dmax + 1)

        def width(d):
            # (area ratio) / (peak ratio): how much broader than the centre.
            # Undefined where the neighbour has no signal at all, so report
            # NaN rather than a zero that would average like a measurement.
            a, p = q.get(c0 + d, 0.0) / q0, pk.get(c0 + d, 0.0) / pk0
            return (a / p) if (p > 0 and pk0 > 0) else float("nan")

        out[view] = {
            "total": tot,
            "share": [q.get(c0 + d, 0.0) / tot for d in ds],
            "peak": [(pk.get(c0 + d, 0.0) / pk0) if pk0 > 0 else 0.0
                     for d in ds],
            # Shift of the neighbour's peak relative to the central channel's,
            # the quantity §9 quotes as +54-61 ns for d=+-1.
            "t_peak_shift_ns": [(tp[c0 + d] - tp0) if (c0 + d) in tp
                                else float("nan") for d in ds],
            "width_rel": [width(d) for d in ds],
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clusters")
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--calib", default=CALIB)
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--n-samp", type=int, default=3200,
                    help="1 ns samples. Must cover drift (~766 ns) + kernel "
                         "window (1000 ns) + the shaped tail (~12 tau = 1.2 us). "
                         "2000 truncates once shaping is on.")
    ap.add_argument("--no-shaper", action="store_true",
                    help="stop at induced charge, before the DREAM front end")
    ap.add_argument("--pzc-residual", type=float, default=0.75,
                    help="un-cancelled fraction beta of the 5 us CSA pole — "
                         "the undershoot mechanism (shaper.py header). A "
                         "hardware unknown: scan {0.5, 0.75, 1.0}; 0 = the "
                         "pre-2026-08-07 no-undershoot behaviour")
    ap.add_argument("--filter-model", choices=("sk", "crrc2"), default="sk",
                    help="sk = manual topology (PZC pole + zeta=0.75 "
                         "Sallen-Key); crrc2 = legacy stand-in for A/B")
    ap.add_argument("--no-ions", action="store_true",
                    help="electrons only (the pre-2026-08-07 behaviour)")
    ap.add_argument("--ion-model", choices=("measured", "analytic"),
                    default="measured",
                    help="measured = the S3 v2 i_elec+i_ion template (needs a "
                         "v2 calib); analytic = the delta + rectangle model")
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--fixed-position", action="store_true",
                    help="do NOT randomise the impact point (see below)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--decoded-out", default=None,
                    help="write sim_decoded.root (T13). Two files are made, "
                         "<prefix>_<XID>.root (X view) and <prefix>_<YID>.root "
                         "(Y), mirroring the two FEUs the real detector reads "
                         "out; the ids come from --feu-ids.")
    ap.add_argument("--feu-ids", type=int, nargs=2, default=(3, 4),
                    metavar=("X_ID", "Y_ID"),
                    help="FEU ids for the X and Y views. These ARE the decoded "
                         "file suffixes, and wft keys its strip map off them "
                         "(wft/io.py), so they must be the target run's own "
                         "ids: det3 (the bench detector this chain mimics) is "
                         "x->FEU 3, y->FEU 4 per its run_config.json. The old "
                         "hardcoded 07/08 is a DIFFERENT physical detector "
                         "(mx17_2), so a det3-configured wft found no files. "
                         "An SPS-run comparison must set these from that run.")
    ap.add_argument("--decoded-tag", default="sim_000",
                    help="the '<date>_<time>_<idx>' field of the decoded file "
                         "name. wft pairs the X and Y files of a subrun by it "
                         "(wft/io.file_tag), so both views must share it — "
                         "see decoded_path().")
    ap.add_argument("--ftst-offset", type=int, default=None,
                    metavar="TICKS",
                    help="fixed X->Y FEU clock offset in 10 ns ftst ticks. "
                         "Default: the value measured on the det3 long run "
                         "(daq.FTST_OFFSET_TICKS_DET3 = 4). Per run, not "
                         "universal — see response/dream/daq.TriggerPhase.")
    ap.add_argument("--noise", default="~/x17/response_sim/dream/noise_det3.json",
                    help="noise spec from response.dream.noise --characterise")
    a = ap.parse_args()

    cf = ClusterFile(a.clusters)
    dig = Digitizer(a.kernel, os.path.expanduser(a.calib), seed=a.seed,
                    with_ions=not a.no_ions, ion_model=a.ion_model)
    shaper = None if a.no_shaper else DreamShaper(
        filter_model=a.filter_model, pzc_residual=a.pzc_residual)
    info = dig.describe()

    print("Stage B over a real ClusterTree\n")
    ci = cf.describe()
    print(f"  clusters  {ci['file']}: {ci['n_events']} events, "
          f"{ci['n_clusters']} clusters, {ci['n_primary_per_event']:.0f} "
          f"primary e-/event, z {ci['z_mm_range'][0]:.2f}..{ci['z_mm_range'][1]:.2f} mm")
    print(f"  kernel    rho_s={info['kernel']['rho_s_MOhm_sq']} MΩ/sq  "
          f"d_k={info['kernel']['d_kapton_um']} µm")
    print(f"  gas       v_drift={info['gas']['v_um_per_ns']:.1f} µm/ns  "
          f"sigma_T(full gap)={info['gas']['sigma_T_um_full_gap']:.0f} µm")
    print(f"  avalanche gain={info['gain_mean']:.0f} sigma0={info['sigma0_um']:.1f} µm"
          f"   transparency={info['mesh_transparency']} [{info['transparency_source']}]")
    if info["with_ions"] and info.get("ion_measured"):
        im = info["ion_measured"]
        io = info["ion"]
        print(f"  ions      {im['f_ion']*100:.1f}% of the induced charge "
              f"[S3 v2 MEASURED], 50/90/99% of the ion charge by "
              f"{im['t_ion_50pct_ns']:.0f}/{im['t_ion_90pct_ns']:.0f}/"
              f"{im['t_ion_99pct_ns']:.0f} ns")
        print(f"            born at <z>={im['z_birth_mean_um']:.1f} µm; "
              f"f_ion cross-check from z-hist {im['f_ion_from_z_hist']:.4f} "
              f"vs {im['f_ion']:.4f} from the currents")
        print(f"            analytic model predicted f_ion="
              f"{io['f_ion']:.3f}, {io['t_ion_transit_ns']:.0f} ns "
              f"({io['t_ion_transit_ns_lo']:.0f}-"
              f"{io['t_ion_transit_ns_hi']:.0f} for mu ±30%)")
    elif info["with_ions"]:
        io = info["ion"]
        print(f"  ions      {io['f_ion']*100:.0f}% of the induced charge "
              f"[ANALYTIC], flat over a {io['t_ion_transit_ns']:.0f} ns transit "
              f"({io['t_ion_transit_ns_lo']:.0f}-{io['t_ion_transit_ns_hi']:.0f} "
              f"for mu +-30%)")
    else:
        print("  ions      DISABLED (electrons only)")
    if shaper:
        sh = shaper.describe()
        print(f"  DREAM     peaking {sh['t_peak_ns_5to100']:.0f} ns (5%->100%, "
              f"code {sh['peaking_code']}), tau {sh['tau_ns']:.1f} ns, "
              f"{sh['mV_per_fC']:.0f} mV/fC")
    else:
        print("  DREAM     DISABLED (charge level)")
    print()

    # RANDOMISE THE IMPACT POINT unless explicitly told not to.
    #
    # The Stage A generator is a PENCIL BEAM at (0, 0), and 512 pads is even, so
    # the active-area centre falls exactly on a pad BOUNDARY. Every muon then
    # lands at the one x where "the channel with the most charge" is a coin flip
    # and the d=+-1 profile is maximally skewed — in the first 60-event run 216
    # events rounded one way and 256 the other, producing a spurious X-view
    # asymmetry (d=-1: 0.308 vs d=+1: 0.093) that is entirely an artifact of the
    # beam position, not of the detector.
    #
    # The measured §9 numbers average over real muons landing anywhere, so the
    # simulation must too. Offsetting each event by a uniform draw over the
    # 31.2 mm superperiod in x and 1.56 mm in y is exact rather than approximate:
    # the whole stack is periodic with those periods, so a shifted track is a
    # physically valid track, and this samples every ESL phase and every
    # sub-pitch position. The proper fix is to randomise the gun in Stage A;
    # this achieves the same sampling without re-running Geant4.
    rng_pos = np.random.default_rng(a.seed + 991)

    # --- T13 decoded-file output ------------------------------------------
    # The detector is 512 x 512 pads read out as 512 X channels and 512 Y
    # channels, i.e. TWO 512-channel FEUs — which is why the data arrives as
    # *_07.root and *_08.root. Stage B already works in absolute pad
    # column/row (digitize.induce), and mx17_m1_map.csv gives
    # channel_num = 0.78 mm * index on both axes, so the map is the identity:
    # X channel = column, Y channel = row.
    #
    # THE REAL WIRING, resolved 2026-08-07 from det3's own run_config.json
    # (`detectors[0].dream_feus` / `dream_feu_orientation`) rather than assumed:
    #
    #     x_1..x_8 -> FEU 3, connectors 1..8      y_1..y_8 -> FEU 4, connectors 1..8
    #     every connector on BOTH views: "inverted"
    #
    # So connector-to-connector is the identity within a view (the structural
    # assumption held), but the channel ORDER inside each connector is
    # REVERSED, on both views — which the identity map here does not apply.
    # This is a pure relabelling: it cannot change any physics and the
    # simulation is self-consistent under it, so every observable computed here
    # is unaffected. It matters only for a channel-by-channel comparison
    # against a real run, i.e. at T14, and it is left explicit rather than
    # silently applied because the exact meaning of "inverted" (per connector
    # vs across the whole FEU) should be read off Mx17StripMap.RunConfig
    # instead of guessed a second time.
    daq = adc_buf = trig = ftst_buf = None
    if a.decoded_out:
        from ..dream.daq import (Daq, N_CHAN, charges_to_fc, write_decoded,
                                 TriggerPhase)
        with open(os.path.expanduser(a.noise)) as fh:
            daq = Daq(json.load(fh), seed=a.seed + 17)
        # ONE trigger phase per event, shared by both FEUs, recorded in ftst
        # with the semantics measured off det3 (daq.TriggerPhase).
        from ..dream.daq import FTST_OFFSET_TICKS_DET3
        trig = TriggerPhase(
            offset_ticks=(FTST_OFFSET_TICKS_DET3 if a.ftst_offset is None
                          else a.ftst_offset),
            sample_period_ns=daq.dt_ns)
        adc_buf = {"X": [], "Y": []}
        ftst_buf = {"X": [], "Y": []}
        print(f"  decoded   writing "
              f"{os.path.basename(decoded_path(a.decoded_out, a.decoded_tag, a.feu_ids[0]))}"
              f" (X, FEU {a.feu_ids[0]}) and "
              f"{os.path.basename(decoded_path(a.decoded_out, a.decoded_tag, a.feu_ids[1]))}"
              f" (Y, FEU {a.feu_ids[1]}), "
              f"{daq.n_sample} samples x {daq.dt_ns:.0f} ns")
        print(f"  trigger   one phase per event, ftst {trig.n_ftst} states x "
              f"{trig.tick:.0f} ns, X->Y offset {trig.offset_ticks} ticks")

    accX, accY, nX, nY = np.zeros(2 * DMAX + 1), np.zeros(2 * DMAX + 1), 0, 0
    pkX, pkY = np.zeros(2 * DMAX + 1), np.zeros(2 * DMAX + 1)
    totX, totY = [], []
    # Peak-time shift and relative width are collected per event and reduced
    # with a MEDIAN, not a mean: §9 quotes a median shift, and width_rel is
    # a ratio that can blow up on events where a neighbour is nearly empty.
    tsX, tsY, wX, wY = [], [], [], []
    t0 = time.time()
    n_done = 0
    for ev, g in cf.events():
        if a.max_events and n_done >= a.max_events:
            break
        ox = oy = 0.0
        if not a.fixed_position:
            ox = rng_pos.uniform(0.0, C.SUPERPERIOD_M * 1e3)
            oy = rng_pos.uniform(0.0, 2 * C.PAD_PITCH_M * 1e3)
        x, y, t, q = dig.transport(cf.x[g] + ox, cf.y[g] + oy,
                                   cf.z[g], cf.t[g], cf.n_e[g])
        if len(x) == 0:
            continue
        cur = dig.induce(x, y, t, q, a.n_samp)
        if shaper is not None:
            # Shape each channel. The budget then integrates the SHAPED
            # waveform, which is what the data analysis sees.
            cur = {k: shaper.apply(v, dt_ns=dig.lut.dt * 1e9)
                   for k, v in cur.items()}
        if daq is not None:
            # Lay the 1 ns shaped channels into a dense 512-wide plane per
            # view, then point-sample onto the 60 ns DAQ grid. The phase is
            # random per EVENT (the sampling clock is uncorrelated with a
            # cosmic) but SHARED by the two views: one trigger, one phase,
            # differing only by the fixed board offset, and recorded in ftst.
            from ..dream.daq import N_CHAN as _NC
            phases = trig.draw(daq.rng)
            for view in ("X", "Y"):
                plane = np.zeros((_NC, a.n_samp), dtype=np.float32)
                for (v, idx), w in cur.items():
                    if v != view or not (0 <= idx < _NC):
                        continue
                    plane[idx, :len(w)] += w[:a.n_samp]
                ph, ftst = phases[view]
                smp, _ph = daq.sample(plane, phase_ns=ph)
                adc_buf[view].append(smp.T)          # (n_sample, n_chan)
                ftst_buf[view].append(ftst)

        b = event_budget(cur, dig.lut.dt)
        if "X" in b:
            accX += np.array(b["X"]["share"]); pkX += np.array(b["X"]["peak"])
            nX += 1; totX.append(b["X"]["total"])
            tsX.append(b["X"]["t_peak_shift_ns"]); wX.append(b["X"]["width_rel"])
        if "Y" in b:
            accY += np.array(b["Y"]["share"]); pkY += np.array(b["Y"]["peak"])
            nY += 1; totY.append(b["Y"]["total"])
            tsY.append(b["Y"]["t_peak_shift_ns"]); wY.append(b["Y"]["width_rel"])
        n_done += 1
        if n_done % 50 == 0:
            el = time.time() - t0
            print(f"    {n_done} events, {el:.0f}s "
                  f"({el/n_done*1e3:.0f} ms/event)", flush=True)

    if not nX or not nY:
        print("no events produced charge — check the kernel and the calib")
        return 1
    accX /= nX
    accY /= nY
    pkX /= nX
    pkY /= nY
    ds = list(range(-DMAX, DMAX + 1))

    print(f"\n  {n_done} events in {time.time()-t0:.0f} s\n")
    print(f"  {'d':>4s}" + "".join(f"{d:>8d}" for d in ds))
    print(f"  {'X':>4s}" + "".join(f"{v:8.3f}" for v in accX))
    print(f"  {'Y':>4s}" + "".join(f"{v:8.3f}" for v in accY))

    print("\n  peak-amplitude ratios  [measured §9: 1.00 / 0.16-0.19 / "
          "0.06-0.08 / 0.03]")
    print(f"  {'X':>4s}" + "".join(f"{v:8.3f}" for v in pkX))
    print(f"  {'Y':>4s}" + "".join(f"{v:8.3f}" for v in pkY))

    with np.errstate(invalid="ignore"):
        tsXm = np.nanmedian(np.array(tsX, dtype=float), axis=0)
        tsYm = np.nanmedian(np.array(tsY, dtype=float), axis=0)
        wXm = np.nanmedian(np.array(wX, dtype=float), axis=0)
        wYm = np.nanmedian(np.array(wY, dtype=float), axis=0)

    print("\n  peak-time shift vs central [ns]  [measured §9: +54-61 for d=±1]")
    print(f"  {'X':>4s}" + "".join(f"{v:8.1f}" for v in tsXm))
    print(f"  {'Y':>4s}" + "".join(f"{v:8.1f}" for v in tsYm))

    print("\n  relative width (area ratio)/(peak ratio)  [measured §9: ~4 at d=±1]")
    print(f"  {'X':>4s}" + "".join(f"{v:8.2f}" for v in wXm))
    print(f"  {'Y':>4s}" + "".join(f"{v:8.2f}" for v in wYm))

    c1X = 0.5 * (accX[DMAX - 1] + accX[DMAX + 1])
    c1Y = 0.5 * (accY[DMAX - 1] + accY[DMAX + 1])
    c2X = 0.5 * (accX[DMAX - 2] + accX[DMAX + 2])
    c2Y = 0.5 * (accY[DMAX - 2] + accY[DMAX + 2])
    xfrac = np.mean(totX) / (np.mean(totX) + np.mean(totY))

    print(f"\n  TRACK-LEVEL ±1 share   X {c1X:.3f}   Y {c1Y:.3f}"
          f"      [measured §9: 0.23-0.28]")
    print(f"  TRACK-LEVEL ±2 share   X {c2X:.3f}   Y {c2Y:.3f}")
    print(f"  X/Y charge balance     {xfrac:.3f} / {1-xfrac:.3f}"
          f"          [measured: 0.49/0.51]")
    # Say what the impact-point sampling ACTUALLY is. --fixed-position means
    # "add no post-hoc offset", which is correct when Stage A's gun already
    # randomised (--beam-spread) and wrong only for a genuine pencil beam. The
    # old wording asserted "pencil beam, on a pad boundary" either way, which
    # mislabels every run over a spread production.
    spread = float(np.ptp(np.mod(np.asarray(
        [cf.x[g].mean() for _, g in cf.events()]) * 1e-3, C.PAD_PITCH_M)))
    if not a.fixed_position:
        src = "randomised in post over the 31.2 mm superperiod"
    elif spread > 0.5 * C.PAD_PITCH_M:
        src = "from Stage A's gun (--beam-spread); no post-hoc offset"
    else:
        src = "FIXED and sub-pitch-degenerate — a pencil beam on a pad boundary"
    print(f"\n  impact point: {src}")
    missing = [n for n, on in (("ion tail", info["with_ions"]),
                               ("DREAM shaping", shaper is not None)) if not on]
    missing += ["ZS", "noise"]
    print(f"\n  still missing: {', '.join(missing)}")

    res = {"n_events": n_done, "d": ds,
           "share_X": accX.tolist(), "share_Y": accY.tolist(),
           "peak_X": pkX.tolist(), "peak_Y": pkY.tolist(),
           "c1_X": float(c1X), "c1_Y": float(c1Y),
           "c2_X": float(c2X), "c2_Y": float(c2Y),
           "x_fraction": float(xfrac),
           "digitizer": info, "clusters": ci,
           "shaper": shaper.describe() if shaper else None}
    if daq is not None:
        from ..dream.daq import charges_to_fc, write_decoded
        res["trigger_phase"] = trig.describe()
        for view, feu in (("X", a.feu_ids[0]), ("Y", a.feu_ids[1])):
            block = np.stack(adc_buf[view])                 # (n_ev, n_samp, ch)
            # Stage B works in elementary charges; the ADC scale is per fC.
            adc = daq.to_adc(charges_to_fc(block), n_ev=len(block))
            path = decoded_path(a.decoded_out, a.decoded_tag, feu)
            write_decoded(path, adc,
                          event_ids=np.arange(len(block), dtype=np.uint64),
                          ftst=np.asarray(ftst_buf[view], dtype=np.uint16))
            occ = float((adc.astype(np.int32)
                         - np.median(adc, axis=(0, 1))).max(axis=1).mean())
            print(f"  wrote {path}  {len(block)} events, "
                  f"mean per-event max excursion {occ:.0f} ADC")

    if a.out:
        with open(a.out, "w") as f:
            json.dump(res, f, indent=1)
        print(f"\n  wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
