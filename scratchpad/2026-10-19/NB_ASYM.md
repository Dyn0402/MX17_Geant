# X left/right neighbour asymmetry and the neighbour timing (2026-10-19)

## Update, later on 10-19: every det3 dataset, the n_TOF beam, the bus solve, displaced feedback

This supersedes two statements below: "det3-only" (other June builds have it too, smaller), and the bus traces as the
candidate for det3's size (the solve gives ~0.01, not 0.1–0.2).

### A. Copy survey (`dump_det.py`, `run_survey.sh`, `survey_table.py`, `survey_plot.py` → `copy_survey.{csv,txt,png}`)

- **Scope:** M3-free, isolated single clusters, 3 files per sub-run.
  - Every decoded det3 sub-run: May/June-6 builds, the June-15 build on 6-22/23/25/26/27/28, every HV and drift scan
    point.
  - n_TOF beam run_79, run_145 and part of run_55, all four chambers.
  - Long runs of the other boards.
- **Estimator:** wiring-independent a = in-phase content of channel p+1 minus channel p−1. It reads ~20 % below the
  M3-vertical value (weekend 0.08 vs 0.10).
- **Closure (`uncopy.py`):** undoing x'_{n+1} = x_{n+1} − a x'_n with one a removes the ±1 asymmetry AND balances ±2
  (6-22: ±1 0.153 → 0.005; ±2 0.092/0.051 → 0.042/0.049). The one-step channel n → n+1 copy model holds.

| det3 build / date | X copy a | notes |
|---|---|---|
| "Bulked 5-5-26" (May 5–6; "Det3 from May beam test" June 6), FEU 3 | **−0.019…+0.005 (≈ 0)** | HV scan 460–530 V, all ≈ 0 |
| "Bulked by Stephan June 15": 6-22/23, FEU 3 | **0.11–0.17** | except 505–525 V, where most peaks rail (0.07 → 0.006) |
| 6-25 (drift 500 V) | 0.097 | |
| 6-27/28, FEU 7 | 0.06–0.11; long run and weekend 0.08 | flat in mesh HV |
| n_TOF run_79 (July 26, Ar/iso 90/10), FEU 3 | **0.053–0.066** (shape fit 0.085, in phase) | beam angle asymmetry ±0.04 (Y views) |
| Y view, every date | ≈ 0 | |

- **Drift field:** the 6-27 drift scan ran 1100 → 100 V in time, and a rose 0.053 → 0.163. The long run that started
  right after at 1000 V gave 0.081, back to the 1000 V value. **The copy grows at low drift field, reversibly.** Over
  weeks at a fixed field it falls (6-22 0.13–0.14 → 6-27 0.08 at 1000 V).
- **Other boards (bench):**
  - det1 (Jan build) 0.00.
  - det2 (Arnaud, June 12) 0.02.
  - det6 (Stephan, June 24) 0.00–0.03.
  - det4 (Stephan batch, June 22) 0.02–0.07.
  - det7 (same batch) 0.06–0.08.
- **Beam controls:** det2 0.06–0.09, within the angle systematic. det6/det7 beam X fail the quality flag (α+1 > 0.6, or
  n−2 ≈ n−1: a position-map/channel problem in those configs, not a copy).
- **Dropped:** the June-17 X of the June-15 board (2 days after bulking). Its central waveform rings with a ~600 ns
  period (`asym_shape_jun17.png`), so the −0.18 there is not a measurement.
- **Reading:** the copy belongs to the June builds, largest on the June-15 board, and **changes with time and with
  drift field**. A fixed board geometry can do neither.

### B. What the board geometry gives (`bus_wpot.py`, `bus_scan/`)

- **Stack thickness:**
  - The gerbers carry no stack-up. "top-resist" is the print area plus labels; the strips are printed after
    fabrication (RESPONSE_SIM_PLAN §3 row 3).
  - The Geant4 stack (`shared/MX17ModuleGeometry.hh`, NEEDED_INPUTS §6): five 26 µm copper layers, 50 µm kapton
    coverlay, ~19 µm glue over the pads, and the CAD 1.70 mm board.
  - That leaves the FR4 as a residual, 264 µm per layer, equal division. That puts L5 ~290 µm and L6 ~580 µm below the
    pad face.
  - A typical 8-layer 1.6 mm board would have 100–250 µm per layer.
- **The solve:** 3D finite-volume prompt weighting potential of one channel: own pads, 0.25 mm vias with 0.5 mm rings,
  390 µm stubs and the one-sided bus, all at 1; everything else grounded, including the other view's buses on the other
  layer. Read at the resistive layer and folded with a cloud (σ 300/450 µm, insensitive).
  - Converged on a 10 µm grid (30 µm: 0.0195 → 20 µm: 0.0114 → 10 µm: 0.0099).

  | geometry | copy |
  |---|---|
  | nominal (L5 bus, 264 µm) | **0.010** |
  | bus on L6 (under the other view's L5 lines) | 0.001 |
  | no bus (control) | 0.002 |
  | L5 at 132 / 66 / 33 / 15 µm | 0.023 / 0.043 / 0.065 / 0.085 |
  | gap 150 / 200 / 250 µm | 0.021 / 0.045 / 0.065 |

- **Reading:** the one-sided bus explains a ~0.01 copy on a normal board, in the view whose buses sit on L5, and ~0 in
  the other view. That matches det1/det2/det6. It cannot reach 0.1–0.2 even touching the pads, and it can't change with
  time or field.
- **The extra on the June builds** therefore needs a dielectric or surface property that changes:
  - Moisture or uncured material from the bulk process (fresh coverlay/photoresist, drying over weeks in the gas flow)
    raises the effective coupling of the buses under the gaps.
  - Or charge-up of the insulating surfaces near the gaps (field dependent).
  - Both are speculative. A moisture origin would tie in with the bench gas reading wet (~1 % H2O, FRESH_LOOK §12).
- **Practical:** det3 data can be corrected per run with `uncopy.py` (one a per sub-run, from the X data themselves).

### C. The far prompt halo is displaced photon feedback (`nb_fb_emul.py`, Stage B `--feedback-spread-um`)

- **Feature:** every board, both views, shows a symmetric in-phase halo at ±3…±5 (~1 % at ±3) that the sim lacks.
- **Model:** a UV photon from an avalanche at the anode reaches the mesh, a height H above, at P(R > r) = H/√(r²+H²).
  Its projection is a Cauchy of half-width H, a 1/k² tail per strip. Stage B already gives feedback bG ≈ 0.45–0.48 of
  the charge, but places it on the parent.
- **Offline emulation** on the sim's mean waveforms (moving f = bG of each strip's signal with that kernel), at
  H = 150 µm = the amplification gap:

  | | X α ±3 / ±4 / ±5 | Y α ±3 / ±4 / ±5 |
  |---|---|---|
  | v3 + feedback | 0.016 / 0.0053 / 0.0026 | 0.020 / 0.0054 / 0.0025 |
  | data | 0.012 / 0.0055 / 0.0034 | 0.026 / 0.0078 / 0.0026 |
  | v3 alone | 0.004 / 0.0005 / 0 | 0.007 / 0 / 0 |

  It moves the ±1 peak only 15–30 ns earlier and adds share, so it is not the ±1/±2 timing fix.
- **Implemented:** `response/digitizer/digitize.py` `add_feedback_spread`, `run.py --feedback-spread-um H
  [--feedback-rmax-mm 10]`. 0 keeps the historical path bit-identical. Test `response/digitizer/test_feedback_spread.py`
  (gain conserved to 1e-4, Cauchy marginal, cap, parent times).
- **Condor 4481160** (`stageBC_w5_fbspread.sub`, pre-registered): v5a/b/c = v3w / v3ws300 / v4c bases + H 150, 5 angles.
  Evaluated by `eval_v5.sh` (yoga). Results: §D.

### E. Next (supersedes the list at the end)

1. Use copy-free X (or `uncopy.py`-corrected det3 data) everywhere X sharing or width is a target. Reconsider v4c's
   σ0 420 µm: v5a matches X without it.
2. ±1/±2 timing on top of v5a (above): the data's ±1 neighbour has less late charge.
3. Y shares (v5a low): a Y-only lever.
4. The det3 copy's origin: what was different about the June-15 bulk (process, drying, coverlay)? A quick check of the
   det3 X D_T excess (+100 µm/√cm, nTof_x17 `x-width-excess` note) on `uncopy.py`-corrected waveforms.
5. If the bulk process is known: a weighting solve with a raised-ε (wet) dielectric under the gaps
   (`bus_wpot.py` has the geometry; it needs an ε option per region).

### D. v5 results (condor 4481160, all 15 done; `ledger_v5*.log`, `xy_extent_v5.log`, `nb_prompt_v5.log`, `nb_side_v5.log`, `nb_width_v5.log`)

Yoga had no AFS token (`kinit` needed), so the points were pulled through the desktop and forwarded:
`SKIP_PULL=1 ./eval_v5.sh`. th00 vs the weekend det3 data (X on the copy-free + side, Y folded):

| | data X | v5a X | v5b X | v5c X | data Y | v5a Y | v5b Y | v5c Y |
|---|---|---|---|---|---|---|---|---|
| α ±3 / ±4 / ±5 | .012/.0055/.0034 | **.011/.0059/.0032** | .013/.0077/.0035 | .019/.011/.0037 | .026/.008/.0026 | .015/.0062/.0032 | .018/.0071/.0038 | .027/.0097/.0049 |
| α ±1 / ±2 | .334/.050 | **.35/.056** | .40/.092 | .42/.148 | .44/.122 | .39/.082 | .46/.098 | .53/.136 |
| width rms, copy-free [mm] | .694 | **.707** | .755 | .817 | .915 | .843 | .872 | .912 |
| ±1 / ±2 peak delay [ns] | +75 / ~+380 | +105..120 / +170 | +75 / +105 | +75 / +50 | +60 / ~+400 | +105 / +330 | +90 / +330 | +45 / +255 |
| data/sim p50 (th00 X / Y) | | **.976 / 1.015** | 1.087 / 1.136 | 1.119 / 1.182 | | | | |

v5a = wet gas, σ0 = calib 94 µm, bG 0.45, H 150 µm.
v5b = v5a with σ0 300 µm.
v5c = v4c (σ0 420 µm, bG 0.48) + H 150 µm.

Ledger, all angles: v5a X 0.93–0.99, Y 1.02–1.15 (inclined Y still dim, §21.9).
Extent: v5a vertical X rms 0.67 vs data 0.66 (`xy_extent` thresholds), 10° 1.47 vs 1.46, 20° 2.75 = data.
v5c's extent is too wide (0.81).

**Against the pre-registration:**
- Far α, all bands pass (X ±5 0.0032 vs a 0.003 edge; Y ±4 0.0062 vs 0.006).
- ±1 shift: 15 ns in X ✔, 45 ns in Y ✘ (band 0–30).
- ±1 share: +0.03 ✔.
- Amplitude: +16 % (0.84 → 0.976) vs a +5–12 % band ✘ (direction right).
- The falsifier did not fire. Displaced feedback reaches the far strips through the chain as the emulation said.

**Reading:**
- **v5a is the first point that matches the copy-free X view with no retuning:** far halo, ±1/±2 shares, vertical and
  inclined widths, amplitude. The σ0 420 µm of v4c is no longer needed for the X width or shares; it had been fitting
  the copy and the halo.
- **Still open:**
  1. The ±1 peak is 30–45 ns late in both views (X +105..120 vs +75, Y +105 vs +60), and ±2 peaks early (+170 vs ~+380
     in X). The data's neighbours carry less late charge at ±1 and more slow charge at ±2. This is the remaining timing
     question, now isolated from the copy and the halo.
  2. Y shares are low in v5a (0.39 vs 0.44; ±2 0.08 vs 0.12). v5b is closer in Y (0.46/0.10) but wider than the data
     in X. A Y-only spread (the resistive strips along y) is the natural lever (§21.9).
- **Candidate next arms:** v5a + Y-side lateral spread (ρ_s or strip line/gap); v5a with H 100/200; v5a + bG refit
  (amplitude), since the spread lowers the peak strip.

Follow-up to `NB_FAR.md` (same day). Question: is the X −1/+1 asymmetry a beat with the resistive strips, and what does
it do to the neighbour timing? Larger statistics come from the 6-28 weekend det3 run (same detector, slot, HV as the
6-27 sample) and the FEU 3/4 runs (6-22, 6-23). Other boards (det2, det4, det6, det7) serve as controls.

## Inputs

- `dump_tagged.py` (yoga) writes per-event ±5 position-ordered waveforms with event id, peak channel, rank, channel ids,
  noise and M3 position/angles, or truth position for sim (`--sim`). Driver `run_dumps.sh`. Outputs are in yoga
  `/media/dylan/data/claude_tmp/tmp/nb/tagged/*.npz`: sat_det3, g_det3_wknd, o22 (= g_det3), o23, g_det2, g_det4,
  g_det6_long, g_det7_long, sim_v3 (v3trw/v3t_th00), sim_v4c.
- Vertical = M3 |ax|, |ay| < 3°. The peak strip must be < 3500 ADC (not railed) and > 200 ADC.
  Weekend: 908 X / 959 Y events (2.7× the 6-27 sample). 6-22: 651 / 725.
- M3 can't resolve the track position within a strip (pointing ~1 mm vs 0.78 mm pitch), so that split is uninformative.
  det3 is rotated 90° to M3: X position ≈ ref_y + 209.1 mm, Y position ≈ −ref_x + 163.0 mm (6-27, z = 702).

## 1. It isn't the beat. It's a one-sided in-phase copy (`nb_asym.py`, `nb_asym_shape.py`, `nb_side.py`)

- **Flat in beat phase.** Peak rank mod 40 (39 ESL : 40 pads, 31.2 mm) shows no modulation. 6-22 sine fit amplitude
  0.020 ± 0.024 on a 0.142 offset. The sim does show a small beat (A2 amplitude 0.028 ± 0.010, Δχ² 7), but its mean is zero.
- It is also flat in connector, channel-in-connector and rank parity.
- **Shape:** the mean difference (−1) − (+1) is the central waveform scaled. A fit d1 = a·c + b·dc/dt leaves a residual
  rms of 0.002–0.010 (from 0.03–0.09), with b ≈ 0.
- **Copy model:** the strip at the lower position picks up a·(its upper neighbour):
  n−1 = s1 + a·c, n+1 = s1 + a·s2, n−2 = s2 + a·s1.
  a from the ±1 difference and a from the ±2 difference agree:

  | run | X FEU | HV mesh/drift | a (±1) | a (±2) |
  |---|---|---|---|---|
  | det3 6-27 | 7 | 490/1000 | 0.085 ± 0.026 | 0.085 ± 0.026 |
  | det3 6-28 weekend | 7 | 490/1000 | 0.145 ± 0.014 | 0.142 ± 0.016 |
  | det3 6-22 | 3 | 495/1000 | 0.196 ± 0.017 | 0.167 ± 0.019 |
  | det3 6-23 | 3 | 495/600 | 0.209 ± 0.017 | 0.127 ± 0.030 |
  | det3 Y, all runs | 8 / 4 | | −0.015 … +0.009 (≈ 0) | |
  | sim v3 / v4c | | | 0.006 / −0.005 | |

- **Direction:** in channel numbers, channel n+1 gets a copy of channel n (`xtalk_map.py`; connector-inverted wiring
  makes that the lower position). No far "routing partner" anywhere on the FEU, so it isn't a fan-out pairing.
- **det3-specific:** per-event interior estimate (`nb_copy_chip_edge.py`): det3 0.09–0.17, det2 0.013 ± 0.011,
  det6 0.015 ± 0.021.
- **Linear:** a is flat with the central amplitude (`nb_copy_vs_amp.py`). 6-23 gives 0.21 ± 0.02 already at 200–700 ADC.
- **Not a tilt:** flat in track angle (`nb_copy_vs_angle.py`). A real lean gives −0.0085 per degree on det2. The det3 copy
  would need ~20°.
- **Not common mode:** far strips on the same chip carry only a flat −0.004 offset (`far_cm.py`).
- **Not in the noise** (`noise_corr.py`, `noise_corr_inject.py`): a copy produced in the electronics would copy noise too.
  Injecting x_n += A·x_{n−1} into the real pre-signal noise raises the channel-distance-1 correlation excess by
  +0.06…0.12 (A = 0.1) and +0.13…0.22 (A = 0.2). det3 X measures +0.03 / −0.03 / −0.02, the same as boards without a
  copy (det2 +0.004, det4 0.000). **So the copy is in the signal only, induced on the detector, not added by the
  FEU/cable/multiplexer.** (A directional regression can't see a copy: it correlates both ways. Checked by injection.)
- **Candidate mechanism, not yet modelled:** on both bus layers every bus line runs midway between two pad rows or
  columns, and every pad joins it on the same side (`gerber_stubs.py`, gerber window dumps). L4 "TrackX" has horizontal
  buses at y = k·0.78 mm with every stub −0.39 mm in y. L3 "TrackY" has vertical buses with every stub +0.39 mm in x.
  So each channel's trace lies under the inter-pad gap on one side and can see charge over that neighbour: a one-sided,
  prompt, linear, signal-only copy. The sim's kernels drive the pads only, with no buses, so they come out symmetric.
  Open: why det3 is ~10× other boards (registration/dielectric?), why det3 moves between runs (0.09–0.21, ~3.5σ between
  the FEU 7 and FEU 3 slots, which also differ in HV), and why only this view.
  Note: `response/common/channel_map.py` read the stub direction as the bus direction ("X layer busses along y"); the
  long bus lines run the other way. The checkerboard result is unaffected.

## 2. Timing with the copy removed (`nb_side.py`, weekend det3; ns vs the central's peak)

The + side in X is clean (to first order a·s2). The folded ±1 used up to now averaged in the copy.

| | data X (+1, clean) | data X (−1) | data Y ±1 | v4c X / Y | v3 X / Y |
|---|---|---|---|---|---|
| ±1 peak delay | **+75 ± 18** | +45 ± 8 | +60…75 ± 12 | +75 / +60 | +135 / +135 |
| ±1 share (peak / central) | **0.298 ± 0.008** (s1) | 0.42 | 0.41 | 0.387 / 0.450 | 0.368 / 0.437 |
| ±1 50 % edge | −21 ± 7 | −16 ± 4 | −1 / −16 ± 5 | +2 / +1 | +22 / +34 |
| ±2 peak | +405 ± 101 | +360 ± 27 | +405 / +390 ± 15 | +90 / +300 | +225 / +330 |
| ±2 prompt α | 0.050 ± 0.003 | 0.085 | 0.12 / 0.13 | 0.098 / 0.10 | 0.082 / 0.09 |
| ±3 prompt α | 0.012 | 0.017 | 0.023 / 0.029 | 0.004 / 0.008 | 0.004 / 0.007 |
| ±4 prompt α | 0.0055 | 0.0058 | 0.006 / 0.0095 | 0.000 | 0.000 |
| ±4 Y 50 % edge | | | +526 / +591 ± 90 | +869 / +902 | +881 / +844 |

1. **±1 delay:** with the copy removed, X is +75 ± 18 ns, which equals v4c. Y +60 also equals v4c. The timing at ±1 is
   closed by v4c in both views. The folded X +45 was partly the copy.
2. **±1 share:** both sims overshare once the copy is removed. X 0.30 data vs 0.39 v4c (+30 %); Y 0.41 vs 0.45 (+10 %).
   v4c's σ0 420 µm was tuned against the folded, copy-inflated X numbers. The same applies to the vertical X width target
   (0.66 mm data), which includes the copy on one side. The share (wants small σ0) and the delay (wants large σ0) now
   pull apart, so σ0 alone can't fit both. The data's ±1 is early *without* extra share, so the sim has too much late
   neighbour charge, rather than too little prompt charge.
3. **±2:** the data's clean X +2 is half the sim's (0.050 vs 0.098) and flat-topped: it rises promptly and peaks late
   (+400). v4c's ±2 is a prompt copy peaking at +90. Y ±2 peaks at +400 in data vs +300 in v4c.
4. **Prompt tail at ±3…±5, both views, every detector:** a roughly exponential in-phase tail (X 0.012 → 0.0055 → 0.0034,
   Y 0.026 → 0.008 → 0.0026) that neither sim has (≤ 0.008 at ±3, 0 beyond). This is the "prompt step" in `NB_FAR.md`.
   It isn't common mode, and it is roughly symmetric (Y + side ~20 % higher, 2σ). At 2–4 mm it is far beyond induction
   through a ~200 µm stack, so an inter-channel coupling ladder is the natural reading. Its noise footprint (~1 % at d = 3)
   is below what `noise_corr.py` can see against the 0.1–0.3 correlated noise.
5. **Y far strips** keep a slow resistive rise in both data and sim, but the data's is larger (±3 0.045 vs 0.031, ±4
   0.019 vs 0.012) and reaches half height ~300 ns earlier at ±4.

## Next

- Use the + side (copy-free) X numbers as the sim targets for ±1/±2 share, the vertical X width and any X neighbour
  timing. Recompute the X width without the − side.
- Model the bus traces in the weighting field (a 3D solve: pads + one-sided buses at L5/L6 under the gaps) to see whether
  the copy size comes out at ~0.01 (other boards) or needs a det3 defect. The PCB stack-up (dielectric thickness to
  L5/L6) is unknown and is the main input.
- Sim: try reducing the late neighbour charge (not adding prompt) at ±1: diffusion (wet), ρ_s, strip line/gap. Add a
  symmetric prompt exponential tail (λ ≈ 0.4–0.5 per strip, ~1 % at ±3) as a readout term and check the ladder.
- Ask whether det3's X cabling/adapter moved with it between slots (6-22/23 vs 6-27/28).
