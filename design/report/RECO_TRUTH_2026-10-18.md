# What the response sim says about our track reconstruction

<p class="meta">Living record, started 2026-10-18 (dylan-Yoga). Source: MX17_Geant <code>design/report/RECO_TRUTH_2026-10-18.md</code>; scripts and logs in <code>scratchpad/2026-10-18/reco_truth/</code>.</p>

## Why

The det3 response simulation (v3 baseline, FRESH_LOOK §21) now matches the cosmic-bench data in amplitude, X/Y charge
sums, inclined strip extent and the central-strip timing. It still misses the ±1 neighbour delay (sim +105..+135 ns,
data +45/+60), the rise (+12 ns) and the inclined-Y peak strip. Unlike the data, it knows the truth: the gun angle
(a pencil beam), the impact point, the drift depth of every electron. Three questions for the cosmic and beam
reconstruction:

1. **Angle closure.** Calibrate wft on the sim exactly as on data (with the true line replacing M3), reconstruct, and
   compare with the gun. Any bias left is a property of the method and applies to data angles too.
2. **Depth-resolved residuals.** The data cluster walks 0.39 → 0.55 mm from the M3 line between mesh and cathode, and
   σ_θ ≈ 1.0–1.1° was called a physics floor on a toy. The sim splits diffusion from reconstruction.
3. **The hit-time ladder.** RECONSTRUCTION_BASIS.md bans hit times because neighbour copies compress the time ladder
   by 20–30 %. Measure the compression with truth, and see how it moves with the neighbour delay.

Ground rules: sim quantities that were *tuned* to data (amplitude, extent, Y wings) validate nothing. Predictions are
written before each run. Angles from wft are in the Geant frame (ANGLED_LADDER §1); the sim drifts at 38.0 µm/ns
(Stage B log at 26.5 mm), not the data bundle's 36.6.

## Status

| item | state |
|---|---|
| Truth sidecars | `run_stageBC_w2.sh` dropped `*_truth.parquet` with the job's work dir. Fixed; v3 and v3+coupling regenerated with truth (condor 4480259, 24 jobs) |
| New angles | Stage A at 5°, 15°, 30° in each view and 10°×10° oblique, made on yoga (seeds 101606–101612) |
| **Wiring bug** | found and fixed for analysis (below). Affects every position-based sim-vs-data comparison at 20° |
| 1. angle closure | **done (sim-calibrated arm)**: falsifier fired, but the data clear wft; the sim is the one off. r06-on-sim arm running |
| 2. depth residuals | **done**: the data's growth with depth is mostly M3's pointing |
| 3. hit-time ladder | **done**: the compression is geometry; the sim's readout adds twice the data's; coupling moves it the wrong way |
| Sim | **wet gas + 0.42 mm prompt anode spread closes the ±1 neighbour delay** (+60 vs data +45/+60; v3 +135) and the vertical X cluster; Y still open. Closure yardstick on it running |

## Found on the way: the sim's channels were not on det3's wiring

Stage B writes pad column k into decoded channel k. det3's run_config marks every connector "inverted", and wft's
strip map (used by every sim-vs-data script, not only wft) reads channel j at 0.78 mm × (64·(j div 64) + 63 − j mod 64).
So inside a connector the sim is read mirrored, and a track that crosses a 64-channel boundary is read as two pieces
~99 mm apart. run.py's comment called the identity a pure relabelling. That is true channel by channel, but not for
anything that reads strip positions.

- The mirror explains ANGLED_LADDER §1's sign puzzle (wft "reproducing the Geant-frame sign"): the mirror and the
  depth convention cancel.
- The sim beam sits on channels ~262–298 at 0–10° (3 % of tracks cross channel 256), but at 20° its centre moves
  to ~320, a boundary, and **more than a quarter of the 20° tracks are torn**.
- Effect on the published comparisons (`xy_extent_m3.py`, v3 baseline): nothing moves at 0° or 10°. **The 10°
  extent behind the 26.5 mm column choice stands.** At 20° the sim's RMS extent drops from 2.90 / 2.88 to
  **2.72 / 2.73 mm** (data 2.75 / 2.71): the 5 % excess at 20° was the wiring.
- Fix: `response/dream/rewire.py` rewrites existing decoded files (verified: the det3 map then reads exactly
  0.78 mm × pad column, and the permutation is its own inverse). New Stage B runs can write the real wiring directly
  with `run.py --wiring connector-inverted` (default stays `identity` so old runs reproduce).
- Not re-checked yet: the amplitude ledgers. They use peak amplitudes and should not care, but a torn 20° track
  also changes which strip is the "peak strip" if the selection is position-based.

## 3. The hit-time ladder

**Method** (`reco_truth/ladder.py`). Per event, in the inclined view: strips ≥ 5σ and ≥ 100 ADC, no saturated strip,
≥ 4 strips. Hit time per strip = 50 % leading edge of its own waveform. Fit time against strip position
(amplitude-weighted least squares). **C = fitted dt/dx ÷ true dt/dx**, true = 1/(v tan θ). C = 1 is an unbiased ladder,
C < 1 a compressed one. Sim: θ = gun angle, v = 38.02 µm/ns. Data (det3 long run, 490 V): θ = each event's own M3 angle,
other view within 3° of vertical, v = 36.6 µm/ns.

The same fit on the **production `combined_hits` times** of the data gives C within 0.01 of the waveform estimator
(10°: 0.841 / 0.818 vs 0.851 / 0.824; `ladder_hits.py`), so this ladder is the hits chain's ladder.

**The new control: an ideal ladder with no readout at all** (`ideal_ladder.py`). Geant4 electrons, drift and diffusion
from the digitizer's own gas table, binned on 0.78 mm strips; strip time = median electron arrival time. No
avalanche, induction, resistive sharing, shaper or sampling. Threshold 7 % of the brightest strip (4 % and 12 % move
C by ±0.01).

Median C, time-on-position fit:

| | 5° | 10° | 15° | 20° | 30° |
|---|---|---|---|---|---|
| data X (M3) | 0.615 | **0.851** | 0.928 | 0.95 (n = 22) | – |
| data Y (M3) | 0.547 | **0.824** | 0.928 | 0.98 (n = 12) | – |
| ideal, no readout (X / Y) | 0.601 / 0.604 | **0.853 / 0.850** | 0.928 / 0.927 | 0.959 / 0.957 | 0.983 / 0.983 |
| sim v3 (X / Y) | | **0.808 / 0.795** | | 0.941 / 0.940 | |
| v3 + coupling 0.05, 15 ns | | 0.771 / 0.768 | | 0.930 / 0.928 | |
| v3 + coupling 0.06, 20 ns | | 0.760 / 0.759 | | 0.924 / 0.923 | |
| v3 at ρ_s 2 MΩ/sq | | 0.808 / 0.762 | | 0.941 / 0.932 | |
| ρ_s 2M + coupling 0.05, 15 ns | | 0.769 / 0.742 | | 0.927 / 0.921 | |

(The position-on-time fit gives the same ordering, all values ~0.03–0.05 higher. Data at 20° is few events
because inclined tracks mostly saturate a strip; 30° has none.)

**What it says:**

1. **Most of the ladder compression is not the readout.** With no readout, the ladder already reads 0.60 at 5°,
   0.85 at 10° and 0.93 at 15°: each strip averages a slice of depth, diffusion moves charge between strips with no
   time correlation (a textbook errors-in-variables dilution of the slope), and the threshold trims the ends.
   The data sit on that curve: the readout adds nothing in X and −0.03 in Y at 10°.
   RECONSTRUCTION_BASIS.md is right that hit times compress the ladder, and right to use the waveforms. Its
   attribution of the compression to the neighbours' delayed copies is mostly wrong for this estimator: the
   neighbours are a few-percent effect on top of the geometry. (Not tested here: the hits chain's own
   per-strip depth assignment with its own t0 and v, the "own t0 and v" row of its table.)
2. **The sim's readout compresses the ladder about twice as much as the data's** (−0.045 / −0.055 vs 0 / −0.03 at 10°).
   That is the sim's late neighbour (+105..+135 ns vs +45/+60) seen through a different observable.
3. **The coupling candidate moves the ladder the wrong way.** (0.05, 15 ns) takes another 0.03–0.04 off C in
   both views; ρ_s 2M takes it off in Y. The data's neighbour share arrives *without* dragging the strips' times
   together; the coupling's prompt copy of the central drags them together. This is an independent vote against
   adopting the coupling as a fitted term (Dylan's open question, HANDOFF 10-18), and for a lateral model that
   makes the neighbour signal earlier without making it a copy of the central (FRESH_LOOK §21.11's reading).
4. Y is more compressed than X in the data (0.824 vs 0.851) and the ideal ladder is symmetric, so the Y readout
   mixes more. That matches wft's kY > 1 on data (r06: 2.9).

**Mechanism checks** (archived August chain, ρ_s 2M, 30 mm, rewired; ideal at 30 mm: 0.872 / 0.965 at 10° / 20°):

| X ladder C | 10° | 20° |
|---|---|---|
| with ions | 0.844 | 0.962 |
| no ions | 0.843 | 0.957 |
| ions 2.45× faster | 0.841 | 0.955 |

- **The ion tail plays no part.** Removing it or speeding it up 2.45× moves C by ≤ 0.007.
- The August chain's readout share is −0.028 against v3's −0.045. Photon feedback (new in v3) is prompt and
  co-located in the model, so it does not change pulse shapes. It does make v3 1.8× brighter, matching the data,
  so more shared-only neighbour strips pass the absolute 100 ADC / 5σ threshold (likely, not tested). The readout share
  would then depend on amplitude through the threshold, and the amplitude-matched v3 is the fair comparison.

## 1. Angle closure: wft calibrated on the sim, scored against the gun

**Method** (`reco_truth/sim_calib.py`, `closure.py`). The data procedure, with the true line in place of M3:
`wft.calibrate`'s ref-pinned fit (cache of waveform windows along the reference corridor, impulse template from bright
inclined strips, kernel hypers fitted with each event's p0 and w pinned to the reference), seeded from det3's
production bundle `calib_bundle_r06` (RC-ladder kernel, c2 = 0.6 c1), v pinned to the sim's own 38.02 µm/ns as
production pins data to 36.6. Training events are spread evenly over the 12 angle points. Then `t13_reco.py`
(seeds from the sim's own waveforms, no truth) on every point, scored against the gun angle and the event's vertex
line. Frame check (`frame_check.py`): the right tan sign beats the flipped one by 4–9× in χ², and vertical views put
the χ² minimum at zero position shift. Inclined views prefer the mesh point moved along the track by ~1.4 mm of depth,
the size of a t0 / depth-origin offset. Rechecked with the sim bundle.

Two arms:
- **sim-calibrated** bundle: is wft unbiased when its calibration is right? That is the method's own closure.
- **data bundle r06 as is**: what the data's calibration does to sim tracks (kernel mismatch plus v 36.6 vs 38.02).

**Pre-registered (2026-10-18, before any closure number):**
- sim-calibrated, inclined 10°: |tan_reco / tan_true − 1| ≤ 1 % in both views; 15–30°: ≤ 2 %. Vertical: |Δθ| ≤ 0.15°.
- σθ at 10° (robust) 0.9–1.3°, the "physics floor" of WAVEFORM_FIRST_THREADING §12.
- Angle constants truth would assign: kw within ±2 % of 1, |w0| ≤ 1 µm/ns (data r06 from M3: kw 1.001 / 0.990,
  w0 −0.04 / −0.17).
- r06 on sim: tan ratio 1.039 (the v ratio) ± 0.03 from the kernel mismatch.
- **Falsifier:** sim-calibrated bias > 3 % at 10° in either view → wft has an intrinsic angle-scale bias that, on data,
  only M3's kw hides. That would transfer to every beam angle (stage 3's k) that is not M3-calibrated.

### Result, sim-calibrated arm (v3 baseline, 1200 events per point; th00 3000)

The truth-pinned calibration on the sim converged to c1 0.050 (at its floor, as on data), τ_s 185 ns, σ_s 71 ns,
**σ_p0 0.034 mm** (at its 0.03 floor) and **Dp 0.0257**. Data r06: σ_p0 0.42 mm, Dp 0.0144. (Below, under "sim".)

Scored against the gun (`closure_table.py`). Flip = the fitted direction reversed (tan sign wrong, p0 moved ~5 mm:
the mirror solution). Core = |tan_reco/tan_true − 1| ≤ 0.25. Residual = median |reco line − true line| at depth z, core events.

| point | view | flip | core tan ratio | core σθ [°] | \|r\| at 3 / 15 / 26 mm |
|---|---|---|---|---|---|
| vertical | X / Y | – | – | 1.10 / 1.26 (Δθ +0.01 / −0.02) | 0.15 / 0.10 / 0.16 ; 0.13 / 0.11 / 0.22 |
| 5° | X / Y | 0.021 / 0.012 | 1.021 / 1.059 | 0.68 / 0.67 | 0.19 / 0.19 / 0.24 |
| 10° | X / Y | 0.010 / 0.004 | **1.045 / 1.066** | 0.88 / 0.80 | 0.23 / 0.26 / 0.34 ; 0.24 / 0.31 / 0.41 |
| 15° | X / Y | 0.018 / 0.020 | 1.009 / 1.028 | 1.05 / 0.88 | 0.34 / 0.38 / 0.46 |
| 20° | X / Y | **0.147 / 0.150** | 0.975 / 1.009 | 1.17 / 1.15 | 0.57 / 0.65 / 0.74 |
| 30° | X / Y | **0.323 / 0.342** | 0.944 / 0.984 | 3.0 / 2.1 | 1.6–2.0 |
| 10° × 10° | X / Y | 0.013 / 0.012 | 1.049 / 1.055 | 0.97 / 0.85 | 0.22 / 0.28 / 0.38 |
| untilted view of an inclined point | | – | – | 0.93–1.25, Δθ ≤ 0.08 | 0.08–0.15 / 0.09–0.11 / 0.14–0.23 |

The same measurement on **det3 data against M3** (`data_flip.py`, production r06 with its M3 angle constants):

| \|tan\| (≈ angle) | 0.06–0.12 (3–7°) | 0.12–0.22 (7–12°) | 0.22–0.32 (12–18°) | 0.32–0.45 (18–24°) |
|---|---|---|---|---|
| flip X / Y | 0.039 / 0.018 | 0.017 / 0.013 | 0.023 / 0.019 | 0.071 / 0.026 (n = 85 / 39) |
| core ratio X / Y | 0.986 / 1.004 | **1.001 / 0.999** | 0.999 / 0.997 | 0.991 / 0.993 |

**Pre-registration:** vertical |Δθ| ≤ 0.15° and σθ 0.9–1.3° hold. **The falsifier fired:** the 10° bias is
+4.5 % / +6.6 % (≤ 1 % predicted, > 3 % falsifies), and the scale is not even linear (+6 % at 10° to −2.5 % at 20°).

**What it means.** The falsifier was written to catch an intrinsic wft bias that M3's kw hides on data. The data
rule that out by themselves: against M3, wft on det3 is linear to ±1 % from 7° to 24° and flips ~2 % of inclined
tracks, while on the sim it swings 7 % between 10° and 20° and flips 15 % at 20°. A bias intrinsic to the method
would show on data too, and it does not. **The data fit wft's forward model better than the sim does.** The sim's
waveforms are what break the method's assumptions, and the obvious suspect is the same one as everywhere else: the
sim's late neighbour (§3, FRESH_LOOK §21). So:

- No angle correction should be carried from the sim to data. The production M3 constants stand.
- Closure becomes a new **yardstick for the sim's lateral model**: a better lateral model must bring the sim's flip
  fraction and linearity towards the data's. Used below.
- 30° (|tan| 0.58) is outside anything det3 saw on the bench (no data events above |tan| 0.45); at 20° the sim flips
  twice as often as the data do (0.15 vs 0.07 ± 0.03).

### The other arm: the data's bundle (r06) on the sim

| | 5° | 10° | 15° | 20° | 30° |
|---|---|---|---|---|---|
| core tan ratio X / Y | 1.081 / 1.052 | **1.082 / 1.074** | 1.058 / 1.061 | 1.031 / 1.056 | 0.978 / 1.014 |
| flip X / Y | 0.022 / 0.014 | 0.011 / 0.007 | 0.009 / 0.009 | **0.027 / 0.019** | 0.117 / 0.103 |

Vertical: Δθ −0.04 / +0.27°, σθ 1.43 / 1.66° (sim-calibrated: 1.10 / 1.26).
- Pre-registered 1.039 ± 0.03 (the v ratio 38.02/36.6): **missed** at 10° (1.08). The kernel mismatch adds +4 %.
- The data's bundle flips **fewer** sim tracks at 20° (2–3 %) than the sim's own calibration does (15 %). That bundle
  has the data's kernel (σ_p0 0.42 mm, slower Dp), so the flips belong to the sim-calibrated kernel. A wft model
  pinned to the sim's narrow-at-the-mesh cloud is the one that loses the direction.

### Refinement: flips belong to the bundle, the nonlinearity to the sim's waveforms

| core tan ratio X / Y, flip at 20° | 10° | 20° | swing 10→20° | flip 20° |
|---|---|---|---|---|
| v3 sim, sim-calibrated | 1.045 / 1.066 | 0.975 / 1.009 | 0.070 / 0.057 | 0.15 / 0.15 |
| v3 sim, data bundle r06 | 1.082 / 1.074 | 1.031 / 1.056 | 0.051 / 0.018 | 0.03 / 0.02 |
| wet + σ0 420 µm sim, sim-calibrated | 1.039 / 1.040 | 0.981 / 1.008 | 0.058 / 0.032 | 0.15 / 0.14 |
| **det3 data vs M3 (r06)** | 1.001 / 0.999 (7–12°) | 0.991 / 0.993 (18–24°) | **0.010 / 0.006** | 0.07 / 0.03 |

- **The flips are a property of the sim-calibrated bundle.** The data's bundle on the same sim flips 2–3 % at 20°, like
  the data; the sim's own truth-pinned calibration (c1 at its floor, σ_s 71–351 ns) flips 15 %. That calibration is
  the weak link, not the waveforms.
- **The nonlinearity is a property of the sim's waveforms.** It appears with either bundle (X swings 5–7 % between 10°
  and 20°; the data swing 1 %), and fixing the neighbour delay (wet + σ0 420 µm) does not remove it.
- Next suspect: the slow pulse. The sim's inclined rise is +12 ns over the data's, with a rigid ~70 ns offset found in
  August. **No-ions diagnostic submitted** (condor 4480264, pre-registered in `stageBC_w3_noions.sub`: swing < 0.02 if
  the ion tail drives it, 0.05–0.08 if not).

## 2. Depth-resolved residuals

The data's |wft line − M3 line| (`data_line_res.py`, same accounting as the sim column; offsets centred):

| depth [mm] | 3 | 9 | 15 | 21 | 26 |
|---|---|---|---|---|---|
| data vs M3, all X / Y | 0.42 / 0.39 | 0.43 / 0.39 | 0.46 / 0.42 | 0.50 / 0.48 | 0.55 / 0.52 |
| data vs M3, head-on X | 0.35 | 0.37 | 0.41 | 0.47 | 0.53 |
| data vs M3, inclined X (\|tan\| 0.14–0.45) | 0.48 | 0.50 | 0.51 | 0.54 | 0.58 |
| sim vs truth, vertical X / Y | 0.15 / 0.13 | 0.11 / 0.10 | 0.10 / 0.11 | 0.13 / 0.17 | 0.16 / 0.22 |
| sim vs truth, 10° X | 0.23 | 0.23 | 0.26 | 0.30 | 0.34 |
| sim vs truth, 15° X | 0.34 | 0.35 | 0.38 | 0.42 | 0.46 |

(The data reproduce RECONSTRUCTION_BASIS's waveform-first row, 0.39 → 0.55, with a different estimator.)

- **On head-on tracks the data's residual is almost all M3.** wft's own error against truth is 0.10–0.22 mm at every
  depth; the data read 0.35–0.53. In quadrature that leaves M3 a pointing error of ~0.32 mm at 3 mm growing to ~0.5 mm
  at 26 mm, i.e. M3's own angle error times the lever arm. **The data table's growth with depth is mostly M3's,
  not the reconstruction's.** That is where the 0.4 mm "M3 pointing floor" came from.
- On inclined tracks the sim's own error is 0.25–0.45 mm and grows with depth (10°–15°), because the angle error
  (σθ 0.9–1.05°) multiplies the depth. For these the sim is the pessimistic case: its angles are worse than the data's (above).
- The "physics floor" σθ ≈ 1.0–1.1° (WAVEFORM_FIRST_THREADING §12, from a toy) is confirmed on the full
  simulation for vertical tracks (1.10 / 1.26°) and is 0.7–1.0° on 5–15° tracks. The data's 1.15 / 1.18° against M3
  includes M3's angle error, so wft on data is at or below the floor.

## Sim: wet gas + a prompt lateral spread at the anode (condor 4480262)

**Why this test.** Three leads point at the same place.
- FRESH_LOOK §12: the bench drift velocity says ~1 % water, but the bench widths look dry
  (D_T 463–507 µm/√cm, the dry Magboltz value). The sim runs the dry table, which is why it drifts at
  38.0 µm/ns instead of the data's 36.6.
- wft's own calibration (item 1): on data the charge cloud is already σ_p0 = 0.42 mm wide at the mesh and grows
  slowly (Dp 0.0144). On the dry v3 sim it starts at 0.03 mm and grows fast (Dp 0.0257).
- FRESH_LOOK §21.11: the data's neighbour carries its share *earlier*, a prompt in-phase share the sim lacks.

**Hypothesis H.** The gas is wet as v says, and the anode has a prompt, charge-conserving lateral spread of
~0.4 mm. Together they mimic dry widths on average while moving charge from late diffusion into an in-phase
neighbour share. Tested with the 1 % H₂O Magboltz table (`design/gas/drift_velocity_Ar_iC4H10_95_5_H2O_1.0.json`,
v 37.2 µm/ns, σ_T over 26.5 mm 506 µm vs dry 785) and Stage B's new `--sigma0-um` (the avalanche spread; the calib
has 94 µm). Arms on the v3 recipe, det3 wiring: wet only (v3w), wet + σ0 300 µm, wet + σ0 420 µm.
Pre-registration in `scratchpad/2026-10-18/stageBC_w3_wetspread.sub`.

**Neighbour shape** (vertical, M3-selected data; `nb_compare.py`, unscaled):

| | ±1 delay X / Y [ns] | ±1 / central X / Y | nb(t = 0) X / Y | nb rms X / Y |
|---|---|---|---|---|
| data | +45 / +60 | 0.368 / 0.408 | 0.219 / 0.225 | – |
| v3 (dry, 94 µm) | +135 / +135 | 0.372 / 0.435 | 0.172 / 0.185 | 0.028 / 0.039 |
| v3 + coupling 0.05, 15 ns (§21.8) | +105 / +105 | 0.361 / 0.422 | 0.226 / 0.236 | 0.018 / 0.032 |
| wet, 94 µm | +120 / +150 | 0.322 / 0.360 | 0.154 / 0.159 | 0.040 / 0.042 |
| wet, σ0 300 µm | +90 / +105 | 0.360 / 0.417 | 0.183 / 0.199 | 0.022 / 0.028 |
| **wet, σ0 420 µm** | **+60 / +60** | 0.382 / 0.461 | 0.202 / 0.242 | **0.018 / 0.038** |

**Widths** (`xy_extent_m3.py`, unsaturated; RMS in mm, n strips, span in mm):

| | vertical X: rms / n / span | vertical Y rms | 10° X rms | 20° X rms |
|---|---|---|---|---|
| data | 0.66 / 4.0 / 2.34 | 0.84 | 1.46 | 2.75 |
| v3 (rewired) | 0.68 / 4.5 / 3.12 | 0.81 | 1.46 | 2.72 |
| wet, 94 µm | 0.52 / 3.0 / 1.56 | 0.72 | 1.40 | 2.71 |
| wet, σ0 300 µm | 0.58 / 4.0 / 2.34 | 0.75 | 1.42 | 2.71 |
| **wet, σ0 420 µm** | **0.65 / 4.0 / 2.34** | 0.79 | 1.45 | 2.72 |

**Amplitude ledger** (data/sim p50, vertical X / Y): v3 0.978 / 1.012; wet 94 µm 0.840 / 0.857 (wet alone makes the
sim 17 % brighter: less diffusion concentrates the charge); wet σ0 420 µm 1.048 / 1.094. Inclined Y (10°): v3 1.128,
σ0 420 µm 1.193; its saturated fraction 0.193 vs data 0.331 (v3 0.218).

**Against the pre-registration:**
1. v 37.2: ✔ (Stage B log).
2. Vertical X width: wet alone 0.52 (band 0.55–0.62, a near miss, narrower than predicted); σ0 420 µm 0.65 (band 0.63–0.72) ✔.
3. 10° extent within ±4 %: ✔ in every arm (1.40–1.45 vs 1.46).
4. nb(0) X ≥ 0.19 with σ0 420 µm: ✔ 0.202. ±1 delay X ≤ +120 ns: ✔ **+60** (and Y +60 = data).
7. Amplitude +3..+12 % (data/sim): ✔ +7 % X, +8 % Y.
5. wft's truth-pinned calibration on the σ0 420 µm sim returns σ_p0 0.25–0.55 mm, Dp 0.008–0.014: **✘ failed.**
   It returns σ_p0 0.031 mm and Dp 0.0256, and puts the extra width into the kernel's time smear (σ_s 351 ns) and
   kY. Wet alone: σ_p0 0.038, Dp 0.0181 (the wet Magboltz value is ~0.010). **wft's σ_p0 and Dp are not clean
   physical readouts**: they trade off against its sharing-kernel parameters. So the data bundle's 0.42 mm is no
   evidence for a prompt spread by itself; that leg of the motivation falls. The neighbour shapes above stand on their own.
6. Ladder C, X at 10° (band 0.80–0.86 for σ0 420 µm): ✔ **0.817**. Wet alone 0.873, σ0 300 µm 0.849.
   The ideal ladder must use the same gas. Ideal wet + 0.42 mm is **identical** to ideal dry with none (0.591 / 0.849 /
   0.926 / 0.958 vs 0.601 / 0.853 / 0.928 / 0.959 at 5/10/15/20°): the ladder only sees the total width, the "mimics
   dry" statement in another observable. Against its own ideal, the sim's readout share at 10° goes from −0.045 (v3)
   to −0.032 (σ0 420 µm); the data's is 0.00 (X) / −0.03 (Y).

**Falsifier not triggered.** One pre-registered band failed (item 5, the wft parameters).

**Reading.**
- **The ±1 neighbour delay, the sim's longest-standing mismatch (×2 since §10.4), closes with a prompt anode spread
  and wet gas.** Neither ρ_s (§21.11) nor the inter-strip coupling (§21.8) did: the coupling took the delay to +105
  and spoiled the ladder (§3). In X every neighbour number now matches the data to 0.01–0.02 and 15 ns, and the
  vertical X cluster has the data's strip count and span, which v3 never had.
- The "widths look dry" puzzle has a reading: wet diffusion plus ~0.4 mm of prompt spread gives the dry total.
  The water is real (v, FRESH_LOOK §12). The prompt spread is a phenomenological term with no mechanism in the
  model yet. Candidates: lateral displacement of the photon-feedback secondaries, which the model places exactly
  on the primary (`add_feedback`, "land where the primary did"). That would also explain why cloud_basics finds the
  data's early width growing with mesh voltage (σ0 0.3 → 0.5 mm over 425 → 505 V), since bG grows with gain.
  The other candidate is a wider avalanche than the 94 µm Garfield calib.
- **Y is not fixed.** Y's ±1/central overshoots (0.461 vs 0.408, v3 already 0.435), the vertical Y width is a little
  narrow (0.79 vs 0.84), and inclined Y gets dimmer (saturates 0.19 vs 0.33). One knob cannot fix Y's core and wings
  together. ρ_s acts only on Y (the stripes run along y), and on both in the same direction. The Y open problem
  (§21.9) stands, with this baseline as the better starting point.
- Amplitude: the σ0 420 µm arm wants bG ≈ 0.48–0.49 instead of 0.45 (≈ 2 % per 0.01, §21.10).
