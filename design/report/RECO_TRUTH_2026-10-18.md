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
| Sim | testing "wet gas + prompt anode spread" (condor 4480262, pre-registered) |

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
