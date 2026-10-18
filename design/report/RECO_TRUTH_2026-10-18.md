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
| 1. angle closure | sim calibration written and frame-checked; runs once the truth points are in |
| 2. depth residuals | after item 1 |
| 3. hit-time ladder | **first answer in** (below) |

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
