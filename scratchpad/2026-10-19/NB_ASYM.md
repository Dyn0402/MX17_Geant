# X left/right neighbour asymmetry and the neighbour timing (2026-10-19)

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
- Sim: try reducing the late neighbour charge (not adding prompt) at ±1: diffusion (wet), ρ_s, stripe line/gap. Add a
  symmetric prompt exponential tail (λ ≈ 0.4–0.5 per strip, ~1 % at ±3) as a readout term and check the ladder.
- Ask whether det3's X cabling/adapter moved with it between slots (6-22/23 vs 6-27/28).
