# det3 X asymmetry: follow-up to the audit (2026-10-21 evening, dylan-MS-7C84)

Continues `audit/AUDIT.md` §4.3. Scripts and JSON in the subdirectories named below; data = raw survey (NB_RAW=~/x17/nb_asym_raw),
sim arms copied from yoga `~/x17/nb_asym_raw/{simgapleak,simfunnel,simxoff}`.

## 1. The even term: the exclusion does not hold (corrected the same evening; `even/`, `loop/`)
**Correction.** The first version of this section (and audit §3.1) said the even term retires the attractor family. It does not:
both compared the sharing at large |a1| with the sharing at delta = 0 *inside the same attractor*, not with an unpulled board.
- **Toy loop** (`loop/toy_loop.py`, the audit toy with a no-pull baseline over a full period of delta): a1(delta) is zero at
  delta = 0 (attractor over the pad centre: cluster narrowest, Δs1 −0.1 … −0.5) and at P/2 (over the boundary: widest,
  Δs1 ≈ 0 … +0.2). Where |a1| ≥ 0.1, Δs1 vs no pull spans negative to ~0, e.g. σ 450 µm, f 0.2: a1 up to 0.29 with Δs1
  −0.23 … −0.02. So "sharing unchanged at a1 0.14" (same PCB, below) is allowed for a pull at offsets between the a1 maximum
  and P/2, and the det2 slope of s1 vs |a1| depends on which zero crossings its cells sit at.
- **Sim arms** (`even/sim_even.py` → `sim_even.json`, 25 arms): E/O = 1.2–2.5 is the within-attractor contrast
  s1(150–350 µm) − s1(< 100 µm), the same kind of comparison as the det2 cell regression. Not an exclusion by itself.
- **Same PCB** (`even/same_pcb.py`): det3 May vs June-15 X sharing at matched Y share equal to ±0.02 over 7 bins while a1 goes
  0 → 0.14. Consistent with no pull, or with a pull sampled around the a1 maximum (board average of a loop ≈ 0).
- **Loop signature in the data** (`loop/det2_crossings_np.py`, residual s1X − F(s1Y) per band, F = board's own binned
  relation): an attractor whose offset sweeps across a board predicts a1 ∝ sin φ and Δs1 ∝ −cos φ (quadrature).
  - det3 June-15 by x: a1 +0.09 → +0.21 → +0.09 while the residual goes +0.08 → 0.01 (at the a1 maximum) → −0.11.
  - det2 by y: residual +0.08 at its rising a1 zero crossing, ≈ 0 near the extrema, −0.04 … −0.08 at the falling crossing.
  - det4 by x: +0.08 / +0.05 / −0.19 with a1 0.03 / 0.11 / 0.03.
  - **But the quiet same-PCB control has residual structure of the same size**: det3 May by x −0.15 / +0.13 in adjacent bands
    with a1 ≈ 0 (bulk-local, not PCB or electronics: June-15 doesn't repeat it). So the loop is suggestive, not evidence.
- **Shape of the transfer** (`even/det2_profile.py`, `even/antisym.py`): per unit a1 on det2, −0.54 / +0.47 at ∓1, −0.08 / +0.04
  at ∓2, ~0 at ∓3; det3 June-15 the same shape. A fixed linear channel operator fails: a nearest-neighbour antisymmetric
  coupling puts ~40 % at ±2, a copy gives n−1 no loss. Position-level mechanisms (where the charge lands) are favoured; an
  attractor is one of them.
- **Reading:** the even term doesn't discriminate at the current level because bulk-local sharing structure (±0.1) is as large
  as the attractor's loop. Registration is neither excluded nor supported by it. The remaining problem for the strip
  attractor is still its size in the full sim (×2–5 short) and its physics (the across-strip pull must be fast without fast
  along-strip spreading, which the Y view forbids); an electrostatic pull (charged gaps focusing avalanche electrons onto the
  strips: prompt, X-only, gain-independent once saturated) is a candidate the sim hasn't tried in that form.

## 2. Time, gain and the session jumps (`time/time_a1.py` → `time_a1_det3.json`)
- On the direct a1 (not the copy estimator, which misreads low-amplitude runs: the 6-23 early runs "0.035/0.057" are
  +0.06 ± 0.05 / +0.16 ± 0.02 on a1, consistent with the night's 0.19), every det3 June-15 session falls on one curve in the
  median Y share: 0.19–0.20 at 0.5, 0.16 at 0.56–0.65, 0.14–0.15 at 0.75, 0.07–0.11 at 1.0.
- Flat with gain at fixed width: 6-23 resist 525 → 465 V (pulse ×6 smaller), a1 0.15–0.21 throughout; 6-22 likewise.
- Flat over hours: 6-22's 8 h scan, 6-27 afternoon → overnight. No charging transient at the survey's resolution (≤ 2 files
  per sub-run; within-run splits agree).
- Charging-up would scale with deposited charge; the data scale with cloud width only.

## 3. 6-19/20 was det2, mislabelled (`s619/`)
- The "det3 FEU 7" 6-19 data (a1 0.03, all 8 chips lit) has det2's readout-chain noise pattern (r = 0.957 with det2 X on
  6-22) and not det3's on the same FEU 7 slot on 6-17 (`fingerprint619.py`). The pattern follows the chain + board, not the
  FEU slot (`fp_validate.py`: same board, other FEU, r ≈ 0; different boards median 0.04).
- Board properties agree (`board_props.py`): FEU 7 lights all 8 chips and reads a1 +0.026 like det2; det3's X readout on
  6-17 lit only chips 0–3, and on 6-19 that 4-chip readout is the FEU 3 file (372 events, a1 +0.077 ± ~0.03, s1 0.75).
- So the one falsifier of "fixed map × cloud width" is a labelling error. Real det3 on 6-19 is low-statistics and ~2σ low.
  The raw_daq_data HV scan + long run could still be decoded for more FEU 3 events.

## 4. A strip pull reproduces the transfer shape at the data's operating point (`loop/toy_shape.py`)
- In the attractor toy, the points where a1 = 0.10–0.20 with |Δs1| ≤ 0.03 vs no pull (the det3 June-15 vs May situation) have
  ±2 odd / ±1 odd = 0.01–0.14 (σ 450 µm: 0.10–0.14) and a small negative ±2 even part (−0.03 … −0.07 per a1). det2: 0.12 and −0.07.
- Representative point (30 % of the charge pulled onto lines 300 µm off the pad centres, σ 450 µm, 200 µm induction): per unit a1
  −0.09 / −0.49 / +0.51 / +0.01 at k = −2 / −1 / +1 / +2, against det2's −0.08 / −0.54 / +0.47 / +0.04; a1 0.16, Δs1 +0.003.
- **Caveat (found next, §5):** this toy pulls a fraction of the WHOLE cloud to the line nearest the cloud centre, a cloud-level
  (non-local) pull that a linear detector can't make. The physical, per-electron version is far weaker unless the electron
  cloud is narrow; see §5. The full resistive-spreading sim's effective pull is much weaker; the
  open question is the physics of a prompt pull that strong (electrostatic focusing by charged gaps is the candidate).

## 5. What a physical pull needs (`loop/toy_local.py`, `toy_local_narrow.py`, `toy_funnel*.py`)
- A linear detector responds electron by electron, so any pull must be per-electron (each electron moved by where it lands).
  Per-electron pulls with σ 300–450 µm Gaussian clouds and 100–200 µm induction give a1 ≤ 0.04–0.07 at |Δs1| ≤ 0.03.
- **Electrostatics of charged gaps:** insulating gaps saturate at σ = ε0E0 ≈ 3.5 nC/cm² (39 kV/cm), after which no field line
  ends on a gap and every gap electron lands on the nearer strip edge. With 540–700 µm strips this funnel alone gives a1 ≈ 0 (each
  gap feeds both its edges); with spreading across a wide strip it widens the cluster (Δs1 +0.2–0.4).
- **The regime that works:** the observed sharing fixes only cloud width ⊗ induction. With a narrow electron cloud
  (σ 100–200 µm) and wide induction (300–350 µm; no-pull s1 0.5–0.67, in the data's range), the saturated funnel onto narrow
  strips gives at |Δs1| ≤ 0.03: 100 µm strips a1 0.13–0.24, 200 µm 0.11–0.22, 300 µm 0.07–0.18, 400 µm 0.03–0.10 (≈ 0 at 540).
  Per-electron pulls onto narrow lines behave the same (f 0.3–0.6: a1 0.1–0.2).
- Why a narrow cloud is plausible: a vertical track's electrons are a depth mixture (30 mm gap, each primary diffused by its own
  depth), whose pitch-scale Fourier content is ~10× that of a Gaussian of the same rms; attachment in the ~1 % H2O gas at low
  drift field removes deep electrons, which is why the Y share narrows at low drift field (0.98 → 0.56, 1100 → 100 V) as a1 rises.
- **Prediction:** the June prints have narrow strips (≲ 300 µm) and wide gaps at a 780 µm pitch; the May prints wider strips.
  A photo decides it. Mechanism is prompt, X-only, gain-independent (saturated), the same over both pad types, bigger for
  narrow clouds.
- **Full-sim test running:** X kernels with 190 / 290 µm strips (condor 4526243), then paired arms w190, w190f1, w290, w290f1
  (`narrow_chain.sh`; funnel width matched with the new `--gap-funnel-width-um`). Note the earlier w690f1 arms funnelled onto
  550 µm edges (the funnel used the production width regardless of the kernel's); fixed by the new option.

## 6. Full-sim narrow-strip test: fails (`paired_odd_narrow.json`, 2026-10-21 19:50)
- Paired arms vs glg0 (36 angle/seed pairs each, condor 4526326; hybrid X kernels 190 / 290 µm strips, 5 MΩ/sq; funnel
  width matched). Odd amplitude O (all / narrow half), ratio to glg0, and X sharing s1 (glg0 0.760):

  | arm      | O all  | ratio | O narrow | ratio | s1    |
  |----------|--------|-------|----------|-------|-------|
  | glg0     | +0.053 | 1     | +0.068   | 1     | 0.760 |
  | w190     | +0.011 | 0.21  | +0.022   | 0.33  | 0.715 |
  | w190f1   | +0.031 | 0.59  | +0.043   | 0.63  | 0.664 |
  | w290     | +0.028 | 0.52  | +0.035   | 0.51  | 0.728 |
  | w290f1   | +0.049 | 0.92  | +0.062   | 0.92  | 0.705 |

  (paired diff errors ±0.004–0.005.) Pass needed ×2.5–5 with no sharing change; every arm is BELOW the 550 µm control,
  and the funnel arms also narrow the sharing (Δs1 −0.05 … −0.10).
- In the full sim the odd amplitude grows with strip width (190 < 290 < 550; funnel adds ~+0.02–0.025 at each width), the
  opposite of the 1-D toy (§5). The sim carries the real depth mixture of a track's electrons (Stage A clusters), so the
  §5 "narrow cloud" regime is not what the sim's clouds look like at the pitch scale, or the toy's induction is not the
  resistive kernel's. Either way the narrow-print prediction is not supported; a strip photo is no longer the decisive
  next step for this branch.
- Y view and the even term were not scored (the odd amplitude already fails).

## 7. Skewed induction: the resistive version is already simulated; ion return paths closed (`skew/`, 2026-10-21 ~21:00)
- **The resistive-strip version of route 9 is not a new test.** The W2 X box is the 31.2 mm superperiod (`kernels.pad_x`
  folds into `C.SUPERPERIOD_M`): 39 strips at 800 µm over 40 pads at 780 µm, so every production and hybrid kernel already
  carries offset resistive strips at every δ with the full time-dependent (Riegler) weighting physics. glg0's O = +0.053 IS
  the skewed-induction + attraction answer for 550 µm strips at 5 MΩ/sq (per-δ: `r5/kernel_delta*.json`).
- **What the sim cannot do, tested in a 2-D toy (`skew/ion_path_toy.py`):** with a saturated funnel the avalanche ions
  drift back along their electron's field line and reach the mesh above where it came from, while the digitizer sends
  them straight up from the strip edge with one shared ion time profile. Integrated charge is path-independent (mesh
  ψ = 0); a1 is an in-phase projection, so the path could matter. Exact layered weighting potentials over the
  superperiod, the sim's DreamShaper, the data's α projection, 20 k clouds per point (O error ≈ ±0.002).
  - ESL transparent: |O| ≤ 0.005 in every case (no funnel, funnel with straight ions, funnel with return paths),
    σ 100–450 µm. (The 200-event test's ±0.05 was noise.)
  - ESL strips floating (`skew/float_psi.py`, zero net charge per strip): the return paths reduce the funnel's |O| by
    5–10 % (σ 100: −0.221 → −0.204; σ 300: −0.106 → −0.094). **Closed: ion paths do not add a skew; they slightly undo
    the funnel's attraction.**
- **The floating-strip limit itself** (`ion_path_toy_float.json`), no funnel / funnel:

  | σ (µm) | O | O funnel | ±2/±1 | s1 |
  |---|---|---|---|---|
  | 100 | −0.160 | −0.221 | 0.00 / 0.02 | 0.37 |
  | 200 | −0.085 | −0.148 | 0.05 / 0.09 | 0.49 |
  | 300 | −0.062 | −0.106 | 0.18 / 0.24 | 0.63 |
  | 450 | −0.025 | −0.052 | 0.76 / 0.71 | 0.88 |

  - It SHRINKS with cloud width (the data's direction), unlike route 9's 3-D `bus_float` + `fwd_parity` note ("growing
    with σ"); the difference is unexplained (2-D full-column X pads here vs the 3-D comb there, and a different estimator).
  - ±2/±1 reaches the data's 1/8–1/4 at σ 250–330 µm.
  - But at the full sim's sharing (s1 ≈ 0.76, σ ≈ 370 µm) O ≈ 0.04, the full sim's own 0.053: **the full sim is already at
    the conducting-strip limit.** At the data's sharing no strip-induction variant reaches ×2.5–5; it would need a cloud
    narrower than the sim's with the same sharing, i.e. wide late spreading after the strip has fixed the centroid (the
    §5 narrow-cloud regime again, which §6 showed the sim does not have).
- Sign: the toy's O is negative in the +x / sim_delta convention; magnitudes are what compare (the data's δ is unknown).

## What this leaves
- The fingerprint (audit §4.1) gains: (2) the transfer is n−1 → n+1 direct, ±2 odd ≈ 12 % of ±1; (5) gain- and
  time-independent at fixed width; no session outlier. The even term is not a usable discriminator at the ±0.1 level.
- Excluded on top of the audit: a fixed linear channel operator (copy or antisymmetric coupling), charging-up.
- Open: what physical, position-level, sum-conserving mechanism moves charge from the n−1 pad to the n+1 pad, X-only (across the strips), with a smooth two-signed geometry per board, June bulks only, and an
  amplitude set by the cloud width.
- Strip-induction routes are now exhausted at the data's sharing (§6, §7): offset resistive strips are in the full sim and
  sit at the conducting-strip limit; ion return paths subtract; a lateral field from strips at different potentials would
  need ~kV between neighbours (E_z ≈ 33 V/µm) to tilt the ions by 0.1, implausible for cosmic-rate currents.
- What every strip route needs and the sim lacks: a cloud whose centroid the strip fixes while it is narrow, with the
  sharing widened afterwards (late spreading), so that s1 is in the data's range while a1 keeps the narrow-cloud size.
  Testable on archive data: compare the data's PROMPT sharing (early samples, before resistive spreading) with the sim's
  at the same in-phase s1. Narrower prompt sharing in the data than in the sim would put the data in that regime.
- Otherwise the remaining classes need the board or the detector: the strip photo (registration, width, gap) and
  the detector tests (collimated scan across a pad, HV ramps and dwell, dry gas).
