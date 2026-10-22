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

## 8. Prompt vs in-phase sharing, data vs sim (`prompt/prompt_share.py`, 2026-10-21 night)
- r(f) = (W₊₁ + W₋₁)/W_c at the central's 25/50/75 % rising crossings, its peak and 50 % falling (ratio of sums, bootstrap
  errors), in bins of in-phase s1. Data `d3jun/d2/d3may_x.npz`, sim glg0 (all 36 arms).
- **All events are dominated by x-inclination.** Data |t(n+1) − t(n−1)| median 134–199 ns vs the sim's 69 ns (the arms are
  th00 + thy only, no x tilt): inclined tracks give one neighbour early and the other late, so the data's neighbour sum is
  wider in time than the central at both ends. Any data/sim timing or sharing comparison must cut on it.
- Near-vertical (`SEL=vert DTN_MAX=80 FW_PCT=25`: |Δt_nb| < 80 ns, central FWHM above its 25th percentile; 380–790 data
  events per board; `prompt_share_vert_80_25.json`), r@25 % vs sim at the same s1:

  | s1 | sim | det3 May | det2 | det3 June-15 |
  |---|---|---|---|---|
  | 0.38 | 0.12 | 0.16 ± 0.02 | 0.28 ± 0.09 | 0.36 ± 0.02 |
  | 0.52 | 0.24 | 0.23 ± 0.02 | 0.23 ± 0.06 | 0.56 ± 0.03 |
  | 0.68 | 0.50 | 0.40 ± 0.03 | 0.34 ± 0.05 | 0.73 ± 0.05 |
  | 0.82 | 0.72 | 0.75 ± 0.08 | 0.58 ± 0.06 | 0.81 ± 0.07 |
  | 0.99 | 0.93 | 0.87 ± 0.06 | 0.62 ± 0.06 | 0.79 ± 0.10 |

  - det3 May (no asymmetry) matches the sim: the method is sound.
  - det2 has the strip routes' regime: prompt sharing narrower than the sim at s1 ≥ 0.68 (by 0.15–0.3), late wider (fall50
    1.12–1.36 vs 0.80–0.85).
  - det3 June-15, the board with the largest asymmetry, has the opposite: prompt sharing wider than the sim by 6–14σ at
    s1 ≤ 0.68, and nearly flat through the pulse (s1 0.37: 0.36 / 0.34 / 0.34 / 0.35 / 0.46 from 25 % rise to 50 % fall),
    with less late content than the sim. Its pulses are also the shortest (FWHM median 318 ns vs 400 det2/May, 598 sim).
    Not neighbour crosstalk on its FEU: injection nearest-neighbour correlation excess 0.08 (det3 May 0.11).
- **Verdict:** the narrow-prompt / late-spreading regime is not common to the asymmetric boards, so it cannot be what
  enables the asymmetry; strip routes that need it fail on det3 June. New fingerprint: the det3 June bulk spreads its X
  charge wide BEFORE the peak (same chamber as det3 May, different bulk), like a fast (low-ρ_s) X spreading; uniform low
  ρ_s failed the Y view (round 8), so this would have to be X-only.
- **No sim arm makes it** (`prompt_share_vert_arms.json`; SIMS=glax10,glax3,glax10w690f1,glr05g0,glr2g0,glw690f1,glg0.1):
  at s1 0.37 / 0.52 det3 June r@25 % = 0.36 / 0.56 vs 0.10–0.17 / 0.14–0.34 across anisotropic strips (σ_xx = 3, 10 σ_yy,
  with 690 µm + funnel), uniform 0.5 and 2 MΩ/sq, 690 µm + funnel and the 10 % gap leak. Every arm's r(f) rises through the
  pulse at low s1; det3 June's is flat.

## 9. det3 June's flat X sharing is response, not coupling; it looks like fast-equalising strips (`prompt/coupling_checks.py`)
Near-vertical cut (|Δt_nb| < 80 ns, FWHM ≥ its 25th pct), `coupling_checks.json`:
- **Rise profile, s1 0.2–0.45** (r = (W₊₁+W₋₁)/W_c at 5/10/15/25/50/75 % rise, peak, 75/50 % fall):
  det3 June X 0.36 / 0.34 / 0.33 / 0.31 / 0.30 / 0.30 / 0.31 / 0.35 / 0.41 (flat from the onset, little late growth);
  det3 May X 0.17 / 0.16 / 0.15 / 0.15 / 0.18 / 0.20 / 0.26 / 0.41 / 0.60; sim glg0 0.09 / 0.10 / … / 0.31 / 0.54 / 0.78.
  At equal in-phase sharing, det3 June's X sharing is all prompt; det3 May's (same chamber, other bulk) and the sim's grow late.
- **s1 floor:** det3 June X p1/p5 = 0.125/0.19, 6 % of events below 0.2 (det3 May 0.04/0.09, 12 %). A uniform pad-to-pad
  coupling is therefore ≤ ~0.06 per neighbour, not the ~0.15 the flat level would need: the flat sharing scales with the
  event, so it is detector response, not crosstalk (and the FEU noise/injection correlations were not elevated, §8).
- **Every det3 June session** with statistics (mesh 465–525 V, drift 100–1000 V, 6-22 / 6-23 / 6-27) has r10 ≈ r25
  ≈ 0.26–0.40 in the s1 0.2–0.45 bin: a stable property of the bulk. Y view: det3 June and det3 May alike (flat-ish early).
- **Interpretation (hypothesis):** a time-flat neighbour shape is what strips that equalise across their width before the
  first samples give: the floating-strip limit of §7 (instant, so r(f) flat by construction). That limit also shrinks with
  cloud width and has data-like ±2/±1, and its size tracks det3 June's per-session survey `a`:
  data (s1 p50 → a) ≈ (0.38 → 0.19–0.21), (0.6 → 0.12–0.15), (0.8 → 0.07–0.10); toy floating (0.37 → 0.16), (0.63 → 0.06),
  (0.88 → 0.03); floating + funnel (0.43 → 0.22), (0.55 → 0.15), (0.68 → 0.11), (0.91 → 0.05). (`a` is the survey's board
  estimator, not the toy's sgn(δ)-weighted O, and the toy is 2-D: suggestive, not a fit.) The sim's 5 MΩ/sq strips
  spread slowly (late-growing r(f), O ≈ 0.05), like det3 May. With audit fact 7 (June builds conduct across their gaps), a
  lower-resistivity June print is a consistent story.
- **The test:** a full-sim arm whose strips equalise across their width in ≪ 60 ns while keeping the Y view: an
  anisotropic sheet far beyond round 8's σ_xx = 10 σ_yy (glax10 still has rising r(f): r25 0.14 at s1 0.38), e.g. σ_xx
  = 100 σ_yy, with and without the funnel. Pass: flat r(f) like det3 June, O ×2.5–5 of glg0, Y view intact.

## 10. Fast-equalising strips fail in the full sim; a strip resistivity gradient does nothing; §9's toy argument was flawed (2026-10-22 ~01:00)
- **σ_xx = 100 σ_yy** (condor 4526340 kernel, 4526745 arms; `paired_odd_ax100.json`, `prompt/prompt_share_vert_ax100.json`):
  glax100 vs glg0 O ×1.04 (diff +0.002 ± 0.002), narrow ×0.99, s1 unchanged; glax100f1 ×1.41 / narrow ×1.71. Identical to
  glax10 (×1.06): across-width equalisation was already fast enough at ×10 and is not what limits the sim. r(f) still rises
  (s1 0.38: 0.14 at 25 % rise → 0.33 peak). Fails size and shape.
- **§9's interpretation was flawed:** the 2-D toys deposit every electron at t = 0, so a static (floating) response is flat in
  time by construction. In the full sim electrons arrive over the drift time, later ones from deeper and wider, and the
  sharing grows through the pulse whatever the strips do.
- **Where the sim's growth comes from, and what the data show instead** (near-vertical, s1 0.2–0.45; r at 5 % rise / 25 % /
  peak / 50 % fall): sim X 0.08 / 0.10 / 0.29 / 0.82, sim Y (th00) 0.05 / 0.10 / 0.25 / 0.89; det3 May X 0.17 / 0.15 / 0.26 /
  0.60, det3 May Y 0.23 / 0.20 / 0.25 / 0.63, det3 June Y 0.20 / 0.22 / 0.29 / 0.57, det3 June X 0.36 / 0.31 / 0.31 / 0.41.
  - General data/sim mismatch, both views, every board: the data's early sharing is flat to the onset, the sim's rises
    2–3× before the peak (the sim's early clouds too narrow, or its width too dependent on arrival time / drift depth).
  - Specific to det3 June X (the asymmetric board and view): late growth suppressed (peak → fall +0.10, vs det3 May X +0.34
    and det3 June's own Y +0.28) on top of a higher early level. Less late X spreading, not faster.
- **Resistivity gradient across the strip (Dylan's idea; `skew/grad_strip_toy.py`, `grad_strip_toy.json`):** 1-D RC
  spreading across a 550 µm strip with σ(x) = σ̄ (1 + g x/(w/2)), insulating gaps, same weighting / shaper / α as the toys,
  12 k clouds per point. g = 0.9 (19:1 across the strip) vs g = 0: |ΔA| ≤ 0.011 ± 0.006, |ΔO| ≤ 0.012, r(f) unchanged, at
  5 MΩ/sq and ×100, σ 150/300, with and without the funnel (toy floor ~0.01–0.02 from the asymmetric δ set). A strip settles
  to one potential whatever its profile, so only the transient acts, and it is < 1/10 of the data. Closed.

## 11. The sim's "too narrow early clouds" was a selection mismatch; the data carry partial-column pulses the sim lacks (2026-10-22)
Scripts and outputs: `../2026-10-22/cloud/` (README in each docstring). Condor 4526808 (cloud arms) and 4526809 (x-tilt arms),
`run_cloud_job.sh`; digitizer options `--diff-t-scale` and `--x-tilt-deg` (default-off, bit-identical: `cg0chk` = glg0 npz).
- **The §8–10 near-vertical cut compared different populations.** FWHM ≥ the set's own 25th percentile keeps short pulses
  in the data: det3 May X |dt_nb| < 30 ns has median FWHM 423 ns vs the vertical sim's 590 ns. In absolute FWHM bins
  (`fwhm_bins.py`) det3 May X matches glg0 at the sim's pulse lengths (s1 0.45–0.7, FWHM 450–550: r05 0.33 ± 0.09 /
  0.39, f50 0.92 / 0.95; 550–650: 0.21 ± 0.06 / 0.33, f50 1.00 / 0.95); early/peak climbs to ~1–2 only below 350–450 ns.
- **The short pulses are partial drift columns, not small tilts** (`short_pulses.py`): det3 May X/Y: onset ~4 samples
  (~240 ns, ~9 mm of drift) later than long pulses, the same 50 % fall, amplitude −25 % (top of the column missing);
  det3 June-15 (6-23, 600 V drift): normal onset, fall 12.6 vs 20.6 samples (bottom missing). Fraction of |dt_nb| < 80:
  det3 May 29 % X / 43 % Y, det3 June 62 / 84 %, sim 1–4 %. Sheared th00 at 1/2/3/5° and Stage A thx05/thx10
  (`mix_score.py`) almost never pass |dt_nb| < 80 (median 120–270 ns) and those that do have LOWER early sharing
  (tilt 2°: r05 0.03); a FWHM-matched tilt mixture does not reproduce the data's profile on any board.
- **Cloud-width arms vs det3 May X in matched-FWHM cells** (`arm_chi2.py`, shape χ² of r05, r25, f50 over the peak, 6 cells
  common to all arms, 18 dof): glg0 23.9, σ0 400 µm 17.5, feedback spread 400 µm 19.0, **dry-gas D_T 86**, no drift
  diffusion + σ0 400 / 250 µm: 57 / 565. The time profile excludes the dry-gas D_T (the board's "widens like dry gas"
  question) and depth-independent clouds; σ0 / feedback spread are unconstrained at this statistics (26–48 events per
  data cell). det3 May Y has too few long pulses for a matched comparison (sim Y pulses are also shorter than X: 521 vs 592 ns).
- **det3 June's reference set is the 6-23 session at drift 600 V** (`d3jun_x/y` = 6-23 long_run); det3 May is 900 V,
  the sim 1000 V. The §8–10 det3 June vs sim timing comparisons mixed drift fields and the short-pulse fraction.
- **det3 June X late growth, per session** (`late_xy.py`, G = r(f50)/r(peak) on each session's longest pulses, X/Y
  within a session): det3 May 0.64–0.81 (≈ 0.72), det3 June 1000 V 0.51–0.68, 600 V 0.47–0.55, 100 V 0.46 ± 0.07, det2
  0.50–0.61, vertical sim 0.56–0.60. det3 May is the high one; the sim sits with the asymmetric boards, so "suppressed
  late X growth" is not a clean det3 June fingerprint (weak drift-field trend within det3 June, ~2σ).
- **Consequence:** compare data and sim timing only at matched absolute FWHM (and drift field). The open item is what
  makes the partial-column pulses (out-of-time tracks, a dead/low-field layer near the mesh or edge, attachment at low
  drift field?), and whether any det3 June X timing feature survives once they are removed.

## What this leaves
- The fingerprint (audit §4.1) gains: (2) the transfer is n−1 → n+1 direct, ±2 odd ≈ 12 % of ±1; (5) gain- and
  time-independent at fixed width; no session outlier. The even term is not a usable discriminator at the ±0.1 level.
- Excluded on top of the audit: a fixed linear channel operator (copy or antisymmetric coupling), charging-up.
- Open: what physical, position-level, sum-conserving mechanism moves charge from the n−1 pad to the n+1 pad, X-only (across the strips), with a smooth two-signed geometry per board, June bulks only, and an
  amplitude set by the cloud width.
- Strip-induction routes are now exhausted at the data's sharing (§6, §7): offset resistive strips are in the full sim and
  sit at the conducting-strip limit; ion return paths subtract; a lateral field from strips at different potentials would
  need ~kV between neighbours (E_z ≈ 33 V/µm) to tilt the ions by 0.1, implausible for cosmic-rate currents.
- Every strip route needs a cloud the strip pins while narrow, widened later; §8 tested it on archive data: det2 has that
  regime, det3 June-15 has the opposite (prompt-wide, flat in time). So no strip route covers both asymmetric boards.
- Fingerprints (§8–10): det3 June X has a high, flat early sharing and suppressed late X growth, stable over all sessions,
  not a coupling. Fast-equalising strips (σ_xx = 100 σ_yy) and strip resistivity gradients are excluded (§10).
- §10's "general sim problem" (data early sharing flat, sim rising) is withdrawn (§11): at matched pulse length det3 May X
  = sim. The data's near-vertical samples are dominated by partial-column pulses (late onset or early end) the sim lacks;
  dry-gas D_T and depth-independent clouds are excluded by the time profile. Identify those pulses next.
- Otherwise the remaining classes need the board or the detector: the strip photo (registration, width, gap) and
  the detector tests (collimated scan across a pad, HV ramps and dwell, dry gas).
