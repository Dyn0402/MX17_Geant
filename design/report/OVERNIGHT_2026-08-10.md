# Overnight worker log — night of 2026-08-09 → 08-10

Running record, appended as work lands so a crash loses nothing. Authority docs
for this session: `RESPONSE_SIM_PLAN.md` §0a, `T14_FREEZE_QUEUE_2026-08-09.md`,
`S3_ION_CLOSEOUT_2026-08-09.md`, `GAS_AND_DRIFT_CAGE_ROADMAP_2026-08-08.md`.

Discipline held to throughout: the frozen T14 default (ρ_s = 2 MΩ/sq · DRY 95/5
· det3 bundle as-analysed) is not touched; everything below the verdict line is
DIAGNOSIS in its own directory. No joint (β, f_ion) fit. No contaminant number
is treated as measured. The drift-cage solve is not started.

---

## 0. Fleet state at 23:10

| host | state |
|---|---|
| laptop | working — all analysis below ran here |
| lxplus | reachable; 1 unrelated job running (`slim_wrapper.sh`, cluster 16705136). No MX17 campaign in the queue. |
| desktop | ⚠️ **BLOCKED — Tailscale SSH re-authentication required.** New sessions get `To authenticate, visit: https://login.tailscale.com/a/l11c9025b3510a7`. Host is up (ping 16 ms, port 22 open). |

**Desktop blocker — needs Dylan, one click.** Tailscale SSH is refusing new
sessions pending an interactive check-in. Consequence for the queue: the T7
slope-hunt (item 3) and the 08-08→09 overnight-chain verification cannot be
driven from this session. A pre-existing ssh channel opened at 22:05 by another
session is still tailing `/media/ucla/mx17_response_sim/slopehunt_chain2.log`,
so **the slope-hunt chain is running** — that channel predates the auth
expiry. Flagged, not worked around.

---

## 1. Queue item 1 — T14 frozen first look: ALREADY RUN

Checked before doing anything, as instructed. The comparison is **not** unrun:
`~/x17/response_sim/stageB_w2/t14_compare/` carries `t14_summary.json`,
`report.html`, and the four `wf_{sim,data}_{x,y}.parquet` legs. The verdict is
written up in `mx17_sim_wft/T14_CAMPAIGN_2026-08-09.md` and published at
<https://dylan-neff.web.cern.ch/notes/t14-sim-vs-data-waveforms.html>.

The T14_FREEZE_QUEUE line calling the comparison "frozen and unrun" is stale —
it describes the state at the moment of the freeze, not the state at 23:00.

### The verdict, det3-native §9 rows, sim vs data side by side

Frozen default, 2500 events per leg per view.

| observable | X sim | X data | X ratio | Y sim | Y data | Y ratio |
|---|---|---|---|---|---|---|
| peak amplitude, median [ADC] | 1429.8 | 2532.3 | **0.565** | 1129.1 | 2139.5 | **0.528** |
| event charge, median [ADC] | 3158.8 | 5316.0 | **0.594** | 3403.8 | 6079.3 | **0.560** |
| q_sum, reco-tight ±3° | 3776.2 | 5885.5 | **0.642** | 2878.2 | 4720.0 | **0.610** |
| strips over threshold, median | 5 | 6 | — | 7 | 8 | — |
| rise 10–90 %, median [ns] | 337.4 | 285.2 | ×1.18 | 307.4 | 261.8 | ×1.17 |
| FWHM, median [ns] | 645.9 | 543.9 | ×1.19 | 552.4 | 396.0 | ×1.39 |
| χ²/dof, median | 19.34 | 83.55 | — | 20.35 | 149.36 | — |
| fractional model residual | 0.034 | 0.027 | — | 0.043 | 0.042 | — |
| saturated fraction (>3500 ADC) | 0.059 | 0.260 | — | 0.038 | 0.110 | — |
| normalised-shape RMS / peak | 0.130 | — | — | 0.145 | — | — |
| seed median amplitude [ADC] | 522 (FEU 7) | 215 | — | 386 (FEU 8) | 184 | — |
| hits per seeded event | 5.54 | 22.26 | — | 7.78 | 28.47 | — |

**Verdict, unchanged and standing: same ballpark, not in agreement.** Amplitude
is the defect — sim is ×0.53–0.64 low on every amplitude-like observable, and
because the data leg clips at 3550 ADC while the sim barely saturates, the true
deficit is *larger* than the ratios show. Shape is close: normalised-waveform
RMS 13–14.5 % of peak, and the same forward model fits both legs at 3–4 %
fractional residual.

Two rows carry caveats already established elsewhere and repeated here so the
table is not read naked: the **Y FWHM ×1.39 contrast does not survive** the
window audit (`hv_slope/HV_SLOPE_2026-08-09.md` §4–5 — it is a θ-window
artifact), and the **data legs are saturation-biased subsamples** (the
reco-quality cut drops railed waveforms one-sidedly, data Y 0.327 → 0.110
against sim acceptance 0.989).

Nothing was re-run to improve this. Per the freeze, it stands.

---

## 2. Queue item 2 — angled ladder: COLLECTED, and the topology question is
## answered NO

The condor campaign (`stageBC_w2_angled.sub`, 5 points) had already completed
and been reduced by another session: `t14_angW_*` (the fixed ±3° window set,
which is the correct one per `ANGLED_LADDER_2026-08-09.md` §5) and the trend
report `t14_ang_trend/report.html`.

**What was missing: the answer.** The report's closing section is titled *"The
rise floor — is it geometry or electronics?"*, lays out the test correctly, and
then stops at the figure captions without stating a verdict. This session
computed it.

### The threshold artifact — a correction to the committed record

Both `ANGLED_LADDER_2026-08-09.md` §4 and the trend report state that *"the
data's sub-200 ns fraction explodes with angle (0.23 → 0.82 → 0.99) while the
sim barely responds (0.020 → 0.048 → 0.066)"*, and read that as the sim failing
to convert inclination into fast rise.

**That reading is an artifact of the 200 ns threshold, which sits below the
sim's own floor.** Re-counting at the 240 ns threshold the rest of the campaign
uses (X view):

| |θ| | sim f<240ns | data f<240ns | sim f<200ns | data f<200ns |
|---|---|---|---|---|---|
| 0° | 0.045 | 0.378 | 0.020 | 0.226 |
| 10° | **0.269** | 0.963 | 0.048 | 0.820 |
| 20° | **0.749** | 0.995 | 0.066 | 0.988 |

The sim's fast fraction goes 0.045 → 0.749 — it responds to inclination
strongly. A counter placed below a distribution registers almost nothing while
that distribution slides past it, which is exactly what the 200 ns cut was
doing.

### The threshold-free test: quantile-by-quantile

`t14_ang_trend/rise_offset_vs_angle.json` (script:
`scratchpad/angle_quantile_offset.py`). Sim − data rise time, in ns, at matched
quantiles of each angle-matched pair:

**X view**

| |θ| | p5 | p10 | p25 | p50 | p75 | p90 |
|---|---|---|---|---|---|---|
| 0° | 90.5 | 94.1 | 81.8 | 40.0 | −33.1 | −31.8 |
| 10° | 76.8 | 82.4 | 78.3 | 75.7 | 71.1 | 66.6 |
| 20° | 79.2 | 79.9 | 72.6 | 69.3 | 69.0 | 72.0 |

**Y view**

| |θ| | p5 | p10 | p25 | p50 | p75 | p90 |
|---|---|---|---|---|---|---|
| 0° | 64.2 | 74.5 | 72.6 | 12.5 | −108.6 | −88.2 |
| 10° | 53.2 | 66.4 | 68.6 | 67.2 | 62.3 | 52.7 |
| 20° | 37.8 | 55.0 | 66.5 | 68.1 | 69.3 | 69.9 |

How much each leg speeds up from 0° to 20°, and the ratio (1.00 = the sim's
angular response matches the data's):

| view | | p5 | p10 | p25 | p50 | p75 | p90 |
|---|---|---|---|---|---|---|---|
| X | sim | −56.1 | −51.1 | −64.2 | −94.2 | −169.1 | −285.5 |
| X | data | −44.8 | −36.9 | −55.0 | −123.5 | −271.1 | −389.2 |
| X | **ratio** | 1.25 | 1.39 | **1.17** | 0.76 | 0.62 | 0.73 |
| Y | sim | −42.6 | −39.4 | −45.4 | −71.3 | −135.3 | −258.1 |
| Y | data | −16.3 | −19.9 | −39.3 | −126.9 | −313.2 | −416.2 |
| Y | **ratio** | 2.61 | 1.98 | **1.16** | 0.56 | 0.43 | 0.62 |

### The verdict

**No — track inclination is not the missing mechanism, and the reason is
sharper than "the sim doesn't respond".** At the two inclined points the sim's
rise distribution is the data's **rigidly displaced by ~65–80 ns at every
quantile**: X 10° spans 76.8 → 66.6 ns of offset from p5 to p90, X 20° spans
79.2 → 72.0. The distributions have the same *width* — at 20° the sim's p5–p90
span is 63 ns against the data's 71 ns — and the sim's response to inclination
at the fast end matches the data's to within 17 % (p25 ratio 1.17 X, 1.16 Y).

So there is no missing fast *population* to generate. There is a **constant
additive delay of ~70 ns on every pulse**, angle-independent, which inclination
cannot remove because inclination is already being converted into rise time at
about the right rate.

This **strengthens the S3 ion closeout rather than competing with it.** A rigid
offset applied uniformly across the distribution is the signature of a fixed
extra time constant in the response chain, and the f_ion dial is measured to
move rise by 96 ns across its range while moving nothing else (S3 closeout, "the
two dials are orthogonal"). The demand is not "add fast pulses" — it is "remove
~70 ns from all of them", which is precisely what f_eff ≈ 0.16 does.

**Item 2 was billed as the last live candidate for the rise-time contradiction.
It is now closed, negative.** The causal list in the S3 closeout stays empty and
the f_eff ≈ 0.2-vs-0.9056 contradiction stands as the result.

### One caveat, and the test that would close it

At **vertical only**, the offset is *not* constant — it runs +90 ns at p5 but
−33 ns at p90 (X), i.e. the data's vertical sample is broader than the sim's in
*both* directions, with a slow tail the sim lacks. The inclined points show no
such thing.

The likely cause is the campaign doc's own open thread #4: `_sel_ids` in
`t14_compare.py` windows the data on **one view's** θ only, so the "vertical"
data leg admits cosmics with unconstrained inclination in the *other* view.
That contaminates a nominally-vertical data sample with genuinely inclined —
hence fast — tracks, and the sim gun is a true pencil beam with no such
contamination.

**Not tested here**: the data-leg events parquet carrying `x_theta_deg` /
`y_theta_deg` is not on the laptop (searched `~/x17`, `/media/dylan/data/x17`);
it was written by another session's `t13_reco` run. The test is one line — cut
the data on `|θ_other| < 3°` before the comparison — and it should be run before
the vertical p75/p90 rows are interpreted as physics. It does **not** touch the
verdict above, which rests on the inclined points where both legs are windowed
in the tilted view.

Products: `t14_ang_trend/fastrise_verdict.json`,
`t14_ang_trend/rise_offset_vs_angle.json`.

---

## 3. Interlude — open thread #4 tested and answered NO

The vertical point's non-rigid offset (§2 caveat) was suspected to be the
one-view θ window admitting other-view inclination. The θ-tagged parquet for the
target run was located (`…/mx17_det3_saturday_scan_6-27-26/long_run_resist_490V
_drift_1000V/mx17_3/wft/events.parquet`, 7093 rows, per-view θ) and intersected
with the vertical legs.

| view | leg | n | offset span |
|---|---|---|---|
| X | published | 2324 | 127.2 ns |
| X | overlap only | 692 | 129.2 ns |
| X | overlap + \|θ_y\| < 3° | 171 | **126.2 ns** |
| Y | published | 1729 | 183.1 ns |
| Y | overlap only | 645 | 185.4 ns |
| Y | overlap + \|θ_x\| < 3° | 145 | **187.3 ns** |

Cutting the other view removes 75 % of the events and moves the span by +1 ns
(X) and −4 ns (Y). **The vertical broadening is real, not a selection
artifact.** The control is overlap-with-cut against overlap-without-cut — both
drawn from the same events — so the subsample is not doing the work.

What stands: at normal incidence the data's rise distribution is genuinely
wider than the sim's at *both* ends, and a pencil-beam monoenergetic gun does
not reproduce the topological diversity of real vertical cosmics. The data's
vertical χ²/dof (83.5 X, 149.4 Y vs sim 19.3/20.4) says the same thing
independently. **The clean place to measure the rise offset is therefore the
inclined points**, and that is what §4 does.

Product: `t14_ang_trend/theta_other_cut.json`.

---

## 4. ⭐ The night's main result — the ion term IS the rise discrepancy,
## measured across the whole distribution at two inclinations

This closes a chain that had been resting on a single statistic.

### Step 1 — f_ion is NOT a rigid delay at vertical

If the sim–data mismatch is a rigid ~70 ns offset (§2), the natural check is
whether the f_ion dial produces a rigid shift. On the vertical DIAGNOSIS legs
already on disk, **it does not** — it *compresses*. Sim rise shift relative to
`noions`, X view:

| leg | f_ion | p5 | p10 | p25 | p50 | p75 | p90 | span |
|---|---|---|---|---|---|---|---|---|
| fion030 | 0.30 | 14.9 | 21.7 | 28.0 | 31.6 | 11.4 | −4.5 | 36 |
| fion050 | 0.50 | 44.9 | 53.8 | 57.2 | 46.7 | 7.7 | −11.0 | 68 |
| fion070 | 0.70 | 81.0 | 86.4 | 80.0 | 55.7 | 0.8 | −20.8 | 107 |
| **default** | **0.9006** | **114.9** | **112.8** | **98.0** | **67.1** | **1.8** | **−28.2** | **143** |

The full ion term pushes p5 back by +115 ns and p90 *forward* by −28 ns. It
imposes a floor; it does not delay. Y behaves the same (+110 / −28).

Taken alone this looks like it refutes the §2 story.

### Step 2 — but at INCLINED incidence it is very nearly a rigid delay

The resolution is geometric and had to be measured. At vertical the sim's rise
distribution spans ~280 ns p5–p90, so the ion floor only bites on the fast
half. At 20° the distribution is already compressed to a ~63 ns span — nearly
everything sits *on* the floor — so there the same dial acts almost uniformly.

The `DIAGNOSIS_noions_thx10/20` legs were already on disk. X view:

**10°**

| leg | p5 | p10 | p25 | p50 | p75 | p90 | mean offset | span |
|---|---|---|---|---|---|---|---|---|
| sim, with ions | 201.3 | 225.0 | 238.9 | 251.4 | 263.9 | 278.1 | | |
| sim, ions removed | 113.0 | 128.5 | 141.2 | 158.6 | 183.8 | 211.6 | | |
| data | 124.6 | 142.6 | 160.6 | 175.8 | 192.8 | 211.4 | | |
| **offset, with ions** | 76.8 | 82.4 | 78.3 | 75.7 | 71.1 | 66.6 | **+75.2** | 15.7 |
| **offset, ions removed** | −10.1 | −13.0 | −18.1 | −15.2 | −6.2 | 5.4 | **−9.5** | 23.5 |
| ion term buys | −88.3 | −96.5 | −97.6 | −92.8 | −80.1 | −66.5 | −87 | 31.2 |

**20°**

| leg | p5 | p10 | p25 | p50 | p75 | p90 | mean offset | span |
|---|---|---|---|---|---|---|---|---|
| sim, with ions | 187.6 | 211.0 | 223.6 | 231.5 | 240.0 | 251.0 | | |
| sim, ions removed | 103.2 | 110.0 | 123.4 | 132.8 | 141.2 | 151.8 | | |
| data | 108.4 | 131.1 | 151.0 | 162.2 | 171.1 | 179.1 | | |
| **offset, with ions** | 79.2 | 79.9 | 72.6 | 69.3 | 69.0 | 72.0 | **+73.7** | 10.9 |
| **offset, ions removed** | −5.2 | −21.5 | −28.4 | −29.6 | −30.0 | −27.1 | **−23.6** | 24.9 |
| ion term buys | −84.4 | −101.0 | −100.2 | −98.7 | −98.8 | −99.2 | −98 | 16.6 |

### The result

At inclined incidence the ion term is a **near-rigid ~90–100 ns delay** (span
17–31 ns across p5–p90), and the data needs about **75 of those ~95 ns
removed**. Deleting the ion term entirely slightly *overshoots* — the sim
becomes 10 ns (10°) to 24 ns (20°) too fast at every quantile.

**So the data's rise is reproduced at f_eff somewhere between 0 and ~0.25, on
six quantiles at two independent inclinations.** That is the S3 closeout's
f_eff ≈ 0.16 demand, but established on the *whole distribution* instead of on
p5 of a single vertical sample — and the vertical p5 statistic was the weakest
leg of that argument, since the vertical data leg is now known to be anomalously
broad (§3).

**The f_eff ≈ 0.2 versus defended 0.9056 contradiction is confirmed and
substantially strengthened.** It is no longer a one-statistic demand curve; it
is a distribution-wide, angle-independent measurement. Nothing here is a fit —
no parameter was adjusted; the legs were already on disk and were simply read at
matched quantiles.

Products: `t14_ang_trend/fion_rigidity.json`,
`t14_ang_trend/noions_angled_closure.json`.

---

## 5. Queue item 6 — Y over-sharing / kY: the strip anisotropy is the SHARING
## asymmetry but NOT the amplitude asymmetry

`RHO_S_SENSITIVITY_2026-08-09.md` closed with *"the Y response model is the more
likely culprit than the ρ_s number"*. This tests the most obvious candidate
mechanism and **refutes it for amplitude**.

### The mechanism proposed, and its falsifier

Reading the Stage B / S1 Y path (`response/solver/kernels.py`,
`response/solver/wpot.py`): the ESL is 550 µm strips on an 800 µm pitch,
patterned in **x** and uniform in **y** (`wpot.py`: *"y is uniform in the model,
so its modes decouple completely"*). So charge spreads freely along y and is
blocked across x by the 250 µm gaps. An X channel is a pad **column** — a comb
running along y, *parallel* to the strips, which therefore keeps charge that
diffuses along them. A Y channel is a pad **row**, *perpendicular* to the
strips, which loses it.

Falsifier stated before running: replacing the patterned sheet with a uniform
sheet of the same ρ_s — the only change — must collapse the asymmetry; if it
does not, the asymmetry lives in the Y code path instead.

### Result — half confirmed, half refuted

Strips vs uniform at ρ_s = 2 MΩ/sq, identical box, grid, pad patterns and
shaper; uniform arm via `solve_uniform_analytic` (the V1 closed form):

| arm | X amp₀ | X rms [strips] | Y amp₀ | Y rms [strips] | Y/X |
|---|---|---|---|---|---|
| strips, prompt | 2.666e-10 | **0.45** | 2.955e-10 | 1.09 | 1.109 |
| uniform, prompt | 2.259e-10 | **1.27** | 2.450e-10 | 1.28 | 1.084 |
| strips, ion-folded | 1.782e-10 | **0.46** | 1.870e-10 | 1.25 | 1.049 |
| uniform, ion-folded | 1.229e-10 | **1.47** | 1.328e-10 | 1.47 | 1.081 |

* **Sharing: confirmed, and it is entirely the strips.** The strip pattern
  compresses X sharing to 0.45–0.46 strips against Y's 1.09–1.25 — a factor
  2.4–2.7 anisotropy — and going uniform makes the two views *identical*
  (1.27/1.28 prompt, 1.47/1.47 ion). The sheet model is strongly anisotropic in
  sharing, exactly as the geometry demands.
* **Amplitude: refuted.** Y/X moves only 1.109 → 1.084 (prompt) and 1.049 →
  1.081 (ion — the *wrong way*). The strip anisotropy does essentially nothing
  to the peak-amplitude ratio.

### What that means for the hunt

The T14 full chain reads sim Y/X = 0.790 on peak amplitude (1129.1/1429.8),
against an unselected detector value of 1.0002. The S1 kernel — the only place
the strip geometry enters — produces Y/X ≈ 1.05–1.11 here and 0.91 in the
archived harness. **Neither is anywhere near 0.79.** So the X/Y amplitude
asymmetry is **not in the electrostatics**; it enters downstream, in Stage B/C's
use of the kernel or in the comparison's per-view selection. That is where item
6 should continue: `kY`, the per-view calibration, drift-diffusion smearing, and
the fact that the two views are read by different FEUs.

⚠️ **Honest discrepancy, not smoothed over.** My Y/X = 1.049 (ion) against
`rho_s_sensitivity`'s 0.908 for W1 rho2M is a 15 % disagreement between two
harnesses. Mine solves fresh at nx = 1560 / ny = 512 with the W1 grounded-gap
boundary (`_pads_at` leaves inter-pad channels at 0) and approximates the ion
transit with a 24-step ramp; theirs uses `CombKernelLUT` at nx = 3120 / ny =
1024 with `apply_ion_transit`. **Only the strips-vs-uniform contrast, which is
internal to a single run, is quoted as a result** — the absolute Y/X from this
script is not a production number.

Products: `design/report/xy_anisotropy_ab_2026-08-10.json`,
`scratchpad/overnight_2026-08-10/xy_anisotropy_ab.py`.

---

## 6. Queue item 7 — wet amp-range gain bracket: SUBMITTED (cluster 16705137)

Runnable despite the desktop being down: lxplus + Kerberos are live (ticket to
08/11 00:12) and this needs neither the desktop nor the field map.

Pre-registration, written before submission and committed first (`287612b`):
`design/report/WET_GAIN_BRACKET_PREREG_2026-08-10.md`.

**The roadmap's premise was wrong, in our favour.** §2 step 1 says the wet
suites are drift-range only and this needs ~2 new Magboltz jobs. **No Magboltz
job is needed** — dry, +0.5 % and +1 % H₂O all already exist at amplification
range on an identical 5 000–60 002 V/cm grid at `Saclay_160m`. The `.gas` files
store **E/p**, so the header's "6.704 … 80.45" is V/(cm·Torr) and reads as a
drift table until multiplied by 745.83 Torr; the genuinely drift-range tables
carry `_drift_` in the filename. Roadmap row to be corrected.

**A confound was found and designed out.** `mm_config.py` gives the dry mixture
`penning: auto` and both wet ones `penning: manual rP = 0.4` — necessarily,
since Garfield has no ternary Ar/iC₄H₁₀/H₂O Penning table. Comparing
dry(auto) against wet(rP 0.4) would have confounded water with a Penning-model
change, and Penning is the very knob the T7 slope hunt is chasing. Penning is
applied at avalanche time (`mm_sim_core.py:69-73`), not baked into the table, so
all three mixtures run at **manual rP = 0.40 on argon**, with dry-at-auto as a
fourth arm carried only to tie back to production.

32 jobs, 490 V, 0.015 cm gap, 8 batches × 200 events per arm:

| arm | mixture | Penning |
|---|---|---|
| A | Ar/iC₄H₁₀ 95/5 | manual rP 0.40 |
| B | Ar/iC₄H₁₀/H₂O 94.5/5/0.5 | manual rP 0.40 |
| C | Ar/iC₄H₁₀/H₂O 94/5/1 | manual rP 0.40 |
| D | Ar/iC₄H₁₀ 95/5 | auto |

Collect with `mm_condor_collect.py` against
`…/garfield_sim/jobs_wetbracket/`. **All water fractions are labelled
fitted-to-data**; this measures a derivative, not an operating point.
