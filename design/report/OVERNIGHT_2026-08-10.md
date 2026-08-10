# Overnight worker log — night of 2026-08-09 → 08-10

Running record, appended as work lands so a crash loses nothing. Authority docs
for this session: `RESPONSE_SIM_PLAN.md` §0a, `T14_FREEZE_QUEUE_2026-08-09.md`,
`S3_ION_CLOSEOUT_2026-08-09.md`, `GAS_AND_DRIFT_CAGE_ROADMAP_2026-08-08.md`.

Discipline held to throughout: the frozen T14 default (ρ_s = 2 MΩ/sq · DRY 95/5
· det3 bundle as-analysed) is not touched; everything below the verdict line is
DIAGNOSIS in its own directory. No joint (β, f_ion) fit. No contaminant number
is treated as measured. The drift-cage solve is not started.

---

## SUMMARY — read this first

**Queue status: items 1, 2, 4, 5, 6 and 7 done or submitted. Only item 3
(T7 slope hunt) is outstanding, blocked on the desktop.**

### One thing needs Dylan

**The desktop needs a Tailscale SSH re-authentication** —
<https://login.tailscale.com/a/l11c9025b3510a7>. The host is up (ping 16 ms,
port 22 open); only the auth check-in is missing. The slope-hunt chain is
running on it unattended and self-merges, so it collects in the morning either
way — but nothing else could be driven there tonight.

### The results, in order of how much they change

1. **The ion term IS the rise discrepancy, now measured across the whole
   distribution** (§4). At 10° and 20° the ion term acts as a near-rigid
   ~95 ns delay and the data needs ~75 of those removed; removing it entirely
   slightly overshoots. So the data's rise is reproduced at **f_eff between 0
   and ~0.25, on six quantiles at two independent inclinations** — where before
   this rested on p5 of a single vertical sample. The f_eff ≈ 0.2 vs defended
   0.9056 contradiction is confirmed and much harder to dislodge. Nothing was
   fitted.

2. **Track inclination is eliminated as the rise explanation** (§2), and the
   committed claim that "the sim barely responds" to inclination is
   **withdrawn** — it was an artifact of a 200 ns threshold sitting below the
   sim's own floor. Corrected in `nTof_x17 ANGLED_LADDER_2026-08-09.md`.

3. **The amplitude ledger closes to a single candidate** (§7, §12). The deficit
   is a **charge** deficit of ×0.63 that is *f_ion-independent*, so it cannot be
   double-counted against the ion thread. It demands a real gain of ~4 × 10⁴
   against the sim's 24 094. Every other row is now dead, controlled or
   separated — including primary ionisation, checked tonight and correct
   (91.2 e⁻/cm, W = 25.97 eV). **The avalanche gain is the only survivor**, and
   the independent 12 σ HV-slope error points at the same defect.

4. **The X/Y asymmetry is an 18 ± 3 % modelling error**, not the factor-0.74 bug
   an earlier section of this very report claimed (§8, corrected in §10). It is
   real, sim-side, and not selection — but modest, and made of two ~10 % pieces.

5. **det3 prefers O₂-like attachment** over attachment-free H₂O-like transport
   (§11) — with the honest caveat that the λ(E) *shape* test does not close, so
   the preference rests only on a decay existing at all.

### Standing open question

**What suppresses the ion-induction term at the readout by ×4–5, when the
charge split, the template, the mesh weighting field, the electronics register,
the gap geometry and the sheet screening all check out individually?** Every
proposed mechanism is eliminated (§9). The one soft spot in that elimination
chain is named there and not papered over.

### Corrections made to the committed record tonight

Three, all dated and loud, per project habit:

* `ANGLED_LADDER_2026-08-09.md` §4 — "the sim barely responds" **withdrawn**
  (threshold artifact); the ladder's missing verdict supplied.
* `GAS_AND_DRIFT_CAGE_ROADMAP_2026-08-08.md` §2 — "~2 new Magboltz jobs"
  **wrong**, zero needed (E/p units trap); and the bracket as specified carried
  a Penning confound.
* This report §8 → §10 — my own factor-0.74 claim **withdrawn**.

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

---

## 7. Queue item 5 — amplitude-deficit ledger

### 7.1 The deficit measured in CHARGE is f_ion-independent — no double-counting

The peak-amplitude deficit and the ion thread do couple, so the ledger has to be
quoted at more than one operating point. Doing that shows they separate cleanly:

| leg | X peak | X q_sum | Y peak | Y q_sum |
|---|---|---|---|---|
| default (f_ion 0.9006) | 0.5646 | 0.6416 | 0.5278 | 0.6098 |
| fion070 | 0.5573 | 0.6265 | 0.5200 | 0.5795 |
| fion050 | 0.5600 | 0.6265 | 0.5183 | 0.5802 |
| fion030 | 0.5809 | 0.6261 | 0.5457 | 0.5826 |
| noions (f_ion 0) | 0.6319 | 0.6330 | 0.5928 | 0.5981 |

**Peak moves ×1.12 across the full f_ion range; q_sum does not move at all**
(X: 0.6261–0.6416, a 2.5 % spread). That is exactly right — f_ion redistributes
charge in time and conserves it, so an integral is blind to it and only the
peak, which is set by shape against the shaper, responds.

**Consequence: quote the deficit as a CHARGE deficit of ×0.63, and it cannot be
double-counted against the ion contradiction.** However f_eff resolves, the
charge deficit survives unchanged. The extra ×1.12 on peak is the ion thread's
and belongs to it.

### 7.2 A signed correction for the data legs' saturation depletion

The reco-quality cut drops railed waveforms, so each data leg is the bottom
(1 − f) of the true amplitude distribution. f follows from saturation fractions
already in the record: detector 0.326 (X) / 0.327 (Y) unselected against the
legs' 0.260 / 0.110. Then the leg's own quantile 0.5/(1 − f) is the true median
— computable read-only from the frozen parquets.

| view | f missing | true median at leg | peak | q_event |
|---|---|---|---|---|
| X | 0.089 | p54.9 | 2532.2 → 2706.1 (**×1.069**) | 5316.0 → 5696.9 (×1.072) |
| Y | 0.244 | p66.1 | 2139.5 → 2569.2 (**×1.201**) | 6079.2 → 7187.1 (×1.182) |

**This correction validates itself.** It is derived only from saturation
fractions, with no reference to the X/Y question — yet it takes the corrected
data peaks to Y/X = 2569.2/2706.1 = **0.949**, recovering the independently
measured unselected value of 1.0002 from a leg that read 0.845. A correction
built for one purpose reproducing a number it was not fitted to is the kind of
check worth stating.

### 7.3 The gain the deficit demands

| view | q_sum ratio | G demanded | saturation-corrected |
|---|---|---|---|
| X | 0.6416 | 37 552 | **40 243** |
| Y | 0.6098 | 39 512 | 46 713 |

against the T7 pooled meshfield sim gain of **24 094** at 490 V.

The X number is the one to quote (Y's saturation correction is the large,
uncertain one). **The deficit demands a real gain near 3.8–4.0 × 10⁴, against
the literature's maximum stable Ar/iso bulk-MM gain of 3–4 × 10⁴** — det3 sits
right at the top of that band, which is consistent with it sparking above 500 V.
So a pure gain explanation is *physically available*, but only just: it puts the
detector at the edge of stable operation and leaves no headroom for any other
factor pulling the same way.

### 7.4 The elimination table

| candidate | leverage | verdict |
|---|---|---|
| **Avalanche gain / α(E)** | needs ×1.56; ×1.52 already measured on the SLOPE | **LIVE, and the only candidate with both the size and the signature.** Data d lnA/dV = 0.449/10 V vs sim 0.296 (≈12 σ). A gain-slope error at the operating point is exactly a gain-scale error. Gated on the T7 slope hunt. |
| ADC scale | 1.2 % | **Control row.** Certified against the datasheet-derived 20.48 ADC/fC (101.1 vs 102.4). Cannot contribute. |
| Electronics gain range | the only discrete ×2 on the data side | **Dead.** 44/44 archived cfgs read `Dream 6/7 = 0xAAAA` = 200 fC = 10 mV/fC. Inference, not read-back (target run archived no cfg), but a one-step range error is excluded. |
| ρ_s / sheet resistivity | ×1.15 (X) over a FACTOR 10 in ρ_s | **Dead by size and signature.** Within the T2b band 2→2.56 MΩ/sq buys ~5 %. A static sheet property also cannot produce a gain-SLOPE error. |
| W2 prompt capture | +25.66 % vs W1, already applied | **Applied, not available twice.** V6/W2 fixed a known-sign boundary error; W1 would be choosing a known-wrong BC. |
| β / PZC residual | 0.6–2.3 % on peak | **Dead.** Freeze-queue measurement across β 0.25→1. |
| f_ion / ion term | ×1.12 on peak, **0 on charge** | **Separated, not eliminated.** Owns the peak-shape part; contributes nothing to the charge deficit (§7.1). |
| Data saturation depletion | ×1.069 (X) | **Sized, signed, and it makes the deficit LARGER.** §7.2. Already folded into the corrected column above. |
| Primary ionisation yield | not sized | **OPEN — the one row with no number.** Geant4's W-value and clusters/mm for Ar/iso 95/5 against literature has not been checked. Needs doing; it multiplies the whole chain. |
| Mesh transparency | 0.955, cross-checked | **Unlikely.** T6's 3D transparency agrees with the avalanche survival 0.9559 independently; no room for ×1.5. |
| Diffusion / time-binning peak dilution | not sized | **OPEN, but bounded** — it cannot touch q_sum, and the deficit is quoted in q_sum. |

**Verdict: the amplitude deficit is a charge deficit of ×0.63, f_ion-
independent, demanding a real gain of ~4 × 10⁴ against the sim's 24 094.**
Every candidate except the avalanche gain is dead, controlled, or cannot touch
an integral. The two rows still owing numbers — primary-ionisation yield and
diffusion dilution — are worth closing, but neither is likely to carry a factor
1.6 on its own, and the gain candidate already has an independent 12 σ signature
pointing at it. **This is the T7 slope hunt's question, and it gates everything
on amplitude.**

Product: `t14_ang_trend/amplitude_ledger.json`.

---

## 8. Queue item 6 continued — the X/Y asymmetry is sim-side, and it is not
## selection either

Before hunting Stage B/C code, selection had to be excluded: the two views' legs
are selected independently and the data's Y leg is cut ~4× harder than its X
leg. Three progressively tighter comparisons, read-only on the frozen parquets —
*paired* keeps only events present in both views' legs, so the two views see an
identical event set:

| selection | n (sim / data) | sim Y/X | data Y/X |
|---|---|---|---|
| as published | 2500 / 2500 | 0.790 | 0.845 |
| paired | 2401 / 304 | **0.782** | 0.937 |
| paired + unsaturated | 2254 / 257 | **0.783** | **0.954** |

**The sim's asymmetry is completely immune to matching** (0.790 → 0.782 →
0.783). **The data's converges on unity** as the bias is removed — 0.845 →
0.954, against the independently measured unselected 1.0002. Two different
routes to the same place, since §7.2's saturation correction gets there too
(0.949) from saturation fractions alone.

So the chain of elimination for the X/Y amplitude asymmetry now reads:

1. **Not the data.** The detector is X/Y symmetric; the leg's 0.845 was
   selection, and it converges to 1.0 under two independent corrections.
2. **Not selection on the sim side.** Immune to pairing and de-saturation.
3. **Not the S1 electrostatics.** The strips-vs-uniform A/B (§5) moves Y/X by
   0.02–0.03 and the wrong way; kernel-level Y/X is ~1.05 (this harness) or 0.91
   (archived), never 0.78.
4. **Therefore it is Stage B/C**, which turns a kernel-level Y/X of ~1.05 into a
   decoded 0.78 — a factor ~0.74 applied to Y relative to X somewhere between
   the kernel and the decoded ADC.

Next places to instrument, in order: per-view charge budgets inside Stage B
(`charge_budget_y` vs `_x` bookkeeping at the LUT level), the per-FEU
configuration carried into Stage C, and drift-diffusion anisotropy. `kY = 1.375`
is worth naming but is *not* obviously the culprit — it lives in the calibration
bundle and `t13_reco` applies it identically to both legs, so it cannot by
itself create a sim-only asymmetry; it would have to be interacting with
something that differs.

Product: `t14_ang_trend/xy_paired_selection.json`.

---

## 9. OPEN QUESTION for the morning — the ion contradiction has no surviving
## mechanism

Recorded as an open physics question rather than pursued further tonight.

**The measurement.** The data's rise is reproduced at f_eff ≈ 0–0.25, now
established distribution-wide at two inclinations (§4). The defended value is
**f_ion = 0.9056**, re-derived through the real woven mesh and gated to 5e-8 on
its own linearity check. That is a factor of four to five.

**Everything proposed has been eliminated, individually and by measurement:**

| candidate | verdict |
|---|---|
| f_ion (charge split) | 0.9056 through the true ψ; shift +0.005 and the wrong way |
| i_ion template shape | validated to 3.8 % at every quantile by an independent reconstruction |
| ion species / mobility | Blanc's law: 4 % slower, wrong direction |
| T10 lateral factorisation | 3.7 ns |
| β / PZC residual | 4 ns across its whole range; rise immune |
| any missing high-pass | falsified as a CLASS by the exchange rate (226 undershoot points needed, budget −6.4) |
| peaking-time register | code 2, 44/44 archived cfgs |
| amplification-gap geometry | 150 µm bulk Micromegas, from the pillar gerber |
| resistive-sheet screening | 8 % where ×4.5 needed |
| **track inclination** | **closed tonight, negative (§2, §4)** |

**The question to put to Dylan:** *what suppresses the ion-induction term at the
readout by a factor 4–5, when the charge split, the template, the mesh weighting
field, the electronics register, the gap geometry and the sheet screening all
check out individually?*

⚠️ **One honest pointer, recorded without pursuing it.** The sheet-screening
sizing that retired the sharpest structural candidate (8 % where ×4.5 was
needed) rests on the assumption that `apply_longitudinal` already carries most
of the sheet dynamics, so the induced-vs-injected distinction is only a ≤10 %
correction. That assumption is the one place in the elimination chain where a
daylight re-derivation could still move a factor — and this project's own record
(Fix1, Fix7, the ρ_s intuitions, the toy that over-predicted amplitude ×6, and
tonight's threshold artifact) is that predictions flip on contact with numbers
more often than is comfortable. It is not a reason to reopen the item tonight;
it is a reason not to call the elimination chain airtight.

---

## 10. ⚠️ CORRECTION to §8 — the X/Y asymmetry is an 18 % effect, not a
## factor 0.74

§8 concluded that "Stage B/C turns a kernel-level Y/X of ~1.05 into a decoded
0.78 — a factor ~0.74 applied to Y". **That comparison was invalid and the
factor is wrong.** It set a point-charge kernel's central-channel amplitude
(`amp0`, no charge cloud, no summation over channels) against a full-chain
per-event peak. Those are not the same quantity and their ratio is not a bug
size.

The right comparison is sim against data in the *same* observable. Paired,
de-saturated, bootstrap 2000×:

| quantity | sim | data | sim/data |
|---|---|---|---|
| peak Y/X | 0.784 ± 0.006 | 0.951 ± 0.028 | **0.825 ± 0.025** |
| q_event Y/X | 1.063 ± 0.005 | 1.186 ± 0.032 | **0.897 ± 0.024** |
| implied width factor | 1.360 | 1.237 | **1.087 ± 0.038** |
| n_over Y/X | 1.400 | 1.333–1.600 | ≈ 1 |

**Charge is not lost, and the effect is not large.** Both legs show the same
qualitative structure — the Y view carries *more* charge than X (q_event Y/X >
1 on both) while showing a *lower* peak, because Y spreads over more channels
(n_over Y/X = 1.4 sim, 1.3–1.6 data). The detector does this too; the sim merely
overdoes it.

The residual decomposes exactly (0.897 × 1.087 = 0.825, closing to three
figures) into two comparable and modest terms:

* **charge partition — 0.897 ± 0.024** (4.3 σ): the sim puts ~10 % less charge
  into Y relative to X than the detector does.
* **extra spreading — 1.087 ± 0.038** (2.3 σ, marginal): the sim spreads Y ~9 %
  more than the detector.

So item 6's target is **an 18 ± 3 % modelling error made of two ~10 % pieces**,
not a cornered factor-0.74 bug. That is a materially different — and much less
dramatic — object than §8 described, and it is worth much less urgency. What
survives from §8 unchanged: the asymmetry is real, it is sim-side, it is not
selection (the sim's Y/X is immune to pairing and de-saturation at 0.790 → 0.782
→ 0.783), and the S1 electrostatics are exonerated.

Caveat on the data leg: n = 257 paired unsaturated events, hence the ±0.028.
The charge-partition term is solid; the spreading term is 2.3 σ and should not
be quoted as established.

Product: `t14_ang_trend/xy_paired_selection.json`.

---

## 11. Queue item 4 — attachment-shape family test on det3

### The test already exists and had already been run for det3

`mx_june_cosmic_qa/17_gap_attachment_test.py` measures the per-strip amplitude
decay length λ against drift depth for each drift-scan point, and
`18_attachment_vs_magboltz.py` compares it to Magboltz 1/η. Both had already
been run on the det3 saturday scan —
`…/mx17_det3_saturday_scan_6-27-26/drift_velocity/mx17_3/gap_attachment_test.csv`
— so this item is the comparison and the verdict, not the measurement.

### The discriminator is sharp

Across the whole 30-mixture `water2d.json` grid, **29 mixtures give η = 0
exactly at 333 V/cm** — every H₂O/N₂ combination, at every isobutane fraction
from 3 % to 8 % and every water fraction from 0.4 % to 1.1 %. Water does not
attach. The single grid mixture containing oxygen
(`iso5.0_h2o0.8_n20.78_o20.21` — the 0.8 % H₂O + ~1 % air row the freeze queue
singled out) gives η = 0.809 /cm, λ = 12.4 mm.

So the family question has a yes/no observable: **any finite amplitude decay
with depth excludes attachment-free H₂O-like transport**, whatever the
concentration.

### det3's numbers

| drift HV | E [V/cm] | measured λ [mm] | air 1 % | air 2 % | O₂ 0.5 % | Ar/CO₂ 90/10 |
|---|---|---|---|---|---|---|
| 500 | 166.7 | 11.12 | 12.80 | 5.84 | 5.64 | ∞ |
| 700 | 233.3 | 8.37 | 16.32 | 7.67 | 7.20 | ∞ |
| 900 | 300.0 | 5.74 | 19.47 | 9.22 | 8.54 | ∞ |
| **1000** | **333.3** | **5.60** | 21.02 | 9.88 | 9.25 | ∞ |
| 1100 | 366.7 | 9.10 | 22.63 | 10.55 | 9.99 | ∞ |

### Verdict — and a caveat that limits it

**det3's waveform data prefers O₂-like attachment over H₂O-like
attachment-free transport.** The decay is finite at every drift field (λ =
5.6–11.1 mm), and no water fraction can produce a finite decay at all. On the
family question as posed, that is the answer.

⚠️ **But the SHAPE half of the method does not close on this dataset, and that
is the part that was supposed to carry the weight.** Every Magboltz λ(E) *rises*
monotonically with field (air 1 %: 12.8 → 22.6 mm over 167 → 367 V/cm). det3's
measured λ *falls* over most of the range (11.1 → 5.6 mm from 167 to 333 V/cm)
and then jumps back to 9.1. That is the opposite trend, and it is not a subtle
mismatch.

Two things drive the caution:

1. The depth scale is `z = v_ridge · t`, and `v_ridge` in this CSV wanders over
   **11.0 – 29.0 µm/ns** across the scan against the bundle's 36.60 µm/ns at
   333 V/cm. A λ in millimetres inherits that scatter directly.
2. λ is not a pure attachment length. Any other depth-dependent amplitude loss
   — threshold effects (the script flags its tail estimator as
   "threshold-steepened"), or transverse diffusion spreading deeper charge over
   more strips — enters the same measurement. Those make the measured λ
   **shorter** than the true attachment length, so the true attachment is
   *weaker* than 5.6 mm implies. That direction is the safe one for the verdict
   above (it stays finite) but it means the concentration cannot be read off.

**So the deliverable, correctly scoped: det3's data is inconsistent with
attachment-free transport and therefore favours an O₂/air-bearing contaminant
family over a pure-water one — but the λ(E) shape does not match the O₂ family
either, so this is weaker than the July 90/10 result was for its run, and it
constrains the FAMILY only by the existence of a decay, not by its shape.**

**Labelling, per the standing discipline:** this is family-constraint inference
from waveform data, not a measurement. It does not establish that oxygen is
present, only that the depth-dependence of the signal is not what an
attachment-free gas produces. No concentration may be quoted from it. The way
out remains an independent species assay on the gas line.

Worth noting how this sits with the rest: the freeze queue observed that
0.8 % H₂O + ~1 % air reproduces the bench v_drift (36.24 vs 36.60) better than
any water-only mixture, and that the air row is "the one that carries real
information". This is an independent observable pointing the same way — the
v_drift argument and the attachment argument now agree on the presence of air,
which is more than either could claim alone. Both remain fitted-to-data.

---

## 12. Ledger row closed — primary ionisation is right, leaving gain alone

The last open multiplicative row in §7.4. Read-only on the frozen Stage A
cluster file (`mx17_muons_500_t0.root`, 500 events, all 500 full-gap crossers,
median path 29.40 mm), normalised per event by its own drift-gap track extent
rather than by an assumed 30 mm.

| quantity | Geant4 median | Geant4 mean | literature (Ar/iso 95/5, MIP, NTP) |
|---|---|---|---|
| total ionisation | **91.2 e⁻/cm** | 107.3 | 90–100 |
| primary clusters | 23.6 /cm | 25.3 | 25–30 |
| energy deposit | 2.37 keV/cm | 2.79 | ~2.4–2.6 |
| implied W | **25.97 eV/pair** | 25.98 | ~26 |
| electrons per cluster | 3.83 | 4.14 | ~3.5 |

Every row lands. The implied W-value of 25.97 eV against a literature ~26 is
the tightest of them and is a genuine internal consistency check, since it is a
ratio of two branches Geant4 fills independently. Clusters/cm sits ~6 % below
the literature band while electrons/cluster sits ~9 % above, so the product —
which is what the chain actually uses — comes out right. Mean above median
throughout is the Landau tail, as expected.

**Verdict: primary ionisation cannot carry the ×1.6. This row is closed.**

### Consequence for the ledger

With diffusion dilution closed by the q_sum argument (an integral cannot be
diluted by spreading) and primary ionisation closed here, **every row of the
amplitude ledger is now either dead, controlled, or separated — except the
avalanche gain.** The ×0.63 charge deficit has exactly one surviving candidate.

That is worth stating plainly because it converges with an independent thread:
the demanded gain is ~4 × 10⁴ against the sim's 24 094 (×1.6), and the data's
gain-vs-HV slope exceeds the sim's by ×1.52 at ≈12 σ. **A single α(E) /
Penning-transfer error at the operating point would produce both**, and neither
number was derived from the other.

### The morning decision this sets up

Two campaigns already in flight are, between them, exactly the test — and
neither was launched for the amplitude thread:

* **T7 slope hunt** (desktop, running unattended, self-merges and shuttles;
  collect in the morning): can any Penning rP reach the data's 0.449/10 V slope?
* **Wet gain bracket** (condor 16705137, submitted tonight): does water move
  gain at all?

**Frame the morning decision as:** if the slope hunt's best rP moves gain toward
×1.6 *while* fixing the slope, the amplitude ledger closes on a single defect
and the gain campaign has its answer. If it fixes the slope at unchanged gain,
then the ledger has no surviving candidate at all and something outside the
current chain decomposition is wrong — which would be the more interesting
outcome, and the one worth being ready for.

Products: `design/report/primary_ionisation_2026-08-10.json`,
`scratchpad/overnight_2026-08-10/primary_ionisation_check.py`.

---

## 13. ⭐ Adversarial re-derivation of the sheet-screening sizing — the
## retirement CONFIRMED, by a stronger argument, with the ladder disputed

§9 named the sheet-screening retirement as the one soft spot in the elimination
chain behind the f_eff contradiction. It is the last eliminated candidate and it
is load-bearing, so it got a second, deliberately different derivation.
Prediction and falsifiers were written into the script header before it ran.

### Method — different on purpose

The original (`nTof_x17 sheet_screening_sizing.py`) tracks a sheet counter-charge
c(k,t) relaxing toward the ion's image exp(−k z(t)) and takes the net as
`img − c`, giving a ρ ladder of **0.81 / 0.86 / 0.92 / 1.02** across ρ_s =
0.5 / 1 / 2 / 5 MΩ/sq — 8 % suppression at the production point where ×4.5 was
needed.

Mine works entirely inside the object the chain already uses — the
time-dependent weighting potential, which *is* the S1 kernel:

* `Ψ_sheet(k, τ) = (S(k)/C(k)) exp(−k²τ/(ρ_s C(k)))`, the repo's own V1
  closed form: the response to unit charge sitting **on** the sheet since τ ago.
* Upward continuation into the gas is exact, because the gas is source-free
  between the sheet and the grounded mesh:
  `Ψ(k, z, τ) = Ψ_sheet(k, τ) · sinh(k(GAP−z)) / sinh(k·GAP)`.
  At k → 0 this is the parallel-plate ramp (GAP−z)/GAP; at large k it tends to
  exp(−kz), *their* image factor — so their expression is the large-k limit of
  mine, which is a real cross-check between the two methods.
* A moving charge is a sequence of dipoles switched on in turn, so the ion's own
  contribution is `Σᵢ [Ψ(zᵢ, T−tᵢ₊₁) − Ψ(zᵢ₊₁, T−tᵢ₊₁)]`. No signal theorem is
  invoked; it is superposition.
* The model's null is the same operator at z = 0.

### The falsifier caught MY error first

The pre-registered falsifier 1 (±0.10 of their 0.92) fired on the first run at
a ratio of 0.068, and negative at 0.5 MΩ/sq. **The bug was mine, not theirs:** I
had included the ion's "appearance" term q·Ψ(z₀,T). The ion is not created from
nothing — it appears with its electron, and that prompt term is the *electron's*,
already accounted as f_e. Including it made my true branch the total
electron+ion signal while the model branch was the ion's arriving charge alone,
so the two were not the same observable. Recorded because the protocol working
on its author is the reason to run it. Falsifier 2 (static-limit f_ion) passed
throughout at 0.8959 vs 0.9000, which is what localised the error to the
accumulation rather than the machinery.

### The result

| ρ_s | true/Q | model/Q | ratio | ÷ bookkeeping | screening-only |
|---|---|---|---|---|---|
| 0.5 M | 0.3044 | 0.5192 | 0.5864 | 0.9968 | 1.0000 |
| 1.0 M | 0.3896 | 0.6643 | 0.5864 | 0.9969 | 1.0001 |
| 2.0 M | 0.4807 | 0.8189 | 0.5869 | 0.9978 | 1.0010 |
| 5.0 M | 0.5917 | 1.0057 | 0.5883 | 1.0002 | 1.0034 |

The ratio is **flat across the whole ladder**, and equals
`f_ion × T_EVAL/T_ION = 0.9000 × 0.6536 = 0.5882` to within 0.3 %.

**That constant is pure bookkeeping — it is the f_ion geometric factor the model
already applies, plus the fraction of the transit inside the window.** Strip it
out and the genuinely ρ_s-dependent part, which is the only part that could be a
resistive-sheet screening effect, is **0.34 % across a factor 10 in ρ_s.**

### Why — and it is structural

The weighting potential in the gas **factorises** as
`Ψ_sheet(k,τ) × cont(k,z)`, and the continuation factor `cont` is
**time-independent** — the gas is source-free and bounded by the grounded mesh,
so the only τ-dependence anywhere is in the sheet's own boundary value. An
in-gas source therefore sees *exactly* the same temporal sheet response as an
on-sheet source, scaled by a purely geometric factor. **The sheet's relaxation
cannot distinguish an induced source from a deposited one.**

So the induced-vs-injected distinction is not a screening effect at all. What it
*is* — a k-dependent low-pass, since `cont(k,z)` suppresses high k, making the
ion's lateral image broader than a point by roughly its own height — is the T10
effect, already measured at "central share 1.0000 → 0.977 at worst". My
derivation reproduces that independently.

### Verdict, and the honest disagreement

**The retirement stands, and this is a stronger argument than the original.**
Not "the effect is only 8 %" but "there is no ρ_s-dependent effect of this kind;
the whole distinction reduces to a geometric factor the model already carries".
×4.5 needs a ratio of ~0.22; the ρ-dependent part offers 0.3 %.

⚠️ **But the two derivations disagree about the ladder, and that is not
resolved.** Theirs swings 26 % across ρ_s (0.81 → 1.02); mine is flat to 0.3 %.
Both agree the effect is far too small to deliver ×4.5 — which is why the
retirement is robust under either — but they cannot both be right about the
ρ-dependence. My reading is that `img − c` double-counts the sheet: the sheet's
response is already inside `Ψ_sheet`, so subtracting a separately-tracked
counter-charge removes it twice. I have not proven that, and I am not claiming
their number is wrong on the strength of my own derivation having already
needed one bug fix tonight. **It is a discrepancy to settle in daylight, not a
correction to make now.**

**Consequence for §9:** the soft spot named there is now much harder. The open
question — what suppresses the ion-induction term by ×4–5 — has no surviving
mechanism, and the candidate that came closest to a structural explanation fails
by two independent derivations rather than one.

Products: `design/report/sheet_screening_rederive.json`,
`scratchpad/overnight_2026-08-10/sheet_screening_rederive.py`.

---

## 14. §0a morning checklist — the 08-08→09 HV-scan chain verified, certs re-run

Doable despite the desktop blocker because the merged product had already been
shuttled: EOS `response_sim/avalanche/aval_calib_meshfield_hvscan.json`
(913 kB, 2026-08-09 10:58). Pulled to
`~/x17/response_sim/avalanche/` and certified on the laptop.

**The chain completed cleanly.** 15 points, schema `aval_calib/3`:
Ar/iC₄H₁₀ 95/5 at 460–530 V (8 points) and 90/10 at 530–590 V (7).

| certification | 95/5 | 90/10 |
|---|---|---|
| gain monotonic in V | PASS | PASS |
| survival in [0.90, 1.00] | PASS (0.9390–0.9594) | PASS (0.9938–0.9989) |
| one distinct field map per voltage | PASS (8/8) | PASS (7/7) |
| Polya θ in [0.3, 3] | PASS (1.087–1.463) | PASS (0.825–1.094) |
| i_ion template non-zero | PASS | PASS |
| f_ion in [0.85, 0.95] | PASS (0.8908–0.9113) | PASS (0.8873–0.9064) |

**All pass.** The "one distinct map per voltage" row is the one that matters
most, because it is the direct guard against the T7 voltage-label incident
recurring — each voltage genuinely used its own `meshfield_vmesh####.txt`.

Three independent cross-checks at the 490 V production point:

* **f_ion = 0.9006** from the template's own current integrals, exactly the
  parallel-plate value — confirming the documented fact that
  `mx17_aval_calib.py` still uses `ComponentConstant` ψ for the *weighting*
  field even in the meshfield branch. The through-mesh 0.9056 (S3 closeout) is
  not in this product, as expected.
* **gain 24 172** against the pooled point's 24 094 — agree to 0.3 %.
* **survival 0.9518** against T6's independent 3D transparency 0.955.

⚠️ **A cert of my own failed first and the bug was mine, again worth recording:**
the i_ion row initially read FAIL because I assumed positive currents. The
calib stores induced currents with a **negative** sign convention (both
`i_elec` and `i_ion`), which is correct — f_ion works because it is a ratio of
two negatives. The cert now tests |Σ|. Nothing was wrong with the data.

### Bonus from the same file — the sim's HV slope, independently

Fitting d ln(gain)/dV over the 95/5 scan gives **0.3106 ± 0.0033 per 10 V**
against the data's 0.4487 ± 0.0093 (`hv_slope/slopes.json`, x/p50_head) — a
ratio of 0.692, i.e. the sim's gain-vs-HV slope is **×1.44 too shallow**. That
is an independent re-measurement of the ×1.52 / ≈12 σ discrepancy the amplitude
ledger rests on (the small difference is fit range and estimator, not physics),
computed here from the raw calib rather than quoted.

---

## 15. Slope-hunt collection harness — built, self-tested, decision rule
## pre-registered

`response/validation/slopehunt_verdict.py`. One command turns
`aval_calib_slopehunt.json` into the verdict table when the desktop chain
finishes:

```bash
python3 -m response.validation.slopehunt_verdict \
    --calib  ~/x17/response_sim/avalanche/aval_calib_slopehunt.json \
    --slopes ~/x17/response_sim/hv_slope/slopes.json
```

**The decision rule is written into the script header, before the data exists**,
with four named outcomes:

* **A** — rP\* fixes the slope *and* moves gain ≥ ×1.4 → **the ledger closes on
  a single defect**; adopt rP\*, re-run Stage B/C.
* **B** — slope fixed at gain < ×1.2 → **no ledger candidate survives** and the
  chain decomposition is suspect. The interesting outcome. Explicitly forbids
  papering over it with a fitted gain factor.
* **C** — no rP reaches the data slope → Penning is not the knob; the α(E)
  error is in the cross sections or the field map.
* **D** — slope fixed but gain overshoots > ×2.2 → the two threads are not one
  defect; report that rather than splitting the difference.

Plus the standing reporting rule, carried over from β and undershoot: **rP\* is
fitted to the slope, so slope agreement is not evidence — the gain is the
out-of-sample prediction and the only part that can confirm anything.**

**Self-tested against the 08-08 HV scan** (`--self-test`), which has no rP axis:
the code path runs end to end, reproduces gain@490 = 24 172 and the sim slope
0.3106, and correctly **abstains** rather than inventing a verdict. The header
notes that an abstain on the *real* slope-hunt output would mean the merge
failed to carry the Penning setting — i.e. the abstain is also a bug detector.

---

## 16. The X/Y spreading term firmed up — both terms are now established

§10 left the extra-spreading term at 2.3 σ, limited by only 257 paired
unsaturated data events. The paired requirement was the binding constraint, and
it can be relaxed without regenerating anything: the vertical-gun DIAGNOSIS
directories each carry their own data leg drawn from the same run, so pooling
them and de-duplicating by `event_id` yields **738 unique paired events (606
unsaturated)**, against 304 (257) from `t14_compare` alone.

Bootstrap 4000×, sim from `t14_compare` (shown in §8 to be selection-immune at
0.790 / 0.782 / 0.783, so it needs no pooling):

| quantity | value | significance |
|---|---|---|
| data peak Y/X | 0.9577 ± 0.0183 | — |
| data q_event Y/X | 1.2048 ± 0.0206 | — |
| sim peak Y/X | 0.7841 ± 0.0063 | — |
| sim q_event Y/X | 1.0627 ± 0.0049 | — |
| **peak Y/X, sim/data** | **0.8190 ± 0.0171** | 10.6 σ from 1 |
| **charge-partition term** | **0.8823 ± 0.0156** | 7.6 σ from 1 |
| **extra-spreading term** | **1.0777 ± 0.0254** | **3.1 σ from 1** |

**Both terms are now established.** The spreading term moves from 2.3 σ to
3.1 σ and its central value is stable (1.087 → 1.078), as is the
charge-partition term (0.897 → 0.882). The decomposition and the ~18 % total
are unchanged; what changes is that neither half can now be dismissed as noise.

⚠️ Caveat on the pooling: the legs come from directories whose θ windows are
derived from each one's own sim distribution, so they are near-identical rather
than identical selections. All are **vertical-gun** points — the tilted ones are
deliberately excluded, since they window on the inclined view and would mix
angular composition into a Y/X ratio. The residual heterogeneity is a
second-order effect on a ratio, and the data's peak Y/X from the pool (0.9577)
sits within 1 σ of the single-directory value (0.951), which is the check that
it did not distort anything.

---

## 17. Queue item 7 COMPLETE — the wet gain bracket landed, 4/4 predictions
## confirmed, and it makes the amplitude problem WORSE

Condor 16705137, 32/32 jobs, 42.4 h CPU, 490 V over a 0.015 cm gap =
32 667 V/cm, 1600 avalanches per arm.

| arm | gain mean | ± | median | vs dry | survival | attached | rel_var |
|---|---|---|---|---|---|---|---|
| A dry 95/5, rP 0.40 | 45 652 | 690 | 40 676 | ×1.0000 | 1.0000 | 0 | 0.373 |
| B +0.5 % H₂O, rP 0.40 | 39 934 | 628 | 34 998 | **×0.8748** (−6.1 σ) | 1.0000 | 2 | 0.396 |
| C +1.0 % H₂O, rP 0.40 | 35 125 | 543 | 31 071 | **×0.7694** (−12.0 σ) | 1.0000 | 0 | 0.382 |
| D dry, Penning auto | 46 392 | 715 | 42 024 | ×1.0162 | 1.0000 | 0 | 0.371 |

### The pre-registration, scored

All four predictions from `WET_GAIN_BRACKET_PREREG_2026-08-10.md` (committed
`287612b` before submission) **confirmed**:

| # | prediction | result |
|---|---|---|
| 1 | gain falls monotonically A > B > C | **CONFIRMED** — 45 652 / 39 934 / 35 125 |
| 2 | 1 % H₂O moves gain 10–40 % | **CONFIRMED** — measured **23.1 %** |
| 3 | survival ≥ 0.90 at 1 % H₂O | **CONFIRMED** — 1.0000, and `n_attached` = 0 |
| 4 | B between A and C, ~half way | **CONFIRMED** — B/A 0.875, C/A 0.769 |

Prediction 3 is worth reading carefully: **survival is exactly 1.0000 in every
arm and the attachment counter is zero.** So the 23 % gain loss is *pure
electron cooling* — water's low-energy vibrational cross sections keeping
electrons below the ionisation threshold — and contains no attachment component
at all. That is consistent with H₂O's 12.62 eV ionisation potential sitting
above both Ar metastables, so water opens no Penning channel while adding a
quencher. It is also the clean complement to §11: water does not attach at
amplification field any more than it does at drift field.

### Two cross-checks, both pass

* **Uniform vs meshfield.** Dry rP 0.40 gives 45 652 in a uniform 32 667 V/cm
  field against the T7 meshfield pooled point's 24 094 — a ratio of **1.895**,
  reproducing §0a's documented **1.85** (uniform-field table vs real map, 32.7
  vs 31.0 kV/cm effective). An independent arm of tonight's work landing on a
  number recorded two days ago from a different campaign.
* **The Penning-model delta is small.** Dry-auto vs dry-rP 0.40 is ×1.0162 —
  1.6 %, far below the 23 % water effect. Designing the confound out was still
  correct (its size was not knowable in advance), and arm D was never
  differenced against a wet arm.

### The roadmap decision — the axis STAYS

Roadmap §2 step 2 said: if water moves gain by less than a few percent, drop the
contaminant axis from the gain campaign entirely. **It moves gain by 23 %.**

> **The contaminant axis stays in the gain campaign**, and the T7 slope hunt's
> working assumption — that composition is fixed while Penning varies — now
> needs its own error bar. A 23 % gain shift is comparable to the effects the
> slope hunt is trying to resolve.

My pre-registered magnitude prediction deliberately contradicted the roadmap's
hoped-for outcome, and it is the prediction that held.

### ⚠️ The consequence for the amplitude ledger — it gets WORSE, not better

This is the part that matters most and it is counter-intuitive, so it is worth
being explicit. The amplitude thread wants **more** gain (§7: demanded ~4 × 10⁴
against the sim's 24 094). Water **reduces** gain.

| | gain at 490 V | factor the deficit demands |
|---|---|---|
| sim as it stands (dry, meshfield) | 24 094 | **×1.66** |
| if the bench gas really held 1 % H₂O | 18 538 | **×2.16** |

So a contaminant hypothesis invoked to explain the slow drift velocity would, if
true, **deepen the amplitude deficit from ×1.66 to ×2.16.** The two
fitted-to-data axes pull against each other: whatever slows the drift also
suppresses the gain, and the gain was already too low.

That is a genuine constraint rather than a curiosity, and it did not exist
before tonight. It also sharpens the slope hunt's outcome B: if Penning fixes
the slope at unchanged gain, and the gas is wet, the surviving discrepancy is
larger than the ledger currently states.

**Labelling, unchanged:** every water fraction here is **fitted-to-data**. This
measures a derivative — *would* water, if present, move gain — not an operating
point, and it is not evidence that the bench gas contained water. **Scope:** one
voltage, so it says nothing about the gain *slope*, which is the live problem.
If this is pursued, the follow-up is a two-voltage version, not a gain
correction.

Products: `design/report/wet_gain_bracket_2026-08-10.json`,
`scratchpad/overnight_2026-08-10/wet_bracket_collect.py`; fragments at
`~/x17/response_sim/avalanche/wetbracket/` and lxplus
`…/garfield_sim/jobs_wetbracket/`.

---

## 18. Triangulating v_drift, the gain deficit and the contaminant axis —
## and why "mutually exclusive" is too strong

§17's consequence can be pushed toward a consistency test between three
previously independent puzzles. It was worth doing, and the result is more
useful than the framing that prompted it, because one step of it is wrong.

**DIAGNOSIS throughout.** This leans on the T14 gain demand, which is itself an
inference conditional on the rest of the chain being right.

### Step 1 — the demanded DETECTOR gain is gas-INVARIANT

The tempting move is to say the wet hypothesis pushes the detector further past
the edge of the literature's stable band. **It does not, and the arithmetic is
unambiguous.** The signal is linear in gain, so changing the sim's assumed gas
rescales `G_sim` and the sim/data ratio together and the inferred detector gain
does not move:

| sim's assumed gas | G_sim at 490 V | sim/data would read | **inferred G_detector** |
|---|---|---|---|
| dry | 24 094 | 0.5987 | **40 243** |
| 1 % H₂O | 18 538 | 0.4606 | **40 243** |

So "the detector sits at the top of the 3–4 × 10⁴ stable band" is a **standing
concern independent of the contaminant axis**, and the wet hypothesis does not
aggravate it. Any framing that has water pushing the *detector* harder is
conflating the required *model error* with the *operating point*.

### Step 2 — what water does move is the required α(E) model error

That does grow, exactly as §17 says: **×1.67 dry → ×2.17 at 1 % H₂O.**

### Step 3 — but the measured slope error generates either one comfortably

This is where "close to mutually exclusive" breaks down. The sim's gain-vs-HV
slope is measured too shallow by **0.1381 ± 0.0099 per 10 V** (§14). A slope
error accumulates into a scale error, so the question is only *where* the sim
and data gain curves would have to agree:

| gas | required model error | curves agree at |
|---|---|---|
| dry | ×1.67 | **V₀ = 453 ± 3 V** |
| 1 % H₂O | ×2.17 | **V₀ = 434 ± 4 V** |

**Both sit inside the 425–530 V scanned range.** So neither required error is
implausible — the measured slope discrepancy produces errors of exactly this
size over entirely ordinary voltage offsets. A single α(E)/Penning defect
remains available in *both* the wet and dry worlds, and the two hypotheses are
**not** close to mutually exclusive.

### What the triangulation does buy

Something weaker but real, and worth stating as such rather than inflated:

* The wet world requires the simulation to be **right about gain only at
  434 V** — 26 V below det3's bench operating range (~460–500 V) and near the
  very bottom of the HV scan. The dry world requires it at 453 V, just under
  that range. **The wet hypothesis is in a worse position, not an excluded
  one.**
* Three previously independent puzzles — the slow v_drift, the gain deficit,
  and the contaminant axis — are now coupled by a measured number rather than
  by argument. Before tonight, the contaminant axis was free to be invoked for
  v_drift with no cost elsewhere. It now carries a cost, quantified: **every
  0.5 % of water invoked to slow the drift adds ~13 % to the gain the
  simulation must be wrong by.**
* ⚠️ Honest limit: V₀ is *derived from* the demanded gain at 490 V, so it is a
  restatement of the same measurement rather than an independent test of it. It
  cannot be checked by looking at 434 V in the existing scan, because the data
  leg carries no absolute gain — only amplitude up to a constant. Making it a
  real test needs an absolute gain calibration on the data side, which the
  HV-slope work already concluded `imon` cannot supply (46 pA of signal against
  1.126 µA of standing current).

### Where this leaves the morning decision

It does **not** pre-empt the slope hunt, and I want to be explicit that it
cannot: it is consistent with outcome A in both the wet and dry worlds. What it
adds is a reading rule for outcome A when it arrives — **if the best rP closes
the gain at ×1.67, that implicitly assumes a dry gas**, and the same rP would
leave a residual ×1.3 if the gas is in fact wet at the level the v_drift
anomaly wants. So outcome A should be reported with the gas assumption stated,
not as an unconditional closure.

---

## 19. Close-out — the desktop wall never dropped, and a caveat for the collect

**Desktop: 72 probes over 6 h (≈ 23:15 → 04:49 UTC), every one refused.** The
Tailscale re-authentication was never completed, so item 3 (T7 slope hunt) is
untouched and remains the single outstanding queue item.

**The slope-hunt product is NOT on EOS.** Checked at 04:49 UTC:
`/eos/experiment/ntof/data/x17/response_sim/avalanche/` holds
`aval_calib_meshfield_hvscan.json`, `…_pooled.json`, `…_diagnosis_grid.json`
and `aval_calib_v2/v3.json`, but **no `aval_calib_slopehunt.json`**. The HV-scan
product was collectible from EOS precisely because its shuttle had run; this
one's has not.

⚠️ **Two possibilities, and they need different actions — check which before
assuming the run failed:**

1. **The chain is still running.** 144 slices is roughly 2.5× the HV scan's 56,
   so this is entirely plausible; it was launched before 22:05.
2. **The chain finished but its EOS shuttle failed.** `shuttle_to_eos.sh`
   pushes over `rsync … lxplus:/eos/…` and **needs a live Kerberos ticket on
   the desktop**. The previous chain's ticket was recorded as valid only
   *through 2026-08-09 11:21*, so by the time this one finished the desktop may
   well have had no valid ticket. The plan already notes the shuttles are
   **warn-and-continue**, i.e. a failed shuttle does **not** fail the chain and
   may leave nothing in the log but a warning.

**In either case the merged JSON should exist on the desktop's local disk** —
the merge step runs before the shuttle.

⚠️ **Correction (2026-08-10, after the re-auth): the merge does NOT write to
`/media/ucla`.** `run_slopehunt_chain2.sh` writes it into the desktop's *repo
checkout* at `~/CLionProjects/MX17_Geant/response/avalanche/aval_calib_slopehunt.json`;
only the raw slices live under `/media/ucla`. Looking in `/media/ucla` for the
merged file — as an earlier draft of this section said to — finds nothing and
looks like a failure. Check both:

```bash
ssh desktop 'tail -40 /media/ucla/mx17_response_sim/slopehunt_chain2.log;
             ls -la ~/CLionProjects/MX17_Geant/response/avalanche/aval_calib_slopehunt.json;
             ls /media/ucla/mx17_response_sim/avalanche/results_slopehunt | wc -l'
```

### What was actually found (2026-08-10, post re-auth)

**The chain is still running** — possibility 1, not a failed shuttle. 88 of 144
slices done after ~11.4 h, four `mx17_aval_calib.py` processes at 100 % CPU,
currently on `penningRP0p65_520V`. Completed arms: `isofrac10` (90/10, auto),
`rP 0.30`, `rP 0.50` — all three voltages each — plus `rP 0.65` at 460 and
490 V. Remaining: `rP 0.65 @ 520 V`, the `rP 0.80` arm, and one further arm,
56 slices in total, so roughly 7 h at the observed rate.

**And the shuttle will fail when it gets there: `klist -s` on the desktop
reports NO Kerberos ticket.** The chain guards the rsync with `klist -s && …
|| echo WARNING`, so it will warn and continue, exactly as §19 anticipated. The
merged file will therefore be correct and local, and simply absent from EOS —
`kinit` on the desktop before the chain finishes, or rsync the merged JSON by
hand afterwards.

**Nothing else is outstanding.** Both watchers have resolved (wet bracket
collected, §17; desktop timed out, here), no jobs are running anywhere, and
both repositories are clean.
