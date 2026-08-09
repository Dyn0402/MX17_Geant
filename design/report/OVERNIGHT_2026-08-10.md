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
