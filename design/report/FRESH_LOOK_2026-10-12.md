# Fresh look at the sim-vs-data waveform mismatch — where we left off (2026-10-12)

The T14 closure campaign was parked on 2026-08-11 to make room for MPGD26
(`DEEP_DIVE_2026-08-11.md`). Nothing on the simulation side has moved since
then; the September and October commits here are syncs. This note sums up the
state in one page and sets out how to restart without re-running the old
campaign's habits.

## 1. The comparison that fails

Target: **det3 cosmic bench, Ar/iC4H10 95/5, 490 V mesh / 1000 V drift over
30 mm**. The default was frozen before the comparison: ρ_s = 2 MΩ/sq, dry 95/5
table, det3 bundle as analysed. The chain runs Geant4 → drift → mesh → Garfield
avalanche (T7 meshfield calib) → S1/W2 induction kernels → ion template → DREAM
shaper → FEU → `sim_decoded_*.root`. The data leg goes through the same `wft`
reconstruction, unmodified.

These numbers include both 2026-08-11 corrections (transparency
double-count, rise-estimator censoring):

| discrepancy | sim | data | size | status |
|---|---|---|---|---|
| **Amplitude** (q_sum, X) | 0.734 of data | — | **×1.46 low**; needs gain ≈ 3.5×10⁴ vs sim 24 094 | open; a pure gain fix is physically possible, but only just |
| **Gain slope** d lnA/dV | 0.311 / 10 V | 0.449 (T14 HV test); **0.419–0.428** over 425–505 V from three independent estimators (Aug 29, see §3) | **×1.4–1.5** | open; the sharpest number |
| **Rise time** 10–90 %, X | 340 ns (vertical) | 299 ns | sim +64–68 ns at 10°/20° | open: matching it needs f_eff ≈ 0.17–0.31 against f_ion = 0.90 (×3–5) |
| FWHM, X | | | ×1.18 wide | same family as rise |
| ±1 neighbour peak delay | +93–142 ns | **+60 ns** (invariant across 6 SPS arms, 2 gases) | ~×2 | points at Stage B sharing |
| X/Y peak ratio | | | sim 18 % low on Y/X (10.6 σ) | Stage B/C |

## 2. What has been excluded (do not re-open without new evidence)

**Slope and amplitude:** Penning transfer (it can give the gain or the slope,
never both; outcome C). Effective gap (closes the gain, ×1.69, but only 9 %
of the slope). T6 field-map shape (mesh vs uniform slope 0.3106 vs 0.3107).
Water and iso fraction (both make it worse). Primary ionisation (measured,
91 e⁻/cm). Mesh transparency (fixed, ×1.1455, now single-counted).

**Rise:** track inclination; t0 and frame; the rise estimator (fixed);
resistive-sheet screening, retired twice by independent derivations (the
sheet's time response factorises out) and a third time by a direct COMSOL
time-domain solve, which finds the ion signal 11–20 % *larger* than the chain
assumes, not smaller (§8); ion mobility within ±30 %
(Garfield-measured transit ~340 ns, a rectangle).

**Self-inflicted defects found on the way.** These are the main lesson. The
old campaign found **seven harness or bookkeeping defects that each moved a
headline number**: the LUT window at 1 µs, a mislabelled voltage scan, a
duplicated transparency, the railed-waveform rise, a bundle overwritten by a
variant run, a cert harness 101 ns out of alignment, and collect grouping
that would merge arms. Every "elimination" in the ledger is only as good as
the code that produced it.

## 3. New since the campaign was parked (found in nTof_x17)

- **The data slope is now solid** (`mpgd26/slides/HANDOFF_hv_window.md`,
  commit d8d1101, 2026-08-29). Over 425–505 V on the 27 June det3 scan, three
  estimators that share nothing but the events agree: peak sample 0.419,
  deconvolved charge 0.4185 ± 0.0023 and the model-free window sum
  0.4280 ± 0.0024 per 10 V. The slope is **curved**, ~0.37 per 10 V over
  425–450 V and ~0.50 near 490–500 V. At the T14 operating point the
  data/sim slope ratio is therefore ~1.6, not 1.44. This removes any hope
  that the slope discrepancy is a data-estimator artifact.
- **The wft kernel was refit (r06, c2 pinned), 2026-08-19**
  (`mx_june_wft/R06_GATE_2026-08-19.md`). The T14 data leg was reconstructed
  on the pre-r06 bundle. That is fine for the waveform-level rise and peak
  comparisons, which do not go through the kernel, but q_sum is a fitted
  quantity, so it depends on the bundle. Re-quote q_sum on r06 before quoting
  ×1.46 again.
- **`quality_ok` is an amplitude cut in disguise**
  (`mx_june_wft/quality_investigation/INVESTIGATION_2026-08-12.md`: a
  pedestal-weighted χ² with exponent 2.0). **T14 applies it to both legs**
  (`mx17_sim_wft/t14_compare.py:197`). The legs differ ×1.46 in amplitude, so
  the same cut keeps different populations on each, which is not
  like-for-like. Re-run the comparison without `quality_ok`, or with the
  gain-normalised shape χ² that replaced it, before trusting any quantile
  offset. (The w→angle constant regression fixed the same day affects
  angles, not waveforms.)

## 4. Proposal for the fresh look

The old approach was to decompose the chain and run one pre-registered test
per hypothesis. It excluded everything inside the model and ended with "the
assumption all the pieces share". A fresh look should **step outside the
chain** and get data to tell us the answer more directly. In order, cheapest
and most decisive first:

0. **Before any of it (minutes):** re-run `t14_compare.py` on the
   TRANSPFIX leg with `quality_ok` dropped (§3). This uses the existing
   parquets and needs no new simulation. If the rise offset moves, the old
   campaign's elimination table is also suspect.
1. **Independent toy forward model (hours, desktop).** Use none of the
   `response/` code: a parallel-plate gap with a prompt electron spike plus a
   rectangular ion current of transit T and fraction f, the manual's DREAM
   filter (PZC plus ζ = 0.75 Sallen-Key, 180 ns 5→100 %), FEU sampling, and a
   straight track drifting at v with diffusion. Scan (T, f, t_peak) and find
   which combinations give the data's rise and FWHM quantiles at 0°, 10° and
   20°. If the toy needs f ≈ 0.2 as well, the problem is physics or inputs,
   not chain code. If the toy matches the data at f = 0.9, there is a chain
   bug, and its location is wherever the chain departs from the toy.
2. **Ask the data for its impulse response (hours).** wft already fits a
   single-cluster temporal template to det3 data: the calib bundle's
   templates, plus `tau_s`, `sigma_s`, c1 and c2. Compare that template
   directly with the sim's (ion ⊗ shaper) impulse response, and run
   `wft/calibrate.py` on the sim's decoded output to get the sim's own
   template. Then compare template to template. That is a like-for-like test
   of the response with tracking, diffusion and selection removed. It also
   puts the ±1 delay (data +60 ns) and X/Y on the same footing.
3. **Garfield literature benchmark for α(E) (days, condor).** This was item 1
   of the post-Prague list. Run the same pipeline as a 128 µm parallel-plate
   bulk Micromegas in Ar/iso 95/5 and compare it with published gain curves.
   If it reproduces the literature slope (~0.45/10 V), look for a
   det3-specific cause. If it under-slopes, the cross sections are at fault,
   and both amplitude and slope become a known Magboltz limitation to be
   calibrated rather than derived.
4. **Bench measurements that settle anchors** (short, needs someone at the
   bench). Fire the FEU pulser once to give the absolute ADC/fC scale (the
   only absolute-amplitude anchor is a datasheet chain). Run a DREAM
   peaking-code scan to test f_ion: at f ≈ 0.9 the rise follows the ion
   rectangle, at f ≈ 0.2 it follows the shaper. Optionally, run an HV ladder
   with known mixtures from the Bronkhorst logs.
5. Only after 1–2: Stage B sharing (±1 delay, X/Y) and the r06 re-quote of
   q_sum.

Rules carried over: pre-register before running; do not tune ρ_s, the gap or
water to close T14 (all three are fitted-to-data axes); check amplitude and
shape separately; certify against code that already passes a closed-form
check.

## 5. Where everything is (checked reachable 2026-10-12 from the desktop)

- Plan and history: `design/RESPONSE_SIM_PLAN.md` §0a, `design/report/`
  (start with `DEEP_DIVE_2026-08-11.md`, `OVERNIGHT_2026-08-10.md`).
- Comparison harness: `~/PycharmProjects/nTof_x17/mx17_sim_wft/`
  (`t14_compare.py`, `t13_reco.py`, `T14_CAMPAIGN_2026-08-09.md`, `hv_slope/`).
- Stage B/C products, T14 comparisons and the data leg: **laptop** (`ssh yoga`)
  `~/x17/response_sim/stageB_w2/`. This includes `w2_rho2M_TRANSPFIX`, the
  DIAGNOSIS legs, `t14_compare_censored/` and `data_leg_satdet3/`. Not on the
  desktop.
- Kernels, avalanche calibs and decoded sim: EOS
  `/eos/experiment/ntof/data/x17/response_sim/` (`s1_w2_ny512/`, `avalanche/`,
  `stageB_w2/`). The laptop cannot reach `/eos/experiment`; go through lxplus.
- Desktop: `/media/ucla/mx17_response_sim/` (avalanche raw, meshfield ladder).
- Trap: `w2_rho2M/default/calib_bundle/` is the v-table **variant**, not the
  frozen default. Use `data_leg_satdet3/calib_bundle/`.

## 6. Step 0 result (2026-10-12): the quality cut hid a third of the amplitude deficit; the rise time is unaffected

**What was run.** `t14_compare.py` was run unmodified on the TRANSPFIX sim leg
against the `data_leg_satdet3` data leg. Both event tables were copied with
`{x,y}_quality_ok` set True everywhere, so the cut is a no-op and nothing else
changes. Outputs are on the laptop under
`~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/`
(`t14_compare_transpfix_noq/`, `run_noq.sh`). The amplitude ledger was
re-derived read-only from the reco tables with
`scratchpad/2026-10-12/quality_cut_amplitude.py`. Rise quantiles come from
`scratchpad/2026-10-12/quality_cut_rise.py`. **Both scripts reproduce the
frozen numbers exactly with the cut in** (q_sum X 0.7344, peak X 0.6462,
rise X 1.106).

**What the cut was doing.** `quality_ok = χ²/dof < 300` with pedestal-only
weights, so χ²/dof ∝ amplitude². On the vertical point the cut removes **27.8 %
(X) and 48.8 % (Y) of data tracks**, but **1.0 % / 0.5 % of sim tracks**,
because the sim is dimmer. In the top q_sum quintile it keeps 17 % (X) and 1 %
(Y) of data. In other words, the frozen comparison removed most of the
brightest data events and almost none of the sim's.

**Amplitude ledger, q_sum median within ±3° of the sim's median angle:**

| selection | data kept X / Y | X ratio | **X deficit** | Y ratio | **Y deficit** |
|---|---|---|---|---|---|
| `quality_ok` (frozen) | 0.72 / 0.51 | 0.7344 ± 0.0075 | ×1.36 | 0.6929 ± 0.0096 | ×1.44 |
| no quality cut | 1 / 1 | 0.6111 ± 0.0063 | **×1.64** | 0.4836 ± 0.0074 | **×2.07** |
| gain-normalised χ² < data p95 | 0.95 / 0.95 | 0.5967 | ×1.68 | 0.4692 | ×2.13 |
| single track candidate, no χ² cut | 0.95 / 0.92 | 0.6180 | ×1.62 | 0.4858 | ×2.06 |

The three cut-free or fairly-cut selections agree with each other. The
waveform-level ratios move the same way: peak X 0.646 → 0.523, Y 0.603 →
0.431. These are lower bounds, because without the cut **44 % (X) / 41 % (Y)
of data peak strips hit the ADC ceiling**, against 11 % / 7 % for the sim.

Consequences:
- **The uncorrected X deficit is ×1.64, not ×1.36**, about 20 % worse. The
  frozen ×1.46 was the saturation-corrected figure, so it now needs redoing
  on the uncut sample. A pure gain fix would need about 4×10⁴, above the
  stable Ar/iso bulk-Micromegas band (3–4×10⁴), so "physically available, only
  just" no longer holds.
- **The Y/X asymmetry grows.** The Y deficit is now 1.26× the X deficit, up
  from 1.06×. The 18 % X/Y term (OVERNIGHT §5) was measured under the same
  cut, which removes 49 % of Y but only 28 % of X on the data side, so it
  needs redoing too.
- **The HV-slope comparison** (`mx17_sim_wft/hv_slope/`) has to be checked
  for the same cut. Its quality-pass fraction collapses with HV, so the cut
  would bias the slope. The 29 Aug threshold-free estimators (0.42–0.43 per
  10 V) do not use the cut, and they agree with the T14 HV test's 0.449,
  which suggests the slope test was not badly affected. Not verified.

**Rise time: unchanged.** X median rise sim/data goes 1.106 → 1.101 and Y
1.162 → 1.175. The sim−data quantile offsets move by ≤7 ns everywhere
(X p5..p90 +103 +102 +83 +31 −57 −69 ns). Railed pulses were already
censored from the rise, so adding them back does not move it. A second
fairness check, keeping each leg's own dimmest 55 % so that both lose the
same top fraction, gives the same offsets (X +108 … −73, Y +94 … +1).
**The rise contradiction is not a selection artifact.** Note the vertical
offsets are not a rigid shift: the sim is slower at low quantiles and
*faster* at p75–p90, i.e. the data's rise distribution is broader. That is
the "anomalously broad vertical data leg" already flagged in
`RISE_SATURATION_CENSOR_2026-08-11.md`.

**Not yet redone with the cut removed:** the angled-closure f_eff numbers
(10°/20°). Those data legs went through the same cut, but given the vertical
result the rise is expected to move only a few ns.

## 7. Step 1 result (2026-10-12): the chain is right, the vertical rise is the wrong observable, and the inclined data pins shaper ⊗ ion only

Scripts: `scratchpad/2026-10-12/` (see its README). The toy is
`toy_pulse.py`, numpy only, with no imports from `response/`. Waveforms are the
real sim and data peak strips with ±3 neighbours, dumped on the laptop with no
quality cut (`waveforms_pm3.npz`, `waveforms_angled_pm3.npz` in
`~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/`). Anchors checked on
the way:

- **60 ns sampling is right.** The run config says `sample_period 60`, and the
  DREAM cfg's `WrClk_Div 6` on the 100 MHz base agrees.
- **The shaper family doesn't matter for a point charge.** CR-RC², CR-RC⁴ and
  the manual's PZC + ζ = 0.75 Sallen-Key, each scaled to 180 ns (5→100 %),
  give a 10–90 % rise of 114–118 ns.

### 7.1 The toy confirms the chain, so the rise mismatch is not a chain bug

At the chain's settings (f_ion 0.908, a 340 ns ion rectangle, 180 ns
shaper), the toy is *slower* than the chain: vertical X rise median 368 vs
339 ns, 20° 265 vs 231 ns. That is expected, because the chain also has the
β = 0.75 PZC residual undershoot and a distribution of ion birth heights.
Neither the toy nor the chain gets near the data at any angle. The chain is
therefore not hiding a bug that makes it slow.

### 7.2 Vertical tracks: the rise mostly measures where the maximum lands

For a vertical track the peak strip sees a ~770 ns flat-topped input, and
"10–90 % of max" depends mostly on where along that top the maximum falls.
In data the 10 % onset is stable (sample ~4–6), while the peak sits 1–15
samples after it. The 10→50 % edge tracks the peak position in both legs
(X data: 76 ns at peak sample 6, 271 ns at sample 18).

**At matched peak delay (peak − 10 % onset), the data's leading edge agrees
with the full ion model and not with the no-ion one:**

| delay 350–450 / 450–550 / 550–700 ns | X edge 10→50 % (ns) |
|---|---|
| data | 159 / 184 / 210 |
| chain (f 0.91) | 154 / 189 / 219 |
| toy f 0.91 | 164 / 184 / 208 |
| toy f 0 | 126 / 141 / 141 |

The vertical rise mismatch is a difference in **delay distribution**: 32 %
(X) / 49 % (Y) of data pulses peak within 350 ns of onset, against 4 % / 15 %
for the chain and 0 % for the toy. So the vertical-point f_eff inference
mixed two things. Removing the ions made the sim peak earlier, which
mimicked the early-peaking data population.

What makes data pulses peak early (data shapes from `wf_sumshape.py`):
- **The peak strip passes late charge to its neighbours faster than in the
  sim, strongly in Y.** Data Y peak strip is 0.37 of peak at +480 ns vs the
  chain's 0.66, while the ±3 sums agree better (0.70 vs 0.82).
- **The summed charge is mildly front-loaded.** In X it is 0.57 vs 0.68 of
  peak at +600 ns, i.e. some loss with drift depth.
- **Drift velocity × gap is consistent.** The summed signal ends at the same
  time in both legs: both reach 0.19 at +840 ns.
- Not modelled: the chain's Y channels carry a −0.2 negative lobe (peak strip
  and neighbours) that the data does not have (data −0.11, sum ~0).
- Not explained: the data has a +0.04–0.08 early foot at −360…−240 ns on the
  peak strip and +0.06–0.16 on the neighbours, absent in the sim. It is also
  why data edges at long peak delays (230–263 ns) are slower than every model.

### 7.3 Inclined tracks: the clean observable, and what it does and doesn't fix

At 10° and 20° each strip sees only ~110 / ~55 ns of drift, so the
waveform is close to the single-electron response. The early-peaking escape
does not apply here (the peak-delay distributions are narrow), and the
mismatch is real:

| rise median (ns) | 0° | 10° | 20° |
|---|---|---|---|
| data | 294 | 180 | 165 |
| chain | 328 | 252 | 231 |
| toy f 0.91, 340 ns, 180 ns shaper | 370 | 279 | 265 |
| toy f 0 | 303 | 176 | 144 |

The averaged 20° data pulse is **narrower than the no-ion toy**. It is 0.50 of
peak 300 ns after the leading 50 % point and 0.02 by 540 ns. The CR-RC² no-ion
toy is at 0.63 and 0.15, the chain at 0.82 and 0.07. A front-loaded ion
current (weighting potential Ψ = (1 − z/g)^p, p ≈ 5–8) fixes the rise (20°:
168 vs 165) but leaves the tail too long.

Fitting the shaper peaking time, shaper family, ion fraction and transit
against the averaged 10° and 20° shapes (`wf_peak_scan.py`, rms residual in
units of peak):

| shaper | f_ion | T_ion | rms 20° / 10° | rise 20° / 10° (data 166 / 181) |
|---|---|---|---|---|
| manual (PZC + SK), 180 ns | 0 | — | **0.021 / 0.023** | 144 / 174 |
| CR-RC², 123 ns | 0.91 | 170 ns | **0.027 / 0.019** | 161 / 181 |
| CR-RC², 180 ns | 0 | — | 0.092 / 0.097 | 144 / 176 |
| manual, 180 ns | 0.91 | 340 ns | 0.183 / 0.172 | 256 / 269 |
| CR-RC², 180 ns | 0.91 | 340 ns | 0.233 / 0.225 | 266 / 278 |

**The inclined data pins only the convolution shaper ⊗ ion current, and that
convolution is about as narrow as a 180 ns shaper by itself.** "Full ions
with a 340 ns transit through a 180 ns shaper" (the chain) is excluded at
any shaper family. But "no ions" is not the only way out: **a ~120 ns shaper
with full ions over a ~170 ns transit fits equally well.** Waveforms alone
cannot separate the two, so the previous campaign's "f_ion suppressed ×3–5"
was one corner of a degenerate family. The others are:

1. **The shaper is faster than 180 ns.** The peaking code was decoded from a
   register word whose word order was *inferred*
   (`response/dream/shaper.py`: "only the reading … makes both sensible").
   Code 1 is 123 ns. One FEU pulser shot, or reading the register back from
   the FEU, settles it.
2. **The ion transit is ~2× shorter than 340 ns.** That would need either a
   *(Superseded by §7.5: ×2 is not plausible, ×1.3–1.5 is.)* Needs either a
   smaller effective gap or higher ion mobility. Mobility ×2 (iC4H10⁺ in Ar
   at ~125 Td, rather than the `IonMobility_Ar+_Ar` table S3 uses), or an
   effective gap ×0.7, both enter as T ∝ g²/μ. Note S3 drives its ions with
   the **Ar⁺ in Ar** table. Ar⁺ suffers resonant charge exchange and is slow;
   the real drifting ion is isobutane-like after charge transfer.
3. A combination of 1 and 2.

Both the ion mobility and the decoded peaking code were taken as inputs the
old campaign did not question. Its exclusion list covered "mobility within
±30 %", but not a wrong ion species.

### 7.4 Next

- **Peaking code.** Read the DREAM register back, or check the shipped cfg's
  word order against the DREAM manual's register map for a *second*
  register whose meaning is unambiguous. If that is not possible, fire the
  FEU pulser once. Either result breaks the degeneracy in §7.3.
- *(Done, see §7.5.)* **Ion mobility.** Rerun S3 at 490 V with an isobutane-ion mobility table
  (or Ar⁺ scaled ×1.5 and ×2 as a bracket) and read T_ion and the 20° shape.
  This is cheap: one voltage, the meshfield map, a few hundred seeds.
- **Amplitude (§6) and rise may connect through the shaper.** With a faster
  shaper and a long ion signal, the ballistic deficit grows, and the peak
  amplitude per unit charge falls. That would make the sim's amplitude
  deficit *worse*, not better, so it does not unify them. Keep the two
  threads separate unless a number says otherwise.
- **The vertical observable** should be retired as an f_ion handle. Use 10°
  and 20°, matched by angle, with no quality cut.

## 7.5 Ion speed alone cannot close it; a 123 ns shaper with a correctly identified ion can (2026-10-12)

**No Garfield rerun needed.** S3 drifts its ions with `AvalancheMC` and no ion
diffusion, so v = K(E)·E along fixed field lines. Scaling K by s at every
field rescales the ion current exactly, I(t) → s·I(st), and the transit
becomes T/s. A uniform scaling of the Ar⁺ table therefore *is* the requested
"Ar⁺ ×1.5 / ×2" bracket. The only thing a rerun would add is a different
*shape* of K(E/N), and over the narrow field range of the gap (~124 Td,
higher only in the funnel) that is second order.

**Scan** (`wf_ion_T_scan.py`, same fit as `wf_peak_scan.py`: f_ion 0.908,
rms of the averaged 10° / 20° peak-strip shape in units of peak; data rise
166 / 181 ns):

| T_ion (ns) | ×K vs Ar⁺ | manual 180 | CR-RC² 180 | manual 123 | CR-RC² 123 |
|---|---|---|---|---|---|
| 340 | 1.0 | 0.183 / 0.172 | 0.233 / 0.225 | 0.133 / 0.115 | 0.164 / 0.151 |
| 260 | 1.3 | 0.109 / 0.107 | 0.169 / 0.169 | **0.037 / 0.032** | 0.069 / 0.069 |
| 227 | 1.5 | 0.085 / 0.086 | 0.148 / 0.150 | **0.036 / 0.030** | 0.037 / 0.042 |
| 200 | 1.7 | 0.067 / 0.070 | 0.133 / 0.137 | 0.056 / 0.046 | **0.020 / 0.024** |
| 170 | 2.0 | 0.051 / 0.056 | 0.120 / 0.125 | 0.080 / 0.066 | 0.027 / 0.019 |
| 110 | 3.1 | 0.029 / 0.035 | 0.102 / 0.108 | 0.116 / 0.099 | 0.059 / 0.041 |
| 80 | 4.3 | 0.024 / 0.028 | 0.097 / 0.103 | 0.129 / 0.110 | 0.071 / 0.050 |

(For reference, the no-ion 180 ns manual shaper scores 0.021 / 0.023.)

- **At the decoded 180 ns shaper, full ions need T ≲ 110 ns**, i.e. an ion
  3–4× faster than Ar⁺ at 124 Td. No plausible ion does that. So with a
  180 ns shaper the ions really must be suppressed, as August concluded.
- **At 123 ns (code 1), full ions fit with T ≈ 200–260 ns**, i.e. K ×1.3–1.7
  over Ar⁺ (manual shape best at ×1.3–1.5, CR-RC² at ×1.7–2). Even at
  123 ns, the Ar⁺ transit (340 ns) still fails (rms 0.12–0.16). The two
  corrections are both needed and point the same way.

**The August "species: 4 % slower" check does not stand.**
`response/meshcell/ion_template_check.py:193` applies Blanc's law to two
*different* ions: Ar⁺ in Ar (95 %) and C8Hn⁺ in iC4H10 (5 %). Blanc's law
combines the *same* ion's mobility in each component gas. In Ar/iC4H10 95/5
the Ar⁺ hands its charge to isobutane within ~ns (IP 10.6 vs 15.8 eV,
k·n_iso ~ 10⁹ s⁻¹), so what drifts is an isobutane-derived ion. Its mobility
*in argon* is the one that matters, and it is not slowed by resonant charge
exchange the way Ar⁺ is. The polarization limit for a mass-58 ion in Ar is
K₀ ≈ 2.2. Measured hydrocarbon-ion mobilities in Ar-rich mixtures are
typically ~1.6–1.9 at low field, against Ar⁺'s 1.21 at the gap's 124 Td.
That is ×1.3–1.6, **exactly the range the 123 ns fit wants**. This is an
estimate, not a measurement: the Coimbra Ar–iC4H10 data (Cortez et al.,
JINST 14 (2019) P04010, 8 Torr, 10–45 Td) is the place to pin the low-field
value, and there is no high-field (~124 Td) measurement for this ion in Ar
that I could find.

**Where this leaves the degeneracy.** It is no longer two equally free
corners. It is:

1. **180 ns shaper:** requires an unexplained ×3–5 suppression of the ion
   signal (the August reading; no mechanism found in August either).
2. **123 ns shaper + isobutane-ion mobility ×1.3–1.5:** each input
   independently physically motivated, no suppression needed.

Option 2 is the more economical explanation, but it rests on the
peaking-code word order. **The shaper readback / FEU pulser is now the single
deciding measurement.** If the shaper is 180 ns, option 1 stands and the
search for an ion-suppression mechanism reopens. If it is 123 ns, set the
chain to code 1, swap the ion table for an isobutane-ion-in-Ar mobility
(×1.3–1.5 bracket until measured), and rerun T14 at 0°/10°/20°.

## 7.6 The shaper is 180 ns: read back from the target run's own cfg (2026-10-12)

**Closed.** The degeneracy in §7.3/§7.5 resolves to the 180 ns side.

- **The target run did archive its cfg.** August's audit said the Saturday
  scan "archived no .cfg", but every subrun has one on EOS:
  `cosmic_bench/june_tests/Run/mx17_det3_saturday_scan_6-27-26/<subrun>/raw_daq_data/CosmicTb_MX17.cfg`,
  plus a copy under `dream_run/`. Of those 50 files, all have identical DREAM register
  lines, including `long_run_resist_490V_drift_1000V` (the T14 data leg; its
  file differs only in `DaqRun Time`).
  Data FEUs 7 and 8: `Feu * Dream * 1 0x081F 0xD023`.
- **Word order, proven from the register itself.** Read low word first
  (0xD023081F), the same register would give 123 ns, but also: CSA
  bypassed (external input = SK filter), positive input polarity with a
  negative-polarity VICM, 5 % discriminator range, HIT off. The FEU 1
  override `0x881F 0xD043` would then leave its peaking time unchanged, even
  though bits 4–7 of the second word are exactly the nibble the cfg comments
  and the commented-out variants edit. That is not a working front end. Read
  high word first (0x081FD023), every field is coherent: negative polarity
  with VICM 1.35 V, CSA in, 17.5 %, ThDAC −127, 5 µs CSA, **180 ns**
  (FEU 1: 50 µs, 283 ns).
- Full decode, sweep over all archived cfgs, and the decoder:
  `design/DREAM_SETTINGS.md`.

**Consequence for the ion question (§7.5 table).** With the 180 ns shaper
the inclined shapes need an ion transit ≲ 110 ns (×3–4 Ar⁺ mobility),
or the ion signal largely absent. No plausible ion species gives that. The
August reading (ion component suppressed, f_eff ≈ 0.2) stands as the only
fit at the true shaper setting. It is **still without a mechanism**, and the
vertical-track argument for it (§7.2) is withdrawn. What is left to test:

- The shaper's *actual* response may differ from its nominal Table 9 value
  (chip spread, PZC trim). This is now the only unmeasured link, and one FEU
  pulser shot through the test capacitor would measure it directly. That is
  a different question from the register setting, which is settled.
- An ion-signal suppression mechanism specific to this **resistive** bulk
  Micromegas: the sheet-screening estimate (≤ 8 %) was the only one sized in
  August.

## 8. Direct COMSOL solves: the induction stage and the true mesh (2026-10-12)

COMSOL 6.4 runs in batch on lxplus (cvmfs, CERN licence). Everything below is
in `response/comsol/` (README there has the full tables, commands and API
notes). The raw outputs are on EOS under
`/eos/experiment/ntof/data/x17/response_sim/comsol/`.

### 8.1 What was tested

The models use no weighting-field theorem. The avalanche is a prescribed
convection current (electrons down onto the ESL, ions up to the mesh over
340 ns). COMSOL solves charge conservation through the real stack: gas, ESL,
50 µm kapton, glue, pads and FR4. The pads are Terminals held at 0 V, and their
currents are the signals. Two cases are run with the same charge:

- `deposit` gives the surface kernel G.
- `pair` is the full truth, electrons plus ions.

The chain's operator, `(1 − f)·G + f·(ion ⊗ G)` with f = ⟨1 − z/g⟩ = 0.901
(`kernel_lut.apply_ion_transit`), is built from COMSOL's own G and compared
with the pair solve after the DREAM shaper (`compare_induce2d.py`). The
comparison therefore tests the operator, not the chain's kernels.

### 8.2 Results

| geometry | central pad: rise, pair / chain | peak, pair ÷ chain | fitted f_eff |
|---|---|---|---|
| 2D cut, ESL strips, 2 MΩ/sq | 227 / 228 ns | 1.03 | 0.88 |
| axisymmetric, uniform sheet, 2 MΩ/sq | 146 / 151 ns | 1.11 | 0.89 |
| axisymmetric, uniform sheet, 0.5 MΩ/sq | 129 / 129 ns | 1.20 | 0.93 |

- **2D cut:** each strip is an infinite floating line, so charge cannot spread
  or drain. The 3 % there sits at that model's numerical floor, which is the
  mismatch of the final charge between the two cases.
- **Axisymmetric:** charge really spreads off the pad. At 2 MΩ/sq only 7 % of
  the deposited charge is still on the central pad at 340 ns; at 0.5 MΩ/sq,
  1 %.

### 8.3 Reading

1. **The resistive layer does not suppress the ion signal.** This holds even
   while the sheet drains the pad during the transit. It is the third
   retirement of screening, now by a full solve and with a sign: the truth is
   larger than the chain, which is the wrong direction for the data
   (f_eff ≈ 0.17–0.31).
2. **The excess is T10's known failure, seen by an independent solver.** The
   fast path freezes the ion's weighting at the k → 0 continuation 1 − z/g.
   The exact lift is sinh(k(g − z))/sinh(kg) (`solver/zextend.py`). At
   pad-scale k it has a steeper gradient near the anode, so the true ion
   current is front-loaded: more of it arrives early, onto the central pad,
   before the sheet spreads.

   August measured this as a 7.55 % shaped residual on W2
   (`T10_SLOWPATH_2026-08-08.md`). It is "a redistribution in time",
   confined to avalanche-sized deposits. The fix (plan §7 step 5, the LUT
   from slow-path templates) was never done.

   COMSOL matches every signature:
   - **Sign:** the peak is higher and the rise faster.
   - **ρ_s trend:** 1.11 → 1.20 as the sheet relaxes faster.
   - **Neighbours:** p1/p2 come out 2–10 % lower relative to p0.
   - **2D:** the effect disappears, because G is flat after landing.
3. **Size of the rise-time effect:** the exact physics makes the rise 0–5 ns
   faster, with the peak time within −8 / +4 ns, against the 64–68 ns the
   data need. It is
   in the right direction for the rise but adds ion amplitude. It is a
   correction worth making (T10 step 5), not the missing mechanism.
4. **Not caused by the 10× slowed electrons** of the default runs: G rebuilt
   from real-velocity (150 µm/ns) electrons gives 1.11 / 1.20.

### 8.4 The true woven mesh

`WovenMesh.java` builds the real plain weave: sinusoidal wires over and under,
9.5 µm radius, the same 1 µm kiss overlap. It is solved on a 67 µm cell whose
walls are the weave's exact mirror planes, so it must be loaded with
**mirror** periodicity. `woven_map.py` exports it on the production grid with
the production flag semantics; all 1.79M flags are cross-checked in numpy,
with zero disagreements. Through the same `gates_check.C` (new `mirror`
argument), ROOT 6.40, 2000 electrons:

| | straight wires (production) | true weave |
|---|---|---|
| amp field, z = −90 µm | 31 014 V/cm | 30 743 V/cm (−0.9 %) |
| electron transparency | 0.960 ± 0.004 | **0.918 ± 0.006** |
| funnel rms | 94.6 µm | 94.3 µm |

The weave costs 4.2 ± 0.7 % of primaries. Adopting it would move the
amplitude deficit from ×1.46 to about ×1.52, slightly worse. It does not touch
the gain slope (the field-map shape is already excluded there, §2) or the
ion question. The map is at `response/meshcell/meshfield_woven_production.txt`
(gitignored) and has not been adopted.

### 8.5 What this changes

- **Induction-side ion suppression is now excluded by a direct solve.** The
  ×3–5 suppression the 180 ns shaper demands (§7.6) cannot come from the
  resistive layer, the dielectric stack or the pad geometry. The remaining
  untested link in that chain is the shaper's *actual* response (the FEU
  pulser).
- **T10 step 5 should be done.** It is a known, now independently confirmed,
  error in the chain. It makes the rise a few ns faster and the central-pad
  ion signal larger.

