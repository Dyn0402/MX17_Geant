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


## 9. Amplitude and gain slope: one mechanism, avalanche feedback (2026-10-13, overnight)

**Short version.** Garfield simulates a pure Townsend avalanche. If each
avalanche can seed further full-gap avalanches, the gain becomes
**G_eff = G/(1 − bG)**. This is the classic feedback form; photons from the
avalanche freeing electrons at the mesh are the likely carrier. It multiplies
the gain **and** the slope d ln G/dV by the **same** factor 1/(1 − bG).
T14's two oldest open rows (§1) are:

- amplitude: data ×1.77–1.82 above the sim at 490 V;
- slope: data ×1.52 steeper.

These are about equal, and the measured slope rises with voltage. Fitted to
the *shape* of the bench gain ladder alone, with no ADC calibration used, the
feedback model predicts the amplitude ratio T14 measured independently. It
also puts the feedback runaway (bG → 1) where the bench starts to discharge.

### 9.1 The fit

The data are the det3 27 June mesh ladder from nTof_x17
`mx17_sim_wft/hv_slope` (`HV_SLOPE_2026-08-09.md`): threshold-free
peak-strip max, M3-selected muons, fixed fiducial, rail-safe points only. A
copy is in `scratchpad/2026-10-13/`. The sim is
`aval_calib_meshfield_hvscan.json`, dry 95/5, 460–530 V, with
ln G = 10.093 + 0.3106 (V − 490)/10. Its local slope is flat; the quadratic
term is +0.003 per (10 V)².

`scratchpad/2026-10-13/feedback_fit.py` fits A(V) = c·G(V)/(1 − b·G(V)) with c
and b free, per view and quantile. Errors are floored at 2 %, the
scan1/scan2 gas drift. Figure: `figs/feedback_fit_2026-10-13.png`.

| view, estimator | points | bG(490) | gain ×1/(1 − bG) at 490 V | χ² feedback / free exp / Garfield shape (ndf) |
|---|---|---|---|---|
| X p50 | 425–490 | 0.52 ± 0.01 | **×2.07** | 32 / 39 / 799 (9 / 9 / 10) |
| X p25 | 425–495 | 0.46 | ×1.85 | 49 / 46 / 936 (10 / 10 / 11) |
| X p10 | 425–500 | 0.44 | **×1.79** | 22 / 78 / 840 (11 / 11 / 12) |
| X p02 | 425–505 | 0.39 | ×1.64 | 22 / 54 / 408 (12 / 12 / 13) |
| Y p50 | 425–490 | 0.49 | ×1.98 | 19 / 41 / 700 (9 / 9 / 10) |
| Y p10 | 425–500 | 0.42 | ×1.72 | 27 / 52 / 652 (11 / 11 / 12) |
| Y p02 | 425–510 | 0.34 | ×1.52 | 31 / 55 / 418 (13 / 13 / 14) |

- **Shape alone predicts ×1.5–2.1, centred near ×1.8. T14 measured ×1.77
  (published) and ×1.82 (rail-bias corrected) on X.** The two numbers share no
  input: one is the curvature of a relative ladder, the other is an absolute
  sim/data ratio through the whole chain and the DREAM calibration.
- **The data's ln A has positive curvature**, +0.011 to +0.014 per (10 V)²
  in a free quadratic: the local slope climbs from ~0.35 at 430 V to ~0.5 at
  485 V. Garfield's is flat.
- **No α(E) can produce that curvature at the real gain.** A Townsend form
  ln G = K·exp(−D/V) with a free scale (3 parameters) fits the curvature only
  with D/V ≈ 4–5, which implies a gas gain of **130–290** at 490 V against the
  detector's ~2–4·10⁴. Penning, gas tables, the iso fraction and the field
  map are all α(E) or geometry changes; this generalises August's "Penning
  gives the gain or the slope, never both" to every α(E).
- **The runaway is where the bench breaks down.** With the fitted bG(490) and
  the Garfield slope, bG reaches 1 at **511–520 V**. The bench discharge
  fraction (`HANDOFF_hv_window.md`) is 11 % at 500 V, 18 % at 510, 26 % at
  515, 45 % at 520 and 54 % at 525; efficiency collapses from 505 V. Feedback
  runaway is a standard breakdown mechanism.
- **Robust to the obvious choices.** The predicted gain factor at 490 V is:
  - ×1.64–2.27 with each pass alone (scan1 or scan2, X/Y, p50/p10), so the
    scan1/scan2 gas drift does not create it;
  - ×1.67–2.74 with Garfield's own slight curvature (quadratic ln G) instead
    of a straight line; the outlier is ×2.74, X p50 scan1 with only 7 points.

  The central values are ×1.8–2.0 throughout.
- **The quantile trend** (bG 0.52 at p50 down to 0.34–0.39 at p02) is not
  understood. It may be rail proximity, since the lower quantiles reach
  higher voltages, or a real event-to-event spread of feedback. Feedback beats
  a free exponential (same number of parameters) in 6 of 7 rows; a free
  quadratic fits better still but has no physics and one more parameter.

### 9.2 Which feedback

- **Ion-impact emission at the mesh** (ions free electrons when they land) is
  **excluded by the pulse shapes.** Its secondaries arrive one ion transit
  (~300–340 ns) late, and bG ≈ 0.45 would put ~30 % of the charge in a delayed
  second pulse. The data are *narrower* than the sim.
- **Photon feedback onto the mesh** is prompt (~1 ns), so it leaves every
  shape observable alone and gives each secondary the full gap gain, as the
  G/(1 − bG) form assumes.

  The magnitude is plausible. The fit needs b ≈ 0.45/24 000 ≈ 2·10⁻⁵
  secondaries per avalanche electron. A rough chain gives ~10⁻⁵ to a few
  10⁻⁵:
  - Ar excitations per ionisation ~0.5;
  - radiative fraction against Penning quenching ~0.1–0.5;
  - transmission through ~135 µm of 5 % isobutane, with absorption length
    ~100–300 µm, ~0.3–0.6;
  - mesh optical coverage 0.49;
  - stainless-steel photoemission yield at ~10 eV, 10⁻³–10⁻²;
  - extraction in gas ~0.5.

  Garfield cannot do this calculation out of the box. Its de-excitation model
  covers Ar levels only, its optical data has n-butane but not isobutane, and
  it does not emit electrons from electrodes.
- **Literature:** over-exponential Micromegas gain curves at *low* isobutane
  fractions, attributed to secondary avalanches from avalanche photons
  ([EUDET memo 2006-10](https://www.eudet.org/e26/e28/e183/e262/eudet-memo-2006-10.pdf);
  [LAPP Micromegas gain study](https://lappweb.in2p3.fr/~chefdevi/Publi/Cnrs/publication_1.pdf)).
  Garfield++ gains below measured ones for MPGDs, with photon feedback named
  as the likely missing process beyond Penning
  ([arXiv:1606.04852](https://arxiv.org/pdf/1606.04852)).

### 9.3 What it explains, and what it does not

- **Explains, with one parameter:**
  - the amplitude deficit (×1.8);
  - the gain slope (×1.5 at 460–490 V; ×1.6 at the operating point, where
    the data slope is ~0.5);
  - its curvature;
  - August's "matching on gain still leaves the slope 1.40× short" (a voltage
    shift cannot do what a gain-dependent factor does);
  - the discharge onset.
- **Does not touch the rise time / ion question.** Secondaries are full,
  normal avalanches with the same f_ion and shape. With §8 excluding
  induction-side suppression, the ion puzzle stands alone.
- **One cross-constraint on it:** the feedback factor reproduces the
  amplitude ratio computed *with full ions* in the sim. If the data's ion
  signal were really ×5 suppressed, the data's peak per unit charge would be
  lower and the needed gain factor larger than the shape fit gives. That
  argues mildly *against* real ion suppression. It is only as good as the
  DREAM ADC/fC scale, which the FEU pulser would pin.

### 9.4 Tests (cheapest first)

1. **Bench ladder at two isobutane fractions** (mixer logged since
   2026-07-07; HV_SLOPE §3 already asks for this). Photon feedback must fall
   steeply with more quencher, so bG(490) shrinks and the 90/10 ladder bends
   much less. A pure α(E) error would not care in this way.
2. **The FEU pulser shot** also fixes the absolute ADC/fC scale, which turns
   the ×1.8 match from "within calibration" into a measurement.
3. **Chain:** adopt G_eff = G/(1 − bG) with b fitted once on the 95/5 ladder,
   and rerun T14 amplitudes. This is a gain calibration, not a shape change.

## 10. The real strips in 3D, and the shaper from the noise (2026-10-13, overnight)

### 10.1 COMSOL 3D on the real strips validates the chain's kernels

`response/comsol/Induce3D.java` puts a point avalanche on the real board:

- the 550 µm ESL strips on an 800 µm pitch, isolated by 250 µm gaps, along
  which the landed charge diffuses;
- the checkerboard pads (680 µm on 780 µm) wired into the real comb channels,
  X columns and Y rows;
- the real stack;
- no weighting theorem.

It uses quarter symmetry about the centre of a strip and of a Y pad; both
mirror planes are exact, including the checkerboard parity. Runs went to lxplus
condor (`response/comsol/condor/`), with outputs in
`/eos/.../response_sim/comsol/condor/`.

`compare_kernels3d.py` reads the chain's own W2 product
(`greens_comb_w2_rho2M_dk50um_g19um.npz`) at the matching source point: column
36, 10 µm off the strip centre. It compares channel by channel.

| | linear elements (228k) | quadratic (135k) | chain W2 |
|---|---|---|---|
| Y0 (the avalanche's row), G at 2 ns | 0.790 | 0.822 | 0.876 |
| Y0 shaped peak, COMSOL/chain | 0.907 | **0.954** | |
| X0 (the crossing column) shaped peak, COMSOL/chain | 0.960 | **0.958** | |
| shaped peak time Y0 / X0 | 156 / 281 ns | 157 / 281 ns | 156 / 281 ns |
| delay X0 − Y0 | +125 ns | **+124 ns** | +125 ns |
| delay Y1 − Y0 | +112 ns | **+107 ns** | +102 ns |
| delay X1 − Y0 | +23 ns | +25 ns | +32 ns |
| late Y0 decay (G at 300 / 1000 ns) | 0.297 / 0.154 | 0.325 / 0.183 | 0.336 / 0.189 |

- **The chain's induction kernels are right.** This is the first check by a
  solver that shares nothing with them. Timing agrees to 0–7 ns on every
  channel that carries signal, and the sheet's spreading along the strip
  (Y0's decay, Y2/Y4 filling at 0.3–1 µs) agrees to a few %. The amplitude
  agrees to ~5 %, and the gap shrinks with refinement (0.907 → 0.954). The
  residual is partly the 10 µm ESL thickness (~2 %) that the chain's
  zero-thickness sheet omits.
- **So the ×2 neighbour-delay discrepancy is not the solver.** Data
  +54–61 ns, sim +93–142 ns. That leaves the inputs, ρ_s first (§10.3).
- **The ion operator holds on the real strips too.** With pair against chain
  operator: rise 186 / 189 ns on Y0, f_eff 0.89, absolute +6 % (quadratic;
  +13.5 % linear). Real strips do not suppress the ion signal either, which
  completes §8.

### 10.2 The shaper's real time scale, measured from the target run's own noise

`scratchpad/2026-10-13/noise_shaper.py` uses the det3 27 June long run, the T14
target, both FEUs (raw decoded_root), quiet traces only. It removes the
per-channel pedestal and the per-64-block common mode, then takes the unbiased
autocorrelation over lags 0–31 × 60 ns. The model is
R(τ) = A·(h⋆h) + B·(h′⋆h′) + C: parallel and series noise through the chain's
own DreamShaper h, stretched in time by s.

| FEU | data R(1..6) | best s | effective peaking | s allowed (χ² < 1.5 min) | series fraction |
|---|---|---|---|---|---|
| 07 (X) | 0.698 0.277 −0.022 −0.154 −0.173 −0.142 | 0.92 | **~166 ns** | 0.90–0.94 | 0.80 |
| 08 (Y) | 0.704 0.289 −0.010 −0.146 −0.167 −0.135 | 0.92 | **~166 ns** | 0.92 | 0.78 |

The model reproduces R(k) to ~0.01 out to lag 11, including the undershoot,
and CR-RC² gives the same s.

**The shaper behaves like its nominal 180 ns to within 8 %, measured in situ.**
It is not the 123 ns that full ions with a ×1.3–1.5 ion mobility would need
(§7.5). This closes the "chip spread / PZC trim" escape hatch of §7.6 from the
data alone, as far as a time scale can. The FEU pulser is still the direct
measurement, and still sets the absolute ADC/fC scale (§9.4).

**What that does to the ion requirement.** The §7.5 scan, rerun at 166 ns
(`scratchpad/2026-10-13/wf_ion_T_scan166.py`, same toy and fit; data rise
166 / 181 ns at 20° / 10.7°):

| 166 ns | shape rms 20° / 10° | rise 20° / 10° |
|---|---|---|
| full ions, T = 340 ns (Ar⁺) | 0.168 / 0.156 | 248 / 263 |
| full ions, T = 200 ns (×1.7) | 0.040 / 0.045 | 188 / 206 |
| full ions, T = 140 ns (×2.4) | **0.010 / 0.017** | 164 / 187 |
| f_ion = 0.2, T = 340 ns | 0.048 / 0.052 | 143 / 173 |
| f_ion = 0, T = 340 ns | 0.028 / 0.017 | 138 / 165 |

At the measured shaper, full ions need T ≈ 140 ns (×2.4 Ar⁺ mobility), down
from ×3–4 at exactly 180 ns. An isobutane-derived ion (×1.3–1.6, §7.5)
covers about half of that factor in mobility. **The rise puzzle is smaller
but still open.**

(The stored `dream/noise.py` autocorrelation "1.00 0.99 0.92 0.81 0.66 0.51" is
a different, per-trace-mean-subtracted estimator on another run, and is not
comparable to these numbers.)

### 10.3 Neighbour delay against ρ_s

The chain's W2 kernels exist at ρ_s = 0.5, 1, 2 and 5 MΩ/sq.
`scratchpad/2026-10-13/delay_vs_rho.py` uses the chain's operator on those
kernels and the 180 ns shaper. It takes Gaussian charge clouds (σ_T) centred
uniformly over a pad cell and computes the event-wise ±1 Y-row peak delay
against the peak row (log: `delay_vs_rho.log`).

| ±1 delay [ns] | 0.5 MΩ/sq | 1 | 2 | 5 |
|---|---|---|---|---|
| σ_T = 200 µm | 36 | 55 | 79 | 119 |
| σ_T = 400 µm | 25 | 34 | 41 | 40 |
| σ_T = 800 µm | 10 | 10 | 9 | 6 |

- The delay is mostly a **cloud-size** observable. Once σ_T ≳ 400 µm the
  neighbours get charge directly and ρ_s hardly matters.
- The data's +60 ns is **invariant across drift fields of 95–243 V/cm and
  two gases** (DEEP_DIVE §5), which change σ_T a lot. In this model that
  happens only if the neighbours see *sheet-transported* charge, the
  small-cloud limit. There, **+60 ns ↔ ρ_s ≈ 1.1 MΩ/sq**, while 2 MΩ/sq
  gives ~80 ns.
- **Caveats:**
  - the +60 ns is det4 at SPS (a different, defective chamber with its own
    ρ_s);
  - the sim's +93–142 ns comes from the full chain with its own depth
    distribution;
  - the toy has no track, only clouds.

  So this is a hint toward a lower ρ_s (or a narrower effective cloud), not a
  measurement.

**But ρ_s cannot fix the neighbour observables on its own: they pull it in
opposite directions.** The Y negative lobe is board physics. COMSOL 3D
reproduces the chain's Y undershoot exactly: deposit kernel shaped min/max
−0.283 vs −0.284 on Y0, and −0.172 vs −0.172 on Y1. With the chain operator,
a point avalanche over a Y pad, and the 180 ns shaper:

| ρ_s | Y0 undershoot | Y1 undershoot | X0 undershoot | delay X0 − Y0 | delay Y1 − Y0 |
|---|---|---|---|---|---|
| 0.5 MΩ/sq | −0.43 | −0.30 | −0.07 | +95 ns | +61 ns |
| 1 | −0.41 | −0.26 | −0.07 | +104 | +77 |
| 2 | −0.39 | −0.21 | −0.06 | +125 | +101 |
| 5 | −0.34 | −0.13 | −0.06 | +184 | +146 |

- The data need a **shorter** ±1 delay, which means *lower* ρ_s.
- They also need a **shallower** Y undershoot (−9.3 % unselected, against
  the sim's −19 % in HV_SLOPE §4–5), which means *higher* ρ_s.
- The X undershoot does not depend on ρ_s.

So the Stage B sharing tension is **not a ρ_s error and not a solver
error**. What remains has to make neighbours arrive earlier *and* drain the
central Y row less, which is beyond the board model as it stands. Two
candidates for a next look:
- capacitive crosstalk between adjacent channels in the readout or
  electronics, which is prompt and bipolar, so it shortens apparent delays
  without touching the drain;
- strips that are not fully isolated, so charge also spreads across them.

Both are untested. The T2b spread measurement put det3 at 1.42–2.56 MΩ/sq. The
  sheet resistance of det3/det4 measured with an ohmmeter on the bus bars
  would settle it.

## 11. Scoreboard after the overnight session (2026-10-13 morning)

| §1 row | status now | how |
|---|---|---|
| **Amplitude** ×1.8 low | **explained, pending test**: photon feedback G/(1 − bG) | §9: the ladder's curvature alone predicts ×1.5–2.1 |
| **Gain slope** ×1.5 shallow, and curved | **explained by the same parameter** | §9: no α(E) gives the curvature at real gain |
| Discharge onset ~510–520 V | consistent | §9: feedback runaway at 511–520 V |
| **Rise / ion fraction** | **open, narrower** | induction-side suppression excluded by full solves (§8, §10.1); the shaper is ~166 ns in situ (§10.2); full ions now need T_ion ≈ 140 ns (×2.4 Ar⁺) instead of ×3–4 |
| ±1 neighbour delay ×2 | **not the solver** (kernels validated, §10.1) and **not ρ_s alone** (it needs ρ_s down while the Y undershoot needs it up, §10.3) | open: crosstalk or inter-strip leakage are untested candidates |
| Y negative lobe, sim deeper than data | **real board physics** (COMSOL reproduces it exactly) | its depth is an input question, tied to the row above |
| X/Y peak ratio | selection artefact on the data side (HV_SLOPE §4) | unchanged |

**What is now excluded, with the method:**
- resistive-sheet screening of the ion signal: 2D, axisymmetric and real-strip
  3D COMSOL time-domain solves;
- errors in the chain's induction kernels: COMSOL 3D;
- any α(E) explanation of the slope (Penning, gas tables, iso, water,
  effective gap, field shape): curvature at real gain;
- ion-impact secondary emission: pulse shapes;
- a DREAM shaper far from nominal: the noise autocorrelation of the target run.

**Next, in order of information per effort:**
1. **Bench ladder at a second isobutane fraction** (the mixer is logged).
   This decides feedback, and fixes the chain's gain model.
2. **FEU pulser shot:** the shaper response directly, plus the absolute
   ADC/fC scale. It turns §9's ×1.8 match into a measurement and closes §10.2
   for good.
3. **Measure det3's bus-to-bus resistance** (ρ_s) with an ohmmeter. It pins
   the T2b band, though §10.3 shows ρ_s alone cannot fix the neighbour rows.
   The same bench session could pulse one channel and record its
   neighbours, which measures crosstalk directly.
4. **Ion mobility of the actual drifting ion in Ar/iC4H10 at ~124 Td** (the
   Coimbra data, Cortez et al. 2019, is at low field). It decides whether the
   ×2.4 of §10.2 is reachable by ion species alone.
5. **Chain:**
   - T10 step 5 (exact z-lift of the ion weighting; a few % and a few ns);
   - adopt G/(1 − bG) once item 1 confirms it;
   - optionally the woven-mesh map (−4 % transparency).

### 10.4 One time scale in the data, but two mechanisms

On vertical tracks (`scratchpad/2026-10-13/xtalk_fit.py`, event-aligned on the
central strip's leading 50 %), the data are faster than the chain by the
**same factor** in two quantities:

| | data | chain | ratio |
|---|---|---|---|
| central peak time, X | 180 ns | 270 ns | 0.67 |
| ±1 neighbour delay, X | +90 ns | +135 ns | 0.67 |
| central peak time, Y | 165 ns | 240 ns | 0.69 |
| ±1 neighbour delay, Y | +75 ns | +120 ns | 0.63 |

- **It is not the sampling clock.** The FEU takes its clock from the trigger
  interface (`Main_Conf_ClkSel TrgIfConClk`, WrClk_Div 6, so 60 ns only if
  that clock is 100 MHz). A 0.67 compression would need ~90 ns samples. Two
  independent anchors exclude that:
  - the noise ACF (§10.2) gives shaper/period = 166/60, so ~65 ns at the
    nominal 180 ns analog shaper;
  - the drift-column time, 757.5 ns over det3's 27.9 mm gap, gives
    36.8 µm/ns, matching wet Magboltz (36.6) only at 60 ns.
- **It is not one physical cause either.** On the validated kernels
  (`delay_vs_rho.py --ion`, ρ_s = 2 MΩ/sq), faster ions or a weaker ion signal
  barely move the ±1 delay:
  - σ_T = 200 µm: 79 / 81 / 83 ns for T_ion 340 / 200 / 140 ns, and 91 / 84
    ns for f_ion 0.2 / 0;
  - σ_T = 400 µm: 40–48 ns throughout.

  Whatever fixes the central rise will not fix the neighbour delay, and ρ_s
  cannot fix the neighbours alone (§10.3). **The rise and the neighbour rows
  are separate open problems**, and the shared 0.67 is either a coincidence
  or a third effect common to both.
- **Crosstalk can't be told apart from extra prompt sharing here.** Fitting
  the data's ±1 average as a·chain_nb + k·data_c + k′·d/dt data_c improves
  the rms ×2–3 (X 0.038 → 0.012, Y 0.059 → 0.027). But it needs
  k ≈ 0.165, a 16 % in-phase copy of the central strip. That is large for
  electronics crosstalk and equivalent to simply more prompt charge sharing
  (a wider effective cloud), so these averages cannot separate the two.
  Pulsing one channel on the bench and recording its neighbours would.


## 12. No more bench time; moisture, and a second gas found in the archive (2026-10-13)

**Constraint from now on (Dylan, 2026-10-13): the detectors are no longer
available.** Every bench item in §9.4, §11 and the HANDOFF (second-iso
ladder, FEU pulser, ohmmeter ρ_s, single-channel crosstalk pulse) can no longer
be done. Only archived data and simulation remain.

### 12.1 What the gas actually was (latest analyses, not the June ones)

The June reading "det3 = ~1 % H₂O + ~1 % air, attachment λ ≈ 16–18 mm"
(`mx_june_cosmic_qa/DET3_WEEKEND_ANALYSIS.md` §4, `RESPONSE_SIM_PLAN.md` §0a
"prefers O₂-like attachment") is **superseded**. nTof_x17
`sps_beam_test_26/analysis/cloud_basics/FINDINGS.md` (2026-10-09/11, on the
`worktree-sps-attachment-transfer` branch, not yet on main) shows that the
June decay came from normalised, threshold-aligned stacking (§10). On
unbiased, trigger-aligned stacks:

| setup | water | O₂ | evidence |
|---|---|---|---|
| bench det3, 27 June (the T14 run) | **0.95 %** | **≲ 2–10 ppm (none)** | 6-field drift scan 35–382 V/cm, plateau flat ±2 % to 2 µs (§21, §25) |
| bench det2/4/7 | 0.6–0.95 % | ≤ 0.005–0.025 % air | single-field fits (§22, §25) |
| SPS det4, Ar/CF₄/iso 88/10/2 | 1.5–1.6 % | **150–230 ppm**, varies over hours | beam-gas fits at 3 fields (§14, §18, §25) |
| SPS det4, Ar/CO₂/iso | ~1.6 % | 270–385 ppm | run_56 (§18) |

Open in that study: the bench widths look *dry* (fitted D_T 463–507 µm/√cm,
matching dry Magboltz), while v(E) says ~1 % water (D_T should then be
~250–300 µm/√cm) (§10, §27).

### 12.2 What moisture does to the chain (T14 ran dry, no attachment)

| stage | effect of ~1 % H₂O (bench) | effect on T14 rows |
|---|---|---|
| drift v | 39.1 → 35.9 µm/ns at 358 V/cm (matches the data's ~35–37) | negligible for shapes; see below |
| D_T, D_L | D_T ×0.63 (Magboltz); data widths look dry | see the toy below |
| attachment | none on the bench (12.1) | none; T14's η = 0 was right for the bench |
| gas gain | ×0.77 at 1 % (wet bracket, cluster 16705137; electron cooling, survival 1.0) | **amplitude deficit grows ×1.8 → ×2.3**; slope untouched at one voltage |
| ion species | H₂O's IP 12.6 eV is above iso's 10.6 eV and its proton affinity below isobutene's, so the drifting ion stays an iso fragment; any water clustering would *slow* it | wrong direction for the rise |
| photon feedback (§9) | VUV absorption length of 1 % H₂O ≈ 4–8 mm and of 150–230 ppm O₂ ≳ 10 cm, against ~0.1–0.3 mm for 5 % iso; under ~5 % change in transmission over the gap | contaminants neither create nor suppress feedback; they cannot be the test |

**The inclined-track rise, with wet transport** (`scratchpad/2026-10-13b/wf_wet_scan.py`,
the §10.2 toy and fit with v, D_T, D_L as parameters, 166 ns shaper; the dry
arm reproduces §10.2 exactly):

| gas arm | full ions, T 340 ns: rise 20° / 10° | T 140 ns | f_ion 0.2 |
|---|---|---|---|
| dry | 248 / 263 | 164 / 187 | 143 / 173 |
| wet 1 % (v and D) | 247 / 258 | 162 / 183 | 140 / 166 |
| wet v only | 249 / 265 | 166 / 190 | 145 / 176 |
| wet D_T only | 246 / 256 | 160 / 179 | 138 / 162 |
| **data** | **166 / 181** | | |

**Moisture moves the rise by 1–7 ns against an 80 ns gap.** It does not
change §10.2's conclusion: full ions still need T_ion ≈ 140 ns. On vertical
tracks the slower wet drift lengthens the pulse, the wrong way for the
common 0.67 factor of §10.4. Bench attachment, the one contaminant effect
that would have pulled vertical peaks earlier, is absent.

**Net:** moisture explains nothing open. It makes the amplitude ledger
~25 % harder (×2.3 needed against §9's shape-predicted ×1.5–2.1). That still
sits inside the DREAM ADC/fC uncertainty, but the pulser shot that would have
pinned it is gone. The SPS O₂ is real but confined to det4 at H4; §27 there
shows it barely moves the sharing constants (≤ 7 % except the unused
prompt/delayed split).

### 12.3 Archaeology: the bench ran other gases in April and May

From `Cosmic_Bench_DAQ_Control` git history; data on EOS
`/eos/experiment/ntof/data/x17/cosmic_bench/det_{0,1,3,4}_old/`; det_3_old is
also on the desktop at `~/x17/cosmic_bench/det_3_old/`. **These are the
pre-June chambers.** The May "det3" (`mx17_1` in its config, "Bulked 5-5-26,
strip resist vertical, ESL no silver paste") is **not** the June det3 of T14
(bulked by Stephan on 15 June).

| run | chamber | gas | ladder | size | analysed before? |
|---|---|---|---|---|---|
| `mx17_det0_He_HV_Scan_4-1-26`, `_4-2-26` | det0 | **He/C₂H₆ 96.5/3.5** | 420–525 V, drift 800 | 65 + 223 GB | efficiency only |
| `mx17_det0_He_stats_4-6-26` | det0 | He/C₂H₆ | overnight | 45 GB | no |
| `mx17_det0_Ar_Iso_HV_Scan_4-7-26` | det0 | Ar/iso 95/5 | 410–520 V, drift 800 | 48 GB | no |
| `mx17_det1_Ar_CF4_HV_Scan_4-25-26` | det1 | Ar/CF₄ (config mislabelled Ar/Iso) | 530, 550 V only | 4.5 GB | no |
| `mx17_det3_HV_Scan_5-5-26` | old det3 | Ar/iso 95/5 | 460–530 V, drift 900 | 37 GB | efficiency only |
| `mx17_det3_ArCF4_HV_Scan_5-6-26` | old det3 | **Ar/CF₄ 90/10** | 530–670 V, drift 900 | 49 GB | efficiency only |
| `mx17_det4_ArIso_HV_Scan_5-7-26` | old det4 | Ar/iso 95/5 | 440–530 V | 51 GB | efficiency only |

Two **same-chamber, two-gas** pairs exist: old det3 (Ar/iso → Ar/CF₄, one day
apart) and det0 (He/ethane → Ar/iso). Neither replaces the second-iso ladder
exactly, but both change the quencher.

### 12.4 First look: old det3, Ar/iso vs Ar/CF₄ — the feedback signature is much stronger without isobutane

`scratchpad/2026-10-13b/`: `extract_may.py` reuses hv_slope's `view_rows()`
unchanged (threshold-free peak-strip max, M3-good rays, so the population is
the same at every voltage). The fiducial was chosen once from both gases'
top-voltage efficiency maps (x −200…100, y −200…50 mm). `ladder_may.py`,
`fit_may.py` (log `fit_may.log`), figure `figs/may_two_gas_ladder_2026-10-13.png`.

Checks first:
- The scans ran top-down in voltage. The Ar/CF₄ 650 V point came 7.5 h
  after the gas change, ran 2 h and is stable (p25 by quarter
  1492/1494/1438/1460), so the gas had settled.
- The Ar/iso 510 V point, the first run after install, drifts +25 % within
  itself and is excluded.
- Points within ~4× of the noise-max floor (~35–40 ADC) are excluded, as
  are railing quantiles (> 0.7 × rail).

| gas, view, median | range | local slopes per 10 V | curvature per (10 V)² | χ² exp / quad / feedback (free slope) | bG at ref | runaway V* |
|---|---|---|---|---|---|---|
| Ar/iso, X | 460–490 | 0.41 0.46 0.47 | +0.016 | 1.5 / 0.13 / 0.25 | 0.26 (490 V) | 526 V |
| Ar/iso, Y | 460–490 | 0.40 0.46 0.48 | +0.021 | 2.4 / 0.10 / 0.25 | 0.32 | 522 V |
| **Ar/CF₄, X** | 590–630 | **0.46 0.48 0.65 0.82** | **+0.063** | **56 / 2.8 / 1.1** | 0.31 (610 V) | **641 V** |
| **Ar/CF₄, Y** | 590–630 | **0.38 0.47 0.65 0.81** | **+0.074** | **84 / 0.6 / 3.8** | 0.38 | **640 V** |

- **Garfield has no curvature in either gas.** Uniform 150 µm gap
  (`garfield_sim/results/*_Saclay_160m.json`): Ar/CF₄ 0.19–0.25 per 10 V,
  flat over 400–530 V; Ar/iso 0.30–0.33 per 10 V. The data slope is ×1.3–1.5
  Garfield's in Ar/iso (as on June det3) and **×2–3.5 in Ar/CF₄**.
  Ar/CF₄'s curvature is 4–5× Ar/iso's on the same chamber.
- **Same chamber, same electronics, only the gas changed:**
  - Ar/iso median 1019 ADC at 470 V, where Garfield gives 24 500;
  - Ar/CF₄ median 811 ADC at 620 V, where Garfield extrapolated (flat
    slope) gives ~1 650.

  Relative to Ar/iso, Ar/CF₄ is ~×12 above Garfield. Primary ionisation
  (~+15 %) and drift attachment (survival 0.93–0.97) do not come close.
- **The runaway lands where the chamber stops behaving.**
  - Ar/CF₄ (V* = 640 V): every quantile flattens above 640 V (p25 1273 → 1466
    → 1514 → 1763 over 640–670), and the mesh current leaves its 0.011 µA
    baseline: 0.035 / 0.064 / 0.11 / 0.18 / 0.32 / 0.50 µA at
    610 / 620 / 630 / 640 / 650 / 670 V.
  - Ar/iso (V* ≈ 520–550 V): the current jumps 0.071 → 0.14 µA between 500
    and 510 V.

  Cosmic avalanches give ~10⁻⁵ µA, so these are self-sustained currents, the
  usual companion of a feedback runaway.
- **The α(E) exclusion of §9.1 holds more strongly here.** For a Townsend
  form, curvature ≤ (s/V)(sV/ln G − 2). With s = 0.065/V at 620 V and gain
  ~10⁴, that gives ≤ 0.025 per (10 V)², against the measured 0.063–0.074.
  CF₄ attachment raises gα above ln G and lowers that bound further.
- **The spread between quantiles grows with voltage in Ar/CF₄.** The p10 slope
  *falls* (0.53 → 0.39 over 610–640 V) while the median's climbs. That is the
  same quantile pattern as June det3 (§9.1, bG 0.52 at p50 against 0.34–0.39
  at p02). A Poisson number of secondaries per avalanche does this: the low
  tail is avalanches without a secondary. It now has a candidate reading.

**Reading.** Isobutane is the gas's UV absorber (absorption length ~0.1–0.3
mm). Ar/CF₄ has none: CF₄ is transparent down to ~120 nm and scintillates
itself. Photon feedback must therefore be much stronger in Ar/CF₄, with an
earlier, sharper runaway. That is what the same chamber shows. Relative to Garfield, Ar/CF₄ is ×2–3.5 in
slope against Ar/iso's ×1.3–1.5. Its curvature is ×4–5 Ar/iso's, and its
amplitude is ~×12 further above Garfield. A gas-table (α(E)) error would
have no reason to scale with the UV absorber. This is the qualitative two-quencher test §9.4 asked
for, found in the archive. **Caveats:**
- it is a different chamber from T14's;
- the Garfield runs are uniform-field (the field-map factor ~1.9 is assumed
  common to both gases);
- the Ar/CF₄ Garfield is extrapolated 100 V past its range on a flat slope;
- it is one ladder per gas.

### 12.5 Next (archive and simulation only)

1. **det0 He/ethane vs Ar/iso** (EOS only, ~300 GB; extraction on lxplus with
   the same `extract_may.py`). Ethane absorbs the He excimer UV well, so the
   prediction is **weaker** curvature than Ar/iso. Together with Ar/CF₄
   (stronger), that brackets the quencher dependence from both sides.
2. **Garfield Ar/CF₄ 90/10 at 560–660 V** through the meshfield maps, so that
   12.4's ×12 is a like-for-like number rather than an extrapolation (condor,
   hours).
3. **old det4 Ar/iso 5-7** (51 GB on EOS): a second chamber's Ar/iso curvature,
   to compare with old and June det3.
4. **Pulse shapes in Ar/CF₄** from the same files. The drifting ion is CF₃⁺ in
   Ar rather than an iso fragment, with a different mobility, and the drift
   is ~×2–3 faster. The shaper and board are unchanged, so the inclined-track
   rise in Ar/CF₄ against the toy is a second, independent handle on the
   T_ion ≈ 140 ns puzzle (§10.2). This is the most direct archive test left for
   the rise.
5. Dropped (no detector): every bench item in §11 "Next" 1–3. The DREAM ADC/fC
   scale stays a datasheet chain.

## 13. Systematic workup: every chamber × gas in the archive (2026-10-13, evening)

Dylan asked for a full workup of det0 (He/ethane against Ar/iso), with the SPS gases as
extra points, to build one consistent story. Scripts, logs and the README are in
`scratchpad/2026-10-13c/`. Figure: `figs/s13_ladders_and_ion_2026-10-13.png`.

### 13.1 Data and method

| set | chamber | gas | scan | notes |
|---|---|---|---|---|
| det3 (June) | det3, bulked 15 June: **the T14 chamber** | Ar/iso 95/5 | 425–510 V, two passes | hv_slope mesh_ladder.csv, as in §9 |
| det3o | old det3, bulked 5 May | Ar/iso 95/5, then **Ar/CF₄ 90/10** a day later | 460–530 / 530–670 V, drift 900 | §12.4 |
| det0 | det0, April | **He/C₂H₆ 96.5/3.5** (4-2 scan; 4-1 repeat), then Ar/iso 95/5 (4-7) | 420–525 / 410–520 V, drift 800 | pulled to the desktop `/home/dylan/det0_old/` (45 GB); NVMe, since the data drive is full |
| det4 at SPS | det4 | Ar/CO₂/iso 95/3/2, Ar/CF₄/iso 88/10/2 | — | 13.5 |

- **Electronics are identical everywhere.** All April/May cfgs (`decode_dream_cfg --sweep`)
  read 180 ns peaking, 60 ns sampling and 200 fC range on data FEUs 3/4, as in June.
- **One extractor for every bench chamber:** `extract_bench.py`.
  - It reads `decoded_root` and `filtered_root`. The latter are uint16[512][32]
    channel-major, verified identical to `decoded_root` transposed.
  - Same observable as hv_slope: the threshold-free peak-strip maximum on M3-good rays.
    Bit-identical to §12.4 on old det3.
  - It also stores the ±2-strip waveforms and the M3 slopes for the shape studies.
- **det0 needs a track window.** Its X FEU has hot channels (ch 50 was the "peak" in 85 % of
  410 V events). With the window (±10 mm + |tan|·30 mm around the M3 prediction,
  alignment fitted from the data) and a noisy-channel mask, the noise floor is ~22–25 ADC.
- **det0's gain varies across the chamber** (its efficiency map follows the gain).
  - Both gases use the same best region: ref_x 50–150, ref_y −200…150.
  - Threshold-free quantiles stay unbiased when the population is the same at every
    voltage.
  - Clean points only: quantile ≥ 5× the noise floor and ≤ 0.7 × rail.
- **Time order is checked for every scan; all scans ran top-down in voltage.**
  - Excluded:
    - det0 Ar/iso 490 V: +40 % within its own 3 h, first point after the He → Ar/iso
      change.
    - det3o Ar/iso 510 V.
    - det0 He 505 V: taken ~6 h out of sequence.
    - Ar/CF₄ above 630 V: the saturation regime.
  - det0's gain also **rose during** its scans: He final 505 V, 24 h after the scan, is
    ×2–5 above the first 505 V. A gain that rises while the voltage steps down **flattens**
    the measured slope, so det0's slopes are lower bounds. A drift linear in time adds no
    curvature.
  - The He 4-1 repeat started at 525 V minutes after the gas change, so it is qualitative
    only.

### 13.2 Gain ladders: one pattern, ordered by gas

`ladders_all.py` → `fits_all.csv` / `fits_all.log`; `width_vs_v.py` → `width_vs_v.log`.

The exponential slope is quoted as a ratio to Garfield's slope for the same gas
(nTof_x17 `garfield_sim/results/*_Saclay_160m.json`, uniform 150 µm gap).
"gfb" is A = c·G/(1 − bG) with Garfield's own slope, two parameters. V* is its runaway.

| chamber, gas | Garfield slope /10 V (no curvature in any gas) | data slope ÷ Garfield | curvature /(10 V)² | gfb χ² vs exp χ² | bG at ref | V* (runaway) | width ln(hi/lo) vs V |
|---|---|---|---|---|---|---|---|
| det3 June, Ar/iso | 0.318 | **×1.26–1.32** | +0.010–0.014 | 9–23 vs 23–56 | 0.32–0.49 (490 V) | 512–526 V | 0.88 → 1.07, mild |
| det3o, Ar/iso | 0.318 | **×1.30–1.45** | ~+0.01 | ~equal (4–5 pts) | 0.18–0.31 (480 V) | 516–533 V | flat |
| det0, Ar/iso | 0.318 | **×1.42–1.51** (lower bound) | +0.02 | 1–6 vs 0.2–16 | 0.35–0.70 (470 V) | 481–503 V | flat |
| det3o, **Ar/CF₄** | 0.208 | **×3.1** | **+0.085** | 13–15 vs 30 | 0.56 (610 V) | **638 V** | **1.40 → 1.76** (X), 1.77 → 2.14 (Y) |
| det0, **He/C₂H₆** | 0.190 | **×2.9–7.2** (lower bound) | **+0.07 to +0.58** | **0.8–4.7 vs 1.3–35** | 0.60–0.83 (500 V) | **510–527 V** | **1.08 → 1.58** (p75/p50) |

- **Every Ar/iso chamber agrees.** Three chambers, bulked in April, May and June, give the
  same ×1.3–1.5 slope excess and mild curvature. The excess is not a property of T14's
  chamber; it is the gas.
- **The excess scales with the gas, not the chamber.** Old det3 gets ×1.4 in Ar/iso and ×3.1
  in Ar/CF₄ a day later. det0 gets ×1.5 in Ar/iso and ×3–7 in He/C₂H₆. Garfield
  (Magboltz, Penning as configured) has *no* curvature in any of these gases.
- **The width signature, which separates mechanisms.** The relative width of the amplitude
  distribution:
  - is flat where bG ≲ 0.5 (all Ar/iso);
  - grows in the last 10–20 V before V* in Ar/CF₄ and He/C₂H₆;
  - shows the same in 50 mm cells, so spatial gain spread does not cause it;
  - is noise-floor subtracted.

  A process that multiplies every avalanche by the same factor cannot do this: an α(E)
  error, gap sag under electrostatic pressure, a field-map shape. Space charge would narrow
  the distribution. A branching process does it: a random number of secondaries per
  avalanche, with mean bG → 1. Higher quantiles also run away first: He V* is 510 V on p90,
  515 V on p75 and 525 V on p50, and June det3 had bG higher on p50 than on p02.
- **The runaway is where each chamber stops behaving.**
  - Ar/CF₄: V* = 638 V, plateau from 640 V, dark current 0.11 → 0.50 µA.
  - He 4-2: V* 510–527 V, with the X median at 2589 ADC by 525 V.
  - He 4-1: V* ~500 V, with a plateau above it as in Ar/CF₄. That scan's gas was still
    changing.
  - Ar/iso: V* ~512–533 V, where the June bench starts to discharge (§9).
  - det0's own mesh current is dominated by a ~2–5 µA ohmic leak, so it is no help here.
- **My §12.5 prediction for He/ethane was wrong.** I expected weaker feedback because ethane
  absorbs the He VUV. It is stronger. He₂* emits at 16–20 eV, so the photons that reach the
  mesh carry a photoemission yield ~10× that at Ar's 9.8 eV, and they also photoionise
  ethane in the gas. The ordering found (Ar/iso < Ar/CF₄ ≈ He/C₂H₆) is consistent with
  photon feedback, but this ordering was not predicted in advance.

### 13.3 Pulse shapes: the Ar/iso ion signal is full but fast

`rise_all.py` → `rise_all.log`; `point_rise.py`; `he_template_fit.py` → `he_template_fit.log`.

**det0 sees almost no drift column.** In both gases:
- the ±2-strip / peak ratio barely grows with track slope (Ar/iso 0.05 → 0.15, He 0.04 → 0.07,
  against old det3's 0.12 → 0.29);
- vertical-track pulses stay above half maximum for only ~300 ns (old det3 Ar/iso: 480 ns).

Its charge comes from a few mm above the mesh, so its pulses are close to the **point
response of shaper ⊗ ion signal**. Why det0 behaves this way (drift gap? cathode?) is
unknown and is a question for Dylan. It makes det0 an unusually clean ion probe.

**The shaper-free test.** He ions are ~5–10× faster than any ion in Ar (Langevin:
K ∝ 1/√(αμ)). So det0's He/C₂H₆ mean pulse is an *empirical* point response: same chamber,
same electronics, same short column. The Ar/iso mean pulse was fitted as the He pulse ⊗
[(1 − f)δ + f·rect(0, T)/T] ⊗ a free Gaussian smear:

| det0, Ar/iso fit to He template | f_ion | T_ion | extra smear | rms |
|---|---|---|---|---|
| vertical tracks, free | **0.96** | **150 ns** | 1 ns | **0.031** (He alone 0.051) |
| ~10° tracks, free | 0.76 | 120 ns | 1 ns | 0.022 |
| fixed f 0.9, T 140 | | | 1–11 ns | 0.023–0.031 |
| fixed f 0.9, **T 340 (Ar⁺, the chain)** | | | | **0.20–0.21** |
| fixed f 0.2, T 340 (suppressed ions, August's reading) | | | | 0.059–0.067 |

**The data prefer a full ion signal with T ≈ 120–150 ns.** Full Ar⁺ ions are excluded ×7
in rms. Suppressed slow ions are 2–3× worse than fast ones. This breaks the f–T degeneracy
that §7.5 and §10.2 could not break with rise times alone.

Independent agreement:
- the toy point-response rise needs T ≈ 130–150 ns for det0's 168–175 ns (`point_rise.log`);
- §10.2 needed T ≈ 140 ns at the in-situ shaper on June det3;
- old det3 and June det3 have the same inclined Ar/iso rises: 180/165 ns against 181/166 ns
  at 10°/20°.

Three chambers and two methods give **one ion transit of ~140 ns in Ar/iso 95/5**.

**What that means physically.** (Correction 2026-10-16: MX17's gap is 150 µm, and the ions see the T6 map's 31.0 kV/cm; the K values in this paragraph assume 128 µm. In the real geometry ×2.4 over Ar⁺ needs K₀ ≈ 2.9; see `reports/Simulating Penning and ion transport.md`, field reconciliation.) The chain uses Ar⁺ in Ar. A 340 ns transit over 128 µm at
490 V needs K ≈ 1.0 cm²/V·s, Ar⁺'s high-field value, which resonant charge exchange holds
down. After ~1 ns the drifting ion in Ar/iso is an isobutane fragment. It has no resonant
channel in Ar, and T ≈ 140 ns needs K ≈ 2.3 cm²/V·s.

That is plausible but unverified. The Coimbra Ar–iC₄H₁₀ measurements are at low E/N
([LCTPC slides](https://indico.cern.ch/event/709670/contributions/3019642/attachments/1670094/2678865/Ion_Mobility_studies_in_Ar-CF4-iC4H10_mixtures_for_the_LCTPC_Collaboration.pdf),
[Cortez, RD51](https://indico.global/event/6159/contributions/50740/attachments/25481/43821/RD51_Presentation_-_Ion_mobility.pdf));
at high E/N they note declustering.

> **Correction (2026-10-13, night): the template fit above overclaims.**
>
> **det0's gap.** det0 has a **6 mm drift gap**: Feb logbook runs 69+, "new 6mm drift
> mounted", never changed after. That explains the short column, and at 800 V the drift
> field is ~1.33 kV/cm.
>
> **He is not a point response.** Magboltz (`magboltz/*.json`, run on lxplus) gives He/C₂H₆
> only ~15 µm/ns there, a 400 ns column, against Ar/iso's ~27.5 µm/ns, a 218 ns column. He
> looks point-like only because it is sparse: ~7 clusters/cm, so ~4 clusters in 6 mm.
>
> **The proper test.** `det0_toy.py` / `det0_toy_scan.py` → `det0_toy_scan.log` model both
> gases: clusters, drift, diffusion, gain spread, strips, sampling, the same selection and
> alignment. The shaper scale is fitted on He, then Ar's (f, T) is scanned:
> - **He calibrates the shaper at tp ≈ 130–150 ns** (best rms 0.028 at tp 130 for 15 µm/ns,
>   0.027 at 140 for 20 µm/ns). That is faster than the nominal 180 ns and than the 166 ns
>   from the noise ACF (§10.2). It is a new tension and unresolved; it could also mean the He
>   column or He ion is mis-modelled.
> - **At tp 130–150, full Ar⁺ ions (f 0.9, T 340 ns) are excluded**: rms 0.17–0.19 against a
>   best of 0.029–0.034.
> - **The best region is T ≈ 140–220 ns with f ≥ 0.5.**
> - **Suppressed slow ions (f 0.2, T 340) are disfavoured but not excluded:** rms
>   0.049–0.059.
> - At tp 150, f 0.2 with T 140–220 fits almost as well as full ions (0.031–0.036).
>
> **So det0 robustly excludes the chain's slow Ar⁺ ion signal, but does not by itself fix
> f.** The f–T degeneracy of §7.5 is only partly broken. "Full but fast" is preferred, not
> proven. The ×1.08–1.31 peak gain below assumes f ≈ 0.9.

**Other gases.** Old det3's Ar/CF₄ inclined rises are 154/149 ns, faster than its Ar/iso
180/165 ns. That is consistent with a fast, non-resonant ion (CF₃⁺) plus a faster drift. It
is not yet fitted, since the Magboltz Ar/CF₄ drift table is running on lxplus.

**Consequence for amplitude** (`peak_vs_tion.py`): at f = 0.91 and fixed charge, T = 140 ns
instead of 340 ns raises the peak by ×1.08 (vertical), ×1.23–1.26 (10.7°) and ×1.27–1.31
(20°), at either shaper. Part of T14's amplitude deficit is the ion speed, not the gain.

### 13.4 What SPS det4 can and cannot add

- **Its gain ladders cannot be used.** det4 is the chamber where 62 % of the area does not
  amplify, in fixed strips. Its run_66 Ar/CF₄/iso resist scan moves the leading-strip
  amplitude only ×1.2 over 780 → 580 V (EXTRACTION_2026-08-05b §11a), and every amplitude is
  ZS-censored. No slope or curvature can be read from it.
- **Its shapes are consistent with 13.3.** In both H4 gases the drifting ion is not Ar⁺:
  Ar/CO₂/iso → iso fragments, Ar/CF₄/iso → CF₃⁺. Both give the same core rise
  (176–179 / 170 ns at matched censoring, TWOGAS §5) and the same +60 ns ±1 delay. Neither
  contradicts fast non-Ar⁺ ions. Neither measures T on its own, because of ZS censoring.
- The SPS O₂ (150–230 ppm) does not move the sharing constants (§12.1).

### 13.5 The story so far

Two separate, gas-ordered effects explain T14's oldest rows. Each now rests on more than one
chamber and more than one gas.

1. **The ion signal is ~2× faster than the chain's Ar⁺ (T ≈ 140–220 ns), probably full.**
   - det0 (6 mm gap, two gases, see the correction in 13.3) excludes slow Ar⁺ at any shaper
     scale He allows and prefers f ≥ 0.5.
   - June det3 needs T ≈ 140 ns at 166 ns shaping, and old det3 has identical rises.
   - It fixes the rise and FWHM rows.
   - It raises the peak per unit charge by ×1.1–1.3.
   - It retires August's "ion signal largely absent" (f_eff ≈ 0.2) reading.
2. **Gain-dependent extra multiplication, G_eff = G/(1 − bG).**
   - Seen on all three Ar/iso chambers (slope ×1.3–1.5) and much stronger in Ar/CF₄ and
     He/C₂H₆ (×3–7).
   - It shows curvature in every case and broadening only near runaway.
   - V* lands where each chamber stops behaving.
   - It explains the slope row, the curvature and the rest of the amplitude.
   - Photon feedback is the leading candidate. α(E) errors and gap sag are excluded by
     curvature and width.

**The amplitude ledger closes within errors.** T14's ×1.8 (X, at 490 V) ≈ ×1.15–1.25 (ion
speed, mixed angles) × ×1.5 (feedback; the §9 shape fit allows ×1.5–2.1). Wet gas
(×0.77 gain, §12.2) would raise the need to ~×2.0, still inside the feedback band. The
DREAM ADC/fC scale is still a datasheet chain.

**Still open:** the ±1 neighbour delay and Y lobe (§10.3–10.4). The 0.67 common time factor
of §10.4 now has one half explained. In the toy (vertical tracks, f 0.91), going from T = 340
to 140 ns shortens the peak delay after the 10 % onset by 96 ns at 166 ns shaping
(453 → 357 ns) and by 84 ns at 180 ns. That is the size of §10.4's central gap (chain 270 ns
against data 180 ns). The neighbour delay was shown in §10.4 not to depend on T_ion, so it
remains a Stage B question.

### 13.6 Next (simulation and archive only)

1. **Chain: replace the Ar⁺ ion mobility with T ≈ 140 ns** (an iso-ion mobility, ×2.4) and
   rerun T14's shapes on the frozen default.
   - Prediction: rise, FWHM and the central peak time move onto the data.
   - Mobility is an exact time rescale (§7.5), so the decoded sim can be re-shaped without
     rerunning Garfield.
2. **Chain: adopt G/(1 − bG)** with b from the June ladder, and re-quote amplitude together
   with item 1.
3. **Old det3 Ar/CF₄ shapes** with the Magboltz Ar/CF₄ drift table (running on lxplus,
   `~/work/mx17_gas_oct13/`). This gives T_ion for CF₃⁺ as a second ion species.
4. **Feedback magnitude per gas from first principles:** a small MC of photon yield ×
   transmission × photoemission (no Garfield support), to see whether Ar/iso < Ar/CF₄ ≈
   He/C₂H₆ comes out unaided.
5. **Ion mobility of iso fragments in Ar at ~120 Td:** literature or Coimbra contact, to turn
   the ×2.4 into an input rather than a fit.
6. Question for Dylan: why does det0 collect only a few mm of drift column?

## 14. Literature check on both mechanisms (2026-10-13, night)

Slides: https://dylan-neff.web.cern.ch/notes/mx17-feedback-and-ions.html (generator
`scratchpad/2026-10-13c/make_lit_deck.py`). Primary texts were read locally; the pdftotext
excerpts are in the job tmp, not kept.

### 14.1 Photon feedback

Neither the mechanism nor the fit form is new.

- **Iguaz et al., microbulk Micromegas, Ar/iC₄H₁₀ 1–35 %, 50 µm gap**
  (arXiv:1110.2641, 1201.3012):
  - "the gain curves for low quantities of isobutane show a deviation from the Rose and
    Korff gain model at high amplification fields. This over-exponential behaviour is due
    to the low quencher concentrations, which cannot avoid the production of secondary
    avalanches by the UV photons generated in the primary avalanche."
  - Below 5 % iso the resolution worsens "as not all the UV photons … are absorbed".
- **Şahin, Kowalski, Veenhof**, JINST 11 P01003 (2016); also NIM A 768 (2014) for Ar–CO₂
  and NIM A 718 (2013) on secondary avalanches:
  - G_T = G/(1 − βG), fitted on top of Magboltz + Penning gains.
  - β falls with quencher fraction and pressure, and is linked to the quencher's
    photo-absorption cross section.
  - **Our model is their standard treatment.**
- **Peskov & Fonte review** (arXiv:0911.0463):
  - "Detectors the cathode of which are not geometrically shielded from the avalanche
    emission, for example MICROMEGAS … can operate at high gains only in quenched gases."
  - Feedback shows as a "gain curve deviation from linear behaviour in log scale" because
    the feedback pulses are prompt (τ ≪ shaping time). That is our "shape untouched,
    charge up".
  - The photocathode quantum efficiency rises steeply with photon energy. This is the
    likely reason He/C₂H₆ (16–20 eV photons) feeds back more than Ar (~10 eV), and why my
    §12.5 prediction failed.
- **Garfield++ MPGD gain shortfalls** (arXiv:1606.04852; MPGD WG4; arXiv:2505.03454) are
  common. Penning fixes much of the gap for THGEMs. Photon feedback is named as a missing
  process: Garfield does not emit electrons from electrodes.

**Not in the literature:**
- a measured β for Ar/iso 95/5 in a 128 µm bulk Micromegas;
- any first-principles β.

What is ours is the size of the effect at 5 % iso (βG ≈ 0.3–0.5 at G ≈ 2·10⁴) and the
same-chamber gas ordering.

### 14.2 Drifting-ion species

The ion that drifts is a hydrocarbon, and its measured mobility predicts the transit.

- **Cortez et al., JINST 14 P04010 (2019)**, Coimbra; Ar–iC₄H₁₀ at 8 Torr, 10–45 Td:
  - Ar⁺ charge-transfers to isobutane. The rate is taken from n-butane,
    1.4×10⁻⁹ cm³/s, i.e. ~1 ns at 1 atm with 5 % iso.
  - The drifting ions are C₄Hₙ⁺, then C₈Hₙ⁺.
  - **K₀ at 95 % Ar: 2.08 (C₄Hₙ⁺), 1.74 (C₈Hₙ⁺)** cm²/V·s, at 15 Td.
- **Garfield `IonMobility_Ar+_Ar.txt`** (what the chain and `funnel_ion_endpoints.py` use)
  gives K₀ = 1.21 at 124 Td. Resonant charge exchange slows Ar⁺ in its parent gas.
- **The chain's 352 ns is Ar⁺ with the 150 µm gap.** T = g²/(KV) at 490 V:

  | ion | 150 µm gap | 128 µm gap |
  |---|---|---|
  | Ar⁺ | 352 ns | 266 ns |
  | C₈Hₙ⁺ | 241 ns | 176 ns |
  | C₄Hₙ⁺ | 201 ns | 147 ns |

  The data band is 140–220 ns.
- **Giomataris (ICFA 1999)**, Ar + 10 % iC₄H₁₀, 100 µm gap: "the electron signal and the
  ion tail are developed within 5 ns and 80 ns"; the ion drift takes place "typically
  within 100 ns". Scaled to 128 µm and 490 V that is ~120–150 ns.
- **Not settled:**
  - hydrocarbon-ion mobility at 130–160 Td;
  - C₄ against C₈ at 1 atm (C₈ is more likely).
- **Gap correction (later the same night):** the MX17 gap is **150 µm, confirmed from the
  pillar gerber** (T14_CAMPAIGN_2026-08-09). So the predictions that apply are the 150 µm
  column: C₄Hₙ⁺ 201 ns (inside the data band), C₈Hₙ⁺ 241 ns (just above it). Giomataris'
  80–100 ns over 100 µm scales as g² to ~180–225 ns.

### 14.3 det0 gap (Feb/May logbooks)

- The Feb logbook records "new 6mm drift mounted" on the Feb detector (det0) from run 69
  on. Dylan says the drift was not changed after the tests.
- The May logbook: "detector 3" had a 30 mm drift at the start of the May test (one later
  entry says 16 mm), and det 4 had 3.6 mm.
- The old det3's 480 ns vertical pulses fit 30 mm.
- The Feb test also ran Ar/CO₂ 70/30 and He/C₂H₆ resist scans on det0 with the 6 mm gap,
  under n_TOF beam. They are more gas points, but beam-environment data.

## 15. Second ion species, and both mechanisms in the chain (2026-10-13, late)

### 15.1 Old det3: Ar/iso and Ar/CF₄ on one chamber

`scratchpad/2026-10-13c/det3o_toy.py` → `det3o_toy.log`, `det3o_toy.csv`. This is the det0 toy
generalised to a 30 mm gap and inclined tracks.
- **Free per gas:** drift velocity v and ion transit T.
- **Shared:** the shaper scale tp.
- **Fitted:** the mean X peak-strip pulse of vertical, ~10° and ~20° tracks (f_ion 0.9).

| gas | best fits (rms over 3 bins) | Ar⁺, T 352 ns | literature T over 150 µm |
|---|---|---|---|
| Ar/iC₄H₁₀ 95/5 (460–490 V) | tp 150–166: T 140–200 ns. Inclined bins alone at tp 166, T 140: 0.016 / 0.019 | 0.14–0.17, excluded | C₄Hₙ⁺ ~207 / C₈Hₙ⁺ ~250 ns at 475 V |
| Ar/CF₄ 90/10 (600–630 V) | tp 166, T 140: 0.027; tp 180, T 80: 0.029; tp 150, T 200: 0.032 | 0.13–0.15, excluded | **CF₃⁺ ~160 ns** (K₀ ≈ 2.11 in 90 % Ar, Cortez et al., LCTPC Nov 2017) |

- **A second ion species, on the same chamber, behaves as its literature mobility
  predicts.** Neither gas is consistent with Ar⁺. Old det3's inclined Ar/iso fits agree with
  June det3's T ≈ 140 ns at tp 166.
- **Caveat on the vertical bin.** The vertical Ar/iso pulses fit poorly (rms 0.09–0.12), with
  v pinned at the grid edge (40 µm/ns). The May logbook calls det 3 "16mm drift" at one
  point (30 mm at another). A shorter gap would explain it. The inclined bins carry the T
  inference; the per-strip time spread there is pitch/(tan·v) and does not depend on the gap.

### 15.2 The chain with both mechanisms

New digitizer options (commits c19081d and the one after):
- `--ion-speed-scale s` time-compresses the measured S3 ion template, i(t) → s·i(st). This
  is exact (§7.5). f_ion is unchanged, and t₉₀ goes 310 → 212 / 177 / 127 ns for
  s = 1.46 / 1.75 / 2.45.
- `--feedback-bg bG` makes each avalanche spawn Poisson(bG·g/Ḡ) secondary Pólya avalanches,
  recursively. The mean is ×1/(1 − bG), verified to 0.3 % at bG 0.3–0.8. The relative
  spread grows 0.82 → 1.43 at bG 0.45.

Condor runs on the W2 ρ₂M angled ladder (th00 / thx10 / thx20), each with a prediction and a
falsifier recorded in the submit file before it ran:
- **cluster 4450898** (`stageBC_w2_ionspeed.sub`): s = 1.46 / 1.75 / 2.45;
- **cluster 4450899** (`stageBC_w2_fullphys.sub`): s = 1.46 / 1.75 with bG 0.45.

The laptop T13/T14 pass is `run_ionspeed_t14.sh`, driven by `orchestrate_t14.sh`, with the
same arguments as `t14_angW_*`. Its s = 1 baseline was reproduced at 20°: sim rise 232 ns
against data 166 ns, peak data/sim ×2.6.

## 16. Chain results: ion speed and photon feedback against the det3 data (2026-10-13, evening)

All 15 T14 points are done: clusters 4450898 and 4450899, laptop T13/T14 pass with the
`t14_angW_*` arguments. Results are in `yoga:~/x17/response_sim/stageB_w2/t14_{ionspeed,fullphys}<s>_<ang>/`;
the table comes from `scratchpad/2026-10-13c/t14_ionspeed_table.py`.

**X view, sim/data** (rise median ns | FWHM median ns | peak-amplitude median, data/sim):

| | vertical | 10° | 20° |
|---|---|---|---|
| s=1 (Ar⁺), no fb | 326/286 · 638/524 · ×1.81 | 251/176 · 411/325 · ×2.24 | 231/162 · 392/302 · ×2.65 |
| s=1.46, no fb | 293/301 · 625/531 · ×1.49 | 213/179 · 353/317 · ×1.68 | 193/166 · 326/297 · ×1.94 |
| s=1.75, no fb | 283/301 · 617/531 · ×1.47 | 199/179 · 338/317 · ×1.62 | 179/165 · 309/296 · ×1.87 |
| s=2.45, no fb | 274/301 · 601/531 · ×1.43 | **183/179 · 322/317** · ×1.55 | **157/166 · 292/296** · ×1.76 |
| s=1.46, bG 0.45 | 304/301 · 647/531 · **×0.81** | 215/179 · 357/317 · **×0.90** | 194/166 · 329/297 · **×1.03** |
| s=1.75, bG 0.45 | 295/301 · 640/532 · **×0.80** | 201/179 · 342/317 · **×0.87** | 180/166 · 312/297 · **×0.98** |

**Y view** peak data/sim: ×1.65–1.77 (s=1) → ×1.25–1.41 (ion speed) → **×0.69–0.76** (with bG 0.45).
Y rise: 13–37 ns too fast for every s ≥ 1.46, with or without feedback.

### 16.1 Ion speed (cluster 4450898)

- At s = 2.45, the inclined X shapes close: rise and FWHM are both within 5 ns at 10° and
  20°. This meets the prediction pre-registered in `stageBC_w2_ionspeed.sub`.
- At the literature hydrocarbon speeds (s = 1.46 / 1.75), inclined X stays 14–34 ns slow.
  A shaper ~8–12 % faster than nominal would cover that (noise ACF: 158 ns on FEU 3, 166 ns
  on FEU 7/8; §10.2).
- Ion speed alone moves the amplitude only ×1.1–1.35.

### 16.2 Both mechanisms (cluster 4450899), checked against the pre-registration

- **Amplitude: closes.** At s = 1.46–1.75, X data/sim is 0.80–0.81 vertical and 0.87–1.03
  inclined. The prediction was 0.9–1.2 for vertical. The falsifier band (0.8–1.4) holds,
  narrowly. With feedback, the 3× angle-dependent amplitude gap of the s=1 chain is gone.
- **Shape: feedback is neutral, as required.** Inclined X moves +1–2 ns against the ion-speed
  twin. The censored vertical rise moves +11–12 ns, which looks like the >10 ns falsifier. It is
  the saturation censoring (§T14 deep-dive): feedback raises the sim's fraction above 3500 ADC
  from 8 % to 35–37 % (data: 26 %). On the **uncensored** rise median, fullphys sits within
  6 ns of its ion-speed twin in every view and angle (X vertical 284 vs 289 at s=1.46 and 273
  vs 279 at s=1.75; Y within 4 ns). So this is a metric artefact, not a bug.
- **New tension: X and Y disagree on the gain.** With bG 0.45, which was fitted to the X
  p25–p50 ladder, Y overshoots: amplitude ×0.69–0.76, saturated fraction 20–28 % in the sim
  against 4 % in the data. Without feedback, Y was already ~15–20 % closer to the data than
  X. Either the X/Y charge split in the chain is off by that much (Stage B/C sharing, cf. the
  X≈Y-at-5M nugget), or bG should be fitted to X and Y together. A rough estimate from the mean scaling
  1/(1 − bG) (not a fit): Y alone wants bG ≈ 0.25, a joint X+Y fit roughly 0.35.

### 16.3 What this leaves

- **Inclined X: a working description.** Hydrocarbon ion speed plus bG ≈ 0.45 matches the
  amplitude to ~10 % and the shape to 14–34 ns. The remaining shape gap matches a shaper at
  the noise-ACF value.
- **Open:**
  - (a) the X/Y amplitude ratio;
  - (b) the Y rise, which is too fast once ions speed up (vertical tracks are a weak rise
    observable, §7.2, but Y is fast at 10°/20° too);
  - (c) whether to adopt shaper 158–166 ns as the default, then run one verdict T14 with
    s = 1.46–1.75, bG fitted to X+Y, and that shaper.

## 17. The FEUs were the *protected* variant: what that can and cannot do (2026-10-13, night)

New information from Dylan: the cosmic-bench FEUs were the spark-**protected** variant.

### 17.1 What the documentation says

- `221217_FeuUsersManual.pdf` (Mandjavidze, Dec 2022):
  - §2, p. 9: the inputs reach the DREAMs "through optional protection circuits". FEUs for
    resistive Micromegas are built without them; FEUs for standard (sparking) Micromegas get
    them at production.
  - No schematic or component values are given.
  - Power register 0x200000, bits 19:4 `ProtFlt` (p. 18): the protection is floating (0) or
    grounded (1), per 32 channels. The target run's cfg sets `Feu_Pwr_PrtFlt 0xFFFF`, so it was
    grounded on every channel (`design/dream_cfg/satdet3_long_run_…cfg:177`).
- The DREAM manual covers only the chip's internal ESD diodes (§8.2).
  - The internal test capacitor (500 fF functional, 4 values in test mode) injects at the CSA
    input. It most likely **bypasses** the board protection, so a pulser shot through it would
    not see the protection network.
- Online sources (CLAS12 MVT papers, DREAM papers, the classic Micromegas reviews) give no
  values either. **Values are being requested from Irakli Mandjavidze.**

### 17.2 Quick test: the network as a filter in front of the shaper

`scratchpad/2026-10-13d/protection_rc_scan.py` → `.log`. This is the §10.2 toy: manual SK
shaper at 166 and 180 ns, full ions (f 0.908), T_ion 340 (Ar⁺) and 233 (C₈Hₙ⁺). It is fitted to
the inclined det3 X average shapes. Data: rise 166 / 181 ns (20° / 10°), **undershoot −0.022 /
−0.023** of peak.

| network | effect on 20° rise (tp 166, T 233; no network: 202 ns) | undershoot | shape rms |
|---|---|---|---|
| high-pass τ = 3 µs | 199 | −0.11 | 0.069 (none: 0.064) |
| high-pass τ = 1 µs | 194 | −0.26 | 0.155 |
| high-pass τ = 200 ns | 176 | −0.61 | 0.363 |
| high-pass τ = 100 ns | 162 | −0.73 | 0.425 |
| low-pass τ = 5–20 ns | 202–204 | −0.009 | 0.065–0.069 |
| low-pass τ = 40 ns | 209 | −0.007 | 0.082 |

- **A high-pass (series C into R to ground) is excluded as the missing mechanism.**
  - A first-order high-pass doesn't remove the ion charge. It adds a negative lobe that
    integrates it away, so speeding the rise costs a large undershoot.
  - Reaching the data's rise needs τ ≈ 100–200 ns and a 60–90 % undershoot. The data has 2 %.
  - Extrapolating the undershoot (≈ 1/τ), the data bounds τ_hp ≳ 15 µs, where the rise moves
    by ~1 ns.
  - The same bound applies to the DREAM's own CSA feedback (R_f C_f = 5 or 50 µs) and to
    anything else that differentiates.
- **A low-pass (series R into strip C)** only slows the pulse. It is harmless for R_s C ≤ 20 ns.
- **Charge division** (C_c in series with C_strip, which differs between X and Y) leaves the
  shape unchanged. It stays open as a candidate for the §16.2 X/Y amplitude disagreement and
  for part of the absolute gain, both at the ~C_strip/C_c level. It needs the values.

**Verdict:** the protection network cannot be what removes the ion signal from the inclined
shapes. The fast-pulse puzzle stays as in §16 (hydrocarbon ion speed plus a shaper at the
noise-ACF value). The network could still matter for amplitudes (X/Y, absolute ADC/fC), so the
component values are worth getting.

## 18. The T14 amplitude ledger is biased by its data-side selection; redone with M3 angles (2026-10-14, night)

Scripts and logs: `scratchpad/2026-10-14/` (README inside). Everything runs on yoga, read-only on
nTof_x17.

### 18.1 Why T14's amplitudes cannot be used to fit bG

T14 (`t14_compare.py`) selects the **data** leg with the det3 reco: quality cut
(χ²/dof < 300, ∝ amplitude²; §6) **and** a ±3° window on the reconstructed angle. The **sim** leg gets
the quality cut only. Both cuts depend on amplitude in the data:

- **Reco angle.** Below |θ| ≈ 4° the reco slope is flagged unreliable, and in **both** legs the reco
  |θ| grows with amplitude: the sim (true θ = 0) reconstructs its brightest events at 4–6°
  (`peak_vs_theta.log`: Y peak median 2266 at reco 0–0.25° vs 2884–3018 at 4–6°). A window on
  the data only, therefore, keeps the data's dim events and drops its bright ones. With the
  window applied to both legs (`t14_symwin.log`), Y moves by +5–10 %, X by ≤ 5 %.
- **Quality cut.** It still removes 21 % (X) / 32 % (Y) of data tracks against 4–7 % of sim. On
  the Y vertical data leg it removes almost every saturated event: **4 % saturated in the T14 data
  leg against 33 % in an unbiased sample**. The §16.2 "Y: 20–28 % saturated in the sim vs 4 % in
  the data" was this cut.

### 18.2 An unbiased ledger: the data by its true (M3) angle

`m3_slopes_longrun.py` adds the M3 reference-track slopes to `hv_slope/peaks_longrun.parquet`
(target run, 490 V, threshold-free peak strip, M3-good rays: 8479 events). The view angle is
+atan(tan_y) for X and −atan(tan_x) for Y (correlation with the det3 reco +0.73 / −0.71).
A fiducial box (ref_x −200…120, |ref_y| < 170 mm) removes rays that miss det3 (beyond ref_x 125,
90 % have peak < 150 ADC). Data = |view angle − a| < 3°, |other projection| < 6°, no det3 reco, no
quality cut. Sim = every simulated event (its true angle is fixed), no cut (`ledger_m3.py`).

Data peak-strip median, true angle vs the T14 selection (`peak_vs_m3.log`, `rise_vs_m3.log`):

| | X 0° | X 10° | X 20° | Y 0° | Y 10° | Y 20° |
|---|---|---|---|---|---|---|
| true angle (M3) | 2932 | 2387 | 2065 (n 70) | 2670 | 2593 | 2097 (n 33) |
| T14-like (reco window + quality) | 2566 | 2221 | 2110 | 1946 | 2154 | 1780 |
| factor | ×1.14 | ×1.07 | ×0.98 | **×1.37** | ×1.20 | ×1.18 |
| saturated, true / T14 leg | 0.38 / 0.26 | 0.30 / 0.24 | 0.23 / 0.24 | 0.33 / **0.04** | 0.33 / 0.16 | 0.27 / 0.15 |

**Data/sim with the unbiased ledger** (`ledger_m3.log`; p50 | p25; sat fraction data/sim):

| sim | X 0° | X 10° | X 20° | Y 0° |
|---|---|---|---|---|
| s=1, no fb | 2.07 \| 1.84 | 2.26 \| 2.05 | 2.58 \| 2.14 | 2.43 \| 2.22 |
| s=1.46, no fb | 1.69 \| 1.51 | 1.72 \| 1.54 | 1.90 \| 1.58 | 1.94 \| 1.79 |
| s=1.75, no fb | 1.66 \| 1.48 | 1.65 \| 1.47 | 1.80 \| 1.50 | 1.90 \| 1.76 |
| **s=1.46, bG 0.45** | **0.91 \| 0.82**; sat .38/.39 | **0.91 \| 0.82**; .30/.25 | **0.98 \| 0.82**; .23/.19 | **1.04 \| 0.95**; .33/.24 |
| s=1.75, bG 0.45 | 0.90 \| 0.80; .38/.41 | 0.87 \| 0.79; .30/.28 | 0.94 \| 0.78; .23/.21 | 1.01 \| 0.93; .33/.26 |

(s=1 inclined Y: 2.77 | 2.37 at 10°, 2.99 | 2.80 at 20°.)

### 18.3 What changes

- **The angle-dependent amplitude gap was mostly the selection.** With hydrocarbon ions, the X p25
  ratio is flat with angle (s=1.46 + fb: 0.82 / 0.82 / 0.82; no fb: 1.51 / 1.54 / 1.58). The median
  keeps a small rise at 20°, where the data have only 70 events.
- **bG ≈ 0.45 stands; the "Y overshoot" of §16.2 reverses.** At s=1.46, bG 0.45, the medians are
  X 0.91–0.98 and Y 1.04, and the saturated fractions match within 0.01–0.05 in X. Without
  feedback the chain needs ×1.66–1.94, i.e. 1/(1−bG) → bG ≈ 0.40–0.48. The p25 ratio being lower
  than the p50 ratio (0.82 vs 0.91 in X) means the **data** distribution is broader on the low side.
  bG hardly changes that (§18.7), so the width cannot split bG from an overall gain factor.
  Gain non-uniformity over the fiducial area is the obvious data-only source.
- **X/Y: the sim's Y/X is ~15 % too low**, the same at every angle: Y/X data/sim 1.14 (fb, p50),
  1.16 (p25), 1.17 (s=1). This is the August charge-partition item (OVERNIGHT §10: sim puts
  ~10 % less charge in Y). It is the opposite sign from §16.2, which was the selection. Charge
  division in the protected FEU input (§17) is still a candidate. It needs C_c vs the X and Y
  strip capacitances.
- **Rise by true angle** (T14 extractor on the M3-selected events): X 10° 187 ns, 20° 169;
  Y 10° 188, 20° 160. The T14 data legs give 176 / 162 / 180 / 160. So the inclined data is
  **7–11 ns slower** than T14 said, and the s=1.46 + fb gap at 10° X shrinks 36 → 28 ns.

### 18.4 The sim's multiplicity excess is a noise-model artefact

T14's "strips ≥ 5σ" is 12 (X) / 13 (Y) in the sim against 6 / 8 in the data. Within ±30 strips of
the peak it is 5 / 8 in both legs. About 10 sim channels (478, 480, 481, 483, 352, 27, 93, 100, …,
the **same numbers on both FEUs**) cross 5σ in ~50 % of events, far from any track
(`hot_channels.log`). The single `noise_det3_satscan_feu07.json` spec is applied to both FEUs, and
those channels' generated noise exceeds the σ stored in the decoded file. Effect: n_over and
q_event in the sim are inflated by ~6 strips × 5σ. The peak amplitude is not affected (signals
are hundreds of ADC). T13 seeding may be affected. **Not fixed.**

The real profile difference is small. X matches. Y in the sim is wider at ±2–3 strips (0.25 / 0.087
of peak vs data 0.165 / 0.048), which is the August "extra spreading" term.

### 18.5 The cause: the common-mode noise is ~50× too large in the sim (fixed)

The hot channels are the **lowest-pedestal** channels (238–262 ADC against a median of 315).
`noise_det3_satscan_feu07.json` has σ_cm = **271 ADC**. That figure is the std (and the
mean-square spectrum) over all events, and ~1 % of target-run events rail whole blocks
(|cm| ≈ 3800). The robust value is 4.5 ADC (MAD), or **5.7 ADC** over the 90 % of events within
10× of that. In the data, the per-block common mode has a 25–75 % range of −4…+3 ADC.

The model then drew a 271 ADC Gaussian common mode on every sim event:
- **11 % of all sim ADC samples sit at 0** (data 0.2 %). Clipped low-pedestal channels become
  fake hits after CNS.
- The rail moves event by event (4095 − ped − cm, i.e. ~3300–3800 after CNS).
- The FeuReader noise reads 8.5–8.9 ADC in the sim against 7.41 in the data.

Fix: `response/dream/noise.py` now takes the common-mode spectrum from quiet events. New specs are
in `design/noise_specs/` and on AFS `calib/*_cmq.json` (σ_cm 5.70 / 4.85 for FEU 07 / 08, and
selftest passes). The older `noise_det3.json` (May run, σ_cm 83) was built by the same estimator
and is probably inflated too. Not rechecked.

**A/B, condor 4451433** (`stageBC_w2_noisecm.sub`): fullphys s=1.46 / bG 0.45 at th00/thx10/thx20
with only the noise spec changed, paired with 4450899. Pre-registered: amplitude moves < 3 %,
inclined shape < 5 ns, n_over drops to the data's.
**Result (`scratchpad/2026-10-14/noisecm_ab.log`): all within the pre-registration.**
- Zero samples 11.4 % → 0.03 %. Hot channels gone: no channel fires far from the track more
  than 2 times in 600 events. FeuReader noise 7.41 ADC, same as the data.
- Ledger p50 moves −1 % (X 0.912 → 0.902, Y 1.038 → 1.024). The X saturated fraction rises
  +0.04 (0.39 → 0.43; data 0.38).
- Inclined rise +3 ns, FWHM +2 ns.

**So the bug inflated multiplicity and q_event (and possibly T13 seeding) but did not bias the
amplitude or shape conclusions above.**

### 18.6 X/Y: it is the Y lateral sharing

The chain splits the charge evenly (`x_fraction` 0.502), but the Y peak strip holds 0.19 of the
view's charge against 0.46 for X (`share_Y` / `share_X` in `stageBC_rho2M.json`). The sim's Y
profile is wider than the data's: mean amp/peak at ±1/±2/±3 is 0.60/0.30/0.11 in the sim against
0.56/0.23/0.09 in the data. X matches. That is a ~10 % lower Y peak, most of the ~15 % Y/X offset
of §18.3. It is the same Y-kernel family as the Y negative lobe and the ±1 rows (§7.2, §10.3–10.4),
and §10.1's COMSOL check did not cover the Y wings. Open.

### 18.7 bG × shaper scan: results (condor 4451431 / 4451432, tp 166 ns, old noise spec)

Amplitude (`ledger_m3`, p50 data/sim, X 0/10/20° | Y 0/10/20°; geometric mean over all six):

| | X | Y | geo-mean p50 / p25 |
|---|---|---|---|
| s 1.46, bG 0.25 | 1.32 1.32 1.42 | 1.51 1.59 1.64 | 1.46 / 1.31 |
| s 1.46, bG 0.35 | 1.14 1.12 1.22 | 1.29 – 1.41 | 1.23 / 1.11 |
| s 1.46, bG 0.45 | 0.96 0.94 1.01 | 1.09 1.13 1.15 | 1.04 / 0.93 |
| s 1.75, bG 0.25 | 1.28 1.26 1.33 | 1.46 1.51 1.54 | 1.39 / 1.26 |
| s 1.75, bG 0.35 | 1.11 1.07 1.15 | 1.25 1.28 1.33 | 1.19 / 1.07 |
| s 1.75, bG 0.45 | 0.94 0.90 0.95 | 1.05 – 1.08 | 0.98 / 0.89 |

- **bG from X and Y together: 0.43–0.46** (p50, log-interpolated). X alone gives ~0.44 and Y alone
  ~0.49. The §16.2 "Y wants 0.25" is gone: that came from the selection bias.
- The saturated fraction in X matches at bG 0.45 (sim 0.35–0.37 vs data 0.38 at 0°). Y stays low
  (0.22 vs 0.33), which is the Y-sharing offset again.
- **The 166 ns shaper lowers the sim peaks by ~8 %** (bG 0.45 twins: geo-mean 0.959 → 1.043 at
  s=1.46). The prediction was 2–6 %, so this is a little outside.

Shape (`shape_all.py`: the T14 extractor on **all** sim events, against the data by true angle; rise | FWHM, sim − data, ns):

| | X 10° | X 20° | Y 10° | Y 20° |
|---|---|---|---|---|
| s 1.46, tp 180 (fullphys, T14 sim leg) | +28 \| +36 | +25 \| +34 | – | – |
| s 1.46, tp 166, bG 0.25 / 0.35 / 0.45 | +21..+22 \| +18..+20 | +19..+20 \| +17..+18 | +18..+19 \| +19..+21 | +28 \| +20..+21 |
| **s 1.75, tp 166, bG 0.25 / 0.35 / 0.45** | **+7..+9 \| +2..+4** | **+5..+6 \| −2..0** | **+5 \| +4** | **+13 \| +3..+4** |

- **Feedback is shape-neutral:** rise and FWHM move ≤ 2 ns across bG at fixed s and tp, as
  pre-registered.
- **Shaper 180 → 166 ns:** −5 to −7 ns rise, −16 ns FWHM (predicted −6..−12 and −10..−20).
- **At s = 1.75 and tp 166 the inclined shapes close to +5..+13 ns in both views.** The §16 "Y
  rise too fast" was the Y vertical point and the selection. At matched inclination Y behaves
  like X.

**→ Verdict 1** (condor 4451436, `stageBC_w2_verdict1.sub`): s 1.75, bG 0.45, tp 166, fixed noise
spec, all five angles. Predictions are in the submit file.

### 18.8 Verdict 1 (condor 4451436): s 1.75, bG 0.45, tp 166 ns, fixed noise — passes

`scratchpad/2026-10-14/verdict1.log`. Data by true M3 angle, sim all events.

| | 0° | 10° | 20° | pre-registered |
|---|---|---|---|---|
| X peak p50 data/sim | 0.93 | **0.89** | 0.95 | 0.90–1.02 (falsifier 0.85–1.10) |
| Y peak p50 data/sim | 1.04 | 1.06 | 1.07 | 1.02–1.15 |
| X saturated, data / sim | 0.38 / 0.41 | 0.30 / 0.29 | 0.23 / 0.23 | within 0.05 |
| rise sim − data [ns] | – | X +10, Y +9 | X +8, Y +16 | +2..+12; Y20 +5..+18 |
| FWHM sim − data [ns] | – | X +4, Y +7 | X +2, Y +6 | ±8 |

- Every prediction holds except X 10° amplitude, which is 0.01 under its band (well inside
  the falsifier).
- **The sim is now ~7 % brighter than the data in X and ~5 % dimmer in Y.** That is the Y-sharing
  offset (§18.6). The overall gain is right to a few %, with no free gain factor: bG came from
  the June ladder shape, before any amplitude fit.
- **Inclined pulses are 8–16 ns slower than data, against 60–70 ns at the August default.** The
  residual is about the size of the ±1 neighbour-delay physics still open (§10.3).
- Noise check (full T14 pass at th00): sim noise is 7.41 ADC, as in the data, and n_over is X 5 vs 6,
  Y 10 vs 8. Both are inside the pre-registered ±2 (Y at the edge, where the sharing is wide).
  Before the fix they were 12 and 13. The remaining T14 passes run on yoga (`~/queue_t14.sh`) for
  continuity.

**What this configuration is.** Three physical changes to the chain: the drifting ion is a
hydrocarbon (C₄Hₙ⁺, ×1.75 Ar⁺ mobility), there is photon feedback at bG 0.45, and the shaper
runs at the noise-measured 166 ns. Plus one bug fix (the common-mode noise). The October fresh
look's question, "why are sim pulses slower and wider and the amplitude ×2 low?", is answered
at the 10 % level on the target run. The ion-speed choice within 1.46–1.75 and the shaper value
are both read from this same run, so the cross-gas test (Ar/CF₄, §16.3 item 4) is the real
out-of-sample check.

### 18.9 Scan bookkeeping

Condor 4451431 (`stageBC_w2_bgtp.sub`: s 1.46/1.75 × bG 0.25/0.35/0.45 × th00/thx10/thx20, all at
the new `--shaper-tp-ns 166`) and 4451432 (`stageBC_w2_bgtp_thy.sub`: the same at thy10/thy20,
added after §18.1). Both have predictions registered in the submit files. T14 and the ledger run
on yoga as the jobs land (`orchestrate_bgtp.sh`).

### 18.10 Y wings: the kernels are right; ρ_s is the lever (2026-10-14, morning)

- **The Y wings pass COMSOL 3D.** `compare_kernels3d.py` on `i3d_dep_2e6_o2` (ρ_s = 2 MΩ/sq), now read on
  all channels Y0..Y5. The shaped peak ratios COMSOL/chain are Y0 0.954, Y1 0.931, Y2 0.945, Y3 0.923
  and Y4 0.950, so the wings sit at the same ~5 % as the centre. Peak times are within 6–38 ns, even at
  1.2 µs. At a given ρ_s the chain's Y lateral sharing is the field solution, and the §18.6 excess is not
  a solver error (`scratchpad/2026-10-14/kernels3d_ywings.log`).
- **ρ_s moves both §18.6 symptoms the same way** (`rho_wings.py`, a shaped point source, no cloud):

| ρ_s [MΩ/sq] | 0.5 | 1 | 2 | 5 |
|---|---|---|---|---|
| Y2/Y0 (next Y channel) | 0.259 | 0.163 | 0.091 | 0.036 |
| Y0/X0 | 0.85 | 1.14 | 1.57 | 2.63 |
| X1/X0 | 0.015 | 0.017 | 0.020 | 0.030 |

  Log-interpolating, ρ_s ≈ 2.5 MΩ/sq raises Y/X by the ~15 % the ledger is missing (§18.3, §18.6). At the
  same point, the next-Y-channel wing drops by ~22 %, against the data's ±2 wing being 23 % below the sim
  (0.23 vs 0.30). X1/X0 hardly moves, which fits X matching already. 2.5 MΩ/sq is at the top edge of
  the T2b band (1.42–2.56). It pulls against §10.3's ±1-delay hint (~1.1 MΩ/sq), which is a
  cloud-size observable and only a hint.
- **Caveats:** this is a point source, so the cloud and track extent are not included, and the profile
  observable is per strip rather than per kernel channel. A higher ρ_s also lowers X0 (−7 % at
  2.5 MΩ/sq), so bG would be refitted. **Test:** a chain point at ρ_s = 5 MΩ/sq (the kernels exist)
  next to verdict 1, which brackets 2.5. Read the Y profile ±1..3 and Y/X from the M3 ledger.

## 19. Irakli's DREAM bench measurements (NSW AC-coupled, TPOT DC-coupled) (2026-10-14, morning)

Irakli Mandjavidze sent two old Saclay bench measurements (to Alexandra, for P2). The files are
Saclay-internal, so they are **not in this public repo**. Originals, text dumps and the fit script are in
`~/x17/irakli_dream_tests_2026-10-14/` on the desktop (dylan-MS-7C84).

- **Protected or not: not stated in either file.** TPOT (resistive MM, "DC-coupled", on-board pulser) is
  almost certainly the unprotected FEU. For NSW, "AC" is more likely an external coupling C than the FEU
  protection. It is still the closest analogue to our protected FEUs: a series C between strip and DREAM.
  He needed an external pulser there because the on-board pulser injects at the DREAM input, after the
  coupling C. This agrees with §17.1: the on-board pulser cannot see our protection network.
- **TPOT, DC, 100 fC range, 180 ns:** gain 18.05 mV/fC at 0 pF (nominal 20, so −10 %), 17.5 at 220 pF,
  16.7 at 470 pF, 14.7 at 1000 pF. ENC ≈ 810 e⁻ + ~10–11 e⁻/pF. Lab conditions with CM removed; his slide
  says not to infer detector capacitance from the noise.
- **NSW, AC, external pulser (6.8 fC/mV, 281–1128 pF total):** gain 16.2 ADC/fC at 281 pF and
  6.7–7.0 at 1104–1128 pF (100 fC range), against 30–35 ADC/fC DC (TPOT). The ratio to DC fits
  k·C_c/(C_c + C_strip) to ±0.01 with **C_c ≈ 510 pF, k ≈ 0.70** (correlated, rough).
  The M1/M2/M4 sheets are NSW strip crosstalk with multiplexing (8–15 % per neighbour at M1) and do not transfer.
- **For MX17:**
  1. The 20.48 ADC/fC datasheet chain is likely ≥ 10 % high at 180 ns before any strip-C or protection
     loss. This assumes the 200 fC range behaves like the 100 fC range, which nobody measured.
  2. If the protection acts like C_c of a few hundred pF, Y/X scales by (C_c + C_Y)/(C_c + C_X).
     The §18.3 15 % Y/X deficit would need C_Y ≈ 75–100 pF below C_X. MX17 X/Y strip capacitances are unknown.
  3. Charge division does not change pulse shape (§17.2 stands).
- **Ask Irakli:** were the NSW numbers from protected FEUs or an external C? Ask for the protection
  schematic/BOM. A P2 repeat should use an external pulser on the detector side of a *protected* FEU.

## 20. Three one-variable tests against verdict 1: ion mobility, ρ_s 5 MΩ/sq, corrected primaries (2026-10-16)

Condor 4472241 (`stageBC_w2_ionK.sub`), 4472242 (`stageBC_w2_rho5M.sub`) and 4472245 (`stageBC_w2_e17.sub`),
evaluated on yoga with `ledger_m3.py`, `shape_all.py` and `profile_tails.py` (`~/ledger_1016_{v1,ionK,r5,e17}.log`,
`~/shape_1016.log`, `~/profile_1016_*.log` on yoga). Predictions were registered in each .sub header before
submission. Ledger values are p50 data/sim, data by true M3 angle. Shape values are sim − data.

### 20.1 Ion mobility (`--ion-model reconstructed --ion-mobility K0=…`): passes; K0 2.05 ≈ verdict 1

| K0 | inclined rise vs verdict 1 [ns] | FWHM vs verdict 1 [ns] | p50 data/sim vs verdict 1 | pre-registered rise / p50 |
|---|---|---|---|---|
| 1.80 | +15..+16 | +18..+20 | +2.9..+7.0 % | +10..+20 / +4..+9 % |
| **2.05** | **+5..+6** | +4..+6 | **+0.9..+2.2 %** | +2..+8 / 0..+4 % |
| 2.25 | −1 | −2..−1 | −0.6..−0.2 % | −6..0 / −3..0 % |

- Falsifier (K0 2.05 > 10 ns or > 6 % from verdict 1) **not triggered**. The X saturated fraction moves
  ≤ 0.03 at every K0. FWHM at K0 2.05 is 1 ns over its 0..+5 band in three of four views. K0 1.80's
  0° p50 (+2.9 %) is just under its band. Everything else lands inside.
- So the reconstructed hydrocarbon ion and the time-compressed Ar⁺ template are equivalent, and
  verdict 1's s = 1.75 reads as **K0 ≈ 2.2** (the 2.25 point reproduces verdict 1 to 1 ns and 0.6 %).
  K0 2.05, the working literature value, is 5–6 ns slower on the inclined rise, which is already
  +8..+16 ns over the data.
- **Proposed default:** `--ion-model reconstructed --ion-mobility K0=2.05`, with 1.80/2.25 as the
  systematic (rise +10 / −6 ns, amplitude ±3–5 %). Not switched yet: it costs 5 ns of inclined-rise
  agreement for a physical parameter, and that trade is Dylan's call.

### 20.2 ρ_s = 5 MΩ/sq: the Y wings close, Y/X does not; **falsifier triggered**

| | verdict 1 (2 MΩ/sq) | ρ_s 5 MΩ/sq | data | pre-registered |
|---|---|---|---|---|
| Y p50 data/sim, 0° | 1.04 | 0.92 | – | 0.75–0.90 |
| X p50 data/sim, 0° | 0.93 | 0.88 | – | 1.00–1.15 |
| sim Y/X change | – | +1.3..+7.2 % | – | +20..+50 % |
| Y ±2 wing (amp/pk, 0°) | 0.244 | 0.172 | 0.165 | 0.4–0.7 × verdict 1 |
| Y ±3 wing (amp/pk, 0°) | 0.089 | 0.040 | 0.049 | |
| X ±1 (amp/pk, 0°) | 0.46 | 0.41 | 0.56 / 0.42 | change < 15 % |
| X 10° rise vs data [ns] | +10 | +9 | | change ≤ 6 ns |

- **The Y wings match the data near 4–5.6 MΩ/sq** (log-interpolating ±3 and ±2), not at the 2.5 MΩ/sq the
  point source suggested (§18.10). At 5 MΩ/sq the ±2 wing is 0.70 × verdict 1, at the top of its predicted
  0.4–0.7 band and within 4 % of the data.
- **Y/X barely moves.** The sim's Y amplitude goes up, but so does X's (X data/sim falls 5–9 % instead
  of rising 7–20 %), so the sim's Y/X rises only 1–7 %. Matching the data's Y/X by ρ_s alone would need
  ~8.5 MΩ/sq at 0° and ≫ 10 MΩ/sq at the inclined angles. The two interpolations disagree by far more
  than the factor 1.5 the falsifier allowed, and Y/X moved < 10 % at three of four angles.
- **Reading:** ρ_s sets the Y lateral spread and the data prefer ~4–5 MΩ/sq for it, but the ~12–20 %
  Y/X amplitude offset is a separate effect. The pre-registered alternative is the next suspect:
  charge division by the protection network and the X/Y strip capacitance (§19, item 2). The X-peak
  prediction failed in sign, so the point-source kernel picture of how X0 responds to ρ_s (§18.10:
  −7 % at 2.5 MΩ/sq) is wrong once the cloud is included.
- The thx20 point is missing from the r5 ledger (`ledger_1016_r5.log` has eight rows); nothing
  above depends on it.

### 20.3 Stage A with the e17c72c primaries: amplitude passes, joint bG lands at 0.47

| bG | X p50 data/sim (mean of 5 angles) | Y | joint (geo. mean) |
|---|---|---|---|
| verdict 1 (old primaries, 0.45) | 0.90 | 1.03 | 0.96 |
| 0.45 | 0.97 | 1.12 | 1.04 |
| 0.50 | 0.87 | 1.01 | 0.94 |
| 0.55 | 0.79 | 0.89 | 0.83 |

- **bG 0.45 vs verdict 1:** p50 data/sim up **+5.2..+10.3 %** (predicted +7..+12 %; falsifier < 4 % or
  > 15 % **not triggered**; the 0° X/Y and X 20° points are 1–2 % under the band). The X saturated
  fraction falls 0.025–0.075 (predicted 0.03–0.08). Rise and FWHM are within 1 ns of verdict 1
  (predicted ±4). The primaries feed the amplitude linearly.
- **The amplitude falls faster with bG than predicted:** −8..−11 % at 0.50 and −16..−22 % at 0.55,
  against −6..−9 % and −13..−18 %. So the joint optimum is **bG ≈ 0.47**, below the pre-registered
  0.49–0.53 (prediction missed). X alone wants ~0.44, Y alone ~0.50. That 15 % split is the §20.2 Y/X
  offset, still present with the corrected primaries.
- With the corrected primaries, verdict 1's bG barely moves (0.45 → 0.47) and the photon-feedback
  picture survives. Refit bG properly once the Y/X offset has a cause, since a joint fit now averages over it.

### 20.4 What changes

1. Ion: the K0 2.05 reconstruction is a drop-in for s = 1.75. Default switch pending (20.1).
2. ρ_s: the Y-wing preference moves from 2.5 to ~4–5 MΩ/sq. That is above the T2b band (1.42–2.56),
   so check the T2b measurement's geometry against what the wing sees. Y/X is **not** ρ_s: next is
   the strip-C / protection charge division (§19), which needs the MX17 X/Y strip capacitances.
3. Primaries: use the e17 clusters from now on, with bG ≈ 0.47 as the interim value.
4. Separately (2026-10-16, night): every PyROOT Garfield run on lxplus/condor used the LCG_109 view's
   Garfield headers against our 927e5c21 library. Class layouts differ (ComponentGrid 688 vs 720 bytes,
   AvalancheMicroscopic 616 vs 680, MediumMagboltz, …), which overruns the heap. All 56 nc1 slices
   died or hung in the guard. Fixed in nTof_x17 `setup_garfield.sh` (faefdac, ROOT_INCLUDE_PATH order).
   The avalanche script is not reproducible per seed (two same-seed runs differ event by event), so the
   audit is statistical: condor 4474927, the August script (baf2873) and tarball with the fixed headers,
   8 × 140 ev at 490 V uniform, against `results_v2` (42.7k, 1120 ev). **Result:** gain 45.0k (θ 1.97,
   σ₀ 33.6 µm) vs 42.8k (θ 1.57, σ₀ 33.7 µm), +5.3 ± 2.6 %. August's own same-setup repro gave 46.4k
   (280 ev), so against the pooled August value (43.5k) it is +3.5 % (~1.4σ). No significant bias: the
   mismatch did not visibly distort the August uniform-field calibration (< ~5 % on gain, σ₀ unchanged).
   Field-map runs were not audited separately.

### 20.5 Null-collision avalanche campaign (condor 4474926): ε × gain moves only +0.7 %

56 slices, `--null-collision`, patched Garfield with the header fix, one 490 V map, pooled as the nc0 file
(`scratchpad/2026-10-16/pool_nc1.py` → EOS `response_sim/avalanche/aval_calib_meshfield_nc1_pooled.json`).

| | nc0 (production) | nc1 | nc1/nc0 | pre-registered |
|---|---|---|---|---|
| survival (= ε_mesh) | 0.9559 ± 0.0026 | **0.9928 ± 0.0011** | ×1.039 | 0.990–0.997 ✓ |
| gain \| survived | 24 094 | **23 368** | ×0.970 (±1.2 % stat) | within ±3 % (edge) |
| θ (Polya) | 1.14 | 1.15 | | |
| σ₀ / t_arrival | 94.6 µm / 6.66 ns | 94.3 µm / 6.59 ns | | |
| **ε × gain (charge per drifted e⁻)** | | | **×1.007** | |

- The electrons that null-collision steps rescue from the wires mostly come with a smaller avalanche,
  so the conditional gain drops ~3 % (2.5σ). The 4-slice desktop test (0.998 ± 0.040) could not see that.
- The chain thins by ε and then draws from the Polya at gain|surv, so its amplitude moves by the product:
  +0.7 %. **No bG refit is needed.** The "mesh transparency ~4 % low" is real, but nearly all of it cancels.
- The chain default is still the nc0 file. Switching to nc1 is a +0.7 % change, worth doing for honesty
  at the next full chain pass. The Mesh card on the site says this (`make_chain_animation.py`, 2026-10-17).

## 21. Y/X is the strip spread, not the electronics; the drift column; neighbours (2026-10-17, desktop)

Yoga was offline, so the evaluation was moved to the desktop. `scratchpad/2026-10-17/m3_slopes_longrun_desktop.py`
rebuilds `hv_slope/peaks_longrun_m3.parquet` from the desktop's det3 decoded files (hv_slope `extract.py`), and
`ledger_m3.py` and `shape_all_desktop.py` then reproduce yoga's verdict-1 ledger and shapes **row for row**.

### 21.1 Pedestal noise: Y is not quieter than X

`noise_det3_satscan_feu0{7,8}_cmq.json` (same run, CM removed): per-channel noise spectrum Y/X = **1.017**
(rms over the 17 bins). The DREAM's ENC is ~810 e⁻ + ~11 e⁻/pF (Irakli, §19), and 7.4 ADC ≈ 2300 e⁻, so most of the
noise depends on the input capacitance. A protection-network charge division large enough for the 15 % Y/X (C_Y
75–100 pF below C_X, §19) would make Y **~25 % quieter**. A 15 % higher FEU 8 gain would make it 15 % noisier.
Neither is seen. This is indirect: the noise model is a bench one.

### 21.2 The X/Y split belongs to the chamber, not the FEUs

`xy_balance_survey.py` → `xy_balance_survey.log`: f = q_X/(q_X+q_Y), the sums of combined-hits amplitudes ≥ 100 ADC,
every June run on the desktop.

| chamber | FEUs | runs | f (median) |
|---|---|---|---|
| **det3** | **7/8** | 6-27 long run, 6-27 overnight | **0.484, 0.484** (107k ev) |
| **det3** | **3/4** | 6-22 long / longer / short | **0.481, 0.491, 0.481–0.494** (130k ev) |
| det2 | 6/8 | 6-22 (four subruns) | 0.502–0.532 |
| det1 | 3/4 | 6-17/6-18 | 0.453–0.480 |
| det4 / det6 / det7 | 6/8, 3/4, 6/8 | 6-24, 6-26 | 0.462 / 0.462–0.495 / 0.428 |

det3 gives the same f on two different FEU pairs, so **FEU gain is excluded** as the Y/X cause. Chambers differ
by ±0.04 in f (±17 % in q_Y/q_X); that is a per-assembly property (PLAN_38 saw det2 vs det3 too).

### 21.3 Sim vs data with one observable: the sums agree, the peak strip does not

`xy_balance_m3.py` → `xy_balance_m3.log`. Both legs go through FeuReader, with strips ≥ 5σ and ≥ 100 ADC. q is the sum of
strip peaks and pk the peak strip. Data are binned by M3 angle (±3° in both projections); "unsat" means no strip
≥ 3500 ADC.

| unsat | f data | f verdict 1 | f ρ_s 5M | pkX/pkY data | v1 | ρ5M | pkY/qY data | v1 | ρ5M | pkX/qX data | v1 | ρ5M |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0° | 0.476 | 0.478 | 0.487 | 1.073 | 1.213 | 1.104 | 0.382 | 0.320 | 0.367 | 0.442 | 0.422 | 0.430 |
| 10° x | 0.573 | 0.599 | 0.609 | 0.915 | 1.034 | 0.990 | 0.365 | 0.318 | 0.367 | 0.244 | 0.219 | 0.231 |
| 10° y | 0.390 | 0.377 | 0.381 | 1.142 | 1.316 | 1.273 | 0.238 | 0.198 | 0.211 | 0.439 | 0.420 | 0.422 |

- **The sim splits the total charge between X and Y correctly** (Δf ≤ 0.03 at every angle). Charge division, strip
  capacitance and FEU gain all scale the *sum*, so none of them is the Y/X offset.
- The offset is in **how each view's charge spreads over its strips**. At 0° ρ_s 5M brings Y's peak fraction from
  0.320 to 0.367 (data 0.382), and pkX/pkY from 13 % off to 3 % off. That agrees with the ledger: Y/X at 0° goes
  1.12 → 1.045 (§20.2).
- What remains is on **inclined** tracks: the inclined view's peak fraction is 5–11 % low in the sim.

### 21.4 Inclined tracks spread over fewer strips in the data: the drift column

`xy_extent_m3.py` → `xy_extent_m3.log`, unsaturated events, charge-weighted RMS across the strips [mm]:

| | X rms data / v1 | Y rms data / v1 / ρ5M | sim/data, inclined view |
|---|---|---|---|
| 0° | 0.66 / 0.73 | 0.84 / 1.11 / 0.87 | – |
| 10° x | 1.46 / 1.66 | | 1.14 |
| 20° x (n = 20) | 2.75 / 3.29 | | 1.20 |
| 10° y | | 1.47 / 1.74 / 1.68 | 1.14–1.18 |
| 20° y (n = 8) | | 2.71 / 3.33 | 1.23 |

At 20° the track term dominates (L·tanθ/√12 = 3.15 mm for 30 mm), and the ratio gives a charge-visible column of
**~25 mm** (~26 mm from the 10° points after removing the vertical width in quadrature). nTof_x17
`mx_june_wft/GAP_STUDY_2026-07-30.md` measured det3's column at **27.9 ± 0.1 ± 1.0 mm** (a dished cathode,
25.7–29.2 mm across the chamber; det2 reads the full 30.5). The June weekend analysis read 23.4 mm. **The chain
uses the CAD 30 mm.**

Consequences, if the column is ~26–28 mm:
- 7–13 % fewer primaries per vertical track, so bG would come down;
- the inclined-view peak fractions and the inclined Y/X rows close;
- the vertical X width (0.73 vs 0.66 mm) moves a little the right way (shorter diffusion path).

The new Stage B option `--drift-gap-mm` drops clusters beyond the cathode and sets E = 1000 V / gap. The default is
30, so earlier runs are unchanged. Smoke test: 27.9 mm keeps 0.9275 of the primaries, and v_drift goes 39.1 → 38.5 µm/ns.
Test points are condor 4477508 (§21.6, pre-registered in `stageBC_w2_gap.sub`).

### 21.5 Neighbour delay on verdict 1: the central strip is fixed, the ±1 is not

`dump_pm3_m3.py` and `2026-10-13/xtalk_fit.py`. Data are vertical by M3 (455 events), sim is th00. Logs are
`xtalk_fit_{v1,r5}_m3.log`.

| | data | verdict 1 | ρ_s 5M | §10.4 (Aug chain) |
|---|---|---|---|---|
| central peak t, X / Y [ns] | 195 / 180 | 210 / 180 | 210 / 195 | 270 / 240 |
| ±1 delay after central, X / Y [ns] | **+45 / +60** | +150 / +120 | +165 / +150 | +135 / +120 |
| ±1/central peak, X / Y | 0.368 / 0.408 | 0.396 / 0.465 | **0.369 / 0.437** | |
| fit: in-phase k, X / Y | | 0.168 / 0.059 | 0.175 / 0.135 | 0.165 / 0.167 |

- The ion speed and the 166 ns shaper fixed the **central** pulse timing (270 → 210 ns against 195).
- The ±1 neighbours have the right **amplitude** at ρ_s 5M but peak **~100 ns late**. Higher ρ_s makes this worse.
  The fit needs a prompt in-phase copy of the central strip (k ≈ 0.13–0.17) plus a derivative term (k′ ≈ 0.15/sample).
  That is the shape of capacitive inter-strip coupling (or of extra prompt sharing). The chain has neither.
  It is the clearest remaining shape row. Candidates: inter-strip/inter-pad capacitance in the readout board
  (the buried X/Y tracks under the pad layer), or a protected-FEU input network that differentiates.
- The data's delay is shorter here (+45/+60) than in §10.4 (+90/+75), because the selection changed from the T14
  reco window to true M3 angles.

### 21.6 Verdict 2 and the drift-gap test (condor 4477507 / 4477508, evaluated on yoga)

All 40 jobs ended with "=== done". `scratchpad/2026-10-17/eval_v2_yoga.sh` is a copy of `eval_v2.sh` without the
condor wait loop, because the desktop was in Windows. Logs: `ledger_v2.log`, `shape_v2.log`, `xy_balance_v2.log`,
`xy_extent_v2.log`, `xtalk_fit_v2r5*_m3.log`. All points use the e17 clusters, bG 0.47, the nc1 calib, tp 166 and `_cmq`.
Ledger values are p50 data/sim, averaged over the 5 angles. Rise/FWHM are sim − data.

**Verdict 2: passes, and the effects add.**

| point | X | Y | joint | pre-registered X / Y | inclined rise X10/X20/Y10/Y20 [ns] |
|---|---|---|---|---|---|
| v2r2s (ρ 2M, s1.75) | 0.922 | 1.064 | 0.991 | 0.91–0.95 / 1.05–1.10 | +10 / +8 / +8 / +16 |
| v2r2k (ρ 2M, K0 2.05) | 0.936 | 1.080 | 1.005 | +1..+3 % | +15 / +13 / +13 / +22 |
| v2r5s (ρ 5M, s1.75) | 0.865 | 0.955 | 0.909 | 0.84–0.90 / 0.92–0.99 | +8 / +7 / +7 / +17 |
| v2r5k (ρ 5M, K0 2.05) | 0.878 | 0.970 | 0.923 | as v2r5s + ion | +13 / +12 / +12 / +22 |

- Every point is inside its band. v2r2s rise/FWHM are within 2 ns of verdict 1 (band ±3).
- The ion shift is +1.5 % X and +1.5–1.6 % Y, with +5..+6 ns rise, at both ρ_s values (and at both gaps).
  So additivity holds to < 0.2 % and 1 ns. **Falsifier not triggered:** the 10-16 one-variable results combine linearly.
- ρ 5M raises the sim Y/X by 4.5 % (band +1..+8 %). Data/sim Y/X goes 1.154 → 1.104, so a **~10 % Y/X offset remains**.
- Y wings: not re-measured (`profile_tails.py` needs T14 outputs, which this chain pass did not make).
  Proxy: the 0° unsaturated pkY/qY is 0.374 vs data 0.382 (r5verdict1: 0.367).

**Drift gap: the column length closes the inclined extent; the falsifier is not triggered; the amplitude prediction missed.**

| inclined-view RMS extent [mm] | data | 30 mm | 27.9 mm | 26.0 mm |
|---|---|---|---|---|
| 10° x | 1.46 | 1.65 | 1.53 (−7 %) | **1.43** (−13 %) |
| 20° x (n = 20) | 2.75 | 3.31 | 3.04 (−8 %) | 2.85 (−14 %) |
| 10° y | 1.47 | 1.66 | 1.55 (−7 %) | **1.46** (−12 %) |
| 20° y (n = 8) | 2.71 | 3.29 | 3.04 (−8 %) | 2.85 (−13 %) |
| vertical X width, 0° | 0.66 | 0.72 | 0.70 (−3 %) | 0.68 (−6 %) |

- The extent moves 7–8 % / 12–14 % (predicted 5–8 / 10–14). The 26.0 point is within 2 % of the data at 10°
  (predicted within 5 %). Interpolating the 10° rows gives a column of **~26.5 mm**. The low-statistics 20° rows
  extrapolate to ~25 mm. GAP_STUDY's 27.9 ± 1 mm is 1–2σ above. Vertical X width −3 / −6 % (predicted 2–5 %).
- Rise moves ≤ 1 ns and FWHM ≤ 1 ns (predicted ≤ 6 ns; the falsifier bound is 10 ns). The ion arm k − s matches 30 mm.
- **Amplitude missed:** p50 data/sim rises +2.2 / +2.9 % (X / Y) at 27.9 mm and +5.7 / +5.4 % at 26.0 mm.
  The prediction was +5..+8 / +10..+15 % (at th00: +3.5–4.7 / +6.8–7.2 %). That is about half of what the primaries
  lose (×0.928 / ×0.867). The removed clusters are the deepest. They are the most diffused, so they put less than
  their share on the peak strip; the pre-registration assumed every primary counted the same.
- Peak fractions (unsat, pk/q of the inclined view, data / 30 / 27.9 / 26.0): 10° x 0.244 / 0.237 / 0.251 / 0.265;
  10° y 0.238 / 0.216 / 0.228 / **0.240**; 20° y 0.193 / 0.150 / 0.163 / 0.170. Each rises by the extent fraction, as
  predicted. Inclined Y closes at 26 mm, inclined X now overshoots by 9 %. The vertical-view fractions rise by 2–5 %.
- **The inclined Y/X does not close.** pkX/pkY is unchanged by the gap (10° y: data 1.142, sim 1.261 / 1.266 / 1.263;
  10° x: 0.915 vs 0.99), because both views' peaks rise together. The Y/X offset is not the column either.
- bG: on the ρ 5M arm, the joint data/sim at bG 0.47 is 0.909 (30 mm), 0.933 (27.9 mm) and 0.959 / 0.974 (26.0 mm,
  s / K0). With §20.3's ~2 %/0.01 slope, the joint optimum is **bG ≈ 0.425 at 30 mm, ≈ 0.45 (s) / 0.46 (K0) at 26 mm**.
  A shorter column needs a *higher* bG than the 30 mm column at the same ρ_s, so §21.4's "bG would come down" had
  the sign wrong.

**Neighbour delay** (`xtalk_fit_v2r5*_m3.log`, th00): the X central peak is now 195 ns, as in the data (210 in verdict 1).
The ±1 delay stays +150 / +120 ns (X / Y) at 30 mm and +135 / +120–135 ns at 26 mm, against the data's +45 / +60.
The ±1 amplitude is right (X 0.363–0.366 vs 0.368, Y 0.423–0.428 vs 0.408). The shorter column trims the chain's late
neighbour tail, so the fitted in-phase k falls (X 0.155 → 0.062). The derivative term stays at k′ ≈ 0.16–0.21/sample.
Inter-strip coupling is still the missing piece.

**What this settles, and what is Dylan's call:**
1. The verdict-2 combination is sound. The chain can move to e17 + nc1 + ρ_s 5M with bG refit.
2. Drift column: the data want **~26–26.5 mm** for det3 (the 10° extent and the inclined Y peak fraction), within
   1–2σ of GAP_STUDY. Proposed: `--drift-gap-mm 26.5` for det3, with 27.9 as the systematic, then refit bG (≈ 0.45–0.46).
3. Ion: unchanged by everything here. K0 2.05 costs +5 ns of inclined rise (which is already +7..+17 ns over the data)
   for +1.5 % amplitude.
4. Still open: the ~10 % Y/X peak offset (neither ρ_s nor the column fixes it; the X/Y sums agree, §21.3), the ±1
   delay (inter-strip C), and the inclined-X peak fraction now overshooting at 26 mm.

### 21.7 bG refit at the det3 column; the v3 baseline (2026-10-18)

Dylan's calls from §21.6: ρ_s 5 MΩ/sq, a 26.5 mm drift column for det3 (27.9 as the systematic), and the
reconstructed ion with K0 2.05 (1.80 / 2.25 as the systematic). Condor 4480068 (`scratchpad/2026-10-18/stageBC_w3_bgrefit.sub`,
pre-registered) scans bG at that column. Logs are in `scratchpad/2026-10-18/*_v3.log`.

| bG | X | Y | joint | Y/X | pre-registered joint |
|---|---|---|---|---|---|
| 0.43 | 0.997 | 1.104 | 1.049 | 1.106 | 1.02–1.07 |
| **0.45** | 0.953 | 1.054 | **1.002** | 1.107 | 0.99–1.02 |
| 0.47 | 0.919 | 1.020 | 0.968 | 1.110 | 0.955–0.975 |

- Every band is hit. **The joint optimum is bG = 0.450** (predicted 0.45–0.46). X alone wants ~0.43 and Y ~0.48: the
  ~10 % Y/X offset again, unchanged by bG (1.106–1.110, band 1.08–1.12).
- bG is a pure amplitude knob here. The 10° extent is 1.46 / 1.48 mm at every bG (data 1.46 / 1.47), and the inclined rise
  and FWHM move ≤ 1 ns. The rise is still +11..+13 ns over the data (Y 20°: +22). **Falsifier not triggered.**
- Inclined Y is the worst view: Y data/sim is 1.13 at 10° y and 20° y, and 0.97–1.03 elsewhere.

**The v3 baseline** is `scripts/condor/recipes/det3_v3.args` (README next to it): ρ5M, e17, nc1, tp 166, 26.5 mm, K0 2.05, bG 0.45.
`run.py`'s defaults stay bare so that old runs reproduce. The s1.75 ion would give bG ≈ 0.443.

### 21.8 Inter-strip coupling: it fixes the neighbour's leading edge, and about half of the delay (2026-10-18)

New Stage B option `--strip-coupling K,TAU_NS` (`response/dream/coupling.py`, test `python -m response.dream.test_coupling`).
It models adjacent channels coupled through a capacitance C_c on input nodes of finite series impedance R. Each neighbour gets
k·i_c + τ·di_c/dt, with τ = R·C_c, and the central strip loses the same (Kirchhoff), so charge is conserved.
It acts on the induced currents before the shaper. The default is off. An offline emulation on v2r5kg260's 60 ns samples
(`scratchpad/2026-10-18/emul_coupling.py`, `.log`) picked the scan points. They ran on condor 4480073 on top of the v3 baseline
(`stageBC_w3_coupling.sub`, pre-registered). `nb_compare.py` → `nb_compare_v3c.log` compares the shapes unscaled;
th00 is vertical by M3.

| (k, τ) | ±1 delay X / Y [ns] | nb(t = 0) X / Y | nb rms X / Y | joint p50 | rise X10/X20/Y10/Y20 | FWHM X10 / Y10 |
|---|---|---|---|---|---|---|
| data | +45 / +60 | 0.219 / 0.225 | – | – | – | – |
| off (v3) | +135 / +135 | 0.166 / 0.181 | 0.029 / 0.037 | 1.002 | +12/+12/+12/+22 | −1 / +14 |
| 0.03, 10 | +120 / +120 | 0.203 / 0.215 | 0.020 / 0.032 | 1.024 | +15/+11/+16/+22 | −3 / +12 |
| **0.05, 15** | **+105 / +105** | **0.226 / 0.236** | **0.018 / 0.032** | 1.037 | +13/+8/+16/+20 | −5 / +9 |
| 0.06, 20 | +90 / +120 | 0.244 / 0.251 | 0.021 / 0.032 | 1.034 | +9/+3/+13/+18 | −9 / +5 |

- **Falsifier not triggered:** the X delay goes below +115 ns, and the neighbour rms improves 39 % in X.
  The delay (X +90..+110, Y +90..+120), leading edge, amplitude (+2..+4 %, band +1..+8), extent (+0..+2 %) and ±1/central (±0.01)
  predictions hold. **Missed:** the Y rms improves only 13 % (band 20–35 %).
- (0.05, 15 ns) is the best point. The leading edge closes in both views. The delay closes by ~30 of the ~90 ns.
  The ±1/central peak ratio stays where it was. A residual in-phase k ≈ 0.05 remains in the X fit: the data neighbour carries
  a prompt copy of the central that a charge-conserving coupling (which also takes it out of the central) does not give.
- Side effects go the right way: the X 20° rise drops 4 ns (12 → 8) and the Y FWHM excess drops 14 → 9 ns.
  The X FWHM goes to −5. At (0.06, 20) X overshoots (−9 FWHM).
- Physically: τ = 15 ns is R·C_c, e.g. 300 Ω × 50 pF. k = 0.05 is C_c over an effective CSA input of ~1 nF. Both are plausible
  for 50 cm strips and a protected input, but neither is measured. **Not in the v3 baseline.** It is a candidate, with its bG refit
  in condor 4480074 (bG 0.46).
- What is left of the delay (~+45..+60 ns over data) is a chain neighbour that stays high too long after its peak. That is the
  resistive-layer spread, so ρ_s and coupling trade against each other here. A joint fit of ρ_s, k and τ on the Y wings and the
  neighbour shape is the next step.

### 21.9 The Y view saturates too rarely on inclined tracks in the sim (2026-10-18)

The saturated fraction per view (ledger `sat d / s`, v3 baseline), at the angle in that view, has no selection behind it:

| | 0° | 10° | 20° |
|---|---|---|---|
| X data / sim | 0.379 / 0.371 | 0.303 / 0.268 | 0.229 / 0.224 |
| Y data / sim | 0.334 / 0.284 | 0.331 / 0.218 | 0.273 / 0.164 |

X matches at every angle. In the data, inclined Y tracks rail as often as vertical ones (0.33), while the sim's fraction falls by
a factor 1.7. So **the sim's Y peak strip is too dim, mostly on inclined Y tracks**. That is the ledger's Y/X (1.17–1.20 at 10°/20° in y).
The unsaturated-only X/Y balance (§21.3) does not see it, because it drops exactly the bright Y events: its Y peak fractions agree
with the data at 26.5 mm. Ledger Y/X at p25, far from the rail, agrees with p50 to within ±0.1 and keeps the inclined-Y excess (1.10–1.25), so this is a real dimness
and not a saturation artefact.

Leads, none tested yet:
- **Resistive strips.** The kernels (`response/solver/kernels.py`) model the resistive layer as strips along y (550/250 µm
  line/gap, insulating gaps) over a checkerboard of pads read out as combs: a Y channel is a pad row, an X channel a pad column,
  both on a 1.56 mm pitch (from the gerbers). All lateral spread on the strips therefore crosses Y rows and stays inside an X
  column. So the strip geometry and ρ_s act on Y far more than on X. ρ_s itself is pinned at 4–5.6 MΩ/sq by the Y wings; the
  printed line/gap has not been measured on a board.
- **Shape.** The data Y central falls faster after its peak than the sim's (−0.08..−0.10 of the peak at +360..+540 ns) and undershoots more
  (−0.125 vs X −0.02). The X tail is slower than the sim's. A per-view shaping difference in the data that the chain lacks
  would change the peak amplitudes per view. `scratchpad/2026-10-18/central_shape_feu_swap.py` was meant to test FEU 7/8 against
  3/4. Its first version is **not usable** (proxy selection, pedestal handling). It needs the 6-22 run's pedestals and an M3 selection.

### 21.10 Systematics of the v3 baseline, and bG with coupling on (2026-10-18)

Condor 4480074 (`scratchpad/2026-10-18/stageBC_w3_sys.sub`, pre-registered) changes one input at a time at bG 0.45.
Logs: `*_v3s.log`. Differences are vs v3k265b45 (joint 1.002; rise X10/X20/Y10/Y20 +12/+12/+12/+22 ns).

| variant | joint data/sim | Δ joint | Δ rise [ns] | Δ FWHM [ns] | other | pre-registered |
|---|---|---|---|---|---|---|
| K0 1.80 | 1.034 | +3.2 % | +10..+11 | +12..+14 | | +2..+5 %, +8..+12, +10..+16 |
| K0 2.25 | 0.985 | −1.7 % | −6..−7 | −7..−8 | | −1..−3 %, −4..−8 |
| gap 27.9 mm | 0.987 | −1.5 % | 0 | ≤ 2 | inclined extent +4..+5 % | −1..−4 %, +5..+8 % |
| coupling (0.05, 15), bG 0.46 | 1.020 | – | as bG 0.45 ± 1 | ± 1 | nb rms X 0.019, Y 0.034 | 1.01–1.03 |

- Every band is hit except the 27.9 mm extent (+4..+5 %, band +5..+8 %), which misses by < 1 %. **Falsifier not triggered:**
  the one-variable differences of §20–21 carry over to the 26.5 mm baseline.
- **v3 systematics as quoted:** ion K0 1.80 / 2.25 → amplitude +3.2 / −1.7 %, inclined rise +11 / −6 ns. Column 27.9 mm →
  amplitude −1.5 %, inclined extent +4..+5 % (the 27.9 point sits 4–5 % above the data's extent; 26.5 sits within ±1 % at 10°).
  In bG terms: K0 1.80 → bG +0.016, K0 2.25 → −0.008, 27.9 mm → −0.007 (slope ~2 %/0.01).
- **With coupling on, the joint optimum is bG ≈ 0.47** (1.037 at 0.45 and 1.020 at 0.46). The coupling takes ~3.5 % of the peak
  amplitude into the neighbours and their bipolar term, and bG puts it back. The shapes are unchanged by the refit.

### 21.11 ρ_s × coupling: the delay and the Y spread pull ρ_s in opposite directions (2026-10-18)

Condor 4480075 (`stageBC_w3_rho2c.sub`, pre-registered): the v3 settings at ρ_s 2 MΩ/sq, without and with coupling (0.05, 15 ns).
Logs: `*_v3r2.log`, `nb_compare_v3r2.log`.

| point | ±1 delay X / Y | ±1/central X / Y | nb rms X / Y | Y rms, vertical Y (data 0.85 mm) | joint | Y/X |
|---|---|---|---|---|---|---|
| data | +45 / +60 | 0.368 / 0.408 | – | 0.85 | – | – |
| ρ5M (v3) | +135 / +135 | 0.361 / 0.423 | 0.029 / 0.037 | 0.82 | 1.002 | 1.107 |
| ρ5M + coupling | +105 / +105 | 0.361 / 0.422 | 0.018 / 0.032 | 0.81 | 1.037 | 1.095 |
| ρ2M | +105 / +105 | 0.385 / 0.458 | 0.034 / 0.051 | 1.05 | 1.091 | 1.161 |
| ρ2M + coupling | +90 / +75 | 0.391 / 0.456 | 0.029 / 0.055 | 1.03 | 1.124 | 1.146 |

- **Falsifier not triggered:** ρ_s does set the neighbour tail (2M is 30 ns shorter than 5M). The delay, ±1/central, leading
  edge and amplitude (+8.9 %, band +4..+9) predictions hold. The inclined extent moves +5 % in Y at 10° y (band ±3 %): a small miss.
- But ρ2M wrecks the Y lateral profile (vertical Y rms 1.05 vs 0.85 mm, ±1/central Y 0.46 vs 0.41, nb rms Y 0.051–0.055)
  and worsens Y/X (1.16). **ρ5M + coupling stays the best compromise.** The delay cannot be traded against ρ_s without breaking the Y wings.
- What is left is a prompt, in-phase neighbour share (the X fit's residual k ≈ 0.05, §21.8) with no loss from the central.
  **Not the insulator:** the August ρ_s × d_k scan (`DESKTOP_RUNS_2026-08-07.md`) has the prompt X ±1 share moving only
  0.205 → 0.212 over d_k 50 → 125 µm. The glue bracket (series 57–82 µm) would move it < 1 %, so a thicker-glue kernel is not worth
  building. **Not diffusion:** the sim's vertical X width is already slightly above the data's (0.69 vs 0.66 mm), and more diffusion
  would widen it. A prompt share that does not widen the cluster and does not drain the central points back at the electronics:
  an input network that couples neighbours with a non-conserving, in-phase term (e.g. a shared reference or common return through the
  protection network). That needs Irakli's schematic. The cheap next test is a non-conserving coupling option (neighbours get k·i_c,
  the central keeps i_c) to see whether k ≈ 0.05 closes the remaining delay.
- **Non-conserving prompt share, tested offline: does not work.** `--strip-coupling K,TAU,KFREE` (a diagnostic third value: neighbours
  gain kf·i_c, the central keeps everything) was emulated on v3c5t15's samples (`scratchpad/2026-10-18/emul_kfree.py`, `.log`).
  kf 0.02–0.08 moves the delay only +105 → +90/+75 ns. It raises ±1/central past the data (X 0.37 → 0.41) and *worsens* the shape
  rms (X 0.018 → 0.035). Not run on condor. **Reading:** the data's neighbour is not the chain's neighbour plus a prompt extra. It carries
  the same share *earlier instead of later*. That wants the slow resistive delivery reduced (higher ρ_s) and a prompt lateral
  delivery added, with the total share held. The Y wings already sit near the data at 5 MΩ/sq, so that is a joint re-fit of the
  lateral model (ρ_s with a prompt spread term), not one more knob on top.
