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
sheet's time response factorises out); ion mobility within ±30 %
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
