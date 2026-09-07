# ⚠️ CORRECTION 2026-08-11 — the rise estimator read railed pulses as fast

**Every rise-time number in `OVERNIGHT_2026-08-10.md` §2 and §4 was measured
with a saturation-biased estimator that speeds up the DATA leg six times harder
than the sim leg. The §4 offsets are re-derived below. The f_eff contradiction
survives — it moves from ≈0.12–0.22 to ≈0.17–0.31 against a defended 0.9056.
This is a numbers correction, not a resolution.**

---

## The defect

`mx17_sim_wft/t14_compare.py::extract_view` measured the 10–90 % rise against
the waveform's own maximum, with no saturation handling. On a railed waveform
that maximum *is* the rail, so the 90 % level sits at 0.9 × rail — a point the
still-rising edge reaches early — and the pulse reads **fast**. The FWHM reads
**wide** for the mirror reason.

That would be harmless if it hit both legs alike. It does not:

| view | data rails ≥ 3500 ADC | sim rails |
|---|---|---|
| X, 10° | 24.1 % | 4.1 % |
| X, 20° | 24.1 % | 3.4 % |
| X, vertical | 26.0 % | 5.9 % |

Six times more often on the data leg — so it is very nearly a pure data-leg
speed-up, **in the direction that makes the sim look slow**, which is the
direction of the discrepancy under investigation.

`wft/model.py` censors at `SAT = 3550` inside the forward fit, so the fitted
reconstruction was never exposed to this; only the waveform-level comparison
was. And this comparison has been caught by the same railed-denominator failure
once already, on the undershoot "cap" that turned out to be a railed-bin
artefact (`mx17_sim_wft/hv_slope/HV_SLOPE_2026-08-09.md:26-34`). Second
instance, same shape.

## The fix

`t14_compare.py` gains `SAT_ADC = 3500` and censors the **shape** observables
(`rise_ns`, `fwhm_ns`) above it, identically on both legs, recording a
`saturated` flag and the uncensored `rise_ns_raw`/`fwhm_ns_raw` alongside.
Amplitude observables (`peak_amp`, `q_event`) are deliberately **not** censored:
the §7 amplitude ledger keeps the railed events and corrects for them
explicitly. The summary now reports `rise_censored_frac` and
`rise_undefined_frac` separately per leg — a deliberate cut and a real
shape statement are different things.

3500 rather than wft's 3550 because the DREAM response bends before the hard
rail; the statistic wants to exclude compression, not only clipping.

The tables below are read-only re-derivations from the **frozen** parquets
(drop `peak_amp ≥ 3500` on both legs, recompute), so they need no re-run:
`scratchpad/2026-08-11/rise_saturation_recensor.py`, machine-readable output
`design/report/rise_saturation_recensor_2026-08-11.json`. The uncensored column
reproduces the committed §4 tables exactly, which is what makes this a
re-derivation rather than a different analysis.

## §4 re-derived — X view, angled closure

Offsets are sim − data at matched quantiles, each sim leg against the data leg
extracted alongside it (the original convention).

### 10°

| | p5 | p10 | p25 | p50 | p75 | p90 | mean | span |
|---|---|---|---|---|---|---|---|---|
| offset with ions — **committed** | 76.8 | 82.4 | 78.3 | 75.7 | 71.1 | 66.6 | **+75.2** | 15.7 |
| offset with ions — **censored** | 63.6 | 70.0 | 71.9 | 71.9 | 67.7 | 62.2 | **+67.9** | 9.7 |
| ions removed — **committed** | −10.1 | −13.0 | −18.1 | −15.2 | −6.2 | 5.4 | **−9.5** | 23.5 |
| ions removed — **censored** | −26.7 | −23.9 | −23.4 | −18.0 | −7.7 | 3.2 | **−16.1** | 29.8 |

### 20°

| | p5 | p10 | p25 | p50 | p75 | p90 | mean | span |
|---|---|---|---|---|---|---|---|---|
| offset with ions — **committed** | 79.2 | 79.9 | 72.6 | 69.3 | 69.0 | 72.0 | **+73.7** | 10.9 |
| offset with ions — **censored** | 46.7 | 62.4 | 67.1 | 66.0 | 67.0 | 70.5 | **+63.3** | 23.9 |
| ions removed — **committed** | −5.2 | −21.5 | −28.4 | −29.6 | −30.0 | −27.1 | **−23.6** | 24.9 |
| ions removed — **censored** | −41.0 | −38.2 | −31.4 | −32.1 | −31.4 | −27.2 | **−33.5** | 13.9 |

### Implied f_eff

Linear interpolation on the mean offset between "ions removed" (f = 0) and
"with ions" (f = 0.9006), i.e. where the offset would vanish:

| | 10° | 20° |
|---|---|---|
| committed | 0.101 | 0.219 |
| **censored** | **0.173** | **0.312** |

## §4 re-derived — vertical f_ion rigidity

Rise shift relative to `noions`, X view. The censor barely moves this table —
which is the expected result, because it is a sim-vs-sim comparison and the sim
legs rail at a common few percent:

| leg | f_ion | p5 | p10 | p25 | p50 | p75 | p90 | span |
|---|---|---|---|---|---|---|---|---|
| fion030 | 0.30 | 18.1 | 21.3 | 27.2 | 28.4 | 10.5 | −6.8 | 35.2 |
| fion050 | 0.50 | 48.5 | 53.5 | 56.4 | 43.1 | 5.9 | −14.5 | 70.9 |
| fion070 | 0.70 | 86.2 | 85.2 | 78.8 | 51.4 | −2.4 | −23.2 | 109.4 |
| **default** | **0.9006** | **119.5** | **112.4** | **97.5** | **63.1** | **−2.0** | **−31.8** | **151.3** |

The §4 Step-1 conclusion is unchanged: the ion term imposes a floor, it does not
delay.

## What changes, and what does not

**Changes.**

- The sim's rise excess at inclined incidence is **~10 ns smaller than
  committed** at both angles (75.2 → 67.9 ns at 10°, 73.7 → 63.3 ns at 20°).
- **Removing the ion term overshoots harder than committed**: the sim becomes
  16 ns (10°) to 34 ns (20°) too fast at the mean, against the committed 10 and
  24 ns. The §4 sentence "deleting the ion term entirely slightly overshoots"
  understates it; it is no longer *slight* at 20°.
- The implied f_eff window moves **up** to ≈0.17–0.31 from ≈0.10–0.22.
- The 20° "with ions" offset span widens from 10.9 to 23.9 ns, so the claim
  that the ion term is a **near-rigid** delay at 20° is weaker than committed —
  the near-rigidity at 10° (span 9.7 ns) is if anything cleaner.

**Does not change.** The contradiction. The data's rise is reproduced at
f_eff ≈ 0.2–0.3 against a defended, measured f_ion = 0.9056 — a factor 3–5,
where it was a factor 4–9. Two independent inclinations, six quantiles each,
nothing fitted. Correcting the estimator moved the number in the direction that
*reduces* the contradiction and did not come close to closing it, which is the
more informative outcome than if it had moved the other way.

## End-to-end on the frozen vertical point, and the check that the fix is inert
## where it should be

The fixed `t14_compare.py` was re-run on the **frozen** default leg
(`w2_rho2M/default` against the sat-scan data leg), output at
`~/x17/response_sim/stageB_w2/t14_compare_censored/`. The frozen
`t14_compare/` is untouched.

| | frozen | censored |
|---|---|---|
| X peak ratio | 0.5646 | **0.5646** |
| X q_event ratio | 0.5942 | **0.5942** |
| X q_sum reco ratio | 0.6416 | **0.6416** |
| Y peak / q_event / q_sum | 0.5278 / 0.5599 / 0.6098 | **identical** |
| X rise median, sim → data | 337.4 → 285.2 | 340.4 → 298.9 |
| **X rise ratio sim/data** | **1.183** | **1.139** |
| Y rise median, sim → data | 307.4 → 261.8 | 309.4 → 267.8 |
| **Y rise ratio sim/data** | **1.174** | **1.156** |
| X FWHM ratio sim/data | 1.188 | 1.184 |

**Every amplitude number reproduces to four decimals.** That is the check that
the censor does what it claims and nothing more: it is confined to the shape
observables, and the §7 ledger — which is quoted in q_sum — is untouched by
this correction and moves only for the transparency fix. It also re-validates
the data leg after it was copied to durable storage (`data_leg_satdet3/`):
same inputs, same answer.

The vertical rise excess falls from ×1.183 to ×1.139 (X) and ×1.174 to ×1.156
(Y). Smaller than the ~10 ns move at inclined incidence in relative terms,
because the vertical distribution is ~280 ns wide and a shift of that size is a
smaller fraction of it.

Censored and undefined fractions, both views:

| | sim | data |
|---|---|---|
| X censored (railed) | 5.9 % | **26.0 %** |
| Y censored (railed) | 3.8 % | **11.0 %** |
| X undefined rise | 0.0 % | **9.0 %** |
| Y undefined rise | 0.0 % | **3.0 %** |

## Two further things this turned up

1. **The vertical DATA leg loses 9.0 % (X) and 3.0 % (Y) of events to an
   undefined rise** (no 10 % or 90 % crossing inside the window) against the sim
   leg's 0.0 % on **both** views. That is a one-sided drop on a comparison that
   is supposed to be like-for-like, and it is a real shape difference rather
   than a cut — the vertical data leg is already known to be anomalously broad
   (§3). It is the same X-heavier pattern as the railing (26.0 % vs 11.0 %),
   which is at least consistent with one cause rather than two. Negligible at
   10° and 20° (0.1–0.5 % both legs), so it does not touch the angled tables
   above. Worth a look on its own.
2. **The 20° data leg is thin**: 811 events uncensored, 615 after the censor,
   against ~2400 on the sim legs. Quantile noise there is not negligible and
   the p5 column in particular should not be over-read.

## Deliberately not done

The estimator is **censored, not corrected**. A proper treatment would fit the
rising edge below the rail and extrapolate the true peak, recovering the railed
quarter of the data instead of discarding it. That is a better estimator and a
bigger change; this is the minimal defensible fix, and it throws away a quarter
of the data leg to get there. If the rise thread needs the statistics back, that
is the next step.
