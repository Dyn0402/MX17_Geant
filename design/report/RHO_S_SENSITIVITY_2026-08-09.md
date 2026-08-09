# ρ_s sensitivity — can a wrong sheet resistivity explain the T14 discrepancies?

**2026-08-09. Verdict: no.** Not for amplitude (wrong size AND wrong
signature), not for rise (wrong direction), and for sharing — the one
observable where the direction is right — the required value is excluded by
the same measurement that constrains sharing. Question asked by Dylan after
the T14 campaign; this note is the quantitative closure.

Companion script: `response/validation/rho_s_sensitivity.py` (numbers
regenerate in ~15 min from the S1 archive). Raw output:
`rho_s_sensitivity_2026-08-09.json` alongside this file. Context: the T14
campaign index (`nTof_x17 mx17_sim_wft/T14_CAMPAIGN_2026-08-09.md`) and the
freeze record (`T14_FREEZE_QUEUE_2026-08-09.md`).

## Why the question is live

ρ_s has never been measured on our modules — it is a scan axis frozen at
2 MΩ/sq before the comparison. The only independent constraint is T2b's
spread measurement (τ_X ≈ 230 ns / τ_Y ≈ 410 ns to spread ~one pad pitch),
which bounds the **product** ρ_s·c′, re-quoted for the production stack as
**ρ_s = 1.42–2.56 MΩ/sq**. And T14's misses all point the way a higher ρ_s
naively pushes: sim amplitude low, sim over-sharing, sim rise slow. So the
hypothesis deserved numbers, not intuition.

## Method

Unit point charge on the ESL at y = 0, x averaged over 8 positions across two
pad pitches → `CombKernelLUT` per-channel induced currents (n_side = 8,
t_max 2.4 µs, dt 2 ns) → `DreamShaper` (180 ns peaking, β = 0.75). Metrics
per kernel: central-channel shaped peak (fraction of deposited charge),
10–90 % rise of the central shaped waveform, and the peak-amplitude sharing
profile across channels (count above 2 %/5 % of max, RMS width in strips).
Run prompt-only (electron delta) and with the analytic ion rectangle
(f_electron = 0.092, 340 ns transit) — the ion-folded rows are the
production-like case and are what is quoted below.

The ρ_s ladder is the **W1 dk50 family** {0.5, 1, 2, 5} MΩ/sq because W2
exists only at rho2M; cross-ρ_s **ratios** within one family are clean, and
the W1-vs-W2 rho2M pair is included to show the boundary/glue change is a
ρ_s-independent offset (~12 % on amplitude), not a scaling.

## Numbers (ion-folded)

| kernel | amp₀ X | amp₀ Y | rise X [ns] | rise Y [ns] | Y n>2% | Y rms [strips] |
|---|---|---|---|---|---|---|
| W1 rho0.5M | 0.205 | 0.148 | 243 | 196 | 11.5 | 2.00 |
| W1 rho1M | 0.208 | 0.170 | 239 | 202 | 8.5 | 1.58 |
| W1 rho2M | 0.217 | 0.197 | 239 | 211 | 5.8 | 1.17 |
| W1 rho5M | 0.235 | 0.238 | 257 | 228 | 5.2 | 0.75 |
| W2 rho2M g19 (production) | 0.199 | 0.174 | 238 | 207 | 8.2 | 1.35 |

X-view sharing is geometric (comb) and flat at ~2 channels throughout.

## The three hypotheses, closed

**1. "Higher R would increase our amplitudes" — right direction, wrong size,
wrong signature.** A factor **10** in ρ_s buys ×1.15 (X) / ×1.61 (Y) on the
shaped peak; within the T2b band, 2 → 2.56 MΩ/sq buys ~5 %. The T14 deficit
is ×1.7–1.9. Decisively: the deficit is **voltage-dependent** (data
d lnA/dV = 0.449/10 V vs sim 0.296, ≈12 σ — `hv_slope/HV_SLOPE_2026-08-09.md`)
and a static sheet property cannot produce a gain-slope error. The amplitude
thread stays with α(E)/Penning (the T7 slope hunt).

**2. "Higher R would make rise times faster" — backwards.** Rise gets
**slower** with ρ_s: ion-folded Y 196 → 228 ns, X 243 → 257 ns over the
factor 10 (the central channel's charge keeps arriving on the sheet
timescale, and slowing the sheet slows it). Data wants faster rise, and the
rise floor was already measured to be the ion term
(`DIAGNOSIS_noions`: remove ions and the rise distribution lands on data at
every quantile; β immune at 4 ns across its full range). Rise belongs to the
f_eff ≈ 0.2-vs-0.9056 contradiction, not to the sheet.

**3. "Higher R would decrease spreading" — yes, and this is the only T14 miss
with ρ_s's signature** (noisefix unmasked sim 11 hits/event vs data 7.4 at
low amplitude; here Y rms 2.0 → 0.75 strips over the factor 10). But matching
the data's sharing needs roughly ×4 in ρ_s ≈ 8 MΩ/sq — far outside the
1.42–2.56 band — and at that value the T2b spread times would be ~2.5× slower
than measured. The sharing-relevant quantity D = 1/(ρ_s·c′) is **exactly what
T2b pins**, so ρ_s cannot fix sharing without contradicting the spread-time
data it would have to explain away first.

## A nugget worth keeping (tension, not a conclusion)

The sim-side X/Y amplitude asymmetry — flagged in the HV-slope audit
(unselected data amplitudes are *identical* X/Y, 2606.5/2607.0 ADC) —
**shrinks monotonically with ρ_s and vanishes at 5 MΩ/sq** in this test:
ion-folded X/Y ratio 1.38 (0.5M) → 1.10 (2M) → 0.99 (5M). One sub-observable
votes for higher ρ_s, in direct tension with T2b. If pursued, the
reconcilable question is not "is ρ_s higher" but **"is c′ wrong (stack), or
is the T2b τ extraction biased"** — only the product is constrained. Given
the deep-Y-undershoot return-path item is also resistive-layer-side, the Y
response *model* is the more likely culprit than the ρ_s *number*.

## Caveats

Point deposit, not a track (no drift-diffusion smear, so absolute sharing
widths here undershoot the full-chain values; the cross-ρ_s ratios are the
result). The Y-channel parity assignment alternates with row offset and is
averaged over deposit x — second order for ratios. β fixed at 0.75; β moves
the peak ≤2.3 % (freeze queue) so it cannot change any conclusion here.
