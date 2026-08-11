# ⚠️ CORRECTION 2026-08-11 — the mesh transparency was applied twice

**Stage B thinned every drifting electron at 0.873 × 0.9559 = 0.8345 where the
physics is a single 0.955. Simulated charge was low by ×0.874 from 2026-08-08,
when the meshfield calibration became production, to today.**

Affects: every Stage B/C product built on `aval_calib_meshfield_pooled.json` —
which is the T14 frozen default, the angled ladder, and all the DIAGNOSIS legs.
Does **not** affect the S3 avalanche calibrations themselves, the field maps,
the S1 kernels, or anything on the data leg.

---

## The defect

Two independent measurements of the *same* quantity were multiplied together.

| | value | what it measures |
|---|---|---|
| `digitize.DEFAULT_TRANSPARENCY` | 0.873 | T6 **1-D** mesh transparency (`mesh_transparency.C`, 2026-08-07) |
| calib `polya.survival` | 0.9559 ± 0.0026 | fraction of seeds launched 180 µm **above the mesh** that produced any avalanche |

The second number is the mesh transparency, measured a second way. The
meshfield S3 calibration seeds electrons above the mesh precisely so that
transparency and funnelling come out *emergent* rather than assumed
(`mx17_aval_calib.py:270-277`, and its own comment at `:474` says so:
"seed electrons above the mesh, so `survival` above already reflects
transparency losses"). A seed absorbed on a wire records gain 0 and lands in
`(gains > 0).mean()`. The residual factor — P(g > 0) for a seed that *did* pass
— is 1.0 to the precision of the uniform-field campaign: 0 of 6400 seeds failed
to multiply across all 56 slices (`DESKTOP_RUNS_2026-08-07.md`).

`digitize.py:349-352` multiplied both into one binomial thinning.

## Why it survived

Three things had to line up, and did.

1. **The runbook's acceptance note deferred the fix and nobody collected it.**
   `FIELD_MAP_RUNBOOK.md:150-159`, accepting the production map on 2026-08-08:
   "The production G7 value is the T6 transparency deliverable and supersedes
   the 1D-model 0.873 (`DEFAULT_TRANSPARENCY` in Stage B) **once accepted**"
   … "supersedes DEFAULT_TRANSPARENCY on Stage B's next touch". Stage B was
   touched repeatedly over the next three days. The constant was never changed.

2. **The agreement was read as a cross-check.** `RESPONSE_SIM_PLAN.md:1341`
   records "survival 0.9559 ≈ T6 transparency 0.955" as a corroboration, and
   `OVERNIGHT_2026-08-10.md` §7.4 puts it in the elimination table as
   > Mesh transparency | 0.955, cross-checked | **Unlikely.** T6's 3D
   > transparency agrees with the avalanche survival 0.9559 independently; no
   > room for ×1.5.
   **That row is wrong and is corrected here.** The two numbers are not an
   independent cross-check; they are one quantity measured twice, and their
   agreement is exactly the signature of the duplication. The row's *verdict*
   still stands on size — a duplicated 0.955 buys ×1.15, not ×1.5 — but its
   reasoning was inverted and it had a real number hiding in it.

3. **The code comment described a calibration that was no longer production.**
   `digitize.py:340-345` stated that the survival factor "IS 1.0 at every
   voltage … so this is exact bookkeeping today, not a correction". True of the
   uniform-field S3 calibs it was written against; false of the meshfield calib
   that replaced them. Read literally it licensed the multiplication.

## The fix

`response/digitizer/digitize.py`, new `split_calib_survival(calib_pt)`. The
calibration's single `survival` number is routed to the factor it actually
measures, keyed on `field_model` — the same key `calib_seed_z0_mm` already uses
for the same class of reason (a meshfield calib's σ₀ already contains a drift
leg that must not be added twice):

| calib | mesh term | P(g>0) term | net thinning |
|---|---|---|---|
| meshfield (production) | `survival` = 0.9559 | 1.0 | **0.9559** (was 0.8345) |
| uniform_field | `DEFAULT_TRANSPARENCY` | `survival` (or 1.0 if absent) | 0.955 (was 0.873) |

Two further changes come with it:

- `DEFAULT_TRANSPARENCY` **0.873 → 0.955**, collecting the supersession the
  runbook asked for on 2026-08-08. This only reaches uniform-field calibs, so
  it does not touch production; it changes any legacy v2/v3 leg by ×1.094.
- In-gap deposits (past the mesh, never drifted) no longer get the meshfield
  `survival` applied. They previously took 0.9559 where the correct factor is
  P(g>0) = 1.0. Worth ~0.4 % of the electrons — real, tiny, and now right.

Regression guard: `response/digitizer/test_transparency_split.py`, which asserts
the net mesh thinning equals the calibration's own number for a meshfield calib,
agrees with T6's G7 0.955, and that no shipped calibration comes back
double-counted. It stubs the LUT so it runs on any host.

```
  production (meshfield)   eps=0.955937 P(g>0)=1.0  net=0.955937   was 0.834533  ->  charge x1.1455
  uniform_field            eps=0.955000 [external, T6 3-D]  net=0.955000
  4 shipped calibrations, none double-counted
```

## Pre-registered consequence for the §7 ledger

Written before the re-run was submitted (condor cluster 11988754, variant
`TRANSPFIX`, `rho2M`, otherwise identical to the frozen default). The fix is a
constant ×1.1455 on surviving electrons, so on a **linear** chain:

| §7 row | frozen | predicted with the fix |
|---|---|---|
| q_sum ratio, X | 0.6416 | **0.735** |
| q_sum ratio, Y | 0.6098 | **0.699** |
| demanded gain, X | 37 552 | **32 780** |
| demanded gain, X, saturation-corrected | 40 243 | **35 130** |
| amplitude deficit | ×1.66 | **×1.46** |
| HV slope d lnA/dV | 0.296/10 V | **unchanged** — a constant factor cannot tilt a slope |

The chain is not perfectly linear: a 14.5 % taller sim waveform rails more
often, and railed events are dropped by the reco-quality cut. So the measured
move should fall **short** of ×1.1455 on peak and track it more closely on
q_sum. A move *larger* than ×1.1455 on either would mean the fix did something
other than what it says, and should be treated as a failure, not a bonus.

## What this does and does not change

- **It does not resolve the amplitude deficit.** ×1.46 still demands a real
  gain near 3.3–3.5 × 10⁴ against the sim's 24 094, still at the top of the
  literature's stable Ar/iso bulk-MM band. The deficit shrinks by 12 %; it does
  not go away.
- **It does not touch the HV slope**, which is the sharper discrepancy
  (0.4487 vs 0.3107 per 10 V) and the one Penning outcome C could not close.
- **It does not touch the rise-time contradiction** at all: charge scaling is
  blind to time structure. See `RISE_SATURATION_CENSOR_2026-08-11.md` for that
  thread's own correction, which is likewise a numbers correction and not a
  resolution.

## Results — 2026-08-11 18:00. **Every pre-registered number landed.**

Condor 11988754 completed 17:14. The cluster log confirms the fix reached the
run: `transparency=0.9559375 [S3 meshfield calib survival …]`,
`aval_survival=1.0` — exactly one mesh term. Reconstructed through the same
`t13_reco.py` driver and compared against the same data leg as the frozen
default; sim and data legs both censored (`RISE_SATURATION_CENSOR_2026-08-11`),
so this is censored-to-censored.

| observable | frozen | **TRANSPFIX** | move | predicted |
|---|---|---|---|---|
| X q_sum reco ratio | 0.6416 | **0.7344** | ×1.1447 | 0.735 |
| Y q_sum reco ratio | 0.6098 | **0.6929** | ×1.1364 | 0.699 |
| X peak ratio | 0.5646 | **0.6462** | ×1.1445 | — |
| Y peak ratio | 0.5278 | **0.6027** | ×1.1421 | — |
| **X demanded gain** | 37 552 | **32 805** | | 32 780 |
| X demanded gain, saturation-corrected | 40 255 | **35 167** | | 35 130 |
| **X deficit** | **×1.67** | **×1.46** | | **×1.46** |
| Y deficit, saturation-corrected | ×1.94 | ×1.71 | | — |
| HV slope | — | unchanged by construction | | unchanged |

The X deficit lands at ×1.460 against a prediction of ×1.46 written before the
job ran. q_sum — the statistic the ledger is quoted in — moved ×1.1447 against
the linear ×1.1455, i.e. 0.07 % short, with the shortfall in the predicted
direction: sim saturation rose from 5.9 % to 8.0 % (X) and 3.8 % to 5.1 % (Y),
and railed events are dropped by the reco-quality cut.

### The pre-registered falsifier was mis-specified, and says so

It read: *"a move LARGER than ×1.1455 on either would mean the fix did
something other than what it says."* **One statistic trips it** — the q_event
median moved ×1.157 (X) and ×1.153 (Y). That is not the fix misbehaving; the
falsifier assumed a single sign and there are two competing effects of
comparable size:

- **Threshold inclusion pushes q_event UP.** `q_event` sums peak amplitude over
  channels above 5 σ, so raising every amplitude 14.55 % promotes channels that
  were below. Measured: mean `n_over` 4.975 → 5.057 (X, ×1.017) and 7.359 →
  7.584 (Y, ×1.031). Real, and absent from `q_sum`, which is fitted charge.
- **Railing pushes it DOWN**, as above.

Their net differs between mean and median of a skewed, rail-truncated
distribution: the q_event *mean* moved ×1.118 (X), *below* linear, while its
median moved ×1.157, above. Recording this rather than quietly dropping it —
the criterion was written too crudely, the numbers are behaving, and which of
those it is should be decided on the mechanism rather than on the fact that the
headline came out right.

### What it changes

**The amplitude deficit is now ×1.46, not ×1.66, and demands a real gain of
~3.5 × 10⁴ against the sim's 24 094.** That is still at the top of the
literature's stable Ar/iso bulk-MM band (3–4 × 10⁴), so §7.3's reading survives
in substance — a pure gain explanation remains *physically available and only
just*. The deficit shrank by 12 %; nothing was resolved.

**Unchanged: the HV slope**, which is the sharper discrepancy (data 0.4487 vs
sim 0.3107 per 10 V) and the one Penning outcome C could not close. A constant
multiplicative factor cannot tilt a slope, and this fix is exactly that.

Incidentally, the sim/data **rise** ratio also fell, X 1.139 → 1.106 and Y
1.156 → 1.162 — a second-order consequence of more sim events railing, not a
statement about the ion term. It is not evidence on the f_eff thread and should
not be read as any.

### Products

- `~/x17/response_sim/stageB_w2/w2_rho2M_TRANSPFIX/default/` — decoded + events
- `~/x17/response_sim/stageB_w2/t14_compare_transpfix/` — the comparison
- `~/x17/response_sim/stageB_w2/t14_compare_censored/` — the frozen leg re-run
  with only the censoring fix, i.e. the correct baseline for the table above
- `~/x17/response_sim/stageB_w2/data_leg_satdet3/` — the data leg, moved out of
  a `/tmp` session scratchpad where it was its own single point of failure

### ⚠️ A trap found on the way

`w2_rho2M/default/calib_bundle/` **is the v-table variant bundle, not the
frozen default's** — the variant ran second (13:29 vs 13:24) and overwrote it.
Reconstructing that directory with the bundle inside it silently reproduces the
variant at v = 39.14 µm/ns instead of the default's 36.60: same event count, no
error, one differing log line. Caught here by reading the log. The correct
bundle survived only inside a `/tmp` scratchpad and is now at
`data_leg_satdet3/calib_bundle/`; a `README_BUNDLE_WARNING.md` sits in the
offending directory. `events_default.meta.json`, not `calib_bundle/`, is the
authority on what the frozen leg used.
