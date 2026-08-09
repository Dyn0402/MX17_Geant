# T14 freeze queue — RESOLVED 2026-08-09

> ## ✅ DEFAULT FROZEN by Dylan, 2026-08-09, before any comparison
>
> **ρ_s = 2 MΩ/sq · DRY 95/5 table · det3 data bundle as-analysed.**
>
> The first comparison runs against this default and **that verdict stands**.
> Iteration afterwards is **pre-declared, not improvised**: ρ_s and the gas/v
> axis are the two variables to vary for agreement and understanding, and
> everything after the default run is **diagnosis, not verdict**. Declaring the
> axes now is what keeps the first look blind and the later looks honest.
>
> Pre-built so later comparisons need no new production:
> * Stage B/C decoded at **all four ρ_s** (the W2 kernels already serve them).
> * A **v-axis reco variant** needing no new simulation at all — v enters at
>   reconstruction, so the same decoded files re-run with a table-v bundle.
> * **Wet gas** is the only axis needing real production, and is
>   **prepared-not-run** (§3 below).

*(The analysis that led to the freeze is kept below, unchanged, because it is
the evidence the decision was made on and the record of when.)*

---

# The two inputs, as they stood before the freeze

**For Dylan. Two decisions, same kind, same discipline.** Both are simulation
*inputs* that the chain does not determine, and for both there is a value that
would obviously make T14 agree better. That is exactly why they have to be
frozen on independent evidence and written down **before** anyone looks at the
comparison — otherwise T14's answer has been laundered into its own input.

Neither is blocking T13, which runs on the data's own calibration bundle
whatever you choose. Both change what T14 means.

---

## 1. ρ_s — the ESL sheet resistivity

**Recommendation: freeze at ρ_s = 2 MΩ/sq.** rho1M sits below the allowed band.

ρ_s has never been measured on our modules; the plan scans {0.5, 1, 2, 5}
MΩ/sq for exactly that reason. The one independent constraint is T2b's charge
spread, which bounds the *product* ρ_s·c′. Re-quoting it for the production
insulator (50 µm kapton + 18.76 µm glue in series, `d_eff` = 70.5 µm) scales it
by c′(75)/c′(70.5) = 0.947:

| stack | band |
|---|---|
| bare 75 µm kapton (as the plan quotes T2b) | 1.5 – 2.7 MΩ/sq |
| bare 50 µm (plan's own re-quote) | 1.1 – 1.9 |
| **production, d_eff = 70.5 µm** | **1.42 – 2.56** |

The mapping is validated against the plan's own numbers: applying it to bare
50 µm reproduces 1.04–1.88 against the 1.1–1.9 the plan states.

⚠️ **The "ρ_s = 1 MΩ/sq nominal" line elsewhere in the plan is a pre-glue
artifact** (2026-08-07, bare d_k = 75 µm) and is not a T14 decision.

**Stakes:** ρ_s alone moved the T10 slow-path residual by 0.65 pp — comparable
to the entire W1→W2 boundary-model change. Stage B/C has been produced at
**both** rho2M and rho1M, so nothing is forced by what happens to exist.

## 2. Gas — dry vs wet Ar/iso 95/5 (plan P1 sub-item)

**No recommendation; this one is genuinely yours.** Both options are defensible
and the distinction between them is the whole point.

Stage B's Magboltz table gives **v_drift = 39.14 µm/ns** at 333 V/cm. The det3
bench bundle for the T14 target run measures **36.60 µm/ns** — the simulation
drifts **6.9 % fast**. The plan already flags this ("+~1 % H2O on the bench —
matters, measured v_drift 36.6 µm/ns is far below dry Magboltz"; water fraction
listed as an open P1 sub-item). This note only makes it a number.

Consequence if left as is: reconstructed depth is scaled ~6.5 % against the
sim's own truth, because reco converts time→depth with the data's v_drift.

- **(a) Keep the dry table**, and quote 6.9 % as an input systematic on every
  depth-dependent comparison. This is what T13 does today.
- **(b) Use a wet (~1 % H₂O) 95/5 table for Stage B.** This is **legitimate
  physics input if and only if** the water fraction comes from the independent
  June bench record — the measured Ar/iso + ~1 % H₂O, det3 drying 3 % → 1 %
  over a week, and the waveform-first result that v(E) matches wet-Magboltz
  95/5. It is **laundering** if the water fraction is chosen because it closes
  the 6.9 %.

The difference between (b)-as-physics and (b)-as-tuning is not in the table
used; it is in what the number was derived from and when it was written down.

---

## 3. Wet gas — an UNCONSTRAINED search axis, not a physics input

> ### ⚠️ Correction, Dylan 2026-08-09 — read this before the rest of §3
>
> **"We have NO measured humidity at any point. We only assume there is humidity
> or some other contamination in the gas because our drift velocity is slow. So
> if dry gas does not agree with data, we will have to search for contaminants
> which match our data unconstrained — not ideal but all we can do."**
>
> This **invalidates the framing used earlier in this note and in §2 above.**
> Both were written as though a "June bench humidity record" existed
> independently of the comparison. It does not. The ~1 % H₂O figure and the
> det3 "dried 3 % → 1 % over a week" history were themselves **inferred from
> v_drift by Magboltz matching** — i.e. derived from the same observable family
> the T14 comparison uses. There is no hygrometer reading behind them.
>
> So the (a)/(b) structure below is wrong where it says a wet table "is physics
> input if the water fraction comes from the June record". Corrected:
>
> * **Default stays dry** — frozen, unchanged, and now on firmer ground: it is
>   the only option not fitted to the observable.
> * **Any contaminant hypothesis is an unconstrained search**, run only after
>   the default comparison, and **labelled fitted-to-data wherever it appears**.
>   It is diagnosis. It cannot become a physics input by being plausible.
> * **The way out is an actual measurement** — a hygrometer on the gas line, or
>   any independent species assay. One line of hardware converts this whole axis
>   from fitted to constrained, and it is worth doing for that reason alone.
> * **One real constraint survives, on FAMILY not concentration.** The July
>   90/10 study discriminated air/O₂ from H₂O by **attachment shape**, not by
>   v_drift alone. That is still inference from waveform data rather than a
>   measurement, but it constrains *which contaminant*, independently of how
>   much. `eta_per_cm` in the wet grid is the handle; record it as a shape
>   constraint, never as a measurement.
>
> The numbers below stand — what changes is what may be concluded from them.

### The mechanics (unchanged)

**No Magboltz run is needed.** The June waveform-first water study already left
the grid, in Stage B's *exact* table schema:
`~/PycharmProjects/nTof_x17/garfield_sim/results/water2d.json` — 30 Ar/iso/H₂O
mixtures, each a list of `{E_Vcm, v_um_per_ns, eta_per_cm, dL_sqrtcm,
dT_sqrtcm}`. Stage B's own table
(`design/gas/drift_velocity_Ar_iC4H10_95_5_Saclay.json`) is the same record type
minus `eta_per_cm`, which the wet grid *adds* — Stage B's `survival()` already
supports an attachment column and currently gets exactly 1 from the dry table.
So building a wet table is an extraction, not a campaign.

**But the wet table does not close the 6.9 %; at the recorded water fraction it
overshoots.** v_drift at 333 V/cm:

| mixture | v [µm/ns] | vs bench 36.60 |
|---|---|---|
| dry 95/5 — Stage B today | 39.14 | **+6.9 %** |
| iso5 + 0.4–0.55 % H₂O | 40.2–41.0 | +10 to +12 % |
| iso5 + **0.95 %** H₂O | **34.81** | **−4.9 %** |
| iso5 + 1.05 % H₂O | 32.99 | −9.9 % |
| iso5 + 0.8 % H₂O + ~1 % air | 36.24 | −1.0 % |

The June record's **~1 % H₂O gives 34.81 µm/ns — 4.9 % too SLOW**, i.e. it
overshoots the correction rather than removing it. The value that would match
is ≈0.8 % H₂O (or 0.8 % with ~1 % air, 36.24).

**Where the line sits, restated after the correction above.** There is no
water fraction that counts as physics input today, because none was measured —
so *every* point on this axis is fitted, including ~1 %. What the table above
shows is that the fit is not even flattering: the previously-assumed ~1 %
overshoots to −4.9 %, and only ≈0.8 % lands on the bench value. That the
"assumed" and "matching" values differ by 0.2 % is not reassurance that the
assumption was nearly right; it is a reminder that a number inferred from
v_drift will always sit near whatever v_drift needs, which is exactly why it
cannot be used to justify v_drift.

The air row is the one that carries real information: 0.8 % H₂O + ~1 % air
reproduces the bench value closely, and air is a physically **distinct
hypothesis** from water with its own attachment signature — oxygen attaches,
water does not. That distinction is testable against `eta_per_cm` /
`attachment_Ar_iso_H2O.json` **without** reference to v_drift, which is what
makes it the only part of this axis capable of constraining anything. The July
90/10 study used exactly that shape argument to exclude air for *that* run;
nothing has tested it for this one, and doing so would narrow the search from
"any contaminant" to a family — still short of a measurement.

**Not run now, by instruction.** What it would take when wanted: extract one
mixture from `water2d.json` into the Stage B table schema, re-run Stage B/C
(~19 min per point), and optionally a wet avalanche point against the
gas-agnostic map ladder.

## First diagnostic result on the v axis — it is NOT the v mismatch

Built as pre-declared (no new simulation: v enters at reconstruction, so the
same decoded files were re-run with the bundle's `v_drift` swapped 36.60 →
39.1424, the dry-table value Stage B actually generated with). Frozen default
vs v-axis variant, 2 980 events each:

| | x fitted | x slope_reliable | y slope_reliable | median χ²/dof (x, y) |
|---|---|---|---|---|
| **default** (bundle v = 36.60) | 98.4 % | 4.6 % | 13.4 % | 19.3, 20.4 |
| variant (v = 39.14, mismatch removed) | 98.4 % | 4.1 % | 10.5 % | 19.3, 20.4 |

**Removing the 6.9 % v mismatch changes fit quality not at all** — χ²/dof is
identical to three figures on both views — and makes slope reliability slightly
*worse*. So whatever drives χ²/dof ≈ 20 and the low `slope_reliable`, **it is
not the drift-velocity mismatch.** That is a useful negative: it removes the
most obvious suspect from the list before the comparison is even run, and it
means the dry-vs-wet gas axis is unlikely to rescue fit quality either, since
its whole effect on reconstruction enters through v.

Two cautions on reading this. It says nothing about whether the *default* is
right — only that this particular knob does not move these particular numbers.
And these are fit-quality diagnostics, not the §9 observables; the comparison
itself remains frozen and unrun.

## What has already been decided, for contrast

These were frozen before looking and are recorded with their evidence, which is
the pattern the two above should follow:

- **Avalanche calib** — pooled 490 V meshfield point, chosen on the grounds
  that all 56 slices used the same field map regardless of their voltage label
  (not on downstream agreement). Reconfirmed 2026-08-08.
- **Drift gap** — 30 mm, decided 2026-08-08, closed.
- **T14 target** — det3 `long_run_resist_490V_drift_1000V`, per P1.

## What is NOT on this list

The W2 kernel itself. It is not a free input: V6 showed W1 grounds copper that
does not exist, and W2 fixes an error of known sign and size (+25.66 % prompt
capture, pre-registered and reproduced). Choosing W1 over W2 would be choosing
a known-wrong boundary condition, not exercising a scan axis.
