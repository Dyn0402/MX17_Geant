# Effective amplification gap, 135 vs 150 µm — PRE-REGISTRATION

**Written 2026-08-11, BEFORE any job was submitted.** Predictions, decision
rule and the falsifier are all below the line and are not to be edited after
the campaign lands; corrections go in a dated section at the bottom.

---

## Why

Nothing pins det3's as-built amplification gap. The only input is the nominal
150 µm pillar height from the gerber; T6's field map implies ~158 µm effective
(`RESPONSE_SIM_PLAN.md:1340`), and det4's beam assessment shows that a varying
mesh height is a real failure mode in this production series. The gap enters
the T14 amplitude deficit directly — a shorter gap at fixed mesh voltage is a
higher field, hence more gain — and it enters the HV slope, which is the
sharper of the two discrepancies:

| | d lnG/dV [per 10 V] |
|---|---|
| det3 data | **0.4487** |
| sim, 150 µm | 0.3107 |

## ⚠️ What is actually being run, and what is blocked

**The meshfield route is BLOCKED. This campaign runs the UNIFORM-FIELD arm.**

`response/meshcell/solve_fieldmap.py` does **not** take the gap as a parameter:
`AMP_GAP = 150.0` is a module constant (line 71) and there is no CLI flag for
it. Three things stand between here and a 135 µm map, and none should be
improvised:

1. **No independent gate at 135 µm.** Acceptance gate G2 compares the solved
   amp-bulk field against `NEBEM_AMP_EZ_AT_490V = 30561 V/cm`, a neBEM
   cross-check computed for the 150 µm geometry. At 135 µm that reference does
   not exist, so the map would either fail G2 or be checked against a number
   rescaled by the very solver under test.
2. **It is not a condor job.** The map solve is gmsh + scikit-fem on the
   desktop, and reading the map back needs the LOCAL Garfield `ComponentGrid`
   3-D region patch (without it transparency reads 1.000 — runbook Decision 5).
   The desktop was behind a Tailscale re-auth wall as of the 08-10 brief.
3. **It is a ladder, not a map.** The per-voltage maps
   (`meshfield_vmesh0460/0490/0520.txt`) came from `overnight_chain.sh`; a
   second gap means a second 3-map ladder, gated.

The uniform-field branch of `mx17_aval_calib.py` **is** genuinely parameterised
by `--gap-um` (`e_field = voltage / gap_cm`, `gap_cm` also sets the weighting
field and the drift region), so it is the branch that can answer this today.

**And the field-map shape has just been exonerated for the slope.** The
slope-hunt's own field-shape A/B arm landed on EOS this morning
(`aval_calib_slopehunt.json`, 03:27):

| arm | 460 V | 490 V | 520 V | d lnG/dV per 10 V |
|---|---|---|---|---|
| meshfield (HV scan, auto) | 9 878 | 24 172 | 61 715 | **0.3054** (0.3106 over all 8 V) |
| uniform field (auto) | 17 062 | 43 395 | 110 085 | **0.3107** |
| mesh / uniform | 0.579 | 0.557 | 0.561 | — |

**The mesh map costs a flat ×0.56 on gain and changes the slope by 0.03 %.**
So for the slope question — which carries the discrimination below — the
uniform-field arm is the right instrument.

The offset is flat in voltage too: ×0.579 / 0.557 / 0.561 across 460–520 V, a
4 % spread over a range across which the gain itself moves ×6.5. So the mesh
map is, to 4 %, a **multiplicative constant in V** — which is a second and
independent reason it cannot contaminate a *slope* measured at either gap, and
is worth stating because the uniform arm's validity is the thing most likely to
be disputed later.

**It does NOT follow that the mesh map factors out of the gap ratio.**
Voltage-independence says the factor is constant along V; it says nothing about
whether it is constant along *d*, and there is a physical reason to expect it is
not — the mesh near-field occupies a fixed ~18.5 µm, which is 14 % of a 135 µm
gap against 12 % of a 150 µm one. The ×0.56 is measured at 150 µm and is
*assumed*, not measured, to carry to 135 µm. That assumption is this campaign's
main caveat and is restated in the decision rule.

## The campaign

Ar/iC4H10 95/5 dry (`Ar_iC4H10_95_5_Saclay_160m.gas`), Penning `auto`, uniform
field, 2 gaps × 3 voltages × 8 seeds = **48 slices**, seeds 81001–81048 all
distinct. The 150 µm arm is **re-run** rather than taken from the slope hunt,
so both arms share a host and a Garfield build (pin `927e5c21`) — the slope
hunt's uniform arm ran on the desktop, and a cross-host leg is exactly the kind
of difference this project has been bitten by. Seeds are not paired across
gaps: the two trajectories diverge at the first collision regardless, so
pairing buys no variance reduction and only risks a filename collision.
`nev` per slice is 230/140/60 at 460/490/520 V for the 150 µm arm (the slope
hunt's own tuning) and 145/90/38 for 135 µm, scaled down by the predicted
×1.57 gain so a slice costs the same wall clock at either gap.

---

# PRE-REGISTERED PREDICTIONS

Derived from the existing 150 µm uniform arm alone, by fitting α(E) = ln G / d
over its three points and re-evaluating at E = V / 135 µm. Nothing here is
fitted to a 135 µm number, because none exists yet.

α(E) from the 150 µm arm: 649.6 / 711.9 / 773.9 cm⁻¹ at 30.67 / 32.67 /
34.67 kV/cm — **linear in E to better than 0.5 %** (quadratic term −0.021).

| observable | 150 µm (measured) | 135 µm (PREDICTED) |
|---|---|---|
| gain at 490 V | 43 395 | **68 000 ± 4 000 (×1.57)** |
| d lnG/dV per 10 V | 0.3107 | **0.309 ± 0.010 — essentially UNCHANGED** |
| ion transit t90 | — | **×0.81** (t ∝ d²/µV, exactly) |

**The slope prediction is the sharp one, and it disagrees with the informal
expectation this scan was proposed under (slope falling to 0.26–0.29).** In a
parallel-plate Townsend picture with lnG = α(E)·d and E = V/d,

    d lnG / dV  =  α'(E)

*exactly* — the explicit gap cancels, and only the E at which α' is evaluated
moves. Since α is linear in E across this range, α' is the same constant and
the slope does not move. Recovering the measured 0.3107 from α' alone is the
consistency check that this is the right decomposition, and it does (0.3107).

## Decision rule

- **Slope stays at 0.31 ± 0.01 and gain rises ×1.5–1.6** → predicted outcome.
  The gap is a *gain* lever and **not** a slope lever, so it cannot close the
  0.4487 vs 0.3107 discrepancy at any gap. Combined with Penning outcome C
  (2026-08-10: Penning delivers gain or slope, never both), that leaves the
  slope with no surviving single-parameter explanation inside the avalanche
  model, and the search moves outside it.
- **Slope falls to 0.26–0.29** → **this is the competing pre-registration**,
  written independently in the MX17_Geant deep-dive session that proposed this
  scan, from an exponential-Townsend extrapolation which carries α″ < 0 where
  the linear fit above carries α″ = 0. The two agree across the measured
  30.7–34.7 kV/cm and diverge only at the 135 µm arm's 36.3 kV/cm, so **the
  scan adjudicates a real sub-question — the curvature of α(E) above 35 kV/cm —
  on top of the headline.** Both are recorded as written; neither is edited
  after the fact. Note the decision content is invariant to which wins: on
  either prediction the gap moves the slope away from the data's 0.4487 or not
  at all, so **the gap cannot close the slope discrepancy under either.**
- **Slope RISES toward 0.4487** → the parallel-plate Townsend framing is wrong,
  and the Penning outcome-C verdict rests on that framing. This is the most
  important possible result and the one to check hardest before believing.
- **Any arm's gain ratio outside ×1.3–1.9** → suspect the run, not the physics;
  check that `gap_um` actually reached the config block of every slice.

## Falsifier for the campaign itself

The 150 µm arm re-run here must reproduce the slope hunt's 150 µm uniform arm
(17 062 / 43 395 / 110 085) to within seed statistics. If it does not, the
lxplus-vs-desktop leg differs and **nothing in this campaign may be quoted**
until that is explained.

## What this cannot do

- It cannot measure the real gap. It measures what the model predicts *if* the
  gap were 135 µm; the as-built gap is still unpinned.
- The mesh map's ×0.56 offset is measured at 150 µm only. An absolute gain at
  135 µm therefore carries an unquantified mesh-geometry systematic, so
  **quote gain RATIOS between the two gaps, not absolute gains.**
- Nothing here is a licence to pick a gap that closes T14. Choosing the gap
  after seeing which value matches the data would launder the comparison into
  its own input — the same rule that governs ρ_s in the T14 freeze queue.

---

## Provenance

- Campaign label on every slice:
  `DIAGNOSIS-GRID / unconstrained-gap-scan / T14-amplitude-and-slope`
- Points: `response/avalanche/mx17_aval_points_gapscan.txt`
- Submit: `response/avalanche/mx17_aval_gapscan.sub` (lxplus condor) —
  **cluster 11988753, 48 jobs, submitted 2026-08-11 16:47 CEST**, after this
  file was written. Smoke slice (460 V, 135 µm, nev = 3) ran first on lxplus926
  and confirmed the gap reaches the output `config` block; ~10 s/event there,
  so slices should land in the 25–40 min band.
- Results land at `~/work/mx17_response/avalanche/results_gapscan/` on AFS;
  merge with `python3 -m response.avalanche.collect <dir> --out ...` and
  **check the merged keys carry `@gap135um`** — if both arms come back under
  one key the grouping fix did not travel with the code.
- Grouping: `collect.py` gained `gap_um` as a campaign axis on 2026-08-11
  (`arm_tags`). Without it the two gap arms share (gas, voltage, Penning)
  exactly and merge into one silently well-formed point — the same failure the
  Penning axis had on 2026-08-10, caught this time before submission.

## Results — 2026-08-11 19:05, on **all 48 slices**

Complete: 48/48, none held, none failed. One slice (135 µm/520 V seed 4, the
largest-avalanche point) was **preempted twice** by the worker pool and landed
on its third attempt at 19:02 — eviction, not a code or memory failure (733 MB
used against 9 000 requested). There is no checkpointing, so each eviction
restarted it from zero, ~48 min per attempt. Worth carrying forward: at this
gain, budget for evictions rather than assuming a slice that stops reporting
has died.

Adjudicated by `response/validation/gapscan_verdict.py`, whose header restates
this decision rule. Product: `gapscan_verdict.json`.

### Campaign falsifier: PASS

| | condor 150 µm | desktop 150 µm | dev |
|---|---|---|---|
| 460 V | 16 868.4 | 17 062.4 | −1.1 % |
| 490 V | 41 841.9 | 43 394.6 | −3.6 % |
| 520 V | 113 712.6 | 110 084.6 | +3.3 % |

Worst 3.6 % against an 11.4 % seed-statistics tolerance. The lxplus and desktop
legs agree; the campaign may be quoted.

### ⭐ The result that outranks both headline observables: ONE α(E) curve

In a uniform field `ln G = α(E)·d` **exactly**, with `E = V/d`, so the two arms
must lie on a single α(E) curve — they sample it over different E ranges, but
it is the same curve. This is a direct test of the parallel-plate Townsend
framing that the entire gain-and-slope argument rests on, across a 15 % change
in `d` and a 26 % change in E:

| gap | V | E [kV/cm] | α [cm⁻¹] |
|---|---|---|---|
| 150 | 460 | 30.67 | 648.9 ± 1.0 |
| 150 | 490 | 32.67 | 709.4 ± 1.3 |
| **135** | 460 | 34.07 | 753.9 ± 1.3 |
| 150 | 520 | 34.67 | 776.1 ± 1.8 |
| **135** | 490 | 36.30 | 827.2 ± 1.6 |
| **135** | 520 | 38.52 | 899.1 ± 2.4 |

Fitting α(E) to the **150 µm arm alone** and extrapolating onto the 135 µm arm:

| E | measured α | predicted | | gain |
|---|---|---|---|---|
| 34.07 | 753.9 | 755.6 | −1.2 σ | ×0.978 |
| 36.30 | 827.2 | 825.5 | +1.1 σ | ×1.023 |
| 38.52 | 899.1 | 895.4 | +1.5 σ | ×1.051 |

**χ²/point = 1.7. The arms lie on one curve, and gains are predicted to 2–5 %
across a gap change never previously tested.** The framing holds. The competing
pre-registration named a rising slope as the signature that it does not; the
slope did not significantly rise, and this test is the stronger version of the
same question.

### The two headline observables — a split decision, honestly

**Slope.** 150 µm **0.3147 ± 0.0048**, 135 µm **0.3273 ± 0.0058**. The matched
move is **+0.0126 ± 0.0075 = 1.7 σ** — consistent with flat, mildly favouring a
small rise. Quote the matched move, not the absolute: it is immune to anything
common to both arms.

- This note's arm (0.309 ± 0.010): measured 0.3273 — **1.4 σ to the nearest
  band edge, 3.2 σ to the band centre. Marginal; consistent on the edge
  reading, uncomfortable on the centre reading.**
- The competing arm (0.26–0.29): **excluded at 6.5 σ to the nearest band edge,
  9.1 σ to the band centre.** The predicted fall did not happen; the slope
  moved the other way.

Both bands are quoted under **both** conventions on purpose. "Excluded at N σ"
conventionally means nearest-edge, and the edge number is the conservative one
for an exclusion — but it is the *flattering* one for this note's own
prediction, so quoting only edges would be picking the convention that suits
whoever is holding the pen. Under either convention the competing arm is dead
and this note's arm survives; the honest summary is that this note's prediction
was correct in direction and about 1–3 σ optimistic in magnitude.

(σ figures come from `gapscan_verdict.py` on the unrounded fit, 0.327280 ±
0.005752. Recomputing them from the rounded 0.3273 ± 0.0058 printed above gives
6.4 / 9.0 — a rounding difference, not a disagreement. The script is the
authority.)

**Gain ratio 135/150.** ×1.560 ± 0.036 (460 V), **×1.691 ± 0.048** (490 V),
×1.643 ± 0.069 (520 V).

- This note's arm (×1.50–1.60): **HIGH**, ~2 σ over across the three voltages.
- The competing arm (×1.6–1.7): **IN.**

**So each pre-registration won one.** Both deviations from this note's numbers
point the same way — α is mildly *convex* above 35 kV/cm rather than exactly
linear, which raises both the 135 µm gain and its slope. That is one cause, not
two, and it is the opposite sign to the α″ < 0 the competing arm assumed, which
is why that arm lost the slope while winning the gain band.

**Ion transit** ×0.826 / ×0.823 / ×0.824 against a predicted ×0.81 — the
d²/µV scaling holds to 2 %.

### A correction to this note's own reasoning

The pre-registration argued the slope must not move because "α is linear in E
to better than 0.5 %". **That claim was weaker than it read.** It rested on a
three-point fit to the 150 µm arm, where a quadratic is exactly determined and
its curvature is therefore unconstrained noise — and indeed extrapolating the
*quadratic* 150-arm fit onto the 135 arm fails at χ²/point = 34 (gains wrong by
up to 26 %). The linear form is supported by the 135 µm arm, which is **new
information from this campaign**, not something known when the prediction was
written. The prediction was right; one of its stated reasons was not yet
earned.

### THE HEADLINE, invariant to all of the above

**The gap cannot close the slope discrepancy.** Data 0.4487 per 10 V; the
135 µm arm reaches 0.3273. Shrinking the gap by 10 % moves the sim **9 % of
the way** to the data and costs a ×1.69 gain change to do it. Both
pre-registrations agreed on this before the run and the data agrees with both:
**the gap is a gain lever and not a slope lever.**

Combined with Penning outcome C (2026-08-10: Penning delivers the gain or the
slope, never both), the slope now has no surviving single-parameter explanation
inside the avalanche model — and the α(E) framing those searches are conducted
in has just been validated across a gap change, so the failure is not the
framing's.

### What this still does not do

It does not measure det3's as-built gap. It measures what the model predicts if
the gap were 135 µm. **A ×1.69 gain lever is large enough to close most of the
amplitude deficit** (×1.46 after the transparency fix), and that is exactly why
the gap must not be chosen by which value closes T14 — the ρ_s rule in the
freeze queue applies unchanged. Pinning the real gap needs a measurement, not a
scan.

Also unchanged: the mesh map's ×0.56 offset is measured at 150 µm and assumed
at 135 µm. The near-field-fraction argument in §"What is actually being run"
signs that caveat — the real-detector gap ratio should sit slightly *below*
the ×1.69 measured here in the uniform field.
