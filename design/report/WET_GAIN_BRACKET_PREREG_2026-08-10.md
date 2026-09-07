# Wet amp-range gain bracket — PRE-REGISTRATION

**Written 2026-08-10 BEFORE the jobs were submitted.** Roadmap
`GAS_AND_DRIFT_CAGE_ROADMAP_2026-08-08.md` §2 step 1. Predictions and
falsifiers below are stated in advance so the result cannot be read backwards
out of whatever comes out.

## ⚠️ Labelling, per the standing discipline

**Every water fraction here is FITTED-TO-DATA, not measured.** There is no
hygrometer reading behind any of it (Dylan, 2026-08-09); the ~1 % figure was
itself inferred from v_drift by Magboltz matching. Nothing in this bracket may
be quoted as evidence that the bench gas contained water.

That is fine for what this actually asks, which is a **sensitivity** question,
not an inference: *if* there were water at the 0.5–1 % level, would it move the
avalanche GAIN enough to matter? The answer is a derivative, and a derivative
does not require the operating point to be established.

## The roadmap's premise is wrong, and it makes this cheaper

The roadmap states the wet suites are drift-range only and that this needs
"~2 new Magboltz jobs before any Garfield". **No Magboltz job is needed.** All
three tables already exist at amplification range and on an identical grid:

| table (`Saclay_160m`, 745.83 Torr) | E grid | present |
|---|---|---|
| `Ar_iC4H10_95_5` | 5 000 – 60 002 V/cm, 20 log points | ✅ |
| `Ar_iC4H10_H2O_94p5_5_0p5` | identical | ✅ (laptop only — shipped to lxplus for this run) |
| `Ar_iC4H10_H2O_94_5_1` | identical | ✅ |

The confusion is a units trap worth recording: the `.gas` files store **E/p**,
so the header's "6.704 … 80.45" is V/(cm·Torr) and multiplies by 745.83 Torr to
5–60 kV/cm. Read as V/cm it looks like a drift-range table, which is how the
roadmap's row came about. The genuinely drift-range tables carry an explicit
`_drift_` in the filename.

## The confound this run exists to avoid

`mm_config.py` gives the dry mixture `penning: auto` and both wet mixtures
`penning: manual, rP = 0.4` — necessarily, because Garfield has no ternary
Ar/iC₄H₁₀/H₂O Penning table. **Comparing dry(auto) against wet(rP = 0.4) would
confound the water effect with a change of Penning model**, and Penning is
precisely the knob the T7 slope hunt is currently chasing.

Penning is applied at avalanche time (`mm_sim_core.py:69-73`,
`EnablePenningTransfer` after `LoadGasFile`), **not** baked into the `.gas`
file, so the tables are Penning-agnostic and the setting is ours to control.
The bracket therefore runs **all three mixtures at manual rP = 0.40 on argon**,
with dry-at-auto carried as a fourth arm purely to tie back to the production
calibration.

## Configuration

490 V across a 0.015 cm gap (the T7 pooled bench point), `Saclay_160m`,
8 batches × 200 events = 1600 avalanches per arm.

| arm | mixture | Penning |
|---|---|---|
| A | Ar/iC₄H₁₀ 95/5 | manual rP 0.40 (ar) |
| B | Ar/iC₄H₁₀/H₂O 94.5/5/0.5 | manual rP 0.40 (ar) |
| C | Ar/iC₄H₁₀/H₂O 94/5/1 | manual rP 0.40 (ar) |
| D | Ar/iC₄H₁₀ 95/5 | auto (tie to production) |

A↔B↔C is the bracket. A↔D measures the Penning-model offset separately, so it
can never be mistaken for a water effect.

This is a **uniform-field** gain scan, which is the right instrument: the
question is a RATIO between mixtures at fixed field and fixed gap, and the mesh
field map is gas-independent by construction (the FEM solve has no gas in it),
so it divides out. Absolute gain is not what is being measured here.

## Predictions, stated in advance

1. **Direction — gain falls monotonically with water, A > B > C.** H₂O's
   ionisation potential is 12.62 eV, *above* both Ar metastables (11.55 /
   11.72 eV), so water opens **no new Penning channel** — it only adds a
   polyatomic with large low-energy vibrational cross sections, which cools
   electrons below the ionisation threshold. It also displaces argon, not
   isobutane (95/5 → 94/5/1). Every term points the same way.
   **Falsifier: gain rises at either water point.**

2. **Magnitude — 1 % H₂O moves mean gain by 10–40 %.**
   **Falsifier: outside that band.** Note this prediction *contradicts* the
   roadmap's hoped-for outcome ("if water moves gain by less than a few %, drop
   the contaminant axis"): I expect the axis has to stay.

3. **Attachment — survival P(g>0) stays ≥ 0.95, unchanged within errors.**
   Water is a poor electron attacher compared with O₂; the attachment story
   belongs to the oxygen hypothesis, not this one.
   **Falsifier: survival falls below 0.90 at 1 % H₂O.**

4. **The 0.5 % point sits between A and C, and roughly half way in log-gain.**
   **Falsifier: non-monotonic, or B outside [C, A].**

## The decision this feeds, and what it does NOT decide

Per roadmap §2 step 2: if water moves gain by less than a few percent, the
contaminant axis drops out of the gain campaign entirely. If it moves gain by
tens of percent, the axis stays — and the T7 slope hunt's assumption that gas
composition is fixed while Penning varies becomes an approximation that needs
its own error bar.

**What this cannot decide**, and must not be read as deciding: whether there is
any water in the bench gas. A gain sensitivity is not a measurement of the
operating point. Separately, gain *scale* is not the live problem — the live
problem is the gain *slope* (data d lnA/dV = 0.449/10 V vs sim 0.296, ≈12 σ),
and a one-voltage bracket says nothing about slope. If this bracket comes back
large, the follow-up is a two-voltage version, not a gain correction.

## Products

Fragments → lxplus `…/garfield_sim/jobs_wetbracket/`, merged with
`mm_condor_collect.py`. Cluster ID recorded in
`design/report/OVERNIGHT_2026-08-10.md` when submitted.
