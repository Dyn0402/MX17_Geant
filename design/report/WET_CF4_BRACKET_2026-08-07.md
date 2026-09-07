# Wet-CF₄ bracket check — 2026-08-07

Audit fix 5 (finding A5). The plan §0a cert row "wet-CF₄ v_drift 74.32 vs 74.7 µm/ns" compares
the table's **dry** column against a published **dry** Magboltz number — dry-vs-dry, and
therefore no evidence at all about the wet curves. The check that matters is the one
`response/params/wet_cf4_drift.py:122-124` names itself:

> if the wet curves do not bracket 13–15 µm/ns the water hypothesis is not reproduced and
> **nothing downstream should be run**.

That result existed nowhere in this repo (the JSON is EOS-only). It does now.

**Source.** `/eos/experiment/ntof/data/x17/response_sim/gas/wet_cf4_drift.json`
(schema `wet_cf4_drift/1`), Ar/CF₄/iC₄H₁₀ 88/10/2 base at 720.8 Torr, 293.15 K, ncoll = 5.
Fetched with `xrdcp root://eospublic.cern.ch/…`. Machine-readable copy of everything below:
`design/report/wet_cf4_bracket_2026-08-07.json`.

**Measured target.** 13–15 µm/ns at 233 V/cm
(`nTof_x17 sps_beam_test_26/analysis/RAW_RUN71_REANALYSIS_2026-08-04.md`).

## Verdict: **PASS** — the wet curves bracket the measurement, P1 is not blocked

Values interpolated to exactly 233 V/cm (the table's own grid straddles it at 231.4 V/cm).
σ_T is over the full 30 mm drift gap — the quantity Stage B consumes and the halo prediction
needs.

| mixture | v [µm/ns] | in 13–15 band | d_T [√cm] | σ_T(30 mm) [µm] | σ_L(30 mm) [µm] | η [/cm] |
|---|---|---|---|---|---|---|
| dry      | 74.66 | no  | 0.02120 | 367.1 | — | 0.0 |
| H₂O 1.0 %| 19.99 | no  | 0.01757 | 304.3 | — | 0.0 |
| **H₂O 1.5 %**| **13.71** | **yes** | 0.01613 | 279.3 | — | 0.0 |
| H₂O 2.0 %| 10.34 | no  | 0.01566 | 271.2 | — | 0.0 |

Wet span 10.34–19.99 µm/ns, which straddles 13–15; H₂O ≈ 1.5 % lands inside it. This is
consistent with the independent nTof_x17 conclusion that the drift deficit is water, and with
the det3 gas history (dried from >3 % to ~1 % over a week).

Note in passing that the interpolated **dry** value, 74.66 µm/ns, reproduces the published
74.7 to 0.05 % — so the old cert row was a correct check of the *wrong thing*, not a wrong
check.

## What this does NOT establish

- **It is not a measurement of the water fraction.** Three wet points bracketing a band is a
  consistency check; 1.5 % is where this table crosses the measurement, not a fitted value.
  The band itself is 13–15 µm/ns, and v is steep in H₂O here (~6 µm/ns per 0.5 %), so the
  implied water is roughly 1.3–1.7 % on this table alone.
- **σ_T moves the wrong way for a "more diffusion" story.** Adding water *reduces* transverse
  diffusion over the gap (367 → 271 µm), so a wet gas does not explain a broader halo by
  diffusion; it explains slower drift.
- **⚠️ Attachment is absent from the table, not confirmed negligible.** `eta_per_cm` is
  identically `0.0` in every one of the 60 field bins of all four mixtures, wet included, over
  the full 40–500 V/cm range. Consequences:
  - Fix 6's survival term e^(−ηz) is a **no-op** against this table. The plumbing is still
    correct to add (and is), but it cannot be validated from this file.
  - The plan's premise that CF₄+H₂O attachment re-weights the depth mixture is **unsupported
    by the only table that could support it**.
  - Most likely cause is the Magboltz `ElectronAttachment` call not populating at ncoll = 5,
    rather than a physics claim that attachment vanishes. **Treat as open**: either re-run the
    suite with attachment verified nonzero somewhere, or state in the plan that P1 proceeds
    with η = 0 as a deliberate, recorded assumption.

## Plan text this replaces

§0a cert table, wet-CF₄ row: it should say the **dry** column reproduces published dry
Magboltz to 0.05 %, and that the wet bracket check passes as recorded here — not that "wet-CF₄
v_drift" is validated by a dry comparison.
