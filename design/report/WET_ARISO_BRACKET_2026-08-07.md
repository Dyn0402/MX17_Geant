# Wet-Ar/iso bracket check — 2026-08-07

Follow-up to the P1 decision (`design/RESPONSE_SIM_PLAN.md` §0a, 2026-08-07): T14 now targets
det3 cosmic-bench data (Ar/iC₄H₁₀ 95/5, 1000 V drift / 490 V mesh, 30 mm nominal gap = 333 V/cm),
not run_71/SPS. The same "the gas probably isn't dry" question that motivated the wet-CF₄ bracket
for run_71 applies to det3's own gas: `~/x17/cosmic_bench/Analysis/fleet_gas_survey.csv` records
**h2o_pct = 1.46** for the exact `long_run_resist_490V_drift_1000V` sub-run behind T12's
pedestals/noise. This is that bracket, run the same way (`response/params/wet_ariso_drift.py`,
mirroring `wet_cf4_drift.py`).

**Source.** Run on lxplus (Garfield pinned master 927e5c21, self-built vs LCG_109),
`/afs/cern.ch/user/d/dneff/work/mx17_response/gas/wet_ariso_drift.py`, ~24 min wall clock, 5
mixtures in parallel. Raw output copied to `response/params/wet_ariso_drift.json` (schema
`wet_ariso_drift/1`), Ar/iC₄H₁₀ 95/5 base at 745.83 Torr (Saclay pressure, matching
`mesh_transparency.C` / the S3 campaign), 293.15 K, ncoll = 5. Machine-readable summary of
everything below: `design/report/wet_ariso_bracket_2026-08-07.json`.

**Gap assumption.** Per user instruction (2026-08-07): use the 30 mm nominal drift gap for this
run. The effective-gap dispute across nTof_x17 docs (19 / 23 / 27.9 / 30 mm — see plan §0a P1) is
a separate, still-open sub-item; sigma_T below is quoted at 30 mm and needs rescaling
(σ_T ∝ √z) once that resolves.

## Verdict: no single bracket — det3 has two *disagreeing* "measured" numbers, not one band

Unlike run_71 (one measured band, 13–15 µm/ns, that the wet curves cleanly bracketed), det3 has
three numbers in circulation for v_drift at this operating point, and they don't agree with each
other (plan §0a P1 / `response/digitizer/gas.py:14-27`):

| candidate | value | status |
|---|---|---|
| dry Magboltz reference | 39.1 µm/ns | this table's own dry column, and the existing `drift_velocity_Ar_iC4H10_95_5_Saclay.json` (39.37 µm/ns near 333 V/cm) |
| plan §1 "measured" | 36.6 µm/ns | provenance not re-verified here |
| micro-TPC (det3, 2026-07-01) | 28.1 ± 0.7 µm/ns | **different gap** — measured over det3's ~19 mm, not 30 mm; not directly comparable until the gap dispute resolves |

Values interpolated to exactly 333.0 V/cm (grid straddles it at 325.9 / 340.1 V/cm). σ_T is over
the assumed 30 mm drift gap.

| mixture | v [µm/ns] | dist to 36.6 | dist to 28.1 | d_T [√cm] | σ_T(30 mm) [µm] | η [/cm] |
|---|---|---|---|---|---|---|
| dry       | 39.15 | 2.55  | 11.05 | 0.04886 | 846 | 0.0 |
| H₂O 0.5 % | 40.62 | 4.02  | 12.52 | 0.03784 | 655 | 0.0 |
| H₂O 1.0 % | 33.91 | 2.69  | 5.81  | 0.02993 | 518 | 0.0 |
| **H₂O 1.46 %** | **25.25** | **11.35** | **2.85** | 0.02385 | 413 | 0.0 |
| H₂O 2.0 % | 17.61 | 18.99 | 10.49 | 0.01934 | 335 | 0.0 |

**The headline number — the actual measured water content for the exact det3 run T12 uses,
1.46 % — predicts v = 25.25 µm/ns.** That is far from the plan §1 value (11.3 µm/ns away, a 31 %
relative miss) and only moderately close to the micro-TPC value (2.8 µm/ns away, a 10 % relative
miss — outside the quoted ±0.7 µm/ns statistical band, but that band doesn't include the gap
confound or this table's own Magboltz uncertainty).

**Reading the bracket the other way** — what water fraction would each candidate imply — sharpens
this: 36.6 µm/ns falls between the 0.5 % and 1.0 % points, implying **~0.8 % water**, well below
the 1.46 % actually measured on that run. 28.1 µm/ns falls between the 1.0 % and 1.46 % points,
implying **~1.3 % water** — close to the measured 1.46 %, despite the different-gap caveat.

**So: the measured water content is markedly more consistent with the micro-TPC number than with
the plan §1 number.** This is not proof — the micro-TPC measurement's gap mismatch (~19 mm vs the
30 mm this table assumes) means it was taken at a different field (≈526 V/cm at 1000 V, not
333 V/cm), and this table has not been evaluated at that field. But at face value, it raises a
real question about which of the two det3 "measured" v_drift numbers is more trustworthy, and it
is a reason to run this same table at the micro-TPC's actual field before trusting either
comparison.

## An unexpected feature, checked and real: v_drift is non-monotonic in water content

Dry → 0.5 % H₂O *raises* v_drift (39.15 → 40.62 µm/ns, +3.8 %) before higher fractions drive it
down sharply (33.91 at 1.0 %, 25.25 at 1.46 %, 17.61 at 2.0 %). This was checked against scan
noise: per-point Magboltz statistical precision at these settings is sub-percent (e.g. dry's own
log reports ±0.03–0.36 % depending on field), so a 3.8 % shift between adjacent mixtures is real,
not noise. The wet-CF₄ bracket did not sample low enough in water (its lowest point was 1.0 %) to
say whether the same feature exists there. Not chased further here; noted because it means "more
water always means slower drift" is not a safe assumption to carry into the gap/field resolution
work.

## What this does NOT establish

- **Not a measurement of det3's water fraction** — 1.46 % is an external number (the gas survey),
  not fit here; this table only says what v_drift *would be* at 1.46 %, at the 30 mm/333 V/cm
  assumption.
- **Does not resolve which "measured" det3 number (36.6 or 28.1) is right** — see above; the honest
  read is that 1.46 % water favors 28.1 over 36.6, but the gap confound on the 28.1 number blocks
  a clean conclusion.
- **⚠️ Attachment is absent from the table, not confirmed negligible** — `eta_per_cm` is
  identically `0.0` in every one of the 60 field bins of all five mixtures, same limitation as
  `WET_CF4_BRACKET_2026-08-07.md` found for the CF₄ table. Most likely the Magboltz
  `ElectronAttachment` call not populating at ncoll = 5. Treat as open.

## Follow-on for §0a P1

This does not close the water-content sub-item — it sharpens it: before T14, either (a) settle the
effective gap first and re-run this table at the resulting field, since the two "measured" v_drift
numbers disagree exactly along the axis (gap/field) that is separately unresolved, or (b) treat
25.25 µm/ns (1.46 % water, 30 mm assumed) as the working T_drift number and flag that both
"measured" comparisons remain approximate. §0a has been updated to point here.
