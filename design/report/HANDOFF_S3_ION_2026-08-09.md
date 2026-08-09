# HANDOFF — S3 ion-signal investigation (rise-time floor follow-up)

2026-08-09 (late). Owner: unassigned — `response/avalanche/` is the
mx17-geant-6b session's domain; this handoff is written so any session can
take it. Coordinator: the T14 session (nTof_x17 side).

## Why this exists — the measurement that triggered it

T14 (first waveform-level sim-vs-data comparison, det3 saturday-scan
`long_run_resist_490V_drift_1000V`, W2 rho2M kernel, published at
`dylan-neff.web.cern.ch/notes/t14-sim-vs-data-waveforms.html`) found the
simulation has a **10–90 % rise-time floor of ~250 ns that the data does not
have**: 38–42 % of data pulses rise in under 240 ns (down to ~155 ns) against
3–5 % of sim pulses. Upper quantiles agree; the sim is missing the entire
fast side. Ruled out by measurement: track angle (matched window, fast risers
are ordinary full-gap tracks), noise model (DIAGNOSIS_noisefix barely moved
it), packet aggregation (packet=False in production), shaper impulse response
(10–90 % = 116 ns at code 2 — data's 155 ns is ABOVE it; the earlier
"model cannot cross" claim was a 5→100 % vs 10–90 % definition error,
withdrawn).

**The discriminator (`DIAGNOSIS_noions`, Stage B `--no-ions` at rho2M)
settled it decisively**: with the ion time-redistribution removed — same
total charge, all prompt — the sim rise distribution lands on the data at
every quantile:

| 10–90 % rise (ns) | p5 | p25 | p50 | p75 | p95 | frac < 240 ns |
|---|---|---|---|---|---|---|
| X sim with ions | 257 | 295 | 337 | 417 | 608 | 3 % |
| X sim no ions | 142 | 197 | 270 | 415 | 634 | **41 %** |
| X data | 150 | 202 | 275 | 409 | 602 | **40 %** |
| Y sim with ions | 240 | 273 | 307 | 380 | 603 | 5 % |
| Y sim no ions | 129 | 168 | 232 | 371 | 617 | **52 %** |
| Y data | 155 | 192 | 243 | 336 | 540 | **49 %** |

Shape RMS of the normalized average waveform collapses from 13.0/14.5 % to
8.6/5.4 % (X/Y). Also: the maximally-prompt no-ions sim still peaks at
×0.63 of data, so the **absolute amplitude deficit is upstream (gain/kernel)
and decoupled from time structure — do not mix the two problems** (Dylan,
2026-08-10: keep amplitude separate; the contaminant grid covers the gas
axis of it).

## The claim under investigation

Data experiences the same ion physics, so the conclusion is NOT "remove
ions". It is: **the modeled slow component — f_ion × ion template × what
survives the β = 0.75 shaper — grossly overstates the real rising edge's
slow content.** Three dials, in decreasing order of prior suspicion:

1. **β / PZC (electronics)** — being scanned NOW (Stage B `--pzc-residual`
   {0.0, 0.15, 0.30, 0.50}, DIAGNOSIS_beta*; 0.75 point = existing
   noisefix). β was pre-declared a fit parameter (`response/dream/shaper.py`:
   undershoot selects β, remaining shape observables are cross-checks; the
   0.6–0.9 prior band is run_71/SPS provenance, NON-constraining for det3).
   Joint target: rise quantiles AND undershoot (sim currently overshoots
   2–3× data: X −9.9 % vs −3.4 %, Y −25.2 % vs −12.0 %, medians). If one β
   fixes both views' rise AND undershoot, S3's job shrinks to validation.
   ⚠️ The peaking-time register (state1 code 2 → 180 ns) is an ASSUMPTION
   from `CosmicTb_MX17.cfg`, not archived with the run — provenance exhibit
   five if wrong; recovery needs CEA DAQ access (Dylan).

2. **f_ion ≈ 0.90 (charge split)** — from the S3 v2 currents
   (`f_ion cross-check 0.9006 vs 0.9006`). Verify the Shockley–Ramo
   weighting treatment on the READOUT electrode through the mesh: mesh
   screening means the readout-induced split can differ strongly from the
   in-gap split. Use the T6 production field map (with the local Garfield
   ComponentGrid 3D region-flag patch — see `setup_garfield` notes), not the
   uniform-field branch.

3. **The i_ion template's time profile** — S3 v2 measured template: 2000
   samples × 0.2 ns; Stage B logs quote 50/90/99 % of ion charge by
   172/310/357 ns, ions born at ⟨z⟩ = 14.9 µm. Cross-checks:
   - Which ion mobility and which ION SPECIES? (Ar+ vs iC4H10+ vs cluster
     ions differ ×~2 in mobility; transit time scales inversely.)
   - Consistency with the analytic model: `--ion-model analytic` vs
     `measured` is a one-flag Stage B A/B (~19 min) nobody has run on W2.
   - Template sanity vs the calib's own `alpha_z_hist` (avalanche depth
     distribution) + amp-gap field: is 172 ns-to-half even kinematically
     right for 50 µm at the S3 field?
   - History: the schema-1 calib shipped `i_elec`/`i_ion` as 2000 ZEROS and
     silently NaN'd the LUT (caught 2026-08-07; digitize.py now hard-fails).
     The v2 arrays are populated — but populated ≠ validated. This subsystem
     has bitten before; treat its provenance with the same suspicion that
     found the noise/FEU/calib/register issues (five provenance failures in
     one day, all the same class: a per-run/per-detector property carried as
     a universal constant).

## Degeneracy warning (the reason S3 work matters even if β fits)

β and the template speed are partially degenerate on the rise observable: a
too-slow template can be compensated by a too-low β. The undershoot breaks
some of it (β drives undershoot; the template mostly doesn't), and the
Y-vs-X asymmetry breaks more (β is view-common; the measured X/Y undershoot
asymmetry −3.4/−12.0 % in DATA is real and must come from elsewhere —
plausibly kY/kernel or view-dependent sharing). S3 should pin the template
and f_ion independently so β is not silently absorbing their errors:
whatever (f_ion, template) S3 defends, β is then fit on top, and the
combination must reproduce rise + undershoot + the average-waveform tail in
BOTH views.

## Concrete work list

1. Re-derive f_ion on the readout electrode with the T6 field map +
   Garfield patch; quote in-gap vs through-mesh splits.
2. Validate the v2 i_ion template: species/mobility used, kinematic
   consistency with alpha_z_hist, and re-measure if the mobility choice is
   wrong.
3. Run the `--ion-model analytic` vs `measured` Stage B A/B on W2 rho2M
   (one condor job, DIAGNOSIS label) — if analytic vs measured differ
   materially on rise quantiles, the template is load-bearing and must be
   right; if not, f_ion/β carry everything.
4. Feed the defended (f_ion, template) back; the T14 session re-runs the
   β fit on top and closes the loop against data.

## Practicalities

- Frozen things: `w2_rho2M/default` decoded + `t14_compare/` = the verdict
  record; never regenerate. All new points: DIAGNOSIS_* naming, own dirs.
- Stage B runs at CERN (laptop cannot read /eos/experiment; LUT builds
  forbidden on the laptop). Scripts: AFS `mx17_s1/` (`run_stageBC_w2.sh`
  pattern), ~19 min/point. Corrected noise spec
  (`noise_det3_satscan_feu07.json`) for all new points.
- Readout on the laptop: `mx17_sim_wft/t13_reco.py` (dual-use seeding,
  σ = 5) then `t14_compare.py` — see the repo's `t14_compare.py` docstring
  and the memory file `mx17-t14-first-comparison` for the sharp edges
  (data-mirror symlink trick, seedstats requirement, estimator caveats).
