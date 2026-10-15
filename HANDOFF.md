# Handoff

## Penning / ion-species review — updated 2026-10-15 (dylan-MS-7C84)

**Resume:** Ion-tail physics review done. Next: explain the T6 map's 31.0 vs 38.3 kV/cm mid-gap field, then swap in a hydrocarbon-ion mobility table.

**Goal:** Work out which ion actually drifts in Ar/iC4H10 95/5, and what r_P and ion mobility the response chain should use. Also find out how far the species alone explains the fast (~140 ns) ion tail in the waveform fits.

**Done:**
- Literature review: `reports/Simulating Penning and ion transport.md`, with notes in `research_notes/Simulating Penning and ion transport/` (commit 5bb6ac2).
- Animated explainer on the notes site, with two tabs (mechanisms / simulating it): https://dylan-neff.web.cern.ch/notes/ion-species-vs-penning.html. Source is `dylan-cern-site/pages/notes/ion-species-vs-penning.html`.
- Told the capsule-alignment session that MX17_Full_Geant (main and ear2-beam @19a6a71) still has the pre-e17c72c ArIso bug. Its gas density is ~7 % high and 95/5 is used as mass fractions, so energy deposited in the gas comes out ~7-11 % high.

**Findings (numbers in the report):**
- Ar+ converts to a hydrocarbon ion (t-C4H9+-like) in ~0.6 ns, a Langevin estimate. `IonMobility_Ar+_Ar` is the wrong species.
- Coimbra measured K0 = 2.26 (98/2) and 2.16 (Ar/CF4/iso 95/3/2) at <= 45 Td. Nothing has been measured at gap fields; working value K0 ~2.05 (1.8-2.25).
- The species switch speeds the transit up by about 1.5-1.85x (340 -> ~185-225 ns). The fits want 140 ns, which is 2.4x.
- The T6 map's mid-gap E_z is 31.0 kV/cm, versus V/d = 38.3 kV/cm (`response/meshcell/ion_template_check.py`). At 31.0 kV/cm, 140 ns needs K0 ~2.75, above every Langevin limit. At 38.3 kV/cm it needs ~2.23, exactly at the ceiling.
- r_P: Garfield hard-codes 0.40 for Ar/iso, taken from a single 10 % / 1 atm point. Use 0.40 dry and ~0.36-0.40 wet, via `EnablePenningTransfer(r, 0, "ar")`. Fit it jointly with feedback beta.

**Next steps:**
1. Find out why the T6 field map's mid-gap field is ~19 % below V/d (map voltage, gap, mesh-plane definition). This decides whether the species alone closes the ion tail.
2. Replace `IonMobility_Ar+_Ar.txt` with a custom flat table at K0 ~2.05, with 1.8 and 2.25 as the systematic, in `response/avalanche/mx17_aval_calib.py` and in the S3 chain.
3. Retract "Blanc gives K0 ~1.16" in `design/report/S3_ION_CLOSEOUT_2026-08-09.md`. It used Ar+'s charge-exchange-limited K. Done correctly, Blanc's law gives ~2.1.
4. Fix the backwards comment at `response/avalanche/mx17_aval_calib.py:174` ("Ar+ dominant ... the faster ion").
5. (Separate repo) Port the GasProperties/PAI fix from e17c72c into MX17_Full_Geant.

**Gotchas / decisions:**
- 153 Td (uniform V/d) and 123.8 Td (field map) are both in circulation. The gas density is fine (1 atm / 293 K); the difference is entirely the field.
- The high-field hydrocarbon-ion mobility is a 3x extrapolation in E/N. Label it as one.

## Primary ionization vs Heed (Stage A fix) — updated 2026-10-15 (dylan-MS-7C84)

**Resume:** Stage A primaries fixed (e17c72c, -11 %); next: regenerate T14 Stage A products, get Alexandra's confirmation.

**Goal:** check Alexandra's finding that Geant4 overpredicts drift-gas primaries versus Heed, and make Stage A as accurate as possible.

**Done:**
- Reproduced her Heed numbers. On a bare slab, option4 is +3.5 % on the mean and +10 % on the median (4 mm Ar/iso); PAI matches Heed to <1 % (e- and mu). About 30 % of her +35 % comes from material around the gas, not from the ionization model. `response/validation/heed_vs_g4/README.md`.
- e17c72c: `GasProperties` (ideal-gas density at 20 C, proper mass fractions, Garfield W/Fano), PAI for e±/mu±/p/pi± in `GasRegion`, fluo+Auger+ignore-cut. mm_sim new/old primaries: 0.893 (100 MeV e-), 0.908 (150 GeV mu-). CPU per event unchanged.
- `response/validation/heed_vs_g4/HANDOFF_ALEXANDRA.md`: claims C1-C6 for her repo's Claude to confirm.

**Next steps:**
1. Regenerate the Stage A products that feed T14 (EOS response_sim, yoga stageB). Until then, they are ~11 % high in primaries (scale by ~0.89 if needed).
2. Send Alexandra the handoff link and collect her confirmation, in particular whether her full-stack ratio falls to ~1.30 with PAI.
3. Optional data check of the absolute primary count, e.g. the most likely cluster charge on MIP tracks.

**Gotchas / decisions:**
- `G4EmParameters::AddPAIModel` is a silent no-op in Geant4 11.3, so PAI goes in via `G4EmConfigurator` in `PhysicsList::ConstructProcess`.
- GasRegion MUST have its own `G4ProductionCuts` (`PhysicsList::SetCuts`). Couples map to regions by cuts pointer; with shared cuts, sub-keV electrons stall in solids and runs are ~10x slower. Repro: `slab.cc` header.
- PAI in the kapton pillars (AmpGas daughters) costs only ~5 %. Early on I wrongly blamed them for the slowdown.
- Heed's sum(E)/sum(n) is not a W (escaping deltas). Neither Heed nor G4 includes the Penning/Jesse effect on W.
- Already on main before this work: the no-ROOT build fails (`RunAction.cc` is missing `<fstream>`), and MT + ROOT output crashes. Use `-t 1` jobs.

**Key files & commands:**
- `include/GasProperties.hh`, `src/PhysicsList.cc`, `src/DetectorConstruction.cc::Construct`, `src/SteppingAction.cc`
- `response/validation/heed_vs_g4/{slab.cc,heed_mip.cc,compare.py}`. Build: Geant4 11.3.2 at `~/Software/geant4/install`, Garfield at `~/Software/garfield` (ROOT 6.40); build mm_sim with `-DROOT_DIR=~/Software/root_6_40/cmake`.

## Irakli DREAM bench data — updated 2026-10-14 (dylan-MS-7C84)

**Resume:** FRESH_LOOK §19: reply to Irakli (protected or external C for NSW? protection BOM?), then get MX17 X/Y strip C.

**Goal:** use Irakli's NSW (AC) and TPOT (DC) DREAM measurements to pin the absolute ADC/fC and a possible X/Y charge-division term.

**Done:**
- Read both attachments. Summary and numbers in FRESH_LOOK §19. TPOT DC gain is 18.05 mV/fC at 0 pF (−10 % vs nominal).
  NSW AC gain fits charge division, C_c ≈ 510 pF, k ≈ 0.70.

**Next steps:**
1. Ask Irakli: NSW "AC" = protected FEU or external coupling C? Protection schematic/BOM?
2. Find MX17 X and Y strip capacitances (fab drawing or a C meter on a spare board), then test the §18.3 Y/X deficit.
3. Lower the sim's ADC/fC by ~10 % (TPOT 0 pF) as a candidate default; check what it does to bG.

**Gotchas / decisions:**
- The files are Saclay-internal and this repo is public, so the raw data stays out of git:
  `~/x17/irakli_dream_tests_2026-10-14/` (desktop only: pptx, xlsx, dumps, `nsw_charge_division_fit.py`).
- The 200 fC range was not measured in the TPOT slides. The NSW sheets have both ranges, about ×1.95 apart.

## Mesh comparison with P2 (Alexandra) — updated 2026-10-14 (dylan-MS-7C84)

**Resume:** rerun G7 / transparency_curve / one S3 HV point with `NULLCOLL=1` (null-collision steps), then fix the 0.955.

**Goal:** compare our mesh/field modelling with github.com/akallitss/P2_Geant_Garfield and her `MeshUnitCellIon_PressureStudy.mph`.

**Done:**
- Note `design/report/MESH_COMPARISON_P2_GEANT_GARFIELD_2026-10-14.md` (side-by-side, what to borrow both ways).
- Test: null-collision steps give G7 = 0.988 ± 0.002 on both woven and straight maps (default stepping 0.928 / 0.947). `gates_check.C` has a `NULLCOLL` env switch. Pushed (`d435c10`).
- Site animation (dylan-neff.web.cern.ch/notes/mx17-response-chain.html) now explains per-event vs ran-once: Geant4 nPrimary → per-electron sampling from Magboltz/Garfield calibrations; Garfield never runs per event (same as P2). Its Mesh stage is amber, citing 0.955 vs 0.988. Source: `design/report/deck/make_chain_animation.py` (`b510371`).

**Next steps:**
1. G7 at 10k electrons with status breakdown; rerun `transparency_curve.csv` with `NULLCOLL=1`.
2. One avalanche HV point of `mx17_aval_calib.py` with null-collision steps (gain/survival may move).
3. Try COMSOL export nodes without a coordinate file for `WovenMesh.java` -> `ComponentComsol`.
4. Ask Alexandra about mesh calendering; send her the surface-charge top BC.
5. Once ε is settled, update the Mesh card in `make_chain_animation.py` (STAGES `mesh`: st/check text, and `EPS_MESH`), regenerate, then `add-note.py ... --slug mx17-response-chain --force --deploy` (command in the script docstring).

**Gotchas:**
- One seed per run, preliminary. Only the `.mph` arrived in ~/Downloads (the second file did not).
- Her repo clone and unzipped XML were in the job tmp dir (gone after the session); re-clone if needed.

## T14 fresh look — updated 2026-10-14 (dylan-MS-7C84)

**Resume:** FRESH_LOOK §18.10: ρ_s ≈ 2.5 MΩ/sq may close Y width + Y/X; next, a chain point at ρ_s 5M (ask Dylan before submitting).

**Morning 2026-10-14 (desktop):**
- The yoga queue finished at 05:00: 34/34 points DONE, no failures, nothing on condor. The harness table
  (`scratchpad/2026-10-14/verdict_harness_table.py`, run on yoga) for verdict 1 gives inclined rise
  +11..+18 ns and FWHM +1..+18 ns. The harness amplitudes (`ratio` = sim/data) are X 0.95–1.13 and
  Y 1.12–1.20. The Y excess there is the §18.1 selection bias; the M3 ledger is the reference.
- **§18.10:** the Y wings pass COMSOL 3D at 2 MΩ/sq (Y0..Y4 peak ratios 0.92–0.95, the same as the
  centre), so the kernels are not the cause of §18.6. ρ_s is the lever. Interpolating from the point
  source, ρ_s ≈ 2.5 MΩ/sq raises Y/X by ~15 % and cuts the next-Y wing by ~22 %, which matches both
  §18.6 gaps, while X1/X0 is ~unchanged. 2.5 is at the top of the T2b band (1.42–2.56). It conflicts
  with §10.3's ~1.1 hint (weak). X0 drops ~7 %, so bG would need a refit.
- **Next 1 replaced:** run one chain point at ρ_s = 5 MΩ/sq (kernels in
  `/media/ucla/x17_resp/s1_w2_ny512/`, and on EOS) with the verdict-1 settings, pre-registered,
  using the `stageBC_w2_*.sub` pattern. Read the Y profile ±1..3 and Y/X from `ledger_m3.py`.
  Dylan has not yet said go.

**Overnight (02:00) state below.**

**Resume:** FRESH_LOOK §18 (all of tonight). Detectors are gone: archive + simulation only.

**Result: verdict 1 passes (§18.8).** Configuration: s = 1.75 (C₄Hₙ⁺), photon feedback bG 0.45, shaper
166 ns (`--shaper-tp-ns`, new), and the fixed noise spec `noise_det3_satscan_feu07_cmq.json`.
- Amplitude p50 data/sim: X 0.93 / 0.89 / 0.95 and Y 1.04 / 1.06 / 1.07 at 0/10/20°. X saturated fraction
  matches to 0.03.
- Inclined shapes: rise +8..+16 ns, FWHM +2..+7 ns (August default: +60..+70 / +90).

**Found tonight:**
1. **T14's amplitude rows are biased by the data-side selection** (§18.1). The quality cut and the
   reco-angle window are applied to the data only, and both correlate with amplitude; Y vertical
   is ×1.37 low. Use `scratchpad/2026-10-14/ledger_m3.py` instead: the data binned by the **true M3
   angle**, no reco and no quality cut, against all sim events. It reverses §16.2: there is no
   "Y overshoot", and bG fits X+Y jointly at 0.43–0.46.
2. **Sim noise bug, fixed** (§18.5). The common mode was 271 ADC instead of ~5, because the estimator
   was owned by railed events. It clipped 11 % of samples at 0 and created hot channels (the sim's
   2× strip multiplicity). The fix is in `response/dream/noise.py`. The A/B run (4451433) shows
   amplitude −1 % and shape +3 ns, so the earlier conclusions are unaffected.
3. **The remaining X/Y offset (~15 %) is the Y lateral sharing in the chain** (§18.6): the Y profile is too
   wide at ±1..3 strips, and X matches. Same family as the Y lobe and the ±1 delay. Open.
4. The bG × shaper scan (4451431/2): feedback is shape-neutral (≤ 2 ns), and 180 → 166 ns gives
   −5..−7 ns rise, −16 ns FWHM and −8 % sim amplitude (§18.7).

**Still running on yoga** (`~/queue_t14.sh`, log `~/queue_t14.log`): the full T14 harness on noisecm →
verdict1 → the bgtp points, for continuity only (~3–4 h). Nothing is left on condor. AFS work is
fine; home quota is untouched.

**Next:**
1. Y lateral sharing (§18.6): compare the chain's Y kernel wings with COMSOL 3D
   (`/eos/.../comsol/condor/i3d_*_o2.csv`, which has the Y channels), and check the routing traces.
2. Irakli's protection values: still wanted for the absolute ADC/fC scale and as a possible X/Y
   charge-division term.
3. Out-of-sample test: Ar/CF₄ (old det3, §15) with the same three changes. It needs an Ar/CF₄
   avalanche calibration, and the CF₃⁺ ion speed instead of C₄Hₙ⁺.
4. Decide whether `t14_compare.py` (nTof_x17) should get a true-angle data selection. It is a different
   repo, so it was not touched tonight.
5. The old `noise_det3.json` (May run, σ_cm 83) was built by the same flawed estimator: regenerate it if
   anything still uses it.

## Previous: T14 fresh look — 2026-10-13 ~21:00 (dylan-Yoga)

**Resume:** FRESH_LOOK §16–17 → get the protected-FEU values from Irakli, then decide new defaults + one verdict T14.
Detectors are gone: archive + simulation only. Slides: https://dylan-neff.web.cern.ch/notes/mx17-feedback-and-ions.html

**This session (laptop):**
- All condor work is finished (4450898 ionspeed, 4450899 fullphys; nothing queued or held; AFS 32 %).
  All 15 T14 points are done on yoga: `~/x17/response_sim/stageB_w2/t14_{ionspeed,fullphys}<s>_<ang>/`.
  Table: `ssh yoga python3 ~/t14_ionspeed_table.py` (a copy is in `scratchpad/2026-10-13c/`).
  The laptop can rsync from lxplus directly; the desktop relay in `orchestrate_t14.sh` is only
  needed when running from the desktop.
- **§16:** the hydrocarbon ion speed s = 1.46–1.75 plus photon feedback bG 0.45 match the X
  amplitude (data/sim 0.80–1.03). Feedback doesn't change the shape: the censored vertical rise
  moves +11 ns only because the saturated fraction rises 8 → 35 %; the uncensored rise moves
  ≤ 6 ns. s = 2.45 closes the inclined X shape on its own.
  **New problem:** Y overshoots with feedback (×0.69–0.76, 20–28 % saturated against 4 % in the
  data). X and Y disagree on the gain by ~15–20 %, and Y's rise is 13–37 ns too fast.
- **§17:** the cosmic-bench FEUs were the spark-**protected** variant (protection built for
  non-resistive Micromegas). The manuals give no circuit values; the `ProtFlt` register
  (0x200000 bits 19:4) was grounded (0xFFFF) in the target run. Toy test
  (`scratchpad/2026-10-13d/protection_rc_scan.py`):
  - a **CR high-pass is excluded**: the data's 2 % undershoot means τ ≳ 15 µs, which moves the
    rise by ~1 ns;
  - a low-pass of ≤ 20 ns is harmless;
  - **charge division** (C_c in series with C_strip) leaves the shape alone but stays open for
    the X/Y amplitude ratio and the absolute ADC/fC scale.

**Next:**
1. Get the protection schematic or BOM from Irakli Mandjavidze (Dylan is asking). With C_c, R_s
   and the diodes known, compute the charge division for the X and Y strip capacitances and
   check it against the §16.2 X/Y amplitude disagreement.
2. Refit bG to X and Y together (back-of-envelope ~0.35; Y alone ~0.25), or fix the X/Y split
   first if (1) explains it.
3. Choose defaults (s = 1.46–1.75, bG from step 2, shaper 158–166 ns) and run one verdict T14
   (the condor submit pattern is `scripts/condor/stageBC_w2_fullphys.sub`; pre-register the
   prediction).
4. Ar/CF₄ chain point as the cross-gas test (needs an Ar/CF₄ avalanche calibration).

**Earlier context (2026-10-12/13, desktop):**

**Goal:** find out why the end-to-end MX17 response simulation (Geant4 → Garfield → DREAM)
gives slower and wider pulses than the det3 cosmic-bench data. The rise-time / f_ion
contradiction of T14 is the problem being worked here.

**Done (2026-10-12, details in `design/report/FRESH_LOOK_2026-10-12.md` §6–7.6):**
- §6: the quality cut hid a third of the amplitude deficit. The rise time is unaffected.
- §7.1: an independent toy model agrees with the chain, so the chain code has no bug.
- §7.2: vertical tracks can't measure the ion fraction. The rise there mostly reflects where
  the peak lands on a flat top. At matched peak delay, data agrees with full ions.
- §7.3: at 10° and 20°, the data pins only the convolution shaper ⊗ ion. The pulse is about as
  narrow as a 180 ns shaper alone.
- §7.5: an ion-mobility rescale is an exact time rescale, so no Garfield rerun is needed. At
  180 ns, full ions need T_ion ≤ 110 ns (×3–4 Ar⁺), which is unphysical. August's Blanc's-law
  "species 4 % slower" check is invalid: it mixed two different ions
  (`response/meshcell/ion_template_check.py:193`).
- §7.6 + `design/DREAM_SETTINGS.md`: **the shaper is 180 ns**, read back from the target run's
  own archived cfg. The word order is proven (high word first). All 3574 archived cfgs were
  swept, and every data FEU in the cosmic, n_TOF and SPS runs is 180 ns. The 123 ns branch is
  closed.

- Published the state of play to the site: 12-slide note
  https://dylan-neff.web.cern.ch/notes/mx17-response-sim-status.html and an animated chain
  walk-through https://dylan-neff.web.cern.ch/notes/mx17-response-chain.html (toy pulse vs the
  real det3 average; the ion toggle shows the rise problem). It supersedes the August
  `t14-sim-vs-data-waveforms` note.

- §8 (`response/comsol/`): **COMSOL direct time-domain solves on lxplus.**
  - The chain's ion operator gets the pulse shape right, but the true ion
    signal is 11–20 % *larger* (axisymmetric, sheet draining the pad). So
    resistive-sheet screening is retired a third time, now by a full solve
    and with a sign.
  - The excess is T10's frozen-1−z/g failure (`T10_SLOWPATH_2026-08-08.md`),
    now confirmed independently.
  - True-weave mesh map: transparency 0.918 ± 0.006, against 0.960 ± 0.004
    for the straight-wire production map through the same macro. Not adopted.
  - **Superseded 2026-10-14** (`design/report/MESH_COMPARISON_P2_GEANT_GARFIELD_2026-10-14.md`):
    with `EnableNullCollisionSteps` (`NULLCOLL=1 gates_check`) both maps give G7 = 0.988 ± 0.002,
    against 0.928 / 0.947 with default stepping. The 4.2 % weave deficit and the 0.955 used for
    `DEFAULT_TRANSPARENCY` are stepping artefacts (one seed each; avalanche campaign not rerun).

- **Overnight 2026-10-13 (FRESH_LOOK §9–11):**
  - §9: **avalanche (photon) feedback**, G_eff = G/(1 − bG). Fitted to the *shape* of the det3
    27 June ladder alone, it predicts the data ×1.5–2.1 above Garfield at 490 V; T14 measured
    ×1.77–1.82. It also gives the slope ×1.5 and its positive curvature, which no α(E) can give
    at real gain (a Townsend fit needs G ≈ 200). The runaway lands at 511–520 V, where the bench
    discharges. Ion-impact feedback is excluded by the shapes. Scripts and figure in
    `scratchpad/2026-10-13/` and `design/report/figs/`.
  - §10.1: **COMSOL 3D on the real strips and comb** (`response/comsol/Induce3D.java`, condor)
    against the chain's W2 kernels: timing exact (neighbour delays within 0–7 ns), amplitude
    within ~5 % and converging. The ion operator holds on real strips too (f_eff 0.89).
  - §10.2: **the shaper from the target run's noise autocorrelation** is ~166 ns (s = 0.92,
    0.88–0.94). At that shaper, full ions need T_ion ≈ 140 ns (×2.4 Ar⁺).
  - §10.3: the ±1 delay is mainly a cloud-size observable. The data's transport invariance
    implies sheet-transported neighbours, which at +60 ns points to ρ_s ~1 MΩ/sq. Only a hint.

**In progress / where it stopped:** nothing half-done. The conclusion is that with the true
180 ns shaper the inclined data needs the ion signal largely absent (f_eff ≈ 0.2). That is
August's reading, and it still has **no mechanism**.

**Next steps** (FRESH_LOOK §11 has the reasoning):
0. **Bench mesh ladder at a second isobutane fraction** (mixer logged since 2026-07-07). Photon
   feedback must weaken with more quencher; this decides §9.
1. Fire the FEU pulser once through the test capacitor (bench). It measures the shaper's actual
   response (§10.2 says ~166 ns from noise) and the absolute ADC/fC scale (§9's ×1.8 match).
   Also measure det3's ρ_s bus to bus (§10.3).
2. **T10 step 5** (plan §7 step 5): build the LUT from the slow-path templates
   (`solver/zextend.py`). This is a known error, now confirmed by COMSOL. It makes the rise a
   few ns faster and the central-pad ion signal larger. Small, but the chain is wrong without it.
3. ~~Validate the chain's kernels G with a full solver~~ **done (§10.1)**: COMSOL 3D on the real
   strips and comb agrees on timing (0–7 ns) and to ~5 % in amplitude. The Y negative lobe
   (item 5) was not in that comparison's window; it is the next thing to read from the 3D runs
   (`/eos/.../response_sim/comsol/condor/i3d_*_o2.csv`).
4. Rerun G7, `transparency_curve.csv` and one S3 HV point with null-collision steps (see above). Woven mesh: rerun S3 gain and `funnel_ion_endpoints.py` in `meshfield_woven_production.txt`,
   then decide whether to adopt it. Its −4.2 % transparency is absorbed into the feedback fit's
   scale if §9 is adopted.
5. Still unexplained: the early "foot" before data pulses, and the −0.2 negative lobe on the
   sim's Y channels (§7.2).
6. Amplitude and gain slope: one candidate mechanism now (§9). Item 0 tests it.

**Gotchas / decisions:**
- **Garfield on the desktop is built against ROOT 6.40**, but `setup_garfield.sh` still selects
  6.30, and libGarfield then fails to load (TGeo ABI mismatch). Source
  `~/Software/root_6_40/bin/thisroot.sh` after it, and pass
  `gSystem->AddIncludePath("-I$GARFIELD_INSTALL/include")` to ACLiC. The MX17 ComponentGrid
  patch is in that build.
- **COMSOL model code cannot read or write files.** Output goes to stdout as `CSV,` lines,
  which `run_lxplus.sh` collects. See `response/comsol/README.md` for this and the other API
  traps.
- Don't use vertical tracks for the ion fraction. Use 10°/20°, matched by angle, with no
  quality cut.
- Archived cfgs live at `<area>/Run|runs/<run>/<subrun>/raw_daq_data/*.cfg` on EOS
  `/eos/experiment/ntof/data/x17/`. A recursive `find` over x17 is too slow; use globs.
- The CEA FEU source (`~/Feu/Firmware5/...`) is on the DAQ PC only. `daq_lxplus` was not
  answering on 2026-10-12.

**Key files & commands:**
- `response/comsol/` — COMSOL models (2D, axisymmetric, 3D strips, woven mesh), condor
  submission, `compare_induce2d.py` (ion operator) and `compare_kernels3d.py` (chain W2 kernels
  vs COMSOL), and a README with the numbers.
- `scratchpad/2026-10-13/` (README inside) — feedback fit, noise-shaper fit, ion-T scan at
  166 ns, ±1 delay vs ρ_s. The data inputs came from the laptop (`yoga`, read-only).
  Outputs are on EOS: `/eos/experiment/ntof/data/x17/response_sim/comsol/`.
- `design/report/FRESH_LOOK_2026-10-12.md` — the running write-up (§6–7.6).
- `design/DREAM_SETTINGS.md`, `design/dream_cfg/` — DREAM settings record, per-run CSV,
  target cfg.
- `scratchpad/2026-10-12/` (README inside) — toy model, waveform analyses,
  `wf_ion_T_scan.py`, `dream_settings_by_run.py`.
- Waveform npz files: laptop `~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/`. Run the
  `wf_*.py` scripts with both npz files in the cwd, using the nTof_x17 `.venv` python (needs
  pandas).
- `python ~/PycharmProjects/nTof_x17_DAQ/dream_scripts/decode_dream_cfg.py [--sweep] <cfg|dir>`
- `design/report/deck/make_sim_status_deck.py`, `make_chain_animation.py`: site generators. Run with
  `--wf-dir <dir with both npz> --out <html>`, then publish with
  `python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py <html> --slug <slug> --force --deploy`.
  Update the LEDGER in the deck script when FRESH_LOOK numbers change.
