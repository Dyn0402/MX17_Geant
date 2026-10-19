# Handoff

## det3 X asymmetry, second round: it lives at the Y pads — updated 2026-10-19 late (dylan-MS-7C84)

**Resume:** deck slides 19–27 (https://dylan-neff.web.cern.ch/notes/det3-x-asymmetry.html), `asym_deck/FINDINGS.md`
"Second round". Mechanism at the Y pad still open; Dylan to pick bench tests (deck "Open").

**Done:**
- Copy only for charge over a Y pad (X/Y matched parity split): det3 June-15 0.22, det4/det7 0.08, det6 0.02, others 0;
  every Y pad alike (no pillar, no stripe-beat pattern). Electronics/cable excluded by this.
- Date + drift-field dependence collapse onto cluster width (Y ±1 share): gas state of the night + field. Gas circulation
  ruled out as the cause (X only, pad-level, flat top-bottom, non-monotonic). Varies along the X-chip axis only.
- Beam drift scans (condor 4502894) of det2/det6/det7: inconclusive (det6 0, det7 too noisy, det2 beam-angle).
- Sims: copper/dielectric/displacement give < 0.02 on the data estimator; floating stripes have the wrong fingerprint.

**Next steps:**
1. Bench: same HV, two gas states; collimated Fe-55 on one Y pad vs X pad; det4/det7 drift scan; ask Stephan about the
   June recipe; microscope det3 Y pads; same PCB/resistive print for det3 May vs June-15?
2. Model candidates not yet built: a partly floating Y pad, a conductive/charged residue at the Y-pad edge on the bus side.
3. Analysis correction: uncopy on Y-pad events only (parity from the Y view).

**Gotchas:** `fwd_parity.py` Y-line placement had a half-pitch bug (fixed; it made ±0.15 per class from symmetric
geometry). Beam needs `raw_asym.py --cm` (flash). Desktop-only derived data: `~/x17/nb_asym_raw/{extra,beam}/`; beam
npz also on EOS `/eos/user/d/dneff/asym_beam/`.

**Key files:** `asym_deck/asym_more.py`, `asym_deck/beam/{beam_drift.py,beamjob.sh,beam.sub,mkjobs.py}`,
`bus_wpot.py`, `bus_float.py`, `fwd_parity.py`, `bus_par/`, `fwd_out/`.

## det3 X asymmetry, raw-ADC deck — updated 2026-10-19 (dylan-MS-7C84)

**Resume:** det3 X copy is real (raw ADC, follows the board); next: the discriminating tests (bus side, cable swap, HV-off pulser).

**Goal:** convince Dylan (sceptical) whether the det3 X left/right neighbour asymmetry is real, from the rawest data up.

**Done:**
- Deck https://dylan-neff.web.cern.ch/notes/det3-x-asymmetry.html (20 slides), built by
  `scratchpad/2026-10-19/asym_deck/make_asym_deck.py`; findings and numbers in `asym_deck/FINDINGS.md`.
- New independent extractor `raw_asym.py` (raw ADC − per-channel pedestal, channel numbers, no map/CMS/clustering) over 16
  detector views + every cosmic-bench sub-run (`run_raw*.sh`, yoga). det3 June-15 X: p+1 bigger in 84 % of events; other 13
  views 50 %. Copy a ≈ 0.19 balances ±1 and ±2; in phase; linear; noise not copied; same FEU 3 clean on det3-May/det6.
- Centre-peaked across the board (chips 2–5 ~0.25, edges ~0.1; 2D map too); rises reversibly at low drift field; flat in mesh HV ≤ 500 V.
- det4/det7 X (Stephan June-22 batch) 0.05–0.06; det1/det2/det6 ≈ 0.

- **By date** (slide 17, added 10-19 evening): centre-chip copy 6-18 0.26 → 6-19/20 0.03 → 6-22 0.16 → 6-24 (600 V) 0.26 →
  6-27/28 0.10 → beam 0.06. Session-stable, jumps between sessions. 6-19 = only run with det3 drift HV sparking (4–9 µA).
  Details in `asym_deck/FINDINGS.md` "By date"; extra npz in desktop `~/x17/nb_asym_raw/extra/`.

**Next steps:**
0. Decode 6-19 HV scan + long_run (raw_daq_data only, EOS june_tests) to see if the copy returns within 6-19/20; ask Dylan
   whether det3 cables/HV were moved on 6-19 (config swapped the boards' HV channels vs 6-17).
1. Pin gerber → channel orientation: does the one-sided bus side predict p+1 receiving the copy? (wrong side ⇒ bus excluded)
2. Ask Dylan/Stephan: cable/adapter swap test; same PCB for det3 May vs June-15?; what differed in the June-15 bulk.
3. HV-off pulser into one det3 X pad (static coupling vs field-dependent bulk effect).
4. Add the hot-channel veto (`asym_lib.hot`) to `../nb_copy_survey.py`/copy_survey.csv; beam rows still need a raw re-check
   (raw_asym's isolation cut rejects every beam event: flash).

**Gotchas / decisions:**
- The old survey's 6-23 "turn-on" (0.035→0.17) was hot channel 492 (97 % of short-run events). Always veto hot channels.
- Junk (flat, no-pulse) events: the pulse-shape cut in `asym_lib.good` is needed.
- `pkill -f <pattern>` over ssh kills the ssh shell itself when the pattern is in the command line (lost run_raw3 once).
- 2D map: X and Y files must be matched by file number + eventId (file lists differ between FEUs).

**Key files & commands:**
- `asym_deck/raw_asym.py`, `run_raw.sh`, `run_raw2.sh`, `run_raw3.sh` — yoga extraction → `/media/dylan/data/claude_tmp/tmp/nb/raw/`
  (copied to desktop `~/x17/nb_asym_raw/`, incl. `survey/`).
- `~/PycharmProjects/nTof_x17/.venv/bin/python make_asym_deck.py` then
  `python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py det3-x-asymmetry.html --slug det3-x-asymmetry --force --deploy`.

## Neighbour asymmetry + timing vs strip distance — updated 2026-10-19 (dylan-MS-7C84)

**Resume:** det3 X copy + displaced-feedback sim (v5a); next: ±1/±2 timing and Y shares on top of v5a (NB_ASYM §E).

Read `scratchpad/2026-10-19/NB_ASYM.md`, its top "Update" sections A–D first (they supersede "det3-only" and
the bus-trace explanation of det3's size). `NB_FAR.md` has the ±1..±5 delay tables.
**Findings:**
- The det3 X −1/+1 asymmetry is a one-sided in-phase copy (channel n+1 gets a × channel n): signal only, not noise, and
  removable with `uncopy.py`.
- It is ≈ 0 on the May det3 build. On the June-15 build it is 0.11–0.17 on 6-22, decaying to 0.08 on 6-27/28 and
  0.06 at n_TOF run_79. It grows at low drift field (reversibly). Other June builds have 0.02–0.08.
- The one-sided bus traces (gerbers) give 0.010 in a 3D solve with the Geant4 stack estimate (264 µm FR4/layer), which
  matches the quiet boards but can't make det3's size or its time/field dependence.
- The far ±3..±5 prompt halo on every board is displaced photon feedback (Cauchy, H ≈ 150 µm = amp gap, f = bG);
  Stage B now has `--feedback-spread-um`. **v5a** (wet, σ0 94, bG 0.45, H 150 µm; condor 4481160) matches the copy-free
  X view with no retuning (halo, ±1/±2 shares, widths, amplitude 0.98/1.02). Open: ±1 peak 30–45 ns late, ±2 early, Y
  shares low.
- Copy-free X: ±1 delay +75 ns = v4c, but share 0.30 vs 0.39; no single σ0 fits shares and delays together.
**Next:** NB_ASYM §E. Gotcha: yoga's lxplus session had no AFS token today, so `kinit` there first, or pull through
the desktop and run `SKIP_PULL=1 ./eval_v5.sh`. Dylan: does anyone know what differed in the June-15 bulk (process, drying,
coverlay), and was det3's X cabling moved between slots? Tagged dumps: yoga `/media/dylan/data/claude_tmp/tmp/nb/{tagged,survey}/`.
v5 points: AFS `mx17_s1/stageBC/v5{a,b,c}_<ang>` and yoga `/media/dylan/data/claude_tmp/v3wet/v5*_<ang>`. Nothing is still running.
`scratchpad/2026-10-17/eval_v2.log` has a stray one-line local change from an old driver; it is deliberately left uncommitted.

## Reco truth study + wet gas / prompt anode spread — updated 2026-10-18 afternoon (dylan-Yoga)

**Resume:** read `design/report/RECO_TRUTH_2026-10-18.md` (live: https://dylan-neff.web.cern.ch/notes/mx17-sim-reco-truth.html,
rebuild + deploy with `scratchpad/2026-10-18/reco_truth/publish.sh`). Dylan asked: use the sim to reinforce the cosmic/SPS
track reconstruction (items 1–3), and keep improving the sim where the path is clear.

**Done (all in the note, scripts in `scratchpad/2026-10-18/reco_truth/`):**
- **Wiring bug:** Stage B wrote pad column k to channel k; det3 connectors are "inverted", so wft's strip map tore sim tracks that
  cross a 64-channel boundary (> 25 % at 20°). `response/dream/rewire.py` fixes existing files; `run.py --wiring connector-inverted`
  (opt-in) for new ones. Only position-based 20° comparisons moved (sim 20° extent now = data). The 26.5 mm choice stands.
- `run_stageBC_w2.sh` now keeps the truth sidecar. v3 / v3+coupling regenerated with truth at 12 angles (condor 4480259;
  byte-identical to the originals). New Stage A angles 5/15/30° and 10°×10° (yoga, seeds 101606–101612).
- **Item 3 ladder:** hit-time compression is mostly geometry (ideal no-readout ladder = data); the sim's readout compresses ~2× the
  data's; the coupling candidate makes it worse (a vote against adopting it).
- **Item 1 closure:** wft calibrated on the sim (truth for M3) is +4.5 % at 10° and nonlinear (flips 15 % at 20°); data vs M3 are
  linear to 1 %. Flips come from the sim-calibrated bundle; the nonlinearity from the sim's waveforms. No correction carried to data.
- **Item 2:** the data's |wft − M3| growth with depth is mostly M3's pointing; wft vs truth is 0.10–0.22 mm on vertical tracks.
- **Sim:** wet gas (1 % H2O, `design/gas/..._H2O_1.0.json` via `MX17_GAS_TABLE`) + `--sigma0-um 420` closes the ±1 neighbour delay
  (X +60, Y +60 vs data +45/+60; v3 +135) and the vertical X cluster (n, span, rms). Y ±1/central overshoots; inclined Y dimmer.

**Finished later the same afternoon:**
- **v4c** (condor 4480263; wet 1 % + σ0 420 µm + bG 0.48, 12 angles, truth, det3 wiring): passes its pre-registration
  (joint 0.987, ±1 delay +75/+60, vertical X 0.64 mm with the data's n and span, 10° extent 1.46). Y/X unchanged (1.12).
  Files on yoga `/media/dylan/data/claude_tmp/v3wet/v4c_<ang>`, AFS `mx17_s1/stageBC/v4c_<ang>`, EOS `response_sim/stageB_w2/v4c_<ang>`.
- **No-ions diagnostic** (condor 4480264): wft's nonlinearity on the sim persists (swing 0.067 X). Not the ion tail. Unexplained.

**Next steps:**
1. Dylan's calls: (a) coupling out? (ladder + delay say yes), (b) adopt the v4 candidate (wet gas + σ0 420 µm + bG 0.48; the bG check passed),
   (c) the physical origin of the prompt spread (lateral displacement of photon-feedback secondaries is the candidate; would explain
   cloud_basics' σ0 growing with mesh voltage). A `--feedback-spread-um` option is the natural next model step.
   The note's last section ("What's next") is the ranked list.
2. Y view: ±1/central too high while the vertical Y width is a bit narrow; inclined-Y dimness worse with σ0. Still open (§21.9).
3. The sim-calibrated wft bundle is the weak link for closure (c1 at floor, σ_s runs away); consider seeding/pinning σ_s.

**Gotchas:** `pkill -f <pattern>` from the Bash tool kills the tool's own shell when the pattern appears in the command line; kill by PID.
Yoga's data disk is at ~98 %; `/media/dylan/data/claude_tmp/v3` (original identity-wired v3, 6.7 GB) and `v3t` (raw truth regen, 2.9 GB)
duplicate the rewired copies and EOS, and can go if space is needed.

## Overnight 10-17/18: v3 baseline, coupling, systematics — updated 2026-10-18 04:30 (dylan-Yoga)

**Resume:** FRESH_LOOK §21.7–21.10 and slides 9–13 of https://dylan-neff.web.cern.ch/notes/mx17-response-sim-state.html.
Dylan to decide: adopt the inter-strip coupling (k 0.05, τ 15 ns, then bG 0.47) as a fitted term, or wait for Irakli's values?

**Done (all pre-registered, every falsifier passed):**
- **v3 baseline** `scripts/condor/recipes/det3_v3.args`: ρ5M, e17, nc1, tp 166, 26.5 mm, ion reconstructed K0 2.05 (Dylan said
  "Ok great" to the switch; `run.py` defaults stay bare so that old runs reproduce), **bG 0.45** (condor 4480068, joint 1.002).
- **Inter-strip coupling** `--strip-coupling K,TAU` (`response/dream/coupling.py` + test; Kirchhoff, conserves charge). Condor 4480073:
  (0.05, 15 ns) closes the neighbour leading edge (0.17 → 0.23, data 0.22), cuts the nb-shape rms 39 % in X and 13 % in Y, and moves the
  ±1 delay from +135 to +105 ns (data +45/+60). It also takes 4 ns off the X 20° rise. With it on, bG → 0.47 (4480074).
- **Systematics** (4480074): K0 1.80 / 2.25 → amplitude +3.2 / −1.7 %, rise +11 / −6 ns; 27.9 mm → −1.5 %, extent +4–5 %.
- **Y/X, sharpened (§21.9):** X saturates as often in the sim as in the data at every angle, but inclined-Y tracks rail 0.33 in the
  data vs 0.16–0.22 in the sim. The sim's Y peak strip is too dim on inclined Y. Strip C / FEU gain are excluded (sums agree).
  The kernels already model the resistive stripes along y. The open lever is the as-built stripe line/gap (never measured).
- Slide note republished (dylan-cern-site 613de81). Evaluation drivers: `scratchpad/2026-10-18/eval_v3*.sh` (yoga).

**Next steps:**
1. Dylan: coupling in or out (above). If in, add `--strip-coupling 0.05,15 --feedback-bg 0.47` to a det3_v3c recipe.
2. Fit ρ_s × (k, τ) jointly on the Y wings + neighbour shape. That needs the Y-wing profile, which needs T14 outputs (not made since 10-16).
3. Inclined-Y peak: measure the printed stripe line/gap on a board (microscope); if it differs, rebuild the kernel at CERN.
4. FEU 7/8 vs 3/4 central pulse shapes: `central_shape_feu_swap.py` is NOT usable yet. The 6-22 run's waveforms plateau at ~0.57
   of the peak after FeuReader's pedestal (mixed 32/37-sample events or ZS?). Sort out that format first, then select by M3 slope.
5. ρ_s × coupling (condor 4480075, §21.11): ρ2M shortens the delay 30 ns but breaks the Y profile (vertical Y rms 1.05 vs 0.85).
   ρ5M + coupling stays best. The remaining prompt neighbour share is NOT the glue (Aug d_k scan: < 1 %) and NOT diffusion
   (vertical width already slightly wide). A non-conserving prompt share (`--strip-coupling K,TAU,KFREE`, diagnostic) was emulated
   offline and makes the shape worse (`emul_kfree.log`). The data neighbour carries the same share earlier instead of later. Next:
   rethink the lateral model (slower sheet delivery and a prompt lateral term with the total held), ideally with Irakli's schematic.
6. Still open: Ar/CF4 out-of-sample; Irakli's questions (protection R/C, inter-channel C).

**Gotchas:** the laptop OOMs building a kernel LUT (`run.py` locally gets killed); run Stage B on condor. `wf_shape.py` (10-12)
loads a file at import; copy `cross()` instead of importing it.

## Verdict 2, drift-gap test, Y/X diagnosis — updated 2026-10-17 late (dylan-Yoga)

**Resume:** §21.6 written. Verdict 2 passes and adds linearly; the data want a ~26–26.5 mm det3 column. Next: Dylan's calls
(ion, ρ_s 5M, `--drift-gap-mm 26.5`), then a bG refit (≈ 0.45–0.46) at that column.

**2026-10-18 00:14 (yoga):** Dylan chose ρ_s 5M and `--drift-gap-mm 26.5` for det3, with 27.9 as the systematic. The bG refit
**condor 4480068** (15 jobs, `scratchpad/2026-10-18/stageBC_w3_bgrefit.sub`, pre-registered) runs K0 2.05 at 26.5 mm with
bG {0.43, 0.45, 0.47} → `v3k265b<bt>_<ang>`. The yoga driver `scratchpad/2026-10-18/eval_v3.sh` waits, pulls and writes
`*_v3.log`; progress is in `eval_v3.log`. To rerun: `setsid nohup bash eval_v3.sh > eval_v3.log 2>&1 < /dev/null &`.
The ion default in the code is still s1.75; the switch to K0 2.05 waits for Dylan's answer.

**Evening update (yoga, desktop was in Windows):** all 40 jobs of 4477507/4477508 finished OK. `eval_v2_yoga.sh` (eval_v2.sh
without the wait loop) pulled them to yoga `/media/dylan/data/claude_tmp/v2` and wrote the `*_v2.log` files. Yoga reaches AFS
fine after `kinit` + `ssh lxplus` (GSSAPI delegation gives the token), so the old "yoga has no AFS token" gotcha is just a
missing ticket. Results → FRESH_LOOK §21.6:
- Verdict 2: every ledger band hit. Ion shift is +1.5 % / +5 ns at both ρ_s values. ρ5M leaves a ~10 % Y/X offset.
- Gap: the 26.0 mm extent matches the data at 10° (1.43/1.46 vs 1.46/1.47 mm). Rise is unchanged. Amplitude moved only
  half the predicted amount (deep clusters are the diffused ones). The inclined pkX/pkY does NOT move, so Y/X is not the column.
- Neighbour ±1 delay unchanged (+135..150 vs +45/+60 ns); X central peak now 195 ns, = data.
- Not done: the Y-wing profile (needs T14 outputs) and the slide-note update (step 4 below).

**Goal:** run the outstanding sim tests: the combined "verdict 2" chain pass, the Y/X cause, the ±1 neighbour delay.

**Done (FRESH_LOOK §21, commits 9eca6f0, a1b638f):**
- **Evaluation now runs on the desktop** (yoga was offline). `m3_slopes_longrun_desktop.py` rebuilt
  `~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet`. `ledger_m3.py` and `shape_all_desktop.py` reproduce yoga's
  verdict-1 ledger and shapes exactly. Sim points are symlinked as `~/x17/response_sim/stageB_w2/w2_<rho>_<kind>/<ang>/decoded_root`
  → `/media/dylan/data/claude_tmp/...` (the desktop `/` is full).
- **Y/X is not the electronics.** With one observable (`xy_balance_m3.py`), the sim's summed X/Y split f matches the data at
  every angle (Δf ≤ 0.03). Only the peak strip is off. det3 gives f 0.484 on FEU 7/8 and 0.481–0.494 on FEU 3/4
  (`xy_balance_survey.log`), so FEU gain is out. Y noise is 1.7 % above X (strip-C division would make it ~25 % quieter).
  ρ_s 5M closes Y/X at 0° to 3 %. The remainder is on inclined tracks.
- **Drift column:** inclined tracks spread over 14–23 % more strips in the sim (`xy_extent_m3.log`), which implies a ~25–27 mm
  column vs the CAD 30. GAP_STUDY (nTof_x17) measured det3 27.9 ± 1 mm. New Stage B option `--drift-gap-mm` (default 30 = unchanged).
- **Neighbour delay on verdict 1** (`xtalk_fit_{v1,r5}_m3.log`): the central strip is fixed (210 vs 195 ns), but the ±1 peaks
  +150/+120 ns after it vs the data's +45/+60. ρ_s 5M matches the ±1 amplitude but makes the delay worse. The fit needs a prompt
  in-phase coupling k ≈ 0.13–0.17 plus a derivative term: inter-strip capacitance, not modelled.

**In progress / where it stopped:**
- Condor running at wrap-up (20 of 40 left at 13:16), outputs to AFS `mx17_s1/stageBC/<tag>/` and EOS `response_sim/stageB_w2/<tag>/`:
  - **4477507** `stageBC_w2_verdict2.sub`: e17 clusters, bG 0.47, nc1 calib, tp 166, × ρ_s {2M, 5M} × ion {s1.75, K0 2.05}
    → `v2r<2|5><s|k>_<ang>`.
  - **4477508** `stageBC_w2_gap.sub`: the v2r5 arm with `--drift-gap-mm` {27.9, 26.0} → `v2r5<s|k>g<279|260>_<ang>`.
  - Predictions and falsifiers are in each .sub header (copies in `scratchpad/2026-10-17/`). **Read those before the numbers.**
- The desktop driver `scratchpad/2026-10-17/eval_v2.sh` was started in the background. It waits for both clusters, pulls,
  then writes `ledger_v2.log`, `shape_v2.log`, `xy_balance_v2.log`, `xy_extent_v2.log`, `xtalk_fit_<k>_m3.log` next to itself,
  with progress in `eval_v2.log`. If it died with the session, rerun it on the desktop:
  `cd scratchpad/2026-10-17 && setsid nohup bash eval_v2.sh > eval_v2.log 2>&1 < /dev/null &` (re-pulls are cheap, ~115 MB/point).

**Next steps:**
1. Check `eval_v2.log` says EVAL DONE (and `eval_v2_failed.txt` is empty). Compare against the pre-registrations and write FRESH_LOOK §21.6:
   verdict 2 (additivity of the 10-16 one-variable results) and the gap test (does the 26–28 mm column close the inclined extent,
   the inclined peak fractions and the inclined Y/X; how far bG must come down).
2. Dylan decides: ion default (K0 2.05 vs s1.75, now also in light of the gap result), ρ_s 5M, the drift column for det3.
3. Neighbour delay: model inter-strip coupling (buried X/Y tracks under the pad layer, or the protected-FEU input).
   Irakli's questions (`scratchpad/2026-10-16/irakli_questions.md`) still matter here; Dylan sends them.
4. Update the slide note (`design/report/deck/make_sim_state_deck.py`): the Y/X slide and the ideas list are superseded by §21.
5. Still open: Ar/CF4 out-of-sample test; ADC/fC −10 %.

**Gotchas / decisions:**
- Condor multi-word args in `queue ... from (...)` get split on spaces. Use separate `ion = ...` + `queue` blocks (see the .sub files).
  One bad cluster (4477506) was removed before it ran.
- The data-side f uses no M3 match in the survey (every cosmic). The sim comparison must be angle-matched (`xy_balance_m3.py`),
  because the f of amplitude sums depends strongly on angle (0.35–0.63).
- `lxstore status` hung on lxplus today (timed out after 30 s). AFS was at 32 % via `fs lq`.

**Key files & commands:**
- `scratchpad/2026-10-17/` — all scripts and logs of this session; `README`-style docstrings at the top of each.
- `response/digitizer/run.py --drift-gap-mm` — the new option.
- `~/PycharmProjects/nTof_x17/.venv/bin/python <script>` — the interpreter for all desktop evaluation scripts.

## Sim to-do sweep: ion K0, rho_s 5M, Stage A e17, mesh G7 — updated 2026-10-17 (dylan-Yoga)

**Resume:** results in FRESH_LOOK §20 + slide note mx17-response-sim-state. Next: Dylan's ion-default call; Y/X via strip C (§19).

**Goal:** work through the outstanding sim list: physical ion mobility, the rho_s lever for Y sharing,
Stage A products with the e17c72c primary fix, mesh transparency with null-collision steps.

**Done:**
- Alexandra **confirmed** the Heed vs Geant4 claims (`HANDOFF_ALEXANDRA.md` marked). Mesh calendering: **unknown**,
  so the uncalendered 37 µm weave stays, and calendering is a +1–2 % field systematic.
- Retracted S3 "K0 ≈ 1.16" (`design/report/S3_ION_CLOSEOUT_2026-08-09.md`); fixed the Ar+ comment in
  `response/avalanche/mx17_aval_calib.py`.
- Synced the tracked `response/` + `shared/` files to the AFS src (`/afs/cern.ch/work/d/dneff/mx17_s1/src`). A smoke run of
  `--ion-model reconstructed --ion-mobility K0=2.05` on lxplus works (90 % of ion charge by 191 ns).
- **Stage A regenerated** at e17c72c+: `clusters_ang_e17/mx17_muons_3k_ang_<ang>_e17_t0.root` on AFS, the same files on EOS
  `response_sim/clusters/`. Recipe: `mm_sim -g ArIso -p muon -e 4000 -n 3000 --beam-spread 31.2 --theta-x/-y <a>
  -s 10160<i> -t 1` (the Aug 9 recipe, read back from the old files; the original seeds weren't recorded). nPrimDrift medians are 0.898–0.915 of the old files.
  The non-angled `mx17_muons_3k_spread_t0.root` was NOT regenerated (th00_e17 is its replacement).
- **Mesh G7 at 10k e-** (`scratchpad/2026-10-16/g7_*.log`): straight map 0.951 → **0.990**, woven 0.928 → **0.990**
  with `NULLCOLL=1`. Every loss is status −5 (on a wire). `gates_check.C` now prints that breakdown. The E_drift curve
  with null-collision steps is in `response/meshcell/transparency_curve_nullcoll.csv`. At 50–200 V/cm the gate's 0.995 ceiling
  "FAILs" on physical ~0.998 values; that is the catastrophe band, not a problem.
- MX17_Full_Geant gas port: already done by another session (81f619b, merged).
- Old `noise_det3.json`: nothing uses it, apart from the fallback defaults (`run_stageBC_w2.sh:79`, `run.py --noise`). Not regenerated.

- **Avalanche nc0 vs nc1** (490 V, 4×140 ev, desktop): survival 0.950 → **0.995**, gain|surv ratio **0.998 ± 0.040**
  (`scratchpad/2026-10-16/aval_nullcoll_490_result.md`; raw `~/x17/aval_nullcoll_2026-10-16/`). So the production calib's 0.9559
  survival is a stepping artefact and the chain's mesh transparency is ~4 % low.
- **lxplus Garfield was unpatched**: the old `garfield-927e5c21.tar.gz` ignores the map's inactive regions (0 % vs 1.64 % on the
  desktop, so electrons cross the wires). Rebuilt as `~/work/garfield_install/lcg109-927e5c21p` + `garfield-927e5c21p.tar.gz`
  (`scratchpad/2026-10-16/build_garfield_927e5c21p_lxplus.sh`). `response/avalanche/check_grid_active.py` now guards each job.

**2026-10-17 session (yoga): results read, nc1 campaign rescued, state-of-play deck**
- **Slide note** https://dylan-neff.web.cern.ch/notes/mx17-response-sim-state.html (13 slides: scorecard, amplitude journey,
  bG split, rise vs K0, Y profile, rho_s vs Y/X, mesh, fixes, open rows, ideas). Generator
  `design/report/deck/make_sim_state_deck.py` reads `scratchpad/2026-10-16/results/` (copies of the yoga ledgers/profiles
  and the nc0/nc1 pooled calibs). Rebuild: `python3 make_sim_state_deck.py --out note.html` then
  `add-note.py note.html --slug mx17-response-sim-state --force --deploy`. Its "Ideas" slide is the ranked to-do list.
- **Relay results → FRESH_LOOK §20** (commit 4d03edf):
  - ionK: passes. K0 2.05 is +5–6 ns rise and +1–2 % p50 vs verdict 1. s = 1.75 corresponds to K0 ≈ 2.2 (2.25 reproduces it to 1 ns).
    **Default NOT switched.** Dylan's call: physical K0 2.05 vs 5 ns of inclined-rise agreement.
  - rho5M: **falsifier triggered.** The Y wings match the data at ~4–5.6 MΩ/sq (not 2.5), but sim Y/X moves only 1–7 %
    (X brightens too). Y/X is not rho_s, so the next suspect is strip C / protection charge division (§19). The thx20 point is missing from r5.
  - e17: passes (+5..+10 % p50 at bG 0.45). The amplitude falls faster with bG than predicted, so the joint optimum is **bG ≈ 0.47**
    (pre-registered 0.49–0.53: missed). X alone wants 0.44, Y alone 0.50.
- **lxplus/condor Garfield was broken by header order** (nTof_x17 faefdac): the LCG_109 view ships Garfield 78fe1bd3 headers
  on ROOT_INCLUDE_PATH, so cling laid out our 927e5c21 classes wrong (ComponentGrid 688 vs 720 B, AvalancheMicroscopic,
  MediumMagboltz…). That causes a heap overrun: all 56 nc1 jobs died (malloc assert) or hung (25 s CPU in 4.5 h) in the grid guard.
  The fix puts `$GARFIELD_INSTALL/include` first in `setup_garfield.sh`. It is deployed to `~/work/mx17_response/{avalanche,avalanche_nc1}` on AFS.
  The AFS clone `~/work/git/nTof_x17` has someone's uncommitted edits, so it was NOT pulled and still has the old script.
  `mx17_response/gas/setup_garfield.sh` is the old version too.
  - **Audit (condor 4474927, §20.4):** the August script and tarball with fixed headers give 45.0k vs August's 42.8k/46.4k at 490 V
    uniform field, which is not significant. Old calibs look OK. Field-map lxplus runs were not audited separately.
  - `mx17_aval_calib.py` is **not reproducible per seed** (two same-seed runs differ), so only statistical A/Bs work.
- **nc1 campaign done** (condor 4474926, 56/56): pooled at EOS `response_sim/avalanche/aval_calib_meshfield_nc1_pooled.json`
  (`scratchpad/2026-10-16/pool_nc1.py`). ε 0.956 → **0.993**, gain|surv ×**0.970**, **ε×gain ×1.007**, so **no bG refit needed**.
  The chain still uses the nc0 file. The Mesh card text on the site is updated (dylan-cern-site 415417e). FRESH_LOOK §20.5.
- The desktop is booted into Windows. Its relay finished before that (ledgers on yoga `~/ledger_1016_*.log`).

**Next steps:**
1. Dylan: switch the ion default to `--ion-model reconstructed --ion-mobility K0=2.05`? (§20.1)
2. Y/X: get the MX17 X/Y strip capacitances and send Irakli the questions (`scratchpad/2026-10-16/irakli_questions.md`, Dylan by email).
3. At the next full chain pass: e17 clusters, bG ≈ 0.47, the nc1 calib file, and rho_s ~4–5 MΩ/sq for the Y wings (check it against the T2b band 1.42–2.56).
4. Still open: the Ar/CF4 out-of-sample test; the true-angle selection in nTof_x17 `t14_compare.py`; ADC/fC −10 %.

<details><summary>2026-10-16 state (superseded)</summary>

**Running (left going at wrap-up, 2026-10-16 ~18:00):**
- **nc1 avalanche campaign on condor** (the 56 slices of `mx17_aval_points.txt`, `--null-collision`, patched Garfield):
  - 4474488: 6 jobs running since 17:02 (470V_s6, 490V_s5, 500V_s2, 500V_s8, 520V_s1, 520V_s8).
  - 4474686: the other 50, resubmitted ~17:55 (`mx17_aval_points_nc1_resub1.txt`) after the first try held: a `| grep` under
    `set -e` made the guard exit silently (fixed), and 520V_s4 went over 6 GB (12.3 GB; request now 10 GB).
  - Job dir `lxplus:~/work/mx17_response/avalanche_nc1`; output `/eos/experiment/ntof/data/x17/response_sim/avalanche/raw_meshfield_nc1_20261016/`
    (expect 56 `aval_meshfield_*.json`, ~55 MB each). The wrapper skips slices already on EOS, so a resubmit is safe.
  - Prediction (pre-registered in `mx17_aval_nc1.sub`): survival 0.990–0.997, gain|surv within ±3 % of 24 094.
- Condor **all finished, exit 0** (outputs on AFS `mx17_s1/stageBC/<name>/`, and on EOS `response_sim/stageB_w2/<name>/`):
  - 4472241 `stageBC_w2_ionK.sub`: verdict 1 with `--ion-model reconstructed --ion-mobility K0=1.80/2.05/2.25` → ionK180/205/225_<ang>
  - 4472242 `stageBC_w2_rho5M.sub`: verdict 1 at rho_s 5M → r5verdict1_<ang>
  - 4472245 `stageBC_w2_e17.sub`: verdict 1 on the e17 clusters, bG 0.45/0.50/0.55 → e17b45/50/55_<ang>
  - Pre-registered predictions + falsifiers are in each .sub header. **Read those before looking at the numbers.**
- **Desktop relay** `~/x17/relay_1016.sh` (copy at `scratchpad/2026-10-16/relay_1016.sh`), log `~/x17/relay_1016.log`.
  It pulls each point from AFS, pushes it to yoga, runs T13/T14 there, and finally runs on yoga: `~/ledger_1016_{v1,r5,ionK,e17}.log`,
  `~/shape_1016.log` and `~/profile_1016_{rho2M_verdict1,rho5M_r5verdict1}_{th00,thy10}.log`. The relay skips finished points,
  so if it dies, rerun `setsid nohup bash ~/x17/relay_1016.sh > ~/x17/relay_1016.log 2>&1 < /dev/null &`.
- **Desktop avalanche slices** `~/x17/aval_nullcoll_2026-10-16/` (`run.sh`, a copy at `scratchpad/2026-10-16/aval_nullcoll_490.sh`):
  490 V meshfield (straight map), seeds 49001–4, 140 ev, default (`_nc0`) vs `--null-collision` (`_nc1`), local
  Garfield (ROOT 6.40). Expect ~30–60 min. Compare `survival` and gain nc1 vs nc0.

**Next steps:**
1. When the relay finishes, compare against verdict 1 (`ledger_1016_v1.log`):
   - ionK: does K0 2.05 ≈ verdict 1 (rise within 10 ns, p50 within 6 %)? If yes, replace `--ion-speed-scale 1.75` with
     `--ion-model reconstructed --ion-mobility K0=2.05` in the default and quote 1.8/2.25 as the systematic.
   - r5verdict1: Y/X and the Y wings (profile logs), then log-interpolate the rho_s that matches each. Do they agree near 2.5?
   - e17b: the bG that best fits X+Y jointly on the corrected primaries (prediction 0.49–0.53).
   - Write the results up as FRESH_LOOK §20.
2. nc1 campaign: once 56 JSONs are on EOS, pool them like `aval_calib_meshfield_pooled.json` (one 490 V point; its note says
   how; `collect.py` merges per-voltage, which is wrong here), check the prediction, refit bG, update the Mesh card
   (`design/report/deck/make_chain_animation.py`, `EPS_MESH`). Held jobs: `condor_q -hold`.
3. Send Irakli the questions in `scratchpad/2026-10-16/irakli_questions.md` (Dylan, by email).
4. Still open: the Ar/CF4 out-of-sample test; the true-angle selection in nTof_x17 `t14_compare.py`; ADC/fC −10 %.

</details>

**Gotchas / decisions:**
- **/media/ucla is failing (2026-10-16)**; Dylan copied it to EOS for recovery. Don't write there. The desktop root and
  `/media/dylan/data` are both <4 GB free; yoga has ~43 GB on `/` and ~68 GB on `/media/dylan/data`.
- Field-map avalanche runs are ~18 s/event per core, so the 56-slice campaign is ~2.7 h on 12 desktop cores: run it on condor.
- Any earlier lxplus/condor run with a field map used the unpatched Garfield; not audited.
- **yoga → lxplus has no AFS token** ("Permission denied" on work/). That is why the yoga-only orchestrators reported every point
  "MISSING". The relay runs on the desktop for that reason. The desktop has only ~4 GB free, so it relays one point at a time.
- `~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh` is broken on the desktop (ROOT 6.30, no Garfield headers,
  no PyROOT). Use `source ~/Software/root_6_40/bin/thisroot.sh; source ~/Software/garfield/install/share/Garfield/setupGarfield.sh`.
- lxplus (LCG_105) has no Garfield Data, so pass `--ion-mobility K0=...`; a bare table name fails there.
- The ionK/rho5M points use the OLD clusters on purpose (one variable at a time against verdict 1).
- Yoga copies of `orchestrate_v4.sh` / `run_ionspeed_t14.sh` take `RHO=` (default rho2M).

## Ion birth heights, reconstructed ion model, field map explained — updated 2026-10-16 (dylan-MS-7C84)

**Resume:** run a T14 chain point with `--ion-model reconstructed --ion-mobility K0=2.05` (vs verdict 1's `--ion-speed-scale 1.75`).

**Goal:** check whether the ion birth distribution is modelled, add a way to use any ion mobility, and explain the T6 map's "low" field.

**Done:**
- Birth heights were already simulated: `alpha_z_hist` in every S3/meshfield calib. They are roughly exponential, with a decay length of 14.7 µm; 90 % of ions are born below 34 µm. The measured template carries the full distribution. It moves ion timing by ≤10 %, so it is not a lever (`scratchpad/2026-10-16/ion_birth_sensitivity.py`).
- New `--ion-model reconstructed --ion-mobility <Garfield table | K0=x>` (`ions.reconstructed_longitudinal`). It keeps the calib's electron spike, ion charge, birth heights and avalanche offset, and rebuilds the ion path in the map's E_z and psi (packed into `response/meshcell/ion_path_profile.json` by `build_ion_path_profile.py`). With the Ar+ table it matches the measured template to 2.8–3.8 %. Test: `python -m response.digitizer.test_ion_reconstructed`.
- **Field map explained, no bug.** "38.3 kV/cm" was 490 V / **128 µm**, a gap MX17 doesn't have; MX17's is 150 µm. The map's 31.0 kV/cm is 5 % under 490/150 because the woven mesh's equivalent plane sits ~8 µm inside the weave (effective gap 158 µm). `response/meshcell/mesh_plane_check.py`. Corrected in the review report, `mx17_inputs.md` and FRESH_LOOK §13.
- In the real geometry, the full ×2.4 over Ar+ would need K0 ≈ 2.9, above Langevin. Hydrocarbon K0 1.8–2.25 gives 1.5–1.85×. Verdict 1 already passes at ×1.75 with feedback and the 166 ns shaper.

**Next steps:**
1. T14 chain point with `--ion-model reconstructed --ion-mobility K0=2.05`, and 1.8 / 2.25 as the systematic. It should land close to verdict 1, since the field is uniform where the ions drift.
2. On lxplus, check that the Garfield Ar+ table resolves (`$GARFIELD_INSTALL/share/Garfield/Data`) before using the default `--ion-mobility`; `K0=` needs no file.
3. Still open from below: retract the S3 "K0 ~1.16" claim and fix the `mx17_aval_calib.py:174` comment.

**Gotchas / decisions:**
- If the real MX17 mesh is calendered (thinner weave), the equivalent plane moves down a few µm and the field goes up ~1–2 %. The map assumes an uncalendered 37 µm weave.
- The reconstruction uses plane-averaged E_z and psi. That is fine in the bulk, where ions spend >90 % of their time, and approximate in the last ~15 µm.

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
1. ~~Find out why the T6 field map's mid-gap field is ~19 % below V/d.~~ Done 2026-10-16: it was a 128 vs 150 µm gap mix-up, and the map is right (see the section above).
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
