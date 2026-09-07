# Audit fix work order — 2026-08-07

Companion to `AUDIT_2026-08-07.md` (the findings) — this file is the implementation spec.
Audience: an implementing agent who has NOT read the audit conversation. Every fix gives the
context, the exact location, the change, an acceptance criterion, and the traps. The
undershoot fix (P6) is ALREADY DONE (`response/dream/shaper.py`, `--pzc-residual`) — do not
redo it.

## Ground rules (read before touching anything)

1. Read `design/RESPONSE_SIM_PLAN.md` §0a and `design/report/AUDIT_2026-08-07.md` first.
2. **Compute policy (plan §10)**: this laptop is orchestration-only. The LUT rebuild (Fix 1),
   any S1 re-solve (Fix 11), and Stage B re-runs go to the **desktop** (`ssh desktop`,
   62 GB); EOS access from the desktop is `rsync … lxplus:/eos/…` with a live Kerberos
   ticket (`klist`). The S3 raw slices (Fix 7) live on EOS and should be processed **on
   lxplus**, not pulled locally (19 GB).
3. **The §9 firewall stands**: nothing here is tuning. Fix 1 will move the halo numbers —
   re-running the ρ_s table afterwards is a *prediction update* and must be recorded as
   such, not iterated against the measured values.
4. **Process rule (plan §0a)**: never re-derive channel indexing ad hoc. Certify new code
   against the solver's own helpers (`kernels.charge_budget_y/_x`, `sum_over_rows`,
   `sum_over_columns`) — four ad-hoc re-derivations were wrong in one day.
5. Several agents share this repo and this machine. `git log` before editing C++; check
   `free -g` before blaming a hang on your own code.
6. After each fix, run the acceptance test named for it, and correct the plan text where the
   fix list says so (the plan is the authority document — plan §0a header).
7. Suggested commit granularity: one commit per numbered fix, message stating what physics
   or bookkeeping changed and what test now guards it.

Order of work: Fixes 2–5 and the C-list are laptop-scale and independent — do them in any
order. Fix 1 (desktop, LUT rebuild + re-runs) after the mechanical fixes so the re-run
picks everything up at once. Fixes 6–7 gate the P1 gas switch, not the current dry chain.

---

## Fix 1 — Rebuild the induction LUT with a window that covers the DAQ frame  [MAJOR]

**Context.** `response/digitizer/kernel_lut.py:47` defaults `t_max_ns=1000.0`. The S1
product is solved out to 10 µs; the DAQ area window is 32 × 60 ns = 1.92 µs (64 × 60 =
3.84 µs for the SPS config). Sheet transport feeds d=±2 on the ~0.6 µs scale and d=±3 past
1 µs (D ≈ 1 m²/s at ρ_s = 2 MΩ/sq), so truncating the kernel at 1 µs guarantees the sim
under-delivers the late halo — the same direction as the headline §9 "missing broad slow
halo" tension. Separately, the longitudinal convolutions truncate at the same edge:
`kernel_lut.py:214-221` (`apply_ion_transit` output cut to `nt`) and `:283` (`y[:, :nt]`
after the FFT convolution) drop any ion charge whose electron-current parent sits in the
last ~340 ns of the window.

**Change.**
- Raise the default `t_max_ns` to **3000** (covers the bench 1.92 µs window + ion transit
  340 ns + margin; make it a documented parameter so an SPS-config run can ask for 4200).
- Pad the two longitudinal convolutions to `nt + len(h_ion)` before truncating, so the
  *stored* window is complete rather than silently losing its own tail. The truncation to
  the output length should happen once, at the end, and be documented as "charge beyond
  t_max is deliberately outside the model".
- `response/digitizer/run.py:122-125`: bump `--n-samp` default so the frame still holds
  drift (~770 ns) + kernel (3000) + shaped tail; 4500 is safe. Update the help text, which
  currently hardcodes the 1000 ns kernel window in its arithmetic.

**Sizing trap.** The LUT build peaked at ~3.7 GB at nt=1000; expect ~3× that. Build on the
**desktop**. Also re-check the plan §0a note that "opening t_max makes it worse at fixed
n_side": that observation was charge dispersing into more channels than the ±`n_side`
window holds — after this fix, verify the audit's `n_side`/`y_window_mm` still capture the
±3 channels the §9 table needs (they should: `n_side=4` ≥ 3), and record the audit numbers
at both `n_side=4` and `8` once.

**Acceptance.**
- `test_lut_vs_solver.py` still passes at its bar (the reference is windowless, so the
  residual should *shrink*).
- `test_charge_audit.py` residual: record before/after. (See Fix 4 — the reference should
  also become per-position.)
- `test_longitudinal.py` passes; add one assertion there: an impulse placed in the LAST
  ion-transit-length of the window must still conserve charge in the padded scheme.
- Then re-run the 400-muon ρ_s table (plan §7 tables) on the desktop and record how the
  d=±2/±3 area shares and the Y-view numbers move. Expectation: halo shares rise. Record
  as a prediction revision in the plan with the audit as the cause.

---

## Fix 2 — `sharing_report`: select the in-gap deposit by ESL phase, not absolute x  [MAJOR]

**Context.** `response/solver/kernels.py:356-368` wants an in-gap deposit at x = 15.6 mm,
but `nearest_column` (`kernels.py:286-291`) minimizes |x_pad − x_want| in absolute x. Near
15.6 mm the 800/780 µm beat puts both candidate pad centres (cols 35/36, x = 15.21/
15.99 mm) at ESL phase ~10 µm from a **strip centre** (production phase convention
`ESL_PHASE_CENTERED_M`: strip occupies phase < 275 µm or > 525 µm; gap centre is 400 µm).
So the plan §3 headline table's "in-gap" row, and `meta.sharing["in-gap"]` in all 12
production products, are a second ON-STRIP deposit differing only in pad-ownership parity.
The G arrays themselves are correct; `validation/tau_g_reinterpretation.py` selects by an
explicit phase cut and is unaffected.

**Change.** Add a phase-based selector: among pad centres (optionally restricted to a
parity), pick the one minimizing circular distance in ESL phase to the target phase
(gap centre = 400 µm under the centered convention). Genuinely-in-gap pad centres exist at
cols 10–21 (x ≈ 0.4–4.3 mm and 26.9–30.8 mm); the reviewer's worked example is col 16,
phase 390 µm. Have `sharing_report` print the selected column, x0, phase, and σ_s(x0) so
this class of mislabel is self-evident next time.

**Acceptance.** For the new in-gap row: `esl_sigma_profile(x0) == 0` (in gap). Physically
expect a qualitative change from the mislabeled row: a gap deposit cannot transport along a
strip, so its late redistribution within the owning view should be suppressed and the
cross-view late feed pattern should differ from on-strip; the near-identical late view
totals (0.366 vs 0.365) that flagged the bug must NOT recur. Then:
- regenerate the sharing block for the nominal production point (a re-run of `run_point`
  is ~6 min on the desktop; do NOT re-solve all 12 points — patch their meta only if cheap,
  otherwise record in the manifest that pre-fix meta carries a mislabeled sharing block);
- replace the plan §3 table's in-gap row (keep the old row struck through with a pointer to
  the audit, per house style).

---

## Fix 3 — One trigger phase per event, real ftst semantics  [MAJOR — blocks T13b/T14 timing]

**Context.** `response/digitizer/run.py:277-284` samples the X and Y planes with two
independent uniform phases (`response/dream/daq.py:107-108` draws inside `sample()`), and
`daq.py:173-174` writes `ftst = 0` for every event. Real data: one trigger, one phase,
recorded in `ftst` (values 0–5, verified nonzero/varying in
`MX17_long_run_datrun_260506_09H28_000_03.root`), and wft *uses* it:
`nTof_x17/wft/calibrate.py:311-336` builds a `dt_xy` correction keyed by `ftst_x − ftst_y`,
applied in `reco.py:270`. The sim as-is injects ~24.5 ns rms of inter-plane jitter that no
real run has and that nothing downstream can correct.

**Change.**
- Draw ONE phase per event (in `run.py`'s event loop, or better: extend `Daq.sample` to
  accept `phase_ns` and have `run.py` pass the same draw to both views).
- Derive and write ftst from that phase with data semantics. The observed 0–5 range over a
  60 ns period suggests 10 ns ticks (6 states); **verify against wft's usage and the FEU
  manual** (`~/x17/Documents/dream/221217_FeuUsersManual.pdf`) before hardcoding — the
  implementer must confirm (a) tick size, (b) whether X and Y FEUs of one event carry the
  same or independently-latched ftst in real data (check a real paired run:
  `ftst_x − ftst_y` distribution width in the det3 long run). Mirror whatever the data does.
- Thread the per-event ftst array through `daq.write_decoded` (`daq.py:152-174` already
  accepts `ftst=`, it is just never passed).

**Acceptance.** A paired-file check on sim output: distribution of (X peak-sample time −
Y peak-sample time) narrows vs before; `ftst` distribution uniform over its states and
IDENTICAL between the two FEU files if that is what data shows. Ideal closure: run wft's
`calibrate` dt_xy machinery on the sim unchanged and confirm it finds the same
phase-vs-ftst structure it finds in data.

---

## Fix 4 — Anchor the real-pad kernel; make the charge audit self-explaining  [MAJOR]

**Context.** The closed-form sum rule certifies only fictitious pitch-sized tiling pads
(`kernels.py:427-429`). The production real-pad capture 0.665023 agrees with the
semi-analytic (0.68/0.78)² × 0.875 = 0.6656 to 0.08 % — but only in a comment
(`test_charge_audit.py:118`), never asserted. `test_lut_vs_solver.py:9-16` claims the LUT
reference is "the same code path whose closed-form sum rule passes" — false:
`charge_budget_y/_x` (`kernels.py:294-316`) is separate indexing from
`sum_over_rows/columns` (`kernels.py:251, 273`). And the charge audit's +1.5 % residual is
an EXCESS the audit's own leak framings (which all predict deficits) cannot produce —
almost certainly the probe positions sampling capture above the grid mean (d=0 capture
spans 0.22–0.35 across columns).

**Change.**
- In `test_charge_audit.py`: assert
  `abs(channel_capture_prompt − (PAD_SIZE/PAD_PITCH)**2 * sum_rule_expect) < 0.005`
  (bar chosen: measured agreement 8e-4; 0.5 % leaves room for fringe fields while still
  catching any indexing-level error).
- Make the audit's reference per-position: for each probe, compute the expected capture at
  THAT (x, y) from `charge_budget_y + charge_budget_x` (total over all channels), instead
  of comparing the probe mean against the grid-mean scalar. Print per-probe residuals and
  their spread. The +1.5 % should collapse to interpolation-level; then tighten `TOL_LEAK`
  from 3 % to 1 %.
- Rewrite the `test_lut_vs_solver.py` docstring to state what it actually certifies: two
  independently-written indexings agree to 1e-4; the closed-form anchor applies to a
  different (tiling-pad) geometry; the real-pad anchor is the assert added above.

**Acceptance.** New assert passes; per-probe audit residual spread printed; a deliberately
broken parity (flip one `% 2` locally, revert after) makes the new per-position audit fail
loudly — verify it does, since the old grid-mean version was blind to redistribution.

---

## Fix 5 — Record the wet-CF₄ bracket check (the one that gates P1)  [MAJOR — data product check]

**Context.** Plan §0a's cert row "wet-CF₄ v_drift 74.32 vs 74.7 µm/ns" compares the dry
control column of the table against a published *dry* Magboltz number — dry-vs-dry.
`response/params/wet_cf4_drift.py:122-124` defines the check that matters: the wet curves
(1.0/1.5/2.0 % H₂O) must bracket the MEASURED 13–15 µm/ns at 233 V/cm, "if the wet curves
do not bracket it … nothing downstream should be run." That result is recorded nowhere.

**Change.** On lxplus (or desktop via rsync): fetch
`/eos/experiment/ntof/data/x17/response_sim/gas/wet_cf4_drift.json`, evaluate v_drift at
233 V/cm for each mixture, and record in the plan (cert table row corrected to say what is
now actually certified) AND in a small committed JSON/text under `design/report/`. Include
transverse diffusion at the same point — it is the input the halo prediction needs.

**Acceptance.** Numbers exist in-repo for all four mixtures at 233 V/cm; the plan row no
longer claims wet validation from a dry comparison. If the bracket FAILS: stop, mark P1
blocked in the plan, do not proceed with any wet-gas Stage B work.

---

## Fix 6 — Attachment through the chain  [MAJOR for P1; no-op dry]

**Context.** Plan §7 step 2 specifies survival e^(−z/λ); nothing implements it.
`response/digitizer/gas.py:53-56` interpolates only `v, dL, dT` (the dry table has no eta
key); `response/digitizer/digitize.py:139-177` (`transport`) has no survival term;
`params/wet_cf4_drift.py:96-103` stores `eta_per_cm` that nothing reads. CF₄+H₂O attaches;
attachment is z-dependent and re-weights the depth mixture that produces the c1-vs-z
diffusion signature — it cannot be folded into a flat gain factor.

**Change.**
- `gas.py`: interpolate `eta_per_cm` when present, return 0.0 when the table lacks it
  (dry tables keep working; log once which branch was taken into provenance).
- `digitize.transport`: survival probability p = exp(−eta_per_cm × z_cm) per electron;
  combine with the mesh-transparency Bernoulli into a single thinning
  (p_surv = ε_mesh × e^(−ηz)) so the statistics stay binomial and the truth accounting
  sees one loss channel with two named factors.

**Acceptance.** With the dry table: bitwise-identical results at fixed seed (η = 0 path).
With a synthetic table carrying constant η: measured survival vs z matches exp(−ηz) over
the 30 mm gap to statistical precision, and c1-vs-z reweighting moves the right way
(deep/high-diffusion clusters suppressed).

---

## Fix 7 — Recover the avalanche survival factor from the S3 raw slices  [MAJOR — feeds P2]

**Context.** Plan §5's "avalanche size cap bites" is false — there is no cap
(`--max-avalanche` defaults to 0/uncapped, `avalanche/mx17_aval_calib.py:96-98, 185-186`);
the falling nev column is the deliberate per-voltage schedule in `mx17_aval_points.txt`.
The REAL defect: `avalanche/collect.py:38` fits the Polya on `g[g > 0]` (conditional on
survival), the per-slice `survival` field present in each raw JSON
(`mx17_aval_calib.py:293`) is dropped by `reduce_file`/`merge`, and
`digitize.py:176` applies that conditional Polya to every mesh-surviving electron. Per-
electron charge is biased high by 1/P(g>0), magnitude unknown from the shipped calib.
Direction feeds P2 (16.5 % ADC saturation).

**Change.**
- On **lxplus** (raw is 19 GB on EOS `response_sim/avalanche/raw/`): extend `collect.py` to
  carry `survival` (per V, with its binomial error) into the merged calib → bump to
  `aval_calib_v3.json`, keep v2 untouched on EOS.
- `digitize.py`: Bernoulli-thin each electron by `survival(V)` BEFORE the Polya draw (this
  is exactly the conditional-distribution decomposition, so the existing conditional Polya
  stays correct as-is). Wire the value from the calib; refuse to run silently if the calib
  lacks the field (explicit warning + survival=1 with provenance tag).
- Correct plan §5's cap sentence (T7 row too).

**Acceptance.** Two routes to mean charge per drifted electron agree:
survival × mean(conditional Polya) vs the direct mean over all seeds including zeros
(computable from raw during the collect run — assert it there). Record survival(490 V) and
the resulting change in the P2 saturation fraction in the plan.

---

## Fix 8 — Name decoded output by the target run's FEU ids  [MAJOR for T14 mechanics]

**Context.** `run.py:374` writes `sim_decoded_07.root` (X) / `_08.root` (Y); the det3
detector the chain mimics is x→FEU **3**, y→FEU **4** (run.py's own comment at 221-229,
det3 `run_config.json`, and the real det3 decoded files `*_03/_04.root`). In wft the file
suffix IS the FEU id (`wft/io.py:39-41`) and keys the strip map; `feus=[7,8]` in
`qa_config.py` is a different physical detector (mx17_2). Left as-is, a det3-configured wft
finds no files; "fixed" by renaming, it inherits the wrong strip map.

**Change.** Add `--feu-ids X_ID Y_ID` (default `3 4` for the det3 bench) controlling both
the filenames and any in-file FEU identification; update `--decoded-out` help text and
`test_daq_wft.py` accordingly. Document in the plan (P5/§8) that the default now matches
det3 and that an SPS-run comparison must set the ids from that run's config.

**Acceptance.** `test_daq_wft.py` passes with the new default names; a wft io smoke read
with a det3-style `feus=[3,4]` config discovers both files.

---

## C-list — smaller but real (one commit each unless trivial)

**C1. Noise scale bias (~3.7 % loud).** `response/dream/noise.py:186-189`: per-channel
traces are generated from the channel-averaged spectrum then scaled by
`sigma_chan / median(sigma_chan)`. The generated trace's rms is the Parseval rms of the
mean spectrum (root of MEAN variance), not the median sigma; the sigma distribution is
right-skewed, so everything is +3.7 %. Fix: divide by the Parseval rms of `amp_chan`
(compute it once from the amplitude spectrum, same convention as `_coloured`). Acceptance:
`noise.selftest()` per-channel residual drops from ~3.5 % to <1 %; tighten its tolerance.

**C2. X-view one-pad-off booking (~1.3 % of avalanches).**
`response/digitizer/digitize.py:195` derives `col` from the true x while
`kernel_lut.py:80-88` defines band offsets relative to the nearest 40 µm LUT sample; every
second pad boundary falls between the two conventions. Fix: derive the column from the SAME
folded LUT sample used for the kernel (one line). Acceptance: extend `test_lut_vs_solver`
with probes at pad-boundary ± 5 µm drawn OFF the LUT grid, comparing against
`charge_budget_*` at the probe's true position: the d-profile must not shift by a whole pad
for any probe.

**C3. Packet-mode Polya variance.** `digitize.py:166-176`: one single-electron Polya draw
scales the whole packet → variance n² Var(g) instead of n Var(g). Fix: draw the packet gain
as `rng.gamma(shape=n_surv*(1+theta), scale=gbar/(1+theta))` (sum of n independent Polyas).
Production `packet=False` is unaffected; fix before profiling §12 item 11. Acceptance: for
n=100 packets, packet-mode charge variance matches per-electron mode to a few %.

**C4. `charge_budget_x` y-fold.** `response/solver/kernels.py:313`: `argmin(|sx.y − y0_m|)`
clamps out-of-box y to the edge. Fold y0 into the 1.56 mm box (mod, mirroring
`charge_budget_y`'s x-fold at `:302`). Acceptance: `charge_budget_x(..., y0_m=5e-3)` equals
`charge_budget_x(..., y0_m=5e-3 % box)` exactly.

**C5. `lut.ix()` seam wrap.** `kernel_lut.py:149-152`: nearest-sample lookup clips at the
superperiod edge; the last half-stride should wrap to x[0] (circular metric mod 31.2 mm).
Also add one probe in `test_lut_vs_solver` that goes THROUGH `ix()` (the current test
bypasses it with argmin). Acceptance: a source at 31.19 mm maps to the x[0] sample.

**C6. ny convergence certificate.** ny=512 (97.5 µm) leaves the PROMPT kernel's pad-edge
shoulder in y wrong by up to ±0.028 absolute (≤5e-4 late; 0.003 after a 150 µm diffusion
smear). One-off, desktop: re-solve the nominal point at ny=1024, difference the prompt
kernels, store the result as `test_ny_grid`-style evidence (mirror `test_time_grid.py`).
If the smeared error is confirmed ≤0.3 %, keep ny=512 with the certificate recorded; if
not, the production grid needs ny=1024 (≈2× cost, still minutes/point).

**C7. Amp-gap/groove clusters.** `digitize.py:158` clips z<0 to 0 and then gives those
clusters full-gap gain AND mesh transparency, contradicting `clusters.py:20-24`. ~0.5 % of
track length, ~9× overweight, prompt and unshared. Decide and implement ONE of: (a) drop
amp-gap clusters with a provenance note (defensible: sub-percent), or (b) partial gain
G_partial = exp(ln(G) · z_birth/gap) with NO mesh transparency for them. Either way make
`clusters.py`'s contract and `digitize.py` agree. Acceptance: a synthetic ClusterTree with
only amp-gap clusters produces the documented behaviour.

**C8. Fano factor — decide.** Plan §7 step 1 specifies F≈0.2; neither stage implements it
(Stage A: `SteppingAction.cc:79-82` floor+Bernoulli; Stage B refuses to re-derive). Either
implement in Stage A (per-cluster n from a truncated Normal(n̄, √(F·n̄)) — cheap, keeps
Stage B untouched) or de-scope with a plan edit stating resolution-type observables are
simulated slightly narrow. Do not leave the plan promising what the code doesn't do.

**C9. Truth block.** Plan §7 step 7 requires per-event truth (true (x,y) incl. the run.py
offsets, t0, per-channel true charge, electrons lost to mesh/survival). Write it as a
sidecar parquet/npz next to `--decoded-out` keyed by eventId. Acceptance: eventId join
between truth and decoded files is total.

**C10. Keep empty events.** `run.py:264` `continue`s zero-electron events before the DAQ
buffer, so decoded files omit triggers real runs record — breaks occupancy comparisons and
eventId alignment. Write noise-only frames for them (they still get pedestal+noise).

**C11. Gas table path.** `gas.py:37` hardcodes `~/PycharmProjects/nTof_x17/...` — breaks
condor. Make it a parameter with env-var override; embed the table's path AND a content
hash in provenance.

**C12. Stale analytic-ion default.** `digitize.py:98` `z_aval_um=5.0` — the S3 v2
measurement says 13.84 µm (mean ion birth height). Update the default and note that
`--ion-model analytic` now reflects v2; the measured template remains production default.

**C13. Portability/robustness nits.** `spread.py:148` `sp.ptp()` → `np.ptp(sp)` (removed
from ndarray in NumPy 2.x). `manifest.py:95-97` formats `None` with `:9.1e` and dies if a
product lacks the sum-rule block → print `FAIL/absent` instead. `wpot.py:247` `_eig_cache`
is dead — delete it, and (optional, free 2×) reuse the eigendecomposition for −k_y (M
depends on k_y² only).

**C14. Stage A tidy-ups** (all MINOR, from the C++ review): `EventAction.cc:35,37` divides
by eV twice in the verbose print (1e6 off; files unaffected). `mm_sim.cc:159-163` vs
`macros/run_default.mac`: supplying the macro runs zero events (no beamOn) while the macro
header claims otherwise — make `main()` issue beamOn after a macro that lacks one, or fix
the macro/header. `RunAction.cc:178-196`: Meta branch addresses point at block-scope locals
— hoist to function scope. Record the primary vertex (x,y) in the EventTree so
`--beam-spread` truth is exact, and stop `run.py` double-randomizing when the Geant4 file
already carries a spread (detect via Meta or a provenance flag; today's default re-offsets
in post and mislabels provenance, statistically benign but wrong bookkeeping). Add the
missing fiducial handling: clamp/flag `col,row` outside 0..511 in `digitize.py:195-196` so
the budget and decoded paths agree on out-of-area delta rays.

**C15. ZS bookkeeping.** Code (`ZsChkSmp=1`) matches det3's config; plan §8.5's "ZsChkSmp=4"
is the July-beam setting — fix the plan line, and parameterize ZsChkSmp from the run config
rather than a constant if ZS is ever switched on. Verify whether firmware keeps pre-samples
before any ZS-mode comparison (unverified in the audit).

---

## Plan-text corrections (no code) — do in one commit

In `design/RESPONSE_SIM_PLAN.md`:
1. §5 / T7: remove the "avalanche size cap bites" explanation (no cap exists; nev follows
   `mx17_aval_points.txt`); note the survival-factor issue and its fix (Fix 7).
2. §8 step 5: "saturation at 3550" → clip is 0..4095 (data reaches 4095); 3550 is wft's
   censoring threshold on pedestal-subtracted W (`wft/model.py:57`).
3. V5: reconcile §0a "certified 7.8e-7" with T5 "V5 still not run" — the script EVALUATES
   exp(−2πng/p); call it an analytic bound, not a test, in both places.
4. Cert table: label the noise round-trip, ADC end-to-end, and ion A/B rows as
   internal-consistency checks (not data anchors), and replace the "98.0 vs 102.4" row —
   that gap was the test contaminating its own pedestal with injections
   (`test_daq_wft.py:55-64`); with clean pedestal events the scale closes at 102.0/102.4.
   Optionally fix the test: use signal-free events for the pedestal estimate, then tighten
   its 10 % tolerance (it would currently pass a +13 % gain error).
5. Per-product sum-rule scope: `kernels.py:492-493` certifies t ≤ ~0.2 ns at reduced grid;
   either extend `times[:6]` to a spread of the full axis (cheap) or say so in T5.
6. P5: the connector-inversion ambiguity is resolved — `Mx17StripMap.apply_orientation`
   reverses per connector (64 ch), unconditionally; the sim's deliberate identity mapping
   stays, but the plan should stop calling it unresolved.
7. V2's description: it tests charge conservation + convergence order (per its docstring),
   not the Galan closed form §3 names.
8. §3 in-gap table row: superseded by Fix 2's regenerated numbers.
9. x/y sign convention (`ActiveAreaFrame.hh:29-38` placeholder): still OPEN — leave open,
   but add the audit's note that a flip inverts checkerboard parity AND the 31.2 mm beat
   phase, and that the beat phase vs absolute channel number at T14 is the in-data check.

## Explicitly NOT in scope here

- T10 slow path (ion lateral shape, P3), T6 field-map export (P4), V6 — tracked in the plan.
- Any comparison against §9 numbers beyond recording how predictions moved (firewall).
- The undershoot (done). β stays a scan parameter until T14.
