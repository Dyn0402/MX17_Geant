# det3 X neighbour asymmetry: ledger of every route tried (2026-10-18 → 2026-10-20)

Read this first; the round-by-round detail follows below it (and in `../NB_ASYM.md`, `../NB_FAR.md`).
Data: raw-ADC survey of every cosmic-bench sub-run (`raw_asym.py`, yoga `/media/dylan/data/claude_tmp/tmp/nb/raw/`).
Sim: v5a recipe on condor (4505153/4505159/4507004/4507242), npz in `~/x17/nb_asym_raw/simxoff/` on yoga.
a1 = α(n+1) − α(n−1) in DAQ order (det3 June-15's excess side positive).

## Current picture (end of 2026-10-20)

**What the asymmetry is:**
- Not a det3 defect. Every June build (det2, det3 June-15, det4, det6, det7) carries an X-only neighbour-asymmetry field.
- The field varies smoothly over each board and takes both signs. The Y view is flat everywhere.
- det2's field comes as bands in y (−0.29…+0.18) that cancel, so det2 looked symmetric. det3 June-15 is positive over
  most of its area (mean +0.14, with a tilted ridge).
- Per board, the geometry is fixed; the amplitude scales with narrower clouds (gas, drift field). That is why det3
  varies between nights.
- It is charge moved within X (n−1 → n+1), prompt, and repeats at the pad pitch (no 31.2 mm beat on any board).

**Leading explanation (partly supported):**
- Resistive-stripe registration: stripe centres offset from the pad columns, with the offset varying across the board
  through print scale and rotation.
- In the sim, a1 vs the true stripe offset reproduces the shape and both signs at ~1/3 of the amplitude.
- This requires a 780 µm-commensurate stripe pitch, which contradicts the 800 µm quote. Dylan is fetching the mask.
- Candidate for the missing strength: the June boards leak charge ACROSS the stripes (late X ±2 content), the sim
  doesn't.

**Open:**
- The size gap (sim ±0.05–0.10 vs data ±0.14–0.29).
- Why only the June builds have it.
- The stripe pitch.

## Routes

Verdict key: ✗ excluded · ~ partial / doesn't fit alone · ✓ supported · ? open.

| # | Hypothesis / route | Test | Result | Verdict | Where |
|---|---|---|---|---|---|
| 1 | Artefact of our code (maps, clustering, CM subtraction) | Independent raw-ADC extractor, DAQ channels, pedestal only | Asymmetry is in the raw ADC: P(n+1>n−1) 0.85 on det3 June-15 X, 0.50 elsewhere | ✗ | top of this file, deck slides 4–8 |
| 2 | Hot channel / selection trap | Hot-channel veto, pulse-shape cut | One evening's "turn-on" was ch 492; the effect survives the veto | ✗ (trap found) | finding 3 above |
| 3 | 39:40 resistive beat (800 vs 780 µm) | a1 vs column mod 40, DAQ then true columns; per-column-means periodogram, permutation null | Flat on every board (det3 sine 0.007 ± 0.004; 1/40 p 0.12–0.86). The sim predicts s1 0.12 / d1 0.07: absent in data | ✗ as a beat | NB_ASYM §1; Third round 1; Fifth round 2 |
| 4 | Electronics copy (FEU, ADC mux, cable) | Noise correlation d=1 vs d=2 with injection; same FEU slot on other boards | Noise not copied; same FEU clean on det3 May/det6 | ✗ | NB_ASYM §1; deck "Noise", "Follows the board" |
| 5 | Connector / M1 card / channel mixing | Linear-mixing fits (any sign, ±1/±2); connector-edge re-extraction (`--edge`) | No mixing fits (χ² 410–750/28); the excess follows the physical neighbour across connectors (+0.21 ± 0.03), DAQ n+1 gets 0 | ✗ | Third round; "Connector edges" |
| 6 | Track angle / tilt | a1 vs M3 angle; vertical-only maps | −0.0085/° on det2: det3 would need ~20°; maps hold for vertical tracks; an x-tilt can't make y-structure | ✗ | NB_ASYM §1; Fifth round 9 |
| 7 | Gas-flow pattern | Map vs inlet→exhaust, X-only, time behaviour | Flat top–bottom, X only, non-monotonic in time | ✗ as cause | deck "Gas idea" |
| 8 | One-sided copper (bus traces under the gap, stubs) | 3D weighting solve (`bus_wpot.py`) at any bus depth, dielectric variants | ≤ 0.01 nominal, 0.085 even touching the pads; fixed handedness | ✗ (also killed by #22) | NB_ASYM §B; deck "Bus traces" |
| 9 | Floating stripes (static conductors) | `bus_float.py` + forward model | Copy in both pad classes, growing with cloud size | ✗ (static version) | Second round 6 |
| 10 | Plain one-step copy x(n+1) += a x(n) | Uncopy closure; parity split with injection; bookkeeping | Uncopy balances ±1/±2 board-wide. But det3 isn't a plain copy: narrow X-pad −0.40, which no copy makes | ~ | NB_ASYM §A; Third round 2–3 |
| 11 | Charge over Y pads coupled into X n+1 ("lives at the Y pads") | X/Y parity split | Split is partly peak-migration selection (injected copy makes 0.20/0.05) | ~ (selection) | Second round; Third round 2 |
| 12 | Rigid X-view offset (X reads the avalanche at x + D) | Static forward fits; full sim condor 4505159 (D 0…−300) | Pre-registered narrow X-pad flip passed (−0.45…−0.51). But zero marginal a1 and fading for wide clouds | ~ needed, not sufficient | Fourth round 1 |
| 13 | Mixed response (1−w)R(x) + wR(x+D) | Exact event-level mixing of offset arms (`raw_asym.py --mix-with`); fingerprint scoring | a1 ≈ 0.00–0.03 for every (D, w): can't make the marginal | ✗ alone | Fifth round 4, 6 |
| 14 | Y→X transfer (floating or cracked-via Y pad, Y-line crosstalk, bus pickup from Y) | ΣX/ΣY by pad class and session | ΣX/ΣY not raised (det3 0.907/0.925 vs det2 0.932/0.944) | ✗ | Fifth round 1 |
| 15 | Extra sensitivity regions on the bus side (gap strip, next-column Y pad or its half, next or skip X line, whole column) | Fast static model (`fwd_region.py`) | All send narrow A_X positive and keep A_Y small (χ² 700–3400/24) | ✗ | Fifth round 7 |
| 16 | Skip transfer n−1 → n+1 (channel level) | Fast model | Right x-conditional trend, wrong class pattern | ✗ | Fifth round 7 |
| 17 | Shift/mix + copy combinations | Fast model, full sim | Best χ² ≈ 240–260/24 (fast), no session fit (full) | ~ | Fourth round 1; Fifth round 7 |
| 18 | Stripe orientation wrong on det3 | Late ±2 content X vs Y per board | det3 June-15 = Y-slow like every June/May board. det0/det1 (Jan/Apr) reversed. Apr-24 "horizontal, broken" shows no spreading | ✗ for det3 | Fourth round 2 |
| 19 | Same PCB for det3 May and June-15 | Dead-line fingerprint | No shared defects; not supported | ? | "Connector edges" 4 |
| 20 | Charge bookkeeping: moved or added? | Integrals per Y-share bin vs det2/det3 May and vs injection | n+1 gains ≈ what n−1 loses; Xc and the sum unchanged; a copy would raise n−1 and Xc | ✓ moved | Fifth round 1 |
| 21 | Timing of the excess | Mean waveforms aligned on the Y central | Excess peaks with the central (+3…+25 ns), not like neighbour charge (+90…+120 for offsets) | ✓ prompt, direct-like | Fifth round 3 |
| 22 | Is it board-uniform? | Per-chip, 16-col, 8×8 maps; connector steps; Y-view maps | Smooth two-signed X fields (det4 −0.16 → +0.19, det6 + → −, det2 bands), no connector steps, Y flat | ✓ X-only field | Fifth round 8–9 |
| 23 | det2 a clean control? | det2 by region; both runs | Positive band = det3's full fingerprint incl. timing; negative bands its mirror; same rows on FEU 3 and 6, amplitude ×2.5 with narrower clouds | ✗ not clean; ✓ fixed geometry | Fifth round 10, 13 |
| 24 | Pitch of whatever causes it | a1 inside det2's positive band in 8-col blocks; 1/40 harmonic | Positive across all 512 columns, no 1/40 (p 0.57): commensurate with 780 µm | ✓ pad-pitch locked | Fifth round 11 |
| 25 | Stripe-registration field (print scale/rotation) | Map fits A sin(φ0 + kx x + ky y); quadrature (sharing ∝ cos φ) | Fits det3/det2 geometry, stable across det3 sessions; sharing quadrature not detected (−0.007 ± 0.009) | ~ geometry yes, quadrature no | Fifth round 12; `quad_test.py` |
| 26 | Can stripe offset produce it (sim)? | a1 vs true δ in the sim (`sim_delta.py`, `sim_delta_fp.py`) | Odd in δ, ±0.05–0.10 (narrow ±0.10–0.24); det3's class pattern at +δ, det2-negative mirror at −δ; ~1/3 of data amplitude | ✓ shape, ~ size | Sixth round 1 |
| 27 | Rotation of the print alone | Physics estimate | Lateral spread along a tilted stripe is µm; only the registration via a 780 pitch makes ~P/θ bands | rotation only via registration | Sixth round 1 |
| 28 | Pin to stripes via reversed-stripe boards (det0, det1) | X and Y 2D maps, det0 He + Ar/iso, det1 all files | No field in either view (det1 weak fine-scale Y hint); May normal boards have none either | ? can't pin | Sixth round 2 |
| 29 | Why only June builds | Along- and across-stripe late content per board | June: along 63–80, across 12–31; May/earlier: 24–37 / 1–9; sim: 64 / −2 | the across-stripe part = gap conduction f ≈ 0.05–0.1 (#34), which does not itself make the asymmetry | Sixth round 3; Seventh 12 |
| 30 | Earlier side threads | ±1/±2 neighbour timing vs v5a; far ±3..±5 halo (displaced photon feedback, H 150 µm) | v5a matches copy-free X; halo explained by feedback spread | separate threads | NB_ASYM §C–E, NB_FAR.md |
| 31 | Static relocation (strip charge spread across/along its stripe; gap charge stays, funnelled, or leaked to neighbour stripes) | Deterministic fast model on the prompt ESL-plane field (`r5/esl_fast.py`) | ≤ ±0.02 for every variant except funnel-all + 1.5 mm along-stripe smear (±0.10). It also misses the sim baseline (0 vs ±0.05), so it is not a calibrator | ~ direction only | Seventh round 1 |
| 32 | Where the sim's a1(δ) comes from | Production X kernel, y-resolved clouds (`r5/kernel_delta.py`) | 0 at t = 0, builds over 100–300 ns of sheet transport; ±0.06 for σ 250 µm next to a strip edge, ±0.015 for σ 450 | ✓ late sheet transport | Seventh round 2 |
| 33 | Is the sim's registration excess prompt like the data's? | Excess waveform vs central, sim \|δ\| ≥ 150 vs data (`r5/delta_timing.py`) | Sim excess leads the central by 70 ns (t50 −90); data by 10 ns (t50 −25); det2/May none | ~ wrong shape by ~1 sample | Seventh round 3 |
| 34 | Cross-gap leak (conducting ESL gaps, σ_gap = f σ_strip) | nx = 1560 X kernels f = 0…0.3 (4520685); full sim with hybrid kernels f = 0/0.03/0.1 (4520810), paired | Kernel: ×1.9 late for σ 250. Full sim: odd amplitude ×1.00 (f 0.03), ×0.90 (f 0.1). But X late ±2 −2.3 → +7.8 → +27.6 ×10⁻³ = the June boards' 12–31 | ✓ explains the June late cross-stripe leak (f ≈ 0.05–0.1); ✗ as the asymmetry amplifier | Seventh round 4, 10, 12 |
| 35 | Prompt funnelling of gap-landing electrons onto the strip edge | Digitizer `--gap-funnel F`, full v5a sim, 72 arms (condor 4520688) paired with the D = 0 seeds 7–12 (`r5/paired_odd.py`) | F = 1: odd amplitude ×1.40 (narrow ×1.54, +0.038 ± 0.003); F = 0.5: ×1.20; sharing unchanged; excess still leads the central by 90 ns | ~ real amplifier, ≤ 1.5×, wrong timing | Seventh round 7 |
| 36 | Strip resistivity ρ_s | Mini kernels at 0.5/1/2/5 MΩ/sq | Registration term saturates at the same size (−0.033, σ 250); low ρ_s only makes it prompt (by 92 ns at 0.5 MΩ) | timing lever, not size | Seventh round 8 |
| 37 | Strip width / gap ratio | nx = 1560 X kernels, insulating gaps, width 410/550/690 µm (effective 420/540/700; condor 4520736) | Odd amplitude (σ 250) ×0.6 at 420, ×1.6 (300 ns) / ×2.4 (1 µs) at 700 | ✓ strong lever; the June print's width is measurable from the mask | Seventh round 11 |

## Next checks (in order)
1. **Cross-gap leak / gap funnelling in the sim:** a digitizer variant giving late charge transfer between adjacent
   stripes (and/or electrons over the kapton gaps moved onto the nearest stripe). Score: the sim's a1(δ) amplitude
   (target ±0.15–0.25 at narrow clouds) and the X late ±2 content (target 12–31 ×10⁻³).
2. **Per-board a1(x, y) calibration map** (`r5/map_fit.py` cells) as the analysis correction, replacing a single
   uncopy a; det2 needs one too.
3. Deck slides for rounds 4–6 (the 2D maps figure `design/figures/asym_maps_2026-10-20.png`, det2 regions, sim a1(δ)).
4. When the mask arrives: stripe pitch and alignment → predicted δ(x, y) per board vs the measured fields.

---

# Raw-ADC check of the X asymmetry (2026-10-19)

Deck: https://dylan-neff.web.cern.ch/notes/det3-x-asymmetry.html, built by `make_asym_deck.py` (observables in
`asym_lib.py`). Rebuild: `~/PycharmProjects/nTof_x17/.venv/bin/python make_asym_deck.py` then
`python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py det3-x-asymmetry.html --slug det3-x-asymmetry --force --deploy`.
Raw survey done (151 detector-views x sub-runs): det3 June-15 X median 0.137 (0.04-0.21, 57 sub-runs), its Y 0.00 +- 0.04;
det3 May X -0.008; det4 X 0.05, det7 X 0.06; det1/2/6 ~0. 2D map (6-24, 4495 events): X blocks 1-2 ~0.25, X block 0 ~0.1.

Extractor `raw_asym.py` (yoga): decoded_root ADC minus per-channel pedestal only. No strip map, no common-mode
subtraction, no clustering, no M3; DAQ channel numbers. Outputs: yoga `/media/dylan/data/claude_tmp/tmp/nb/raw/`,
copied to desktop `~/x17/nb_asym_raw/`. Batch 3 (`run_raw3.sh`, every cosmic sub-run in ../survey_jobs.tsv, both views)
was started on yoga at wrap time and writes `raw/survey/` + `index.tsv`. Beam files give 0 events with this selection (flash).

Good-pulse cut used below: 200 < A_c < 3500, central peak sample 4..22, |mean(samples 0-2)| < 0.1 A_c.
P = fraction of events with A(p+1) > A(p-1); a1 = in-phase (p+1 minus p-1) projection on the central pulse.

| sample | P | a1 | per chip (0..7) |
|---|---|---|---|
| det3 June-15 X, 6-24 long (FEU 3) | 0.85 | +0.177 | .08 .10 .22 .29 .28 .22 .17 .10 |
| det3 June-15 X, 6-23 long (FEU 3) | 0.62 | +0.158 | .07 .15 .19 .21 .18 .16 .12 .09 |
| det3 June-15 X, 6-28 weekend (FEU 7) | 0.56 | +0.092 | .04 .11 .10 .09 .13 .08 .10 .02 |
| det3 June-15 X, 6-27 drift 100 V | 0.74 | +0.159 | -.03 .06 .09 .33 .23 .18 .14 .08 |
| det3 June-15 X, 6-27 drift 1100 V | 0.54 | +0.064 | |
| det3 June-15 Y | 0.49 | -0.008 | flat ~0 |
| det3 May build X (same FEU 3) | 0.49 | -0.012 | flat ~0 |
| det1 X / det2 X / det6 X (FEU 3) | 0.50/0.48/0.50 | +.004/-.009/+.012 | flat |
| det4 X / det7 X (Stephan June-22 batch) | 0.56/0.53 | +0.054/+0.034 | |

Findings:
1. The asymmetry is in the raw ADC: on det3 June-15 X the p+1 neighbour is bigger in 85 % of events; every
   other board/view through the same code sits at 50 %. (R-L)/c is 0.20-0.23 at every sample: an in-phase copy.
2. It is NOT uniform across the board: high in chips 2-5 (board centre), ~0.1 at the edges, same shape on 6-23, 6-24, 6-27.
   A PCB trace geometry would be uniform. Next: 2D map (match X/Y by eventId, d3jun_x + d3jun_y).
3. TRAP: the survey's low 6-23 short/longer values (0.035 -> 0.057, looked like a turn-on over the evening) are a hot
   channel (ch 492, 3259/4493 events of the short run). With a hot-channel veto the longer run has the long run's chip profile.
   The survey (copy_survey.csv) needs a hot-channel veto before any time-dependence claim.
4. Noise (quiet events, common mode removed): corr(d=1)-corr(d=2) is -0.03 on det3 June X, same as det3 Y and May X;
   injecting a 0.15 copy gives +0.13. Not copied in noise -> not the FEU chain after the input.
5. Open: whether the one-sided bus side predicts the observed direction (gerber -> map orientation unpinned), and
   which data view sits on L5 vs L6 (bus solve gives 0.010 on L5, 0.001 on L6). May and June-15 det3 may be different PCBs.

Original deck plan (slide-note, `make_asym_deck.py` to write): setup cartoon; suspects map; raw event; 3-step processing with
struck steps; coin-flip bars; L/c vs R/c hists; mean waveforms; copy shape; copy model + uncopy; hot-channel trap;
same-FEU grid; noise; board map; amplitude; survey; time/drift; bus geometry; verdicts; what it doesn't rule out.

## By date (slide 17 "By date", added 2026-10-19 evening)
Board-centre chips 2-5 only (the copy is centre-peaked; several 6-22 high-mesh points have 70-90 % of hits on chip 7, so
the board average there is low for geometric reasons). Session means, det3 June-15 X:
| session | FEU | drift | a (chips 2-5) |
|---|---|---|---|
| 6-18 long_run (EOS, new; only chips 0-3 read) | 7 | 1000 | 0.26 |
| 6-19/20 longer_run (EOS, new; 7458 ev, flat over 8 files, no centre peak) | 7 | 1000 | 0.03 |
| 6-22/23 | 3 | 1000 | 0.16 |
| 6-23/24 | 3 | 600 | 0.26 |
| 6-25 | 3 | 500 | 0.17 |
| 6-27/28 (resist HV 3:3 -> 3:4) | 7 | 1000 | 0.09-0.10 (drift scan to 100 V: 0.21) |
| 7-26 beam run_79 (old method, whole board, 90/10) | 3 | 700 | 0.05-0.07 |
Stable within a session, jumps between sessions; not a monotonic decay; FEU slot can't explain it (6-18 and 6-19 both FEU 7).
6-19 is the only run with det3's drift channel (0:7) misbehaving: 4-9 uA spikes, sag to 890-960 V (others flat 0.34 uA).
The 6-19 config swapped HV channels between the two boards vs 6-17; the other board that night (FEU 3) reads 0.08 +- 0.03
(372 noisy events), so no clean label swap, but a cabling mix-up isn't excluded. 6-26 weekend scan is empty on EOS;
6-19 HV scan + long_run exist only as raw_daq_data (could be decoded to extend 6-19/20).
Extra extractions: ~/x17/nb_asym_raw/extra/{d617_long_x,d619_longer_x,d619_longer_feu3_x}.npz (desktop only, derived).

## Second round (10-19 evening): the copy lives at the Y pads (deck slides "Y pads" .. "Simulations")
Numbers from `asym_more.py` → `asym_more.json`, `beam/beam_drift.py` → `beam_drift.json`, `../fwd_parity.py` → `../fwd_out/`.
1. **Parity split.** X/Y matched by (file index, eventId); (p_x + p_y) % 2 = pad type under the charge; the Y-pad class is
   the one with the smaller Y ±1 share and A_x/A_y. Copy (X-pad / Y-pad): det3 June-15 0.020 / **0.221**; det4 0.000 /
   **0.078**; det7 0.018 / **0.078**; det6 −0.004 / 0.017; det3 May, det1, det2 ≈ 0 / ≈ 0. Same phenomenon on all three.
   Flat mod 3/4/6/12 (pillars 4.68 mm = 6 pads excluded), no 40-ch period (stripe beat excluded), X channel parity alone flat.
   For Y-pad events it is a redistribution between the two X neighbours (α+1 ≈ 0.43, α−1 ≈ 0.18 of the central).
2. **One curve.** For Y-pad events, excess on n+1 / Y signal vs the Y ±1 share falls on one curve for all sessions and
   drift voltages (0.32 at Ysh 0.47 → 0.12 at 0.9). Date + drift-field dependence = cluster width (gas of the night via
   diffusion/drift velocity, and the field). Within a run per-event width matters little → width is a proxy for the gas state.
   At fixed 1000 V pulse widths are 1400 (6-22), 1140 (6-18), 840 (6-27), 600 ns (6-19): the drift velocity differed
   between nights. 6-19 (0.03) stays an outlier on the width curve (X share 0.84 → ~0.13 expected).
3. **Board axes.** Varies along the X-chip axis only (6-24: 0.55 → 0.15), flat along the Y-chip axis (= top-bottom with
   Dylan's orientation, the inlet→exhaust direction). Local pulse width does not follow it → not local gas.
4. **Gas circulation: not the cause** (flow speed, X-only, pad-level, flat top-bottom, non-monotonic in time). Gas
   composition as a modulator: yes.
5. **Beam drift scans** (condor 4502894, runs 58/61/64 90/10 etc., raw_asym.py --cm): det6 ≈ 0 at all fields; det7
   unmeasurable (noisy X FEU, Y fails); det2 0.03–0.12 (≈ 0 on bench: beam angle); det3 0.05–0.14 flat in 90/10.
   No bench drift scan of any board but det3.
6. **Simulations.** `bus_wpot.py` (now with parity fold, --eps-gap/--eps-kap), `bus_float.py` (stripes as floating
   conductors), `fwd_parity.py` (uniform clouds, data's estimator). Copper at any bus depth, wet gaps, kapton ε, a 60 µm
   displacement over Y pads: |a| < 0.02 in both classes. The old 0.010 was a kernel asymmetry that peak selection absorbs.
   Floating stripes offset ±200 µm: ±0.06–0.15 but equal in both classes and growing with σ: wrong fingerprint.
7. Electronics/cable excluded: they cannot know the pad type under the charge; copy linear.
Open: what at a Y pad couples ~30 % of its charge signal into the n+1 X line on Stephan's June bulks.

## Third round (2026-10-20): stripes, a selection caveat, and the bus side (deck slides "Stripes" .. "Inside", "Bus traces")
Scripts and JSON in `scratchpad/2026-10-20/`; event cache `~/x17/nb_asym_raw/esl_cache.npz` (desktop, derived; `esl_cache.py`).
1. **Resistive stripes: no mark.** With TRUE pad columns (col = 64·(ch//64) + 63 − ch%64; the old mod-40 check used DAQ
   channels, which mixes phases between connectors), Y-pad copy vs column mod 40: sine amp 0.007 ± 0.004 on 0.221; X-pad
   0.006 ± 0.009; Y share 0.007 ± 0.005; det3 May / det2 flat too. A floating-stripe coupling (fwd_out/f_off±200) would
   swing ±0.10 over the beat. So no mechanism tied to how the stripes sit over the pads/bus gap.
   Per-event "landed on stripe vs kapton gap" has no handle with cosmics: cloud σ ≈ 450 µm (fwd model matches Ysh 0.65)
   and the peak pad fixes x only to ±390 µm ≈ one 800 µm period → harmonic attenuated to ~0.002. Needs a point source.
2. **Caveat on the parity split** (`inject_parity.py`): injecting a plain channel copy x(n+1) += a·x(n) into copy-free
   boards and running the same analysis gives det2 + 0.3 → Y-pad 0.20 / X-pad 0.05 (det1 0.18/0.10, det3 May 0.16/0.09).
   Peak migration: the copy lifts n+1 above n, parity flips, the event lands in the X-pad class with a negative a1.
   The old control (symmetric fwd model → 0) did not test this. So the split does not by itself show pad-locality or
   exclude electronics (the noise test still excludes electronics). det4/det7 are consistent with a plain copy a ≈ 0.15–0.2.
3. **det3 June-15 is not a plain copy** (`uncopy_parity.py`, closure on injected controls is exact): undoing a = 0.33
   flattens the Y-pad class in every Y-share bin, but the narrow-cluster X-pad bins stay at −0.43/−0.32/−0.19 for any
   a in 0..0.35. Data X-pad by Ysh: −0.40 → +0.20; a plain copy gives ≈ 0 → +0.15. Two-copy fit (a_Y, a_X;
   `fit_two_copy.py`): a_Y 0.30–0.35, a_X 0.05–0.15 but χ² 230–380/21. Migration (`migrate.py`): 13.4 % of Ysh < 0.4
   clusters on det3 have the X peak on the X-pad parity (77 % with the true column at n−1); det2/det1 4–4.5 % at a ≤ 0.3,
   13 % only at a ≈ 0.8 (but det3 May baseline is 20 %, so not decisive alone).
4. **Forward model** (`fwd_variants.py`, solved weighting fields): an X-view sensitivity offset of ~200 µm toward n+1
   (relative to Y) reproduces the X-pad sign flip (−0.39 → +0.21); offset + copy 0.1 also lifts the Y-pad class. The
   model's pad-type purity vs Ysh is lower than the data's (0.53–0.65 vs 0.87 narrowest bin): shape only, no χ².
5. **Inside the Y pad** (`ypad_probe.py`): copy flat in η_y (y inside the Y pad, copy-free view): 0.21–0.25 everywhere,
   so not a via/stub/defect at one edge. Excess waveform peaks with the Y central (720 ns), before the X central (840) and
   the neighbours' own signal (900); fit 0.17·X_c + 0.07·Y_c leaves 7 %. Prompt: induced as the charge arrives.
6. **Bus side, checked:** column buses (DAQ X) = `L3-TrackY.gbr` = physical L5 (shallow), row buses on L6 → copy in X only
   ✓. Gerber: column buses fan out at the top, row buses on the right; detector (pad side): X connectors bottom, Y right →
   gerber = pad-side view mirrored in y, x kept. Stubs +0.39 mm in x → column c's bus in the gap toward c+1, beside the Y
   pads of column c+1 → Y-pad charge reaches the bus of the line one column lower = DAQ n+1 (connector reversal) ✓.
   Relies on ActiveAreaFrame/Mx17StripMap conventions. Copper solve: 0.010 nominal; det3 needs ~20×, det4/det7 ~6×.
Reading now: a one-sided reach of each X line toward its own bus, right side, right view, far too strong on the June
bulks; the date/field dependence is the cluster width acting on a fixed one-sided kernel (6-19 still an outlier).

## Readout question + per-session fits (2026-10-20, later; slides "Readout?", "Sessions")
1. **Not the one-sided readout (M1/connectors/cables/FEU).** All X lines are read at the same edge: not handed in x;
   passive neighbour coupling is reciprocal. Any readout effect = fixed linear channel mixing. Data:
   - profile across the board (Y-pad, 16-col bins, true columns) is smooth through all 7 connector edges: step rms at
     edges 0.026 vs 0.041 between bins inside connectors (`m1_profile.py`);
   - flat along the X line (by Y chip, distance from the readout edge): 0.21–0.27;
   - general mixing x'(n) = x(n) + Σ m_k x(n+k), k = ±1, ±2, either sign, injected into det2/det3 May/det1 and run through
     det3's analysis (`mix_fit.py`, `mix_reach2.py`; scored on class asym + ±1 sums by Ysh): best χ² 410–750/28, narrow
     X-pad bins −0.10…+0.23 vs data −0.40/−0.29. Only m₋₁ ≈ 1 gets there, and that wrecks the shares. Structural reason:
     det3 acts like an X–Y registration shift, which an X-only map can't make.
   - Not testable: DAQ-neighbour pairs across a connector edge (extractor drops chip-edge peaks). det4/det7 (plain-copy-like)
     rest on the profile + noise test only.
2. **Per-session forward fits** (`fwd_session_fit.py`, Y-share sextiles, 12 class asymmetries, 0.02 floor): offset + copy
   beats pure copy everywhere (χ² 44 vs 211; 214 vs 535; 8 vs 59; 87 vs 196; all 59 vs 365). Copy 0.10–0.15 in every
   session; offset −250 (6-23/24, 600 V; 6-27 drift ≤ 500 V), −200 (6-22/23), −100 µm (6-27/28 1000 V): it tracks the
   session's median Y share (0.52 → 1.01). With quantile bins the absolute width is not used, so this is the old width
   curve reappearing as an offset; deciding "one fixed geometry" needs a model with the data's Y-share scale (model
   medians 0.37–0.49: no resistive spreading / shaping in the static fields).

## Connector edges, boards, X-offset sim (2026-10-20 evening; slides "Connector edge", "Boards")
1. **Connector-edge test** (`raw_asym.py --edge` on yoga → `~/x17/nb_asym_raw/edge/`, `edge_test.py`): peaks at a
   connector's lowest column (col % 64 == 0), whose physical col−1 is on the other connector/cable/DREAM chip.
   det3 June-15 6-23 (FEU 3): physical α(col−1) − α(col+1) = +0.206 ± 0.025 there vs +0.21…+0.23 inside; DAQ n+1 at that
   boundary (next chip, 127 columns away) α = 0.00. The copy follows the pads. 6-28 (FEU 7): +0.03 ± 0.03 vs +0.09 inside
   (2σ, not decisive). Controls at col%64==0: det2 −0.01, det3 May +0.02, det4 +0.07. At col%64==63 every board reads
   +0.04…+0.14 (cross-connector neighbour reads low on all boards): use the ==0 case.
2. **run_config.json** descriptions (all local + EOS configs): May det3 "Bulked 5-5-26. Strip resist vertical. ESL no
   silver paste"; det1 April "Bulked 4-24-26. Strip resist horizontal … Resist strips broken"; beam mx17_A = det3 "Bulked
   by Stephan June 15", mx17_B = det2 "Bulked by Arnaud June 12. Giant pillars on parts of the detector", mx17_C = det6
   (Stephan June 24, board D), mx17_D = det7 (Stephan June 22 batch, board B); det4 = June 22 batch, board C. Stephan
   (via Dylan, 10-20): nothing different about the June bulks.
3. **Stripe direction from data** (`stripe_orient.py`, ±2 late content): Y > X on every board except det1 (Jan build,
   reversed); June boards' Y ±2 peaks 240–420 ns late. det3 June-15 looks like det2/det6/det7: orientation does not
   single it out.
4. **PCB fingerprint** (`pcb_fingerprint.py`): det3 June-15 Y row 23 dead on both FEU 4 and FEU 8 (follows the board);
   X 112 dead on all FEU-3 boards (electronics). det3 May has a different set (Y 52/182/299, X clusters 151–163, 335–357,
   462–507) and row 23 alive: no shared board-level defect. Same PCB not supported (bulk-level defects can't be excluded).
5. **Full-sim X-view offset**: digitizer `--x-view-offset-um` (X band reads the avalanche at x + D). Pilot condor
   4505153 (v5a recipe, th00, D = 0/−100/−200/−300): D = 0 is bit-identical to v5a_th00 (raw windows equal). Only ~190
   matched isolated events per arm → too few. Statistics run condor 4505159: 4 offsets × 12 angle files × 2 seeds
   (`run_xoff_job.sh`, decoded ROOT to EOS `response_sim/stageB_w2/xoff/`, npz to AFS `simxoff_npz/`); evaluate with
   `sim_xoff_eval.py` against 6-27/28 1000 V (the v5a target) and all det3 June-15.

## Fourth round (2026-10-20, dylan-Yoga): X-offset full sim read out; stripe direction on every pre-June board
1. **X-offset full sim (condor 4505159, 96/96 done)**, `sim_xoff_eval.py` (now skips the empty steep-angle arms,
   env `XOFF_ARMS`, `XOFF_COPIES`, `XOFF_SESSIONS`, `XOFF_OUT`) + `sim_xoff_rescore.py` (adds the Y-pad fraction fY per
   Y-share bin to the χ²; the eval's χ² drops bins where the sim has no events of a class, so arms with no narrow X-pad
   events never paid for the data's −0.40). Arms restricted to th00 + thy* (the thx/txy arms are one-signed in x and
   give a −0.04 baseline a1 at zero offset). Absolute Y-share bins 0/.4/.55/.7/.85/1/1.3. Outputs
   `sim_xoff_eval_{sym,sess}.json`, `sim_xoff_rescore_sess.txt`.
   - **Pre-registered test passed:** narrow X-pad class at −300 µm = −0.49/−0.45 (falsifier was ≥ −0.1). At −200:
     −0.51…−0.60/−0.44. At 0 µm there are no narrow X-pad events at any copy.
   - **The offset is needed:** it is the only ingredient that makes the narrow X-pad flip and the pad-type migration.
     Narrow-bin Y-pad fraction: data 0.87/0.74/0.67; offset −200: 0.91/0.78/0.77; plain copy 0.94/0.90/0.89.
     fY χ² 16–27 with the offset vs 80–163 for any plain copy.
   - **A rigid offset is not enough:** its asymmetry fades as clusters widen (−200 µm, Y-pad class +0.27 → −0.02 over
     the Y-share bins), but the data's does not. On 6-23/24 the Y-pad class is +0.29 … +0.36 flat out to the widest
     bin, and the X-pad class goes +0.13 … +0.28 there. That part is copy-like (a fixed fraction).
   - Adding a copy on top overshoots the narrow Y-pad bins (+0.41 vs +0.26). No (offset, copy) on the grid fits any
     session. Best totals: 6-22/23 −100 µm + 0.20 (χ² 74, ~20 bins); 6-23/24 −100 + 0.20 (195); 6-27/28 plain 0.25 ≈
     −100 + 0.10 (36–40).
   - **Reading:** det3 behaves like a shift at narrow widths and a copy at wide widths. Both fit a one-sided
     pad-level coupling: X line n+1 senses a fixed fraction w of the charge over the Y pads on its bus side. That
     looks like a shift when the cloud sits on one Y pad (n+1 can overtake n, so the peak migrates) and like a copy
     for wide clouds. Next model: add that cross-kernel to the digitizer (X(n+1) += w × induced signal of bus-side Y
     pads), not a rigid offset.
2. **Stripe direction, every board** (`stripe_orient2.py` → `stripe_orient2.{txt,json}`). Charge spreads along the
   stripes, so the slow view shows the direction. Metric: ±2 late content (600–960 ns after the central peak),
   bootstrap errors. New extractions on yoga `/media/dylan/data/claude_tmp/tmp/nb/raw/orient/`:
   - det0 4-1 He/ethane, 505/520 V.
   - det4 May build ("Bulked 5-6-26", 5-7 Ar/CF4).
   - The April-24 bulk ("Strip resist horizontal. Waves on mesh. Resist strips broken", EOS
     `cosmic_bench/det_1_old/mx17_det1_Ar_CF4_HV_Scan_4-25-26`, 550 V only). It has no decoded_root, so it was read
     from the M3-filtered arrays with the new `raw_asym.py --array`.

   | build | n X/Y | ±2 late X / Y (×10⁻³) | slow view |
   |---|---|---|---|
   | det0 (4-1-26, He/ethane) | 4092/4335 | 23.6 / −0.8 | **X (−50 σ)** |
   | det1 (Jan) | 5455/5066 | 37.0 / 9.4 | **X (−41 σ)** |
   | Apr-24 bulk "horizontal" | 256/825 (60 ADC floor) | 4.0±3.3 / 6.2±2.0 | none (+0.6 σ) |
   | det3 May "vertical" | 15123/16913 | 8.5 / 31.0 | Y (+75 σ) |
   | det4 May (5-6 bulk) | 23069/23059 | 0.8 / 31.9 | Y (+183 σ) |
   | det2, det3 June-15, det4, det6, det7 | | | Y (+9 … +130 σ) |

   - "Vertical" (det3 May) = Y-slow. det0 (April 1) and det1 (January) are X-slow: their stripes run the other way.
     Both sit under EOS `det_0_old`, so they are probably one PCB (not proven). Every May and June build runs the same way, so none of det2/3/4/6/7 is reversed.
   - The April-24 re-bulk shows no slow spreading in either view, which fits its "resist strips broken" note. Its
     cosmics are faint at 550 V (most below 200 ADC; σ_noise ≈ 10 ADC). It also has a fixed pulsing spot at X 487–494 /
     Y 48–52, which the hot-channel veto removes. Its direction can't be read from data; the config says horizontal.
   - None of this touches the det3 X asymmetry: det1 (reversed) has copy ≈ 0, and det3 June-15's direction is the
     common one.

## Fifth round (2026-10-20 day, dylan-Yoga): bookkeeping, stripe pitch, timing, mixed responses
Scripts in `scratchpad/2026-10-20/r5/` (`r5lib.py` = matched X/Y events per build from the raw survey).
1. **Charge bookkeeping** (`probe_charge.py`, `probe_charge_bins.py`, `charge_sessions.py`). Integrals of positive
   samples, per Y-share bin, Y-pad events, relative to ΣY:
   - Against det2 at matched Y share, det3 June-15's n+1 gains about what its n−1 loses (Ysh 0.55–1.0:
     +0.061/−0.052, +0.067/−0.056, +0.048/−0.041), and Xc and the ±1 sum are unchanged.
   - An injected plain copy raises n−1 and Xc too, because each channel gets its lower neighbour's share.
   - ΣX/ΣY is not raised on det3 (Y-pad 0.907, X-pad 0.925; det2 0.932/0.944) and doesn't rise with the copy
     between sessions on the same FEUs (6-22/23 vs 6-23/24).
   - **Reading:** the excess is charge moved inside X, not added, and not taken from Y. That excludes any Y→X
     transfer: bus pickup, a floating or cracked-via Y pad, and Y-line crosstalk.
   - Caveat: the gas changes the X-vs-Y width relation between sessions (Xc at fixed Y share 0.38 → 0.50), so matched
     Y share is not a perfect control.
2. **Stripe pitch from the data** (`periodogram.py`, `beat_sig.py`, `beat_col.py`). Per-event statistics overstate the
   significance. Because of per-channel gain and sharing offsets the noise unit is the column, so the test uses
   per-column means, detrended at σ = 30 columns, against a column-permutation null.
   - **No board shows a stripe beat in X** at 1/40 or anywhere in 0.015–0.04 (det3 June-15: s1 p = 0.12, d1 p = 0.86).
     Y-view controls look the same.
   - The asymmetry has a 5-column (and 2.5) line on several boards (det2 too): common to the readout, not the stripes.
   - So the stripes leave no measurable imprint at pad scale. A stripe mechanism strong enough to make 0.2 on det3
     would make a large beat on every 800 µm board, unless the det3 June-15 print alone had a pitch locked to 780 µm.
3. **Timing of the excess** (`timing.py`, `run_timing.py`). Y-pad events, aligned on the Y central:

   | sample | excess peak | excess t50 | Xc peak / t50 |
   |---|---|---|---|
   | 6-22/23 | +21 ns | −151 | +18 / −181 |
   | 6-23/24 | +25 | −130 | 0 / −165 |
   | 6-27/28 | +3 | −165 | +22 / −199 |
   | plain copy (det2 + 0.25) | +8 | −168 | +10 / −175 |
   | sim offset −200 / −300 µm | +120 / +88 | −121 / −124 | +21 / −265 |

   The data's excess peaks with the central (a direct-like signal), not late like the neighbour charge that a
   ≤ 300 µm offset moves. Its rise is ~30 ns behind the central's.
4. **Mixed responses**: X = (1 − w)·X(D = 0) + w·X(D). The condor 4505159 arms share seeds: Y is bit-identical and X
   differs only in signal, so the mix is exact. Closure is bit-exact (w = 0 → xo0, w = 1 → xom200).
   - `raw_asym.py --mix-with/--mix-w`, `run_mix.sh` (D = −100/−200/−300 × w 0.2–0.8), `eval_mix.py` + `fingerprint.py`
     (one score: A_Y, A_X absolute; fY and the m1/p1 bookkeeping as effects, det3 − det2 vs model − sim D0).
   - The effects form is needed because the sim's X central fraction is flat in Y share (~0.5) while every board's
     falls 0.5 → 0.32: a baseline sim/data difference unrelated to the asymmetry.
   - Condor **4507004**: D = −450/−600/−780 µm, th00 + thy*, seeds 5/6 (36 jobs), to reach a one-pitch mix.
5. **x-position conditional** (`xpos_probe.py`). For Y-pad events, E[α(n+1) | α(n−1)] at matched Y share:
   - Copy-free boards: it falls with α(n−1) (0.37 → 0.24), as geometry requires.
   - det3: it stays flat (~0.45–0.5). The excess is ≈ 0 when the charge sits on the n+1 side and grows to ≈ 0.2
     on the n−1 side. A plain copy is flat (+0.24).
   - So n+1 takes a share of what n−1 should get.
   - The full sim has too few events per (Y share, u) bin to trace this (`eval_xpos.py`).
6. **Mixed responses in the full sim** (`eval_mix.py`; D −100/−200/−300 × w 0.2–0.8). Every mixture gives a1 ≈ 0.00–0.03
   and a coin flip ≈ 0.5; data a1 0.14, coin flip 0.63–0.81. A position-level effect alone can't make the X marginal.
   The sim's narrow X-pad bins are empty, so the class test is weak there. On the bins it has, copy 0.2 scores best
   on 6-22/23, and copy 0.1 or mixtures on 6-27/28.
7. **Fast static model** (`fwd_xpos.py`, `fwd_region.py`). Every candidate is scored on all fingerprints at once
   (Y-share quantiles: the static model's sharing scale differs).
   - Extra sensitivity regions on the bus side (gap strip, Y pad of the next column, its near half, the next or
     skip X line, the whole next column), each with its donor, all fail. They send the narrow X-pad asymmetry
     positive and leave A_Y small.
   - Only displacement-type models (shift, or mix with D ≈ −400…−600) reproduce A_Y and the narrow A_X. Best χ²
     ≈ 240 on 24 bins, so no static model fits.
   - The Y-pad fraction is ~0.10 low on det3 at every width (0.551 vs 0.63); an injected copy doesn't change it.

### The asymmetry is a smooth, two-signed, X-only field over every board (2026-10-20 afternoon)
8. **Across-board profiles** (`chip_profile.py`, `block_profile.py`, true column order):
   - det4 changes sign: −0.12/−0.16 in chips 0–1, +0.13 in chips 4–5, smooth through chips 2–4.
   - det6 changes sign too: +0.04 → −0.04.
   - det7 rises from 0.01 to 0.12.
   - det3 is a smooth bump: 0.09 → 0.18 (chip 3) → 0.07.
   - No step at any connector boundary (det3: every step ≤ 0.02 ± 0.03), so not wiring.
   - A fixed-layout coupling (bus, stubs, vias) has one sign everywhere: **excluded**.
9. **2-D maps** (`map2d.py`, `map_fit.py`, `map_y.py`; 8 × 8 cells, matched X column × Y row). Constant-model χ²:

   | board | X | Y |
   |---|---|---|
   | det3 | 332/63 | 74/63 |
   | det2 | **736/61** | 59/61 |
   | det6 | 47/25 | 37/25 |
   | det7 | 21/8 | 3.5/8 |
   | det3 May | 78/56 | 76/56 |
   | det1 | 34/27 | 46/27 |

   - **The Y view has no structure anywhere: the field is X-only.** That excludes vector misregistrations (via
     drill offsets, inner-layer scale/rotation); only something running along y and acting across x can do it.
   - det3: a tilted ridge, max at x-quarter 1 at the bottom and x-quarter 2 at the top. Its two FEU-3 sessions give
     the same geometry and phase (−123° / −153°), with amplitude 0.048 vs 0.124.
   - **det2, the old "copy-free control", carries a local asymmetry as large as det3's**:
     - By y quarter: −0.06 / +0.14 / +0.01 / −0.06; cells up to ±0.18 at 8σ.
     - It is in both det2 runs on two FEU slots (6-17 FEU 3: −0.087/+0.166/+0.010/−0.078; 6-22 FEU 6:
       −0.008/+0.079/−0.003/−0.033).
     - It averages to zero over the board.
   - det1 and det3 May are quiet.
   - Track angle can't do it: the tagged dumps' M3 angles are < 1° per cell, and an x-tilt can only follow x, while
     det2's field follows y.
10. **det2's positive band is a little det3** (`det2_regions.py`, rows 128–255):
    - a1 +0.142, coin flip 0.714.
    - A_Y +0.16…+0.27; narrow → wide A_X −0.18 → +0.12; n−1 depleted, n+1 enriched.
    - Excess timing +27 / −126 ns (det3 +20 / −145).
    - Its rows 0–127 and 384–511 are the mirror image (a1 −0.058, coin flip 0.395, A_Y −0.10…−0.17, narrow A_X
      +0.31); the middle rows are flat.
    - **The det3 mechanism is generic; det3 June-15 is a board where its sign is positive over most of the area.**
11. **Pitch-locked.** Inside det2's positive band, a1 is positive in every 8-column block across all 512 columns,
    with no 1/40 harmonic (p = 0.57).
    - An 800 µm stripe pattern would flip it every 40 columns (31.2 mm), and det2's band means would average to
      zero. So the structure repeats at the pad pitch (780 µm) to ≲ 0.1 %.
    - No board shows the beat the sim predicts (sim s1 0.12, d1 0.07 at 1/40; data ≤ 0.03 / ≤ 0.011,
      `beat_col.py`, `beat_sim.py`).
12. **Working hypothesis: resistive stripes printed at the pad pitch, with per-board registration.**
    - Charge evens out across a stripe's width (fast), which pulls narrow clouds toward the stripe centre. If the
      centres sit at an offset δ(x, y) from the pad columns, the X view sees n+1 ≠ n−1 ∝ sin(2πδ/P).
    - It explains: X only; charge moved, not added; prompt; larger for narrow clouds (gas, drift field); a smooth
      two-signed field from the print's scale and rotation errors (det3 ~0.07 % scale + ~0.25 mrad; det2 ~2 mrad
      rotation); and the absence of the 31.2 mm beat on every board.
    - Against it so far:
      - The stripe pitch was given as 800 µm (memory, 2026-08-06, "borrowed from another detector").
      - The sim's stripe pull at its strongest phase is ±0.065 vs data ±0.17–0.2.
      - The predicted X-sharing modulation (∝ cos φ) isn't detected (`quad_test.py`: det3 −0.007 ± 0.009); only
        an in-phase anti-correlation (−0.031 ± 0.005).
    - Pending: the sim's fingerprint per stripe phase (condor 4507242, 60 jobs).
13. **Same bands in both det2 runs** (6-17 FEU 3, 6-22 FEU 6), by 64-row band:
    - 6-17: −0.29 / +0.02 / +0.18 / +0.15 / +0.03 / −0.01 / −0.05 / −0.13
    - 6-22: −0.07 / +0.02 / +0.07 / +0.09 / +0.03 / −0.01 / −0.02 / −0.08
    - The geometry is fixed. The amplitude is ~2.5× larger in the run with narrower clouds (median Y share 0.35 vs
      1.07). det3's session-to-session swings fit the same picture: fixed geometry, width-set amplitude.
14. **Sim, phase-resolved (condor 4507242, 60 jobs, `sim_phase_fp.py`): the stripe hypothesis is NOT reproduced by
    the sim's stripe physics.**
    - Ranking phases on seeds 5–10 and measuring on held-out seeds 11–16, the "most positive" phases give
      a1 +0.016 (ranked +0.11…+0.17: mostly noise); the negative set −0.057. So the sim's stripe pull is ≲ 0.03–0.05,
      against ±0.2 in data, and its class pattern isn't det3's.
    - The sim never fills phases 12–19 (the gun covers ~50 columns).
    - Either the real stripes act more strongly or differently than the S1 kernels (transverse equalisation, charge
      on the kapton gaps, a different print), or the cause is another pad-pitch structure running along y.
15. Pure large offsets (condor 4507004, −450/−600/−780 µm) swap the pad-class labels (−450, −600) or look like no
    offset (−780). They add nothing on their own; the decoded-level mixtures were not built, since the stripe
    picture replaced them.

**State at the end of 2026-10-20 (dylan-Yoga):**
- Established:
  - The det3 X asymmetry is one instance of an X-only, smooth, two-signed field that every board with a resistive
    print carries in some form (det2 as strongly as det3; det1 and det3 May quietly).
  - Its geometry is fixed per board, its amplitude is set by the cloud width (gas, drift field), it repeats at the
    pad pitch along x, it is charge moved within X, and it is prompt.
  - Excluded: copper, bus and vias, readout/connectors, Y→X transfer, track angle, and an 800 µm beat.
- Working hypothesis: the resistive-stripe print, at the pad pitch, registered differently on each board.
- Not yet shown: a model that produces ±0.2 (the sim's stripes give ≲ 0.05).
- Question for Dylan/Stephan: the stripe screen's actual pitch, and how it's aligned to the board.
- Figure: `design/figures/asym_maps_2026-10-20.png` (`r5/fig_maps.py`).

## Sixth round (2026-10-20 evening, dylan-Yoga): can stripe registration make it; the reversed-stripe boards
Scripts in `scratchpad/2026-10-20/r5/`.
1. **Sim asymmetry vs the true stripe offset** (`sim_delta.py`, `sim_delta_fp.py`; 36k X-only and 3.7k matched events,
   xo0 th00 + thy*, seeds 5–16). δ is the distance from the peak column's centre (PAD_ORIGIN + c·780) to the nearest
   stripe centre (stripes at k·800, `solver/kernels.py`). It changes 20 µm per column, so locally it acts as a
   pitch-locked offset.
   - a1(δ) is odd and monotonic: −0.16 (−325 µm, next to a gap) … −0.05 (−225) … 0 … +0.04 (+75) … +0.085 (+225).
     The narrow-cloud half is ~2× larger (±0.10, −0.24 at −325). X sharing is lowest at δ ≈ 0 (0.63) and highest
     near ±300 µm (0.85–1.1).
   - At δ = +150…+350 the sim reproduces det3's class pattern:
     - narrow → wide A_X: −0.19 → +0.12/+0.13 (det3 −0.12 → +0.20);
     - n−1 below n+1 for narrow Y-pad events;
     - δ = −350…−150 gives the mirror (A_X +0.40 → −0.17, like det2's negative bands).
   - **Stripe registration does generate this asymmetry, with the right sign structure and width dependence, at
     ~1/3 of the strength:** A_Y +0.02–0.07 vs det3 +0.25; matched a1 ±0.035 vs ±0.14.
   - The earlier held-out phase test (`sim_phase_fp.py`) failed because ranking noisy phase means on a few seeds
     selects noise.
   - A rotation alone can't make it (a 2 mrad tilt moves charge sideways by µm). As registration, a rotation of a
     780 µm print gives moiré bands ~P/θ ≈ 390 mm apart, running nearly perpendicular to the stripes, i.e. along x:
     det2's horizontal bands. With an 800 µm pitch a rotation only shortens the 31 mm moiré.
2. **Reversed-stripe boards** (`orient_maps.py`; det0 4-1 He all 17 sub-runs, det0 4-7 Ar/iso via `--array` on lxplus
   (9 sub-runs), det1 Jan all 5 files, det4 May build). Field strength = noise-subtracted cell rms:

   | board | stripes | across-stripe view field | X χ²/ndf (8×8) | Y χ²/ndf (8×8) |
   |---|---|---|---|---|
   | det0 He (10.8k) | along x | Y 0.004 | 50/61 | 49/61 |
   | det0 Ar/iso (9.4k) | along x | Y ~0.016 | 80/61 | 72/61 |
   | det1 all files (3.0k) | along x | Y 0.036 | 53/50 | **103/50** (p 1.6e-5; 4×4 p 0.11) |
   | det4 May (16k) | along y | X 0.015 | 76/58 | 77/58 |
   | det3 May (4.6k) | along y | X 0.020 | 78/56 | 76/56 |
   | June boards | along y | X 0.05–0.15 | | flat |

   - **These boards can't pin the field to the stripes.** The reversed boards show at most a weak, fine-grained Y
     structure (det1), but the May boards with normal orientation show no X field either. The field belongs to the
     June builds.
   - 1D profiles along each view's own coordinate are dominated by single-block channel defects, so they are not
     usable for this.
3. **What distinguishes the June builds** (stripe_orient metric, late ±2 content ×10⁻³):

   | | along the stripes (Y) | across the stripes (X) |
   |---|---|---|
   | June builds | 63–80 | 12–31 |
   | May / earlier | 24–37 | 1–9 |
   | sim (v5a) | 64 (+360 ns), June-like | −2: no cross-stripe leak at all |

   - The sim's stripes are June-like along the stripes, so stripe conductivity is not why its pull is 1/3.
   - The June boards also leak charge late ACROSS the kapton gaps between stripes, which the sim lacks. That's a
     concrete candidate for the missing physics (gap surface conduction, plausibly humidity-dependent with ~1 % H2O
     in the gas). It isn't yet shown to produce asymmetry.
4. Next:
   - A digitizer/kernel variant with a cross-gap leak (and/or electrons funnelled from gaps onto stripes), checked on
     the sim's a1(δ) amplitude.
   - A per-board a1(x, y) calibration map for analysis.
   - The stripe mask from Dylan (pitch, alignment).

## Seventh round (2026-10-20 afternoon, dylan-Yoga): cross-gap leak and funnelling, in the solver itself
Scripts in `scratchpad/2026-10-20/r5/` (`esl_fast.py`, `kernel_delta.py`, `delta_timing.py`), condor files in
`scratchpad/2026-10-21/`.
1. **Static fast model** (`esl_fast.py`, deterministic clouds on the prompt X-line field at the ESL plane): moving strip
   charge across or along its stripe, funnelling gap charge, or leaking 10–30 % to neighbour stripes all give
   |a1| ≤ 0.02 at every offset. The one large case is all charge on stripes and spread ~1.5 mm along them (±0.10), but
   the same model gives 0 for the sim's own case (insulating gaps), where the full sim gives ±0.05. So a static
   picture can't stand in for the kernel; it only says that charge locked to the stripes carries the registration and
   charge sitting in the gaps dilutes it.
2. **The production kernel itself** (`kernel_delta.py` on a 9-time slice of G_X from
   `greens_comb_w2_rho5M_dk50um_g19um.npz`, sliced on lxplus; slice in `/eos/user/d/dneff/tmp/gx/`): with y-resolved
   gaussian clouds and the peak-column selection, A = (G[c+1] − G[c−1])/G[c] vs δ is zero at t = 0 and builds up over
   100–300 ns. It is ±0.04–0.06 for σ = 250 µm next to the strip edges (|δ| 270–330) and ±0.015 for σ = 450 µm. The
   registration asymmetry in the sim is a sheet-transport effect.
3. **Timing** (`delta_timing.py`, same estimator on sim and data, X-only, aligned on the central peak):

   | sample | n | excess amp | excess peak / t50 vs central (ns) |
   |---|---|---|---|
   | sim \|δ\| ≥ 150 (sign-folded) | 21.9k | 0.069 | −71 / −90 |
   | data det3 June-15 X | 42.8k | 0.142 | −10 / −25 |
   | data det2 X, det3 May X | | ≤ 0.005 | — |

   The sim's registration excess leads the central by about one sample; the data's has nearly the central's shape.
   So registration as the sim models it is not prompt enough. That isn't decisive (one sample, and the sim's pulses
   rise ~70 ns slower overall), but a stronger mechanism has to be prompter, not later.
4. **Conducting gaps in the W2 solver.** New, default-off: `wpot.esl_sigma_profile(gap_frac=)`,
   `W2Solver(esl_gap_frac=)`, `w2_production.py family --nx --gap-frac` and `combine-x` (X-view products only).
   Production is bit-identical at the defaults. A local mini pipeline at nx = 390 (80 µm; the pad edges fall off-grid,
   so per-δ values carry ±0.03 lattice artefacts, but f = 0 and f > 0 share them):

   | 5 MΩ/sq, σ 250 µm | odd a1 at 300 ns | at 1 µs | late ±2 across stripes (×10⁻³, 1 µs) |
   |---|---|---|---|
   | f = 0 | −0.029 | −0.033 | 1.0 |
   | f = 0.03 | −0.029 | −0.043 | 9.8 (σ 450: 13 → 34) |
   | f = 0.03 at 0.5 MΩ/sq | −0.047 | −0.021 | 252 |

   A leaky gap reproduces the June boards' late across-stripe content. It does not make the prompt asymmetry 3×
   stronger: what it adds comes late, which is the wrong timing for the data.
   Running: condor 4520685, nx = 1560 (20 µm, exact pad edges), f = 0, 0.01, 0.03, 0.1, 0.3, one job per f (both
   X families + combine-x); products → `/eos/user/d/dneff/tmp/gx/gap<f>_nx1560/gx_w2_rho*M_*.npz`. Score with
   `python3 r5/kernel_delta.py <npz>…`.
5. Reading: the candidate for the missing strength has to be prompt and stripe-locked. Prompt electron funnelling onto
   the stripes (the gaps' kapton surface charged up, so arriving electrons are steered onto the stripes) fits that
   better than a slow leak. The static model says funnelling only pays once the charge also spreads along the stripe.
   That is the next thing to put in the digitizer: move each avalanche that lands in a gap onto the nearer stripe
   edge before the kernel lookup (the existing kernel, no re-solve).
6. **Analysis correction: X-centroid bias maps** (`r5/corr_maps.py` → `r5/corr_maps.json`; next step 2). The 3-strip X
   centroid is pulled by dx = −(α₊₁ − α₋₁)/(1 + α₊₁ + α₋₁)·P along the true column axis. Per 8×8 cell (64 columns ×
   64 rows), minus the symmetric reference (det3 May + det1, −0.7 µm):

   | build | events | cells | mean dx (µm) | cell rms (µm) | χ²(dx = 0) |
   |---|---|---|---|---|---|
   | det3 June-15 | 14.6k | 64 | −74.5 | 76.6 | 3997/64 |
   | det2 | 5.6k | 62 | −1.6 | 68.7 | 1090/62 |
   | det4 | 0.9k | 13 | −31.7 | 43.5 | 46/13 |
   | det6 | 1.4k | 26 | −7.4 | 31.6 | 56/26 |
   | det7 | 0.7k | 9 | −23.1 | 48.2 | 27/9 |
   | det3 May (control) | 4.6k | 57 | −0.1 | 22.5 | 83/57 |
   | det1 (control) | 1.6k | 28 | +4.2 | 20.8 | 34/28 |

   det3 sessions as a scale on the all-session map: 6-22/23 1000 V 0.96 (χ² 47/60), 6-23/24 600 V 1.57 (126/61),
   6-27/28 1000 V 0.44 (47/59). So the correction is one fixed map per board times a per-run amplitude (from the
   run's own board-mean a1, or its cloud width). det2's bias averages away but is ±70 µm locally. det4/6/7 need more
   statistics than the bench has for a cell map.
7. **Gap funnelling in the full sim** (condor 4520688: v5a recipe + `--gap-funnel` 1 and 0.5, th00 + thy05…30,
   seeds 7–12, npz in yoga `~/x17/nb_asym_raw/simfunnel/`, decoded ROOT on EOS `response_sim/stageB_w2/funnel/`).
   The funnel has its own rng stream, so each arm is paired event-for-event with the D = 0 arm of the same seed.
   Odd amplitude O = ⟨sgn δ · a1⟩ over |δ| ≥ 50 µm (`paired_odd.py`, 36 pairs each):

   | arm | O all | O narrow half | X sharing |
   |---|---|---|---|
   | D = 0 (4507242) | +0.056 | +0.071 | 0.760 |
   | funnel 0.5 | +0.067 (×1.20, +0.011 ± 0.003) | +0.086 (×1.21) | 0.766 |
   | funnel 1.0 | +0.078 (×1.40, +0.022 ± 0.004) | +0.108 (×1.54, +0.038 ± 0.003) | 0.765 |

   Timing (`delta_timing.py fu1 …`): excess 0.109 (|δ| ≥ 150), peak −91 ns, t50 −131 ns vs the central; the
   unfunnelled sim −71/−90, the data −10/−25. Funnelling strengthens the registration asymmetry by up to ~1.5× and
   makes it no prompter, so funnelling alone doesn't close the 2.5–5× gap or the timing.
8. **ρ_s** (mini kernels, insulating gaps): the registration term saturates at the same value at every ρ_s (odd
   −0.033 for σ 250 µm at 0.5…5 MΩ/sq). Lower ρ_s only gets there sooner (by 92 ns at 0.5 MΩ, by ~500 ns at 5 MΩ).
   So a lower-resistivity June print would make the registration excess prompt without making it larger.
9. Strip-width scan: see item 11.
10. **Gap-conductivity scan at nx = 1560** (condor 4520685, products `/eos/user/d/dneff/tmp/gx/gap<f>_nx1560/`,
    `r5/kernel_delta.py` → `r5/kernel_delta_gap1560.json`, log `r5/kd1560.log`). The f = 0 grid control reproduces the
    production kernel (−0.024 at 300 ns, same as the 10 µm slice). Odd amplitude ⟨sgn δ · A⟩, σ = 250 µm:

    | | 300 ns | 1 µs | late ±2 at 1 µs (×10⁻³, kernel metric) |
    |---|---|---|---|
    | f = 0 | −0.024 | −0.024 | 0.7 |
    | f = 0.01 | −0.026 | −0.033 | 2.3 |
    | f = 0.03 | −0.028 | −0.042 | 7.5 |
    | f = 0.1 | −0.030 | −0.046 | 38 |
    | f = 0.3 | −0.023 | −0.021 | 129 |
    | 0.5 MΩ, f = 0 | −0.024 | −0.025 | 0.7 |
    | 0.5 MΩ, f = 0.03 | −0.053 | −0.027 | 230 |
    | 0.5 MΩ, f = 0.3 | −0.003 | −0.001 | 777 |

    For σ = 450 µm the odd term stays ≤ 0.004 throughout. A moderate leak roughly doubles the registration asymmetry
    for narrow clouds: it drains the gap charge onto the neighbouring strips, so more of the charge ends up
    registered. Too much leak makes the sheet uniform and the asymmetry goes. At 5 MΩ the doubling arrives late; at
    0.5 MΩ it is there by 300 ns. Doubling from the leak and ~1.4× from funnelling together reach the 2.5–5× the data
    need. This is the first combination that gets there in size.
    Full-sim test running: condor 4520810 (144 jobs). Hybrid kernels (`2026-10-21/build_hybrid.py`: production Y
    kernels + the nx = 1560 X kernels upsampled to 10 µm) for f = 0 (grid control), 0.03, 0.1, and 0.1 + `--gap-funnel
    1`, all at 5 MΩ, th00 + thy*, seeds 7–12, npz → `simgapleak_npz`. Score: `paired_odd.py glg0 … glg0.03 … glg0.1 …
    glg0.1f1 …` and `delta_timing.py glg0.1 …`; also the late ±2 across-stripe content with `stripe_orient` metric
    against data (June 12–31, May 1–9). The Y view keeps insulating-gap kernels, so this tests X only.
11. **Strip width** (condor 4520736, insulating gaps, products `/eos/user/d/dneff/tmp/gx/gap0_w<w>_nx1560/`,
    `r5/kd_width.log`, `r5/kernel_delta_width.json`). Odd amplitude, σ = 250 µm, 5 MΩ/sq:

    | strip width (gap) | 300 ns | 1 µs | late ±2 (σ 450, ×10⁻³) |
    |---|---|---|---|
    | 420 µm (380) | −0.016 | −0.014 | 4.9 |
    | 540 µm (260), production | −0.024 | −0.024 | 10.3 |
    | 700 µm (100) | −0.037 | −0.056 | 25.7 |

    Wider stripes (narrower gaps) raise the registration asymmetry strongly, and the late ±2 grows with them. That is
    the opposite of "more gap charge, more asymmetry". It's the registered (strip) charge that carries it, so a fuller
    print registers more of the cloud. If the June print's stripes are wider than 550 µm, that alone gives up to
    ~2.4× at 5 MΩ. **For the mask:** read the stripe width as well as the pitch and alignment. A June-vs-May
    difference in width or fill is the cleanest candidate yet for "why only the June builds".
12. **Full sim with leaky gaps** (condor 4520810: hybrid kernels, production Y + nx = 1560 X upsampled; npz in yoga
    `~/x17/nb_asym_raw/simgapleak/`, decoded ROOT on EOS `response_sim/stageB_w2/gapleak/`). Paired on the same 36
    (angle, seed) sets (`paired_odd.py` → `paired_odd_gapleak.json`; `sim_late2.py`):

    | arm | O all (vs grid control) | O narrow | X sharing | X late ±2 (×10⁻³) |
    |---|---|---|---|---|
    | glg0: f = 0 grid control (vs xo0: ×0.95) | +0.053 | +0.068 | 0.760 | −2.3 |
    | glg0.03 | ×1.00 | ×1.02 | 0.784 | +7.8 |
    | glg0.1 | ×0.90 (−0.005 ± 0.002) | ×1.04 | 0.828 | +27.6 |
    | glg0.1f1 (+ funnel 1) | ×1.29 | ×1.49 | 0.844 | +27.9 |
    | data: June builds / May and earlier | | | | 12–31 / 1–9 |

    - A gap conductance of 5–10 % of the stripes' reproduces the June boards' late across-stripe content with the
      data's own metric. That closes the "June boards leak across the stripes, the sim doesn't" item (#29).
    - In the full sim it adds nothing to the asymmetry. The kernel's ×1.9 is late, and the in-phase peak estimator
      doesn't see it; the extra sharing it brings dilutes a1 instead.
    - With leak plus funnel, the funnel's ×1.3–1.5 is all there is.
    - So the June builds differ in their gap surface (conducting), but that difference is not the asymmetry
      mechanism. What remains are size levers acting at the peak: strip width (item 11) and prompt funnelling.
    - Running: strip-width full-sim arms (condor 4520870: hybrid w410 / w690, insulating gaps, 36 pairs each, paired
      with glg0). Score with `paired_odd.py glg0 ~/x17/nb_asym_raw/simgapleak w410 <dir> w690 <dir>`
      (npz land in `simgapleak_npz` as `glw410_*` / `glw690_*`).
