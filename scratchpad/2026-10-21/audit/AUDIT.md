# Audit: the det3 X neighbour asymmetry (2026-10-21, fresh eyes)

Scope: everything in `HANDOFF.md` (top two sections), `scratchpad/2026-10-19/asym_deck/FINDINGS.md` (routes 1–40,
rounds 1–8), `NB_ASYM.md`, `NB_FAR.md`, the r5 scripts/JSON, and the git log 10-18 → 10-21. Plus new checks
from this audit (below, files in this directory) and the main session's Y-view-by-class / route-14 rerun
(`scratchpad/2026-10-21/yview/yview_class.json`).

Notation: a1 = α(n+1) − α(n−1) in DAQ order, α = in-phase content of a neighbour relative to the central;
s1X = α(+1) + α(−1) (the X ±1 sharing), s1Y = the same in Y ("Y share"). Y-pad / X-pad class = parity of
(X peak + Y peak), i.e. the pad type under the charge.

> **Addendum 2026-10-22 (strip film photo, FOLLOWUP §14–15).** The print film of the resistive strips was photographed
> (Dylan, macro shots with a steel ruler): pitch **760.7 µm** (±1 % systematic), fill 0.694 (spec 550 / 800 = 0.6875).
> Step 6 below named the pitch as the decisive test (registration needs a 780 µm, pad-locked print). If the film is the
> production pattern (being confirmed), registration is excluded: with any pitch ≠ 780 µm the strip offset recurs every
> ~30 mm, so it can reach ≤ 0.025 per 64-column cell (data 0.07–0.20) and must leave a ~30 mm ripple in s1 and a1 that the
> data do not show (≥ 5–10× below the sim). The sim's strip imprint is too strong on every board. See FOLLOWUP §15.

## 0. Summary

- The observation is solid and well characterised as data. It is a fixed, smooth, two-signed, X-only field of
  a1 on every June-built board, with an amplitude that changes between sessions.
- The leading explanation, resistive-strip registration, is **weaker than the ledger says**. Registration
  pulls the cloud toward pad-pitch attractors, and any attractor also changes the ±1 *sum* (the sharing) along
  with a1. On det2, the best board for this test because its field is two-signed, the sharing does not follow
  |a1| at all.
  - With all events (no pad-class labels) s1X − s1Y is −0.08 … −0.01 in every det2 band while a1 runs
    −0.21 → +0.14 → −0.11. The per-cell slope is Δs1/|a1| = −0.19 ± 0.06.
  - With the (migrating) Y-pad labels: within ±0.05 where a1 is −0.33 … +0.23, slope +0.11 ± 0.05.
  - The sim's own registration predicts ×20 that: Δs1 ≈ 2.2·|a1|.
  - A generic strip-quantisation toy predicts ≳ 0.1.
  - The data's charge is *displaced*, not *redistributed*. Route 20 already said "n+1 gains what n−1 loses,
    the ±1 sum unchanged", and the 2026-10-20 quadrature test found no cos φ term. Neither was connected to
    the even term that registration requires.
- The sim's "1/3 of the size" and the levers that bring it to ×1.8 (width, funnel, anisotropy) are all inside
  the same attractor family. They inherit the even-term problem, so closing the size gap with more levers
  inside that family is unlikely to converge. That, more than lack of effort, is why the week felt like
  wheel-spinning: the sim was being tuned on one projection of the fingerprint (odd a1 and its timing), while
  a projection that discriminates (the even term) was never scored.
- The "lives at the Y pads" story (route 11) is overturned by the main session's X-independent label: a1 is the
  same over X and Y pads (`scratchpad/2026-10-21/m3class/`). Every fingerprint built on X-peak parity classes
  (A_Y / A_X by class, narrow X-pad flips) measures the effect convolved with peak migration, not its location.
  That voids the class arguments used for and against several models (routes 9, 12, 13, 15–17, 26).
- The same label problem hit the width story. Night to night, a1 tracks the median cloud width (0.18 / 0.14 /
  0.10 at median Y share 0.52 / 0.77 / 0.93). Within a night there is no per-event width trend once labels are
  dropped, because the per-event Y share mostly says which pad type the charge sits on. The old "one curve"
  was partly migration.
- Several "controls" were not clean (det2, det6) or rested on an unvalidated fast model. Re-graded in §2. The
  ones that matter most are re-established by within-board tests: electronics, Y→X transfer, gap conduction as
  a local cause. The fast-model exclusions (routes 15–17, 31) carry little weight.
- Way forward (§4): make the even term a mandatory score, test a skewed-response family (the X weighting field
  leaning toward one neighbour, e.g. steered by offset strips, with no charge concentration), run the charging/time tests the session data already allow, and get
  the physical print: a photo or the mask of a June board in two regions of opposite sign. That photo is the
  one measurement that can settle registration either way.

## 1. The observation (data only)

All numbers from the raw-ADC survey (pedestal-subtracted ADC, DAQ channels; no strip map, CM, clustering or M3),
matched X/Y events with the edge cut (`r5lib.ev`), recomputed here in `audit_data.py` → `audit_data.json`.

| build (bench) | events | a1 X | a1 X, Y-pad / X-pad class | a1 Y | X field (8×8 χ²/ndf, const) |
|---|---|---|---|---|---|
| det3 June-15 | 14 550 | **+0.137 ± 0.003** (coin 0.71) | +0.221 / +0.020 | −0.003 ± 0.003 | 332/63 |
| det2 (Arnaud, June 12) | 5 567 | −0.001 ± 0.005 | −0.002 / +0.002 | −0.010 ± 0.004 | **736/61** (bands ±0.2–0.3) |
| det4 (June 22) | 914 | +0.042 ± 0.013 | +0.078 / +0.000 | +0.014 ± 0.011 | sign change across x |
| det6 (June 24) | 1 432 | +0.012 ± 0.007 | +0.017 / −0.004 | −0.016 ± 0.006 | 47/25, + → − across x |
| det7 (June 22) | 689 | +0.052 ± 0.015 | +0.078 / +0.018 | +0.012 ± 0.012 | 21/8 |
| det3 May (control) | 4 572 | +0.004 ± 0.005 | −0.004 / +0.020 | −0.010 ± 0.005 | 78/56 |
| det1 Jan (strips along x) | 1 646 | −0.005 ± 0.007 | −0.010 / +0.012 | +0.001 ± 0.009 | 34/27 |

Facts that are data, not interpretation:
1. **In the raw ADC.** On det3 June-15, P(n+1 > n−1) = 0.71 (matched sample) and 0.85 on the 6-24 long run.
   Every other board and view through the same code is at 0.49–0.58.
2. **X only, in both pad classes of Y.** The main session found the Y view flat in each class on every board
   (det3 June-15 a1Y −0.002 / −0.003 ± 0.004; Y maps 30/15 and 17/15). A shared a1Y ≈ −0.01 appears on quiet
   and June boards alike, so it is instrumental.
3. **A smooth two-signed field, fixed per board.**
   - Bands by 64 rows/columns (`bands.py`):
     - det2 by y (all events): −0.21, +0.02, +0.14, +0.14, +0.03, −0.01, −0.04, −0.11 (Y-pad labels: −0.33 … +0.23).
     - det4 by x: −0.17 … +0.25.
     - det6 by x: +0.04 … −0.04.
     - det3 by x: +0.12 … +0.29, flat in y.
   - The same det2 band geometry appears in its two runs on two FEU slots. det3's two FEU-3 sessions give
     the same map phase.
   - No step at any connector boundary.
4. **Same over X pads and Y pads** (main session, `m3class/ysh_class.py`).
   - The X-peak parity labels (det3 Y-pad 0.221 vs X-pad 0.020; det2 Y-pad map 971/15 vs X-pad 15/15) use the X
     peak, which the asymmetry itself moves, so they can't locate it.
   - With an X-independent label, the Y-view share s, det3 June-15's a1X is flat in s: +0.17/+0.19 (s < 0.3),
     +0.11…+0.14 (0.5–1.0), +0.15/+0.18 (s > 1.25). On det3 May the Y-pad fraction falls 0.80 → 0.20 over that
     range, so an effect over Y pads only would fall 3–4×.
   - Within one s bin the old labels split it +0.26 / −0.48: pure migration.
   - det2's positive band: +0.17…+0.19 at small s, +0.10/+0.03 (±0.04) at s > 1.1 (weak stats).
   - M3 can't give an X-independent column either (pointing 0.9–2.0 mm vs 0.78 mm pads).
5. **Charge moved, not added.**
   - Within board (det3 June-15, Y-pad): X+1 0.216, X−1 0.116 of ΣY, against 0.160 / 0.163 on det3 May.
     n−1 is depleted by about as much as n+1 gains.
   - ΣX/ΣY is *lower* in the Y-pad class (0.907) than in the X-pad class (0.925), not higher by about 0.10 as
     a Y→X transfer would make it.
   - Caveat: the cross-board numbers (det3 May) carry gas/width differences.
6. **The ±1 sum does not follow the field (new, §3.1).** Per det2 band (all events), s1X − s1Y = −0.05, −0.06,
   −0.08, −0.07, −0.03, −0.04, −0.04, −0.01 while a1 goes −0.21, +0.02, +0.14, +0.14, +0.03, −0.01, −0.04,
   −0.11. Per-cell fit (8×8): ds1X/d|a1| = −0.19 ± 0.06 (all events), +0.11 ± 0.05 (Y-pad labels).
7. **Pad-pitch locked.** Inside det2's positive band, a1 > 0 in every 8-column block across all 512 columns. No
   1/40 (31.2 mm) harmonic on any board (det3 s1 p = 0.12; sim at 800 µm predicts 0.12 amplitude, data ≤ 0.03).
8. **Amplitude set per session, by cloud width.**
   - det3 June-15 sessions (centre chips): 0.26 (6-18), **0.03 (6-19/20)**, 0.16 (6-22/23), 0.26 (6-23/24,
     drift 600 V), 0.17 (6-25), 0.09–0.10 (6-27/28), 0.05–0.07 (n_TOF run_79, July).
   - A fixed map × per-session scale fits det3: 0.96 / 1.57 / 0.44, χ² 47/60, 126/61, 47/59.
   - Night to night it tracks the cloud width: median Y share 0.52 / 0.77 / 0.93 with a1 +0.18 / +0.14 / +0.10
     (6-23/24, 6-22/23, 6-27/28; `sess_width.py`, all events). det2's two runs differ ×2.5 at median Y share
     0.35 vs 1.07.
   - **Within a night it does not** (all events): the per-event Y share mostly labels the pad type, and a1 is
     the same over both. The older "one curve" (Y-pad-labelled events vs Y share, all sessions) used the
     migrating labels and overstated a per-event width dependence [sel]. Within 6-23/24, a1 is +0.08 for the
     shortest-pulse quartile and +0.29 for the longest.
   - The drift scan 1100 → 100 V raised it 0.05 → 0.16, reversibly.
   - **6-19/20 is off that curve** (0.03 where ≈ 0.13 is expected). It is also the only det3 night with drift
     current spikes (4–9 µA) and a sagging drift HV.
9. **Within a board the width trend is not obvious.** On det3 June-15 the Y share falls 0.89 → 0.48 from chip
   0 to chip 7 while a1 peaks in the middle (0.29). The board's own width gradient is a confounder for the
   map's shape (pooled sessions with different illumination add to it).
10. **Timing depends on pulse width.** The excess (n+1 − n−1) relative to the central is:
    - shortest-pulse quartile (near-vertical): late, +52 ns (fit delay) / +49 ns (peak);
    - middle quartiles: −4 and −52 ns;
    - longest-pulse quartile: early, −54 / −29 ns.
    The all-angle average looks prompt (−10 ns). Pulse width is a mix of track angle, drift speed and the
    sample's composition.
11. **Only the June builds.**
    - det3 May, det4 May, det1 and det0 are quiet (X field rms ≤ 0.02).
    - det1 has a weak, fine-grained Y hint (Y 103/50 on 8×8, p 0.11 on 4×4).
12. **June builds also conduct across the strips.** The late (600–960 ns) ±2 content in the across-strip view
    is 12–31 ×10⁻³ on June boards vs 1–9 on May and earlier. **But it is flat inside each board** (§3.3),
    whatever the local a1.
13. **The production position carries it.** wft x residual vs M3, per cell, K = −0.77 ± 0.13 times the
    raw-ADC centroid bias (det2, det3 runs).

## 2. Every hypothesis, re-graded

Key: **X** excluded · **W** weakened / doesn't fit alone · **S** supported · **U** untested or not properly
tested. Flags: [ctl] contaminated control · [avg] a board average hides cancelling structure · [sel] selection
effect · [sim] rests on an unvalidated model · [rev] verdict reversed during the week.

| # | Hypothesis | Ledger | Audit | Why |
|---|---|---|---|---|
| 1 | Code artefact | ✗ | **X** | Independent raw extractor; the same code gives 0.50 on every other board/view. Solid. |
| 2 | Hot channel | ✗ (trap found) | **X** | Veto applied; time-dependence claims made before it were withdrawn [rev]. |
| 3 | 39:40 strip beat | ✗ | **X** | Per-column periodogram with permutation null, all boards. Solid as a data result. It also says no 800 µm strip pattern leaves an imprint on any board, May included, while the sim says it should (s1 0.12 at 1/40). So either every board has a pad-pitch print, or the sim overstates what strips do. |
| 4 | Electronics after the input | ✗ | **X** | Rests on noise-copy injection on det3's own channels, the connector-edge test (det3's own data, 6-23: physical neighbour +0.21 ± 0.03, DAQ neighbour 0.00) and the sign change. "Same FEU clean on det6" is [ctl] (det6 carries a field), but det3 May on FEU 3 is clean and the other tests don't need a control. Solid. |
| 5 | Connector / M1 / linear mixing | ✗ | **X** | Edge test is within-board (6-28 only 2σ). Mixing fits injected into det2 [ctl], det3 May and det1; det3 May / det1 alone carry it. Solid. |
| 6 | Track angle | ✗ | **X** | Slope measured on det2 [ctl, avg], but the vertical-only maps and an x-tilt's inability to make y-bands don't need det2. Solid. |
| 7 | Gas flow | ✗ | **X** as cause | Flat top–bottom, X only. Gas as a width modulator stays in. |
| 8 | Bus traces / one-sided copper | ✗ | **X** | The 3D solve depends on an unknown stack-up, but any fixed layout has one sign; the two-signed field kills it. Solid. [rev] (was the candidate on 10-19). |
| 9 | Floating strips (static) | ✗ | **U** (reopened) [sim][sel] | Rejected because the copy came out "equal in both pad classes" (now what the data do, fact 4) and "growing with cloud size" (the data shrink). Static conductor model only. A strip that steers the induced signal without concentrating the charge is the one strip variant that could pass the even-term test (§4.2), so this deserves a proper redo. |
| 10 | Plain one-step copy | ~ | **X** | Within-board bookkeeping (n−1 depleted, det3 Y-pad 0.116 vs 0.160–0.163 on det3 May), ΣX/ΣY not raised (the narrow X-pad −0.40 is a migration fingerprint of a non-copy, supporting only). The original bookkeeping was vs det2 [ctl, avg]; the within-board and det3-May versions hold. [rev] (the "copy" framing ran 10-19 to 10-20, and `uncopy.py` and the "+ side clean" sim targets came from it, see §3.4). |
| 11 | Lives at the Y pads | ~ (selection) | **X** [sel][rev] | Overturned by the X-independent label (fact 4): a1X is flat in the Y-view share, so the effect is the same over X and Y pads; the class split is migration. Excludes anything tied to the Y pads (vias, Y-pad copper, the bus beside the Y pads). |
| 12 | Rigid X-view offset | ~ | **X** | A rigid X-vs-Y offset is translation-invariant in X: zero marginal a1 by construction. It was "needed" only for the narrow X-pad flip, which is a migration fingerprint [sel]. |
| 13 | Mixed responses | ✗ alone | **W** [ctl][sim][sel] | Scored as (det3 − det2) vs (model − sim D0). det2's board average carries class structure (its negative bands have narrow A_X +0.31, positive −0.18), so the "effects" baseline is contaminated. The marginal-a1 argument (≈ 0 for every mixture) doesn't need det2 and stands. |
| 14 | Y→X transfer | ✗ | **X** (re-established) | Original vs det2 [ctl], and cross-board ΣX/ΣY spans 0.86–1.01 (gain ratios), so it was weak. Rerun within board (main session): Y-pad 0.907 < X-pad 0.925, where transfer predicts +0.10. Excluded regardless of det2; and since the effect is the same over X pads (fact 4), a transfer at the Y pads is excluded independently. |
| 15 | Extra bus-side sensitivity regions | ✗ | **U** [sim] | Fast static model, which misses the full sim's own baseline (route 31). Not evidence either way. |
| 16 | Skip transfer n−1 → n+1 | ✗ | **U** [sim] | Same fast model. Notably the data's x-conditional (route 5.5: n+1 takes a share of what n−1 should get) is exactly a skip-like, sum-conserving transfer. Deserves a proper test (§4). |
| 17 | Shift/mix + copy | ~ | **W** [sim][ctl][sel] | Fast-model χ² and det2-based effects. |
| 18 | det3 strips rotated | ✗ for det3 | **X** | Strip direction from the slow view, every board. Solid. |
| 19 | Same PCB May/June | ? | **U** | Dead-line fingerprint doesn't match; bulk-level can't be excluded. |
| 20 | Moved or added | ✓ moved | **S** | Original vs det2 [ctl]; holds within board and vs det3 May (fact 5). |
| 21 | Excess is prompt | ✓ | **W** [sel][rev] | "Prompt" is the all-angle average of a late (short pulses) and an early (long pulses) excess. Replaced by #39. |
| 22 | X-only smooth two-signed field | ✓ | **S** | Data, many boards, Y flat in both classes. Solid. Map shapes partly confounded by the cloud-width gradient and session mix (fact 9). |
| 23 | det2 a clean control | ✗ | **X** | Correctly found [rev]. Every earlier use of det2 as "copy-free" is flagged here. |
| 24 | Pad-pitch locked | ✓ | **S** | Data. Holds for whatever the cause is, not only strips. |
| 25 | Registration geometry (A sin(φ0 + k·r)) | ~ | **W** | The geometry fit works, but the predicted cos φ sharing term is absent (−0.007 ± 0.009). Strengthened against by §3.1. |
| 26 | Sim strip offset makes it | ✓ shape, ~ size | **W** [sim][sel] | Odd term at ~1/3 of the size, but the sim's even term (s1 0.63 → 0.77–0.88) is ~20× what det2 allows (§3.1). The "det3 class pattern" it reproduces is migration-convolved, so it is a size check, not a shape check. |
| 27 | Rotation alone | physics | **S** (as stated) | Only via registration of a pad-pitch print. |
| 28 | Pin to strips with reversed boards | ? | **U** | det0/det1 quiet; can't pin. |
| 29 | Why June: across-strip late content | conduction | **S** as a June marker; **X** as the local cause | New (§3.3): flat inside det2 (8–13 across its bands) and det3 (27–32) whatever a1 is. Across boards: det6 has the most (31), with a weak field; det2 the least (12), with a strong one. |
| 30 | ±1/±2 timing, far halo | separate | unaffected, but see §3.4 | v5a tuning used "copy-free + side" numbers; under "moved" that side is depleted, not clean. |
| 31 | Static relocation (fast model) | ~ | **U** [sim] | Self-declared non-calibrator. |
| 32 | Sim a1(δ) is sheet transport | ✓ | **S** about the sim | Says how the sim makes it, not how the board does. |
| 33 | Sim excess timing | ~ | **W** [sim][sel] | Different pulse lengths in sim and data; no matched-width comparison for short pulses. |
| 34 | Cross-gap leak | June marker ✓, amplifier ✗ | **S / X** | The marker part is data; "not the amplifier" is now also data (§3.3), not only [sim]. |
| 35 | Gap funnelling | ~ ×1.4 | **W** [sim] | Inside the attractor family: raises the odd and even terms together. |
| 36 | Low ρ_s | ✗ (Y view) | **X** | The Y/X kept ratio (0.12 / 0.001 vs 0.9) is a gross effect; robust even with an imperfect sim. [rev] (was a size+timing lever for a few hours). |
| 37 | Strip width | ✓ lever | **W** [sim] | Same family; s1 rises with it in the sim (0.760 → 0.808). |
| 38 | Anisotropic strips (A = 10) | ✓ best model | **W/U** [sim] | One tuned parameter. Y view not computed (production Y kernels kept), so "keeps Y" is untested. Inside the attractor family. Its timing pattern is the one real success. |
| 39 | Excess timing vs pulse width | ✓ | **S** | Data. Pulse width is a mixed variable (angle, drift speed, session). |
| 40 | Production wft carries it | ✓ | **S** | Data; the correction is mechanism-independent (map × scale), which is its strength. |

Verdicts reversed during the week:
- "det3 only" → "every June board".
- "bus traces" → excluded.
- "one-sided copy" → "charge moved, not a copy".
- "copy lives at the Y pads" → "partly selection".
- "strip print not reproduced by the sim" (round 5) → "reproduced at 1/3" (round 6).
- "det2 copy-free" → "det2 carries the field".
- "decays over weeks" → "cloud width".
- "prompt" → "late or early by pulse width".
- Gap leak: "first combination that reaches the size" (kernel) → ×0.9 (full sim).
- ρ_s: lever → closed by the Y view.

Each reversal was a correct correction; the pattern is that kernel- or fast-model-level claims kept failing in
the full sim, and board-average claims kept failing in maps.

## 3. New checks in this audit

### 3.1 The even term (`share_vs_a1.py`, `bands.py`, `toy_attractor.py`, `toy_strips.py`)
- Any mechanism that pulls the cloud toward pad-pitch-periodic points offset by δ from the pad centres
  (registration, funnelling, anisotropic spreading, wider strips) makes a1 odd in δ *and* changes the ±1 sum
  s1 evenly in δ.
  - The sim's own table (`r5/sim_delta.json`): s1 0.63 at δ = 0, 0.70–0.77 at |δ| = 75–225, 0.80–1.14 near
    the gaps, while a1 is ±0.05. Fit: s1 = 0.64 + 2.2|a1|.
  - Strip-quantisation toy (1-D: cloud σ, strips of width w at the pad pitch, the charge landing on a strip
    spread across it, pads read through box ⊗ Gaussian): Δs1/a1 ≈ 0.4–4 depending on σ, w and the induction
    width. At the δ needed for a1 ≥ 0.1, Δs1 ≈ +0.05 … +0.25.
  - A weak generic pull (20 % of the charge moved fully to a line) is the most forgiving case: a1 0.2 with
    Δs1 ≈ +0.05.
  - Hard ceiling: with σ = 450 µm clouds, even full equalisation across 700 µm strips gives |a1| ≤ 0.06.
- Data: per det2 band, s1X tracks s1Y to ±0.04 with all events (±0.05 with Y-pad labels) while a1 swings by 0.35 (0.55). Per-cell regression
  s1X = c0 + c1·s1Y + k|a1| + k'·a1:
  - det2 Y-pad: k = +0.11 ± 0.05 (all events −0.19 ± 0.06);
  - det3 June-15: k + k' ≈ −0.3 (mostly one sign there, so k and k' are degenerate);
  - controls det3 May / det1: k consistent with 0 at |a1| rms 0.03.
- Reading: the board's mechanism moves charge from n−1 to n+1 with the sum unchanged (a displacement), not by
  concentrating it (an attractor). That is consistent with route 20 and with the null quadrature test, and it
  is the strongest data argument against registration as modelled. It is not airtight: a response that is skewed rather than attracting (§4.2) could escape.

### 3.2 Y view by class and route 14 without det2 (main session)
Incorporated in §1 (facts 2, 4, 5) and §2 (#11, #14). The later X-independent label (`m3class/`) replaces the
class-split reading of fact 4. The Y view is flat in both classes to ~0.005 on det3 June-15, so
a cancellation between classes in Y is excluded: the effect is a property of the across-strip (X) view only.

### 3.3 Gap conduction is a board-wide marker, not the local cause (`late_vs_a1.py`)

| board | by | a1 per band | X ±2 late (×10⁻³) per band | corr(\|a1\|, late) |
|---|---|---|---|---|
| det2 | y | −0.21, +0.01, +0.15, +0.16, +0.03, −0.01, −0.04, −0.11 | 9.7, 12.3, 13.0, 9.6, 11.7, 11.4, 12.9, 8.1 (±1.6–2.3) | −0.48 |
| det3 June-15 | x | +0.09 … +0.22 … +0.09 | 26.9 … 32.4 (±1.5) | −0.35 |
| det6 | x | +0.05 … −0.08 | 23 … 42 (±7–9) | −0.08 |
| det3 May | x | ≈ 0 | 3.7 … 11.6 | |

The conducting gaps are uniform over each June board; the field is not. Whatever the bulk process did to the gaps,
it isn't what sets the local sign or size of a1.

### 3.4 Copy-model leftovers that were never revisited
- NB_ASYM §2 treated the det3 X "+ side" (position order) as copy-free, giving ±1 share 0.30 against 0.42 on the
  other side, and concluded that v4c overshares and that σ0 can't fit share and delay together.
- Under "charge moved" (route 20) that side is *depleted*, not clean. The symmetric value (~0.36) is the right
  target, and v5a's 0.35 already matches it.
- Sim tuning targets for the X ±1 share and width should come from det3 May / det1, never from det3 June-15
  or det2.
- `uncopy.py` corrections assume the copy model and should not be used; `wft.xasym` (map × scale) replaces them.

## 4. The fingerprint, what was missed, and the way forward

### 4.1 Fingerprint any explanation must reproduce
1. X view (across the strips) only; Y flat in both pad classes.
2. Odd transfer n−1 → n+1 with the ±1 sum and ΣX unchanged (no even term: Δs1 ≲ 0.03 at |a1| 0.2).
3. A smooth field over the board, two-signed, pad-pitch locked (no 31.2 mm beat), with fixed geometry per board.
4. The same over X pads and Y pads (X-independent label); X-peak parity splits show only migration.
5. Amplitude larger on nights with narrower clouds and at low drift field (reversibly), not tied to the
   per-event Y share within a night; larger for longer pulses within a night; one session (6-19/20) off the
   night-to-night trend.
6. Excess late (+50 ns) for short pulses, early (−50 ns) for long ones.
7. June builds only. The June builds also conduct across their gaps, uniformly.
8. Same geometry on different FEU slots; nothing at connector edges.

### 4.2 What was missed
- **The even term** (§3.1). It was available from day 6 (quad_test, route 20 bookkeeping), but never used as a
  score for the sim levers. Every lever tried (width, funnel, leak, ρ_s, anisotropy) was scored on the odd
  amplitude and the timing only.
- **Mechanisms that skew the X response without concentrating the charge.** Uniform illumination plus a
  translation-invariant, mirror-symmetric response gives a1 = 0. Two ways out:
  - an attractor (registration as simulated), which brings the even term with it;
  - a response that is not mirror-symmetric: X pad n's weighting field leans toward n+1, for example because a
    partly conducting, offset strip steers the induced signal.
  The second kind shifts the centroid with the sum unchanged at first order, and it is the same over X and Y
  pads. Only its static version was ever modelled (route 9, `bus_float.py`), and that was rejected on a class
  criterion that is now void. The time-dependent version (strips as resistive electrodes in the weighting
  field, which is what the W2 kernels do for the pads) was never checked for an odd, sum-conserving term.
- **Time-dependent / charging explanations were dropped too early.** The facts that look like a surface state
  were all reassigned to "cloud width":
  - the drift-field dependence (reversible);
  - the session jumps, with 6-19 off the width curve on the night of drift-HV trouble;
  - June-only, with uniformly conducting gaps on those boards.
  The width curve is a correlation within the det3 sessions. No test was made of a1 vs time since HV-on or
  after a ramp or trip, or of a1 vs rate.
- **6-19/20 was never resolved.** Under "fixed map × width" it is a falsifier unless it is a different board or
  a different board state. Its raw_daq_data (HV scan + long run) can be decoded.
- **Within-board width gradients** (det3: Y share 0.89 → 0.48 across x). Maps pooled over sessions mix
  amplitude and illumination. The fixed-map claim should be redone per session in width-matched cells.
- **The [ctl] items above.** Fortunately the important ones survive when redone within a board.

### 4.3 Ranked next steps

Doable on archive data now:
1. **Score the even term everywhere.** Add Δs1 vs a1 (per δ bin in the sim; per cell or band in data, s1X − s1Y)
   to `paired_odd.py` and to the deck's lever table. Any sim arm whose even term exceeds ≈ 0.1 × its odd
   amplitude per unit |a1| fails. That will likely retire the whole lever family. One afternoon.
2. **Time and charging tests.**
   - a1 (centre chips, hot-channel veto) vs time within each long run, in 15–30 min bins.
   - The same around every HV change: 6-27/28 resist 3:3 → 3:4, the drift scan's step edges, 6-19's spikes.
   - Beam vs bench at matched width (rate).
   - Charging predicts transients with minutes-to-hours time constants after a change; pure width predicts
     a step that follows the gas.
3. **Decode and map 6-19/20** (raw_daq_data): does it show det3's geometry at a low scale (fixed map, state
   change) or nothing at all (a different board or cabling)?
4. **Per-session, width-matched maps** for det3 (6-22/23 vs 6-23/24 vs 6-27/28). Is the geometry really fixed,
   or does the ridge move?
5. **Test a skewed-response family.**
   - First a cheap check: in the existing sim arms, measure the even term per δ bin separately from the odd one.
     How much of the sim's a1 is attraction and how much is induction skew from the offset strips?
   - Then redo route 9 with resistive (not floating) offset strips in the weighting field.
   - Score on the marginal a1, the even term, the width curve and the timing, and never on X-peak parity classes.

Needs the board, the print spec or the detector:
6. **Photograph or measure the strips** on a June board (det2 is ideal: rows 0–63 and 128–255 have opposite
   sign) and on det3 May. Registration predicts opposite strip-centre offsets (±100–300 µm) in those two regions
   and a 780 µm pitch. Seeing no offset difference kills it outright; seeing it revives it, and then the even
   term has to be explained. Also read strip width and gap from the photo or mask.
7. **Detector tests (when access returns):**
   - a collimated X-ray / ⁵⁵Fe scan across one pad in a region of known sign: a1 vs sub-pad position tells an
     attractor (position-dependent, s1 compressed) from a displacement (constant shift);
   - a1 vs resist/mesh HV ramps and dwell time (charging);
   - a1 after flushing dry gas for a day (humidity, given the bench's ~1 % H₂O).
8. Keep `wft.xasym` as the analysis answer regardless of mechanism: it is data-driven (K = −0.77 ± 0.13) and
   only needs a per-run scale.

### 4.4 Honest status of the leading hypothesis
Strip registration explains (1), (3) and (7) naturally. It needs a 780 µm pitch against an 800 µm quote. In the sim
it gives ≤ 1/2 of the size, and only with a stack of tuned levers. Its necessary companion, a sharing change of
the same order as the asymmetry, is absent in the data. It should be held as one candidate among a small family
(registration as an attractor, a skewed response, surface charging), with the even-term score and the strip photo as the
tests that separate them, not as the working answer.

**2026-10-22:** the strip photo came in (FOLLOWUP §14–15): 761 µm, not 780. If the film is the production pattern, both
registration readings (attractor and skewed response) are excluded, and of this family only surface charging (or another
non-strip, pad-locked mechanism) remains.

## Files (this directory)
- `share_vs_a1.py` → `share_vs_a1_nc8.json`: per-cell s1X vs |a1|, a1 (data) and the sim reference slope.
- `bands.py` → `bands.json`: a1, s1X, s1Y, a1Y, α±1 per 64-row/column band.
- `toy_attractor.py` → `toy_attractor.json`, `toy_strips.py` → `toy_strips.json`, `toy_strips_prf.txt`: even vs odd
  terms for generic attractors and strip quantisation.
- `late_vs_a1.py` → `late_vs_a1.json`: across-strip late content vs a1 per band.
- `audit_data.py` → `audit_data.json`: survey, maps, width curves (all events) for the deck; `sessions.py` adds the
  det3 session table.
- `sess_width.py` → `sess_width.json`: per-night a1 vs median Y share / pulse width, all events.
- `toy_curves.py` → `toy_curves.json`: toy curves for the deck.
- `make_audit_deck.py`: the new deck (published at the old slug).
