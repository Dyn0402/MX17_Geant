# W2 night shift, 2026-08-08 → 09

**Status: COMPLETE.** Handoff executed:
`design/report/HANDOFF_W2_OVERNIGHT_2026-08-08.md`. Handoff §§1–5 done; §6
(Stage B/C + T13, a stretch goal) **not reached** — see §8. Commits: `e7759e0`,
`f187a93`, `9cc1db6`, `99b6e4d`, `1d9257c`, `8611934`, `937c184`, `12d4d1f`,
`a56c87a`, `96c5962`, `5838690`.

---

## Verdict

**The W2 grid is produced, complete and accepted.** All 66 family slabs landed
(04:12), the combine ran clean (05:46), and **all four ρ_s products pass every
pre-registered acceptance bar**. The absolute-amplitude defect that V6 exposed
in W1 — every kernel ~27 % low — is fixed: prompt capture moves 0.670026 →
**0.841977**, i.e. **+25.66 %**, against V6's static +27.2 % reduced by the
known ny=512 grid term (predicted +25.5 %).

**T10 still FAILS on W2, so plan §7 step 5 remains sequenced.** Like-for-like
at ρ_s = 2 MΩ/sq, the worst DREAM-shaped residual moves **8.26 % → 7.55 %**
against a 2 % bar. W2 improves it by 0.71 pp (−8.6 % relative) and does not
come close to certifying the fast path. The plan's note that V6's re-solve
might overturn the T10 verdict is now answered: **it does not.**

Attribution, stated with its caveat rather than as a clean result: 0.71 pp is
above the ~0.45 % ny=512 grid bound (audit C6) and well above the ~0.1 % stack
term measured last night, so some of it is genuinely the boundary model. But
**ρ_s alone moves the same number by 0.65 pp** within W2 (rho1M gives 8.20 %,
rho2M 7.55 %) — comparable to the W1→W2 move itself. So the improvement is not
a sharp discriminator, and the honest statement is that W2 changes the number
by about as much as ρ_s does, while the FAIL verdict is robust to both.

c1 moved the wrong way — W1 +28.5 % → W2 +29.2 % (rho2M), +34.6 % (rho1M) —
confirming last night's finding that c1 is the softer comparator and the shaped
residual is the one to quote.

The W2 caching cert returned **0.0187 / 0.0188**, exactly the harness artifact
predicted in §4 before the products existed. It says nothing about W2 and must
not be read as either a pass or a regression until the harness is fixed.

## Headline numbers

| quantity | pre-registered | measured | verdict |
|---|---|---|---|
| `channel_capture_prompt` | **0.841977** | **0.841977** (all 4 ρ) | ✅ exact |
| `channel_capture_late` | flat, ≈ prompt | 0.845369–0.845381 (+0.40 %) | ✅ no decay |
| `x_fraction_prompt` | **0.5 to 1e-4** | **0.50000002** | ✅ to 2e-8 |
| `x_fraction_late` | record only | **0.5** to 1e-8 | see below |
| ratio vs W1 | 1.20–1.30 | **1.2566** | ✅ |
| T10 slow path, W2 rho2M | vs 8.26 % (W1) | **7.55 %** | ❌ FAIL (bar 2 %) |
| LUT caching cert, W2 | expected ~1.87 % (harness artifact) | **0.0187 / 0.0188** | as predicted |

**The capture number is a real cross-check, not a tautology.** 0.841977 was
pre-registered from an independent **full-grid CG** solve; the measured value
comes from the **family-assembled** kernels via a completely separate path
(66 Bloch families → ifft2 → row/column sums). They agree to the 6 digits
printed. Together with the `gd_rank` identity (§6) that is a second
full-scale certification of the family assembly, at an N the solver's own
battery cannot reach.

**Two departures from expectation, both benign and both recorded:**

- **Late capture rises 0.40 % rather than staying flat.** W1's prompt and late
  agree to 4e-9; W2's late is 0.845377 against a prompt of 0.841977. It is a
  *rise*, not the decay the handoff warned about, and it is the same 0.40 % at
  every ρ_s. Not chased.
- **`x_fraction_late` is 0.5 too**, to 1e-8. A late-time X/Y asymmetry was
  expected on the grounds that the 800 µm ESL period breaks the 780 µm
  translation symmetry once the sheet conducts; it does not appear in this
  observable. The W1 baseline behaves the same way (0.4999999999750454), so
  this is not a W2 property — the expectation was simply wrong.

**Sharing** (ρ_s = 1 MΩ/sq, vs the W1 `s1_ny1024` product):

| | W1 | W2 | ratio |
|---|---|---|---|
| X `view_total_prompt` | 0.877482 | 0.880258 | 1.0032 |
| Y `share_prompt` d=0 | 0.4834 | 0.4980 | — |
| Y `share_prompt` d=±1 | 0.2585 | 0.2452 | — |
| X `tau_1e` (on-strip) | 5.15e-8 s | 3.49e-8 s | 0.68 |

The X peak moves **+0.32 %**, against the ~+0.2 % the handoff predicted from
V6's `pad_split` — right direction, right order. Relaxation is ~30 % faster in
W2, which is expected once the inter-pad channel is no longer pinned to ground.
The Y view's `view_total_prompt` ratio is 2.52, but on a quantity of 0.0007 →
0.0018: that is the checkerboard null, where the Y row through an on-strip
deposit owns no pad, so it is a large ratio on a negligible number rather than
a discrepancy.

The capture and x_fraction predictions were **pre-registered before any product
existed** (2026-08-08 23:55), from the solver author's independent full-grid CG
on the production box: two view unions at 0.420988 each. Recorded in
`scripts/condor/` tooling and the acceptance checker at the time of writing.

---

## 1. Two changes to the handoff's plan, both forced

### 1.1 The combine cannot run on the laptop

Handoff §3 gives a laptop `xrdcp` recipe for pulling the 66 slabs. It cannot
work: `/eos/experiment` redirects from `eosuser.cern.ch` to
**`eosexperiment.cern.ch`, which does not resolve off-site** (`getent hosts`
returns nothing). xrootd reports this as the unhelpful `[FATAL] Invalid
address`, including for a literal IPv4 target, which is why it reads like a
client bug. ssh and https to CERN are unaffected; it is that one name.

So combine runs at CERN as a condor job (`scripts/condor/run_w2_combine.sh`,
`w2_combine.sub`): pull the slabs to local scratch, combine, push the products
back to EOS, keep an AFS copy for rsync.

### 1.2 The recertification also moved to CERN

Independently of the above, `RESPONSE_SIM_PLAN`'s machine-roles section forbids
LUT builds on the laptop — the build peaks ~3.7 GB, the machine is shared, and
the documented failure mode is a **silent OOM kill that reads as a hang**. The
laptop had ~8 GB free with 2 GB already swapped. Condor gives it an
uncontended, requested 16 GB.

### 1.3 T10 is run on rho2M, not rho1M

The handoff says "re-run t10_slowpath against the W2 product" without naming a
ρ_s. The published 8.26 % baseline used `greens_comb_rho2M_dk50um_g19um.npz`
(`design/report/t10/t10_prod_ny1024.json`), so rho2M is the only choice that
changes one thing at a time; rho1M would move ρ_s and boundary model together.
rho1M is carried as a sensitivity point. The baseline's own calib is used
(md5 `ece7ccd5…`, shipped to AFS as `calib/aval_calib_t10baseline.json`) —
none of the calib files already on AFS matched it, so a default run would
**not** have been like-for-like.

### 1.4 `test_time_grid.py` cannot be re-run on W2

It requires a 4× denser (nt=240) re-solve of the *same* boundary model — a
second 66-job fleet. No dense product exists. The W1 certification (<0.5 %)
carries: the time axis is `kernels.log_times(60)`, byte-identical between W1
and W2, and W2 changes only the spatial boundary condition. Documented rather
than skipped silently or faked. See §4 — this carry-over became less
comfortable during the night.

---

## 2. The control: the environment reproduces the published baseline exactly

Before trusting any W2 number, the whole cert chain was dress-rehearsed against
the *W1* product that produced the published T10 result, so that a difference
could be attributed to physics rather than to LCG/AFS/condor.

| | published (laptop, 08-07/08) | rehearsal (LCG_105, condor) |
|---|---|---|
| worst DREAM-shaped residual | 8.26 % | **8.26 %** |
| raw current | — | 9.41 % |
| channel charge | ≤0.09 % | **0.09 %** |
| worst c1 shift | +28.5 % | **+28.5 %** |
| ladder (aval / 2 mm / 10 mm / 30 mm) | 8.26 / 2.66 / 0.75 / 0.20 % | **identical** |
| in-gap Y c1 | 0.380 → 0.365 | **0.380 → 0.365** |
| k=0-only self-check | — | 4.14e-15 |

The environment question is closed.

---

## 3. Three environment traps, measured not guessed

1. **LCG_105's `setup.sh` exports `BASE`** (→ the gcc release dir). A variable
   named `BASE` set before the source and used after silently resolves into
   CVMFS; this killed the first rehearsal instantly. Blast radius measured
   rather than patched on suspicion: of `BASE SRC WORK OUT SLABS OUTD CALIB
   PROD EOSDIR EOSBASE RHOTAG`, **`BASE` is the only name LCG_105 clobbers**.
   That mattered — `run_w2_combine.sh` consumes `WORK`, `SLABS` and `OUT` after
   its own source, and would have failed at ~03:00 with the whole fleet's
   output already on EOS.
2. **Condor logs are stale until job end.** CERN removed `stream_output`, so
   `.out`/`.err` appear only when the job finishes — and a resubmitted job
   reuses the filenames, so while it runs the files still hold the *previous*
   attempt's content. Mine still showed the already-fixed `BASE` traceback.
   Check `ClusterId`/`JobStatus` before believing a log.
3. **Condor liveness fields lag by ~13 min** — see §5.

---

## 4. The LUT caching cert does not reproduce its published 1e-4

Discovered by the rehearsal, and **independent of W2**. On the same W1 product
that the plan certifies at 1e-4, `test_lut_vs_solver` gives:

    worst per-channel residual 0.0187 (1.87 %) against a 2 % tol → PASS, barely
    col 32 X  max rel 0.0001     ← this is the published 1e-4
    col 32 Y  max rel 0.0187     ← Y carries the entire discrepancy

Two hypotheses were tested during the night.

**Hypothesis A — the kapton+glue stack. FALSIFIED.** The `_g19um` products were
rebuilt 08-08 14:52–14:54, after the 08-07 cert, so a glue-induced change was
the natural guess. But the *no-glue* `dk75um` product gives **the same 0.0187**
(X again 0.0001) to the 4 dp the log prints, on two products whose kernels
genuinely differ (ref 0.0733/0.0924 vs 0.0745/0.0941). The residual is
product-**independent**.

**Hypothesis B — `y_stride` decimation. ALSO FALSIFIED.** Commit `e36a39c`
"LUT: decimate y like x" landed 2026-08-08 14:00, *after* the certification,
giving ny=1024 products `y_stride = 2` while ny=512 keeps stride 1 — and both
products above are ny=1024. But the identical cert on an **ny=512** product
(`s1/greens_comb_rho2M_dk75um.npz`) also returns **0.0187 / 0.0001**. The
residual is grid-independent as well as product-independent.

### Root cause: the test compares two different times

Found in `test_lut_vs_solver.py` lines 90–93 and 152–153:

```python
t_ref = float(lut.t[-1])                       # 3000 ns, the LUT's last sample
it    = int(np.argmin(np.abs(t_src - t_ref)))  # nearest LOG-grid sample …
t_ref = float(t_src[it])                       # … = 3101.2 ns — t_ref is reassigned
kt    = int(np.argmin(np.abs(lut.t - t_ref)))  # clamps back to 3000 ns
...
r = np.array([ref[dd][it] for dd in ds])       # solver at 3101 ns
g = np.array([got[dd] for dd in ds])           # LUT integrated to 3000 ns
```

`lut.t` is a **uniform 1 ns** grid; `t_src` is the **61-point log** grid. Snapping
`t_ref` onto the log grid can land *beyond* the LUT's coverage, and `kt` then
clamps silently. The comparison is solver-at-3101 ns against LUT-to-3000 ns —
a **101 ns (3.4 %) misalignment**.

This explains every observation, including the ones that killed both hypotheses:

- **Y-specific**: the Y d=0 charge is still *decaying* at ~3 µs (resistive-sheet
  relaxation), so 101 ns of missing time is visible; X is prompt-dominated and
  flat by then, hence 1e-4.
- **Sign**: the LUT value (earlier time) reads *higher* than the reference —
  0.0924 vs 0.0907 — which is what a decaying Y channel requires.
- **Product- and grid-independent**: it is an artifact of the two time axes, not
  of kernel content.

**And it explains the published 1e-4 rather than contradicting it.** The nearest
log samples are 961.7 ns and 3101.2 ns:

| `T_MAX_NS_DEFAULT` | nearest log sample | relative to LUT end | result |
|---|---|---|---|
| 1000 (pre-Fix 1) | 961.7 ns — **below** | LUT has a real sample at 962 ns | match to 0.3 ns → **1e-4** |
| 3000 (current) | 3101.2 ns — **above** | `kt` clamps to 3000 ns | 101 ns gap → **1.87 %** |

So **Audit A1 / Fix 1** — which correctly raised `t_max` from 1000 to 3000 ns
because the DAQ integrates that long — flipped the nearest log point from just
*below* the LUT's coverage to just *above* it, exposing a latent clamp. The
plan's `1e-4` was true when written and was invalidated as a side effect of an
unrelated, correct fix.

**The fast path is not degraded.** This is a harness defect, not a caching
defect, and the fix is a one-liner (choose the source index at or below
`lut.t[-1]`, or interpolate the reference to the LUT's last sample). Not applied
tonight: it changes a certification and belongs in daylight with the W2 products
as its target.

**⚠ Reading note for the W2 cert.** Because the cause is the shared time axis,
**W2 will show the same ~1.87 %**, and that is not a W2 defect. Equally, the
W2 products are ny=512 → stride 1, so a *near-1e-4* result would not have been
evidence that W2 fixed anything either. The W2 caching cert is uninformative
about W2 until the harness is fixed; judge W2 on the slow path and the product
acceptance instead.

**Free sensitivity result: T10 is stack-insensitive.** The same diagnostic job
also ran the slow path on `dk75um`:

| | shaped | raw | charge | c1 |
|---|---|---|---|---|
| W1 rho2M dk50um_g19um (glue) | 8.26 % | 9.41 % | 0.09 % | +28.5 % |
| W1 rho2M dk75um (no glue) | 8.34 % | 9.39 % | 0.09 % | +24.1 % |

Across a real change of insulator stack the headline moves ~1 % relative. With
the grid term bounded at ~0.45 % (audit C6) and the stack term now measured,
a substantial W2 move is attributable to the boundary model by elimination.
Quote the shaped residual as the robust comparator; c1 is the softer of the two.

**Not fixed tonight, deliberately**: it still passes its 2 % tol, the cause is
not yet established, and any re-tune wants the W2 products as its target in
daylight. The plan's `caching ✅ 1e-4` row is stale for current production
regardless of W2 — flagged here, amendment left for the morning.

---

## 5. Fleet operations

**Timing.** All 66 slabs landed by **04:12**, 8.1 GB total. The handoff
estimated 3–3.5 h per family; the *fast* nodes did it in **56 min**, because
the estimate came from a benchmark on the oversubscribed lxplus *login* node.
But the fleet is enormously heterogeneous:

| `wall_s` | min | p25 | median | p75 | max |
|---|---|---|---|---|---|
| seconds | 3376 | 4676 | 5492 | 8296 | **20724** |

That is a **6.1× spread at identical N and identical `gd_rank`** — pure node
heterogeneity, since every job solves the same size problem and the CPU-second
totals scale with the wall time (slower cores need proportionally more
CPU-seconds for the same flops). 487 core-hours total; 122 h if run serially,
compressed into a ~6 h span by the fan-out.

The tail dominates the schedule: the slowest five were `x001` (5.8 h), `y008`
(3.7 h), `y011`, `y004`, `y001` (~3.2–3.3 h each). `x001` alone held the
completion gate for **two hours** after 65/66 were down — worth knowing when
sizing the ny=1024 follow-up, where the same tail would be doubled. The
practical lesson for a hard-gated fleet is that the median tells you almost
nothing; plan against p95.

**A wrong diagnosis, made cheap by a reversible action.** `y000` appeared
wedged: `RemoteUserCpu` *and* `RemoteSysCpu` frozen at 24105.0/2567.0 across
four reads spanning 13 min (two of them 100 s apart) while wall climbed past
7900 s. A running `dsyevd` must accumulate user CPU, so this looked decisive.
**It was wrong** — at 23:58 the counters jumped to 27291.0/2872.0, i.e. +3186 s
of CPU across 960 s of wall (3.3 of 4 threads), meaning the job had been
computing continuously and condor's ClassAd simply had not updated for ~13 min.

The action taken was to **race a duplicate rather than `condor_rm`**, on the
grounds that a kill is irreversible and could not be undone if the original was
minutes from done, while a duplicate costs one slot and cannot lose (same code,
same inputs, idempotent `xrdcp -f`). That choice turned a wrong diagnosis into
one idle slot instead of ~2.2 h of destroyed progress.

**Final wedge count: 0 of 66.** An earlier "1-in-66 wedge rate" claim is
withdrawn. The standing lesson is not that jobs wedge, but that **neither
`MemoryUsage` nor `RemoteUserCpu` is trustworthy below ~20 min of sampling**;
`MemoryUsage` is worse than useless here, reading exactly 24415 for running,
wedged-looking and successfully completed jobs alike. Prefer a non-condor
signal (an actual output artefact) before any irreversible action.

---

## 6. Two full-scale certifications that were not planned

**`gd_rank` = 6200 is an exact geometric identity.** Every landed slab reports
it. In cells, the snapped pad is **67 x-cells (10 µm) × 7 y-cells (97.5 µm)**
out of **78 × 8** per pitch, so

    gap fraction = 1 − (67·7)/(78·8) = 155/624 = 0.2483974358974359
    N × gap fraction = 24960 × 155/624 = 6200   exactly (integer, verified in
                                                rational arithmetic)

The Gd pseudo-inverse rank equals the number of gap cells in the family-reduced
cell space — the theoretical rank of `R diag(D) R` for a hard-mask projection.
This certifies at production N, across 34+ × 24 960-mode eigendecompositions,
that (a) every job built the identical snapped mask geometry (**670 × 682.5 µm**
is the discrete gap, if anyone asks) and (b) `pinv_rtol = 1e-10` cuts at exactly
the geometric rank, with zero spurious kept or dropped modes. This was the one
production-scale numerical unknown left open in the solver battery, which the
laptop cannot reach; the slab survey closed it as a by-product.

**`x_fraction_prompt` is pinned at 0.5 by a symmetry, not merely "≈0.5".**
Translating by one pad pitch in x maps every X-owned pad onto a Y-owned pad
while leaving the metal/gap mask invariant (the mask is pad-lattice periodic;
only ownership parity flips), and capture is translation-invariant. So the
110 µm vs 97.5 µm gap anisotropy — which is real, and follows from the same
cell counts as `gd_rank` — **cannot** reach this observable at t=0. This also
explains W1's measured 0.49999999658 more fundamentally than the
fractional-drive argument did. Genuine X/Y asymmetry can appear at *late* time,
where the 800 µm ESL period breaks the 780 µm translation symmetry once the
sheet conducts (W1's `x_fraction_late` ≈ 0.49/0.51).

Consequence: the earlier guidance to tolerate a ~1 % x_fraction deviation was
**withdrawn before the products landed**, and the acceptance bar tightened from
0.48–0.52 to 0.5 ± 1e-4. A 1 % deviation is now a real inconsistency.

**The symmetry is visible in the production slabs themselves, at machine
precision.** The per-drive full-grid captures recorded in the slab metas
decompose the two views in completely different bases — 40 X phase drives
against 64 Y row-combs — and they agree:

    X view: 40 × 0.010524708864218921 = 0.420988355
    Y view: 64 × 0.006577943040138035 = 0.420988355
    |X − Y| = 7.7e-14
    total = 0.841976709            (pre-registered: 0.841977)

So the view balance is not merely close to 0.5, it is 0.5 to ~1e-13 in the
production data, from two unrelated drive decompositions. Measured before the
combine ran.

### A predictor that was right for reasons unavailable at the time

An early attempt to forecast the capture summed each slab's `prompt_capture`
over families. The arithmetic looked right — 128 × the per-drive 0.006577943 =
0.8419767, inside the expected band. That agreement is what prompted reading
`solve_family`, which showed `prompt_capture` is `v0.mean()` from
`solve_prompt_cg(vpad)`: **full-grid and `ifam`-independent**, hence a per-box
per-drive constant every family recomputes. The predictor was deleted rather
than kept with a caveat.

It was later shown to be arithmetically correct after all: the metal decomposes
into **128 identical row-combs** (64 per view), each a pure translate of the
others under the same symmetry as above, so superposition is exact and
0.006577943 × 128 = 0.841977. Deleting it was still right at the time — a
predictor that agrees with the expected answer for reasons one cannot state is
a coincidence not yet disproved, not evidence. The justification arrived
afterwards, from a measurement.

---

## 7. What this does not establish

- **It does not certify the digitizer's fast path.** T10 still fails at 7.55 %
  against a 2 % bar. A W2 kernel is a better kernel; it is not a substitute for
  §7 step 5.
- **It does not make the W2 caching cert meaningful.** The 1.87 % reading is a
  harness artifact (§4) and is not evidence either way about these products.
  No caching claim should be quoted until the harness is fixed and re-run.
- **It does not validate absolute amplitude against data.** Everything here is
  internal consistency — sum rules, symmetries, rank identities, and agreement
  between two code paths. The claim that W2 is *right* rests on V6's physics
  argument, not on a measurement. T14 remains the only test that can fail this
  against reality, and it has deliberately not been run.
- **It does not make `s1_w2_ny512` strictly better than `s1_ny1024`.** W2 fixes
  a 27 % absolute error and re-introduces the ~0.45 % ny=512 pad-edge shoulder
  term. For shape-sensitive shallow-deposit work the finer W1 grid is still the
  better instrument, and the W2-at-ny=1024 upgrade is not done.
- **The T10 improvement is not a sharp result.** 8.26 % → 7.55 % is larger than
  the grid and stack terms, but ρ_s alone moves the same number by 0.65 pp. The
  robust statement is the FAIL, not the size of the gain.
- **The 0.40 % late-capture rise is recorded, not explained.** It is uniform
  across ρ_s and absent in W1, and nothing here says why.
- **Nothing downstream has been rebuilt on W2.** The digitizer LUT, Stage B/C
  and T13 still run on W1 products; §6 of the handoff was not reached (see §8).

## 8. Deferred to daylight

- **Handoff §6 (Stage B/C regeneration + T13) was NOT reached.** The fleet's
  straggler tail consumed the headroom that the fast families had bought:
  `x001` alone held the completion gate for two hours after 65 of 66 slabs were
  down, moving combine from ~02:00 to 04:13 and the certs to ~05:50. Since
  nothing downstream has been rebuilt on W2, the digitizer LUT, Stage B/C and
  T13 all still run on W1 products — i.e. on kernels ~27 % low in absolute
  amplitude. That is the first thing to fix, and it is now unblocked.
  Note also that the W2 caching cert cannot certify the rebuilt LUT until the
  §4 harness defect is fixed, so the two follow-ups are coupled: fix the
  harness, then rebuild and re-cert.

- `_gather()` re-reads each slab per (drive, ρ): the X stage alone does ~640 GB
  of repeated decompression against ~50 GB for Y. Accepted tonight rather than
  optimised, because the edit risks a silent indexing slip in `slab[di, ri]`
  and would need the identity battery re-run unattended. **Deferred to the
  solver author**, with battery test 8 as the gate.
- **The LUT cert harness fix (§4) — owned by the solver author**, who will take
  the nearest-at-or-below index (or interpolate the reference onto `lut.t[-1]`),
  re-cert against the W2 products, and amend the plan's T10 row. The story there
  is now complete: the `1e-4` was **true when measured**, was invalidated the
  same week by an unrelated and **correct** fix (Audit A1/Fix 1 raising `t_max`
  1000 → 3000 ns), was surfaced by the W1 dress rehearsal, and had its mechanism
  pinned by three successive falsifications (glue stack, `y_stride`, then the
  time-axis snap).
  Being folded in as a **guard, not just a fix**: "snap-to-nearest across two
  mismatched time axes" is the same bug family as the s-vs-ns accident that
  voided the 2026-08-08 19:12 T10 run. This repo now has two *measured*
  instances of axes that silently disagree, which is enough to justify a
  standing check rather than a one-off correction.
- A W2 grid at ny=1024 — doubles the family to 49 920 modes (~20 GB matrices)
  and needs the y-parity/mirror-family reduction first.
- `test_time_grid` on W2 (§1.4), which §4 makes more interesting than it was.
