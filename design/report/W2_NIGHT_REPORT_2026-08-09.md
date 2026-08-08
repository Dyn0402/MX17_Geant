# W2 night shift, 2026-08-08 → 09

**Status: IN PROGRESS — written as the night ran. Sections marked PENDING were
still running when this line was last saved.** Handoff being executed:
`design/report/HANDOFF_W2_OVERNIGHT_2026-08-08.md`.

---

## Verdict

PENDING (combine + cert).

## Headline numbers

| quantity | pre-registered | measured | verdict |
|---|---|---|---|
| `channel_capture_prompt` | **0.841977** | PENDING | |
| `channel_capture_late` | flat, ≈ prompt | PENDING | |
| `x_fraction_prompt` | **0.5 to 1e-4** | PENDING | |
| `x_fraction_late` | record only | PENDING | |
| T10 slow path, W2 rho2M | vs 8.26 % (W1) | PENDING | |
| LUT caching cert, W2 | vs the W1 story below | PENDING | |

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

**Hypothesis B — `y_stride` decimation.** Commit `e36a39c` "LUT: decimate y like
x" landed 2026-08-08 14:00, *after* the certification, because the ny=1024 LUT
hit ~10 GB and was OOM-killed. `kernel_lut.py` sets
`y_stride = max(1, round(tgt/dy_um))`, so an ny=512 product keeps **stride 1**
(previous behaviour) while ny=1024 gets **stride 2** — and both products tested
above are ny=1024. Status: PENDING (cluster 13353186 runs the identical cert on
an ny=512 product).

**⚠ Reading note for the W2 cert, which is the easiest misattribution
available.** The W2 products are **ny=512 → stride 1 → no decimation**. If the
W2 cert returns near 1e-4 that is **not** evidence that W2 fixed the fast path;
it may be nothing but stride 1 vs stride 2. W2 must be compared against
whichever W1 story the ny=512 test selects.

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

**Timing.** The handoff estimated 3–3.5 h per family; measured `wall_s` is
~3400 s (**~56 min**), because the estimate came from a benchmark on the
oversubscribed lxplus *login* node. Wall times span 3376–7327 s — a **2.2×
spread at identical N and identical `gd_rank`**, i.e. node heterogeneity
(AVX2 vs AVX-512).

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

PENDING — to be completed with the results.

## 8. Deferred to daylight

- `_gather()` re-reads each slab per (drive, ρ): the X stage alone does ~640 GB
  of repeated decompression against ~50 GB for Y. Accepted tonight rather than
  optimised, because the edit risks a silent indexing slip in `slab[di, ri]`
  and would need the identity battery re-run unattended. **Deferred to the
  solver author**, with battery test 8 as the gate.
- The LUT caching-cert investigation (§4) and the plan's stale `1e-4` row.
- A W2 grid at ny=1024 — doubles the family to 49 920 modes (~20 GB matrices)
  and needs the y-parity/mirror-family reduction first.
- `test_time_grid` on W2 (§1.4), which §4 makes more interesting than it was.
