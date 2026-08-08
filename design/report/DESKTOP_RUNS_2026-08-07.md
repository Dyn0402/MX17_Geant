# Audit fixes — the runs that needed real compute (2026-08-07)

Companion to `AUDIT_FIXES_2026-08-07.md`. Everything here was run on the **desktop**
(16 cores, 62 GB) or streamed from **EOS**, per plan §10. Results only; the code changes are
in the commits each section names.

**Connecting to lxplus: use the `lxplus` ALIAS, never the hostname.** `ssh lxplus` works;
`ssh dneff@lxplus.cern.ch` is refused from both machines even with a live Kerberos ticket.
The alias in `~/.ssh/config` sets `GSSAPITrustDns yes`, and that is the whole difference —
`lxplus.cern.ch` is a round-robin, so without DNS trust the service principal does not match
the `lxplusNNN` node you actually land on and every auth method fails. Several hours of this
report's EOS work was done via `xrdcp root://eospublic.cern.ch/…` before that was understood;
that route also works and needs only a ticket (`eosuser` and `eosntof` give
"Invalid address").

---

## Fix 1 — the LUT window. The finding is NOT what the audit predicted.

A1 said: truncating the kernel at 1 µs starves the late halo, "biasing the halo area shares
LOW, **the same direction as the headline §9 tension**", and predicted "halo shares rise"
after the rebuild.

**They fall.** Measured on 400 real muons (`mx17_muons_3k_spread_t0.root`), ρ_s = 1 MΩ/sq,
d_k = 75 µm, S1 re-solved at production resolution:

| config | c1_X | c2_X | c1_Y | c2_Y | X/Y balance |
|---|---|---|---|---|---|
| n_side = 4, t_max = 1000 (PRE-FIX) | 0.2073 | 0.0431 | 0.1516 | 0.1099 | 0.514 / 0.486 |
| n_side = 4, t_max = 3000 | 0.2074 | 0.0431 | 0.1282 | 0.0963 | 0.581 / 0.419 |
| n_side = 8, t_max = 1000 | 0.2073 | 0.0431 | 0.1451 | 0.1053 | 0.503 / 0.497 |
| **n_side = 8, t_max = 3000 (POST-FIX)** | **0.2074** | **0.0431** | **0.0981** | **0.0737** | **0.515 / 0.485** |
| *measured §9* | *0.23–0.28* | | *0.23–0.28* | | *0.49 / 0.51* |

Three things fall out, and only the third was anticipated:

1. **The X view is completely insensitive to the window** — c1_X moves by +0.1 % over a 3×
   change in t_max, and the peak profile is identical to four decimals. Whatever the X-view
   halo tension is, the LUT window is not its cause.
2. **The Y view moves a lot, and DOWNWARD.** c1_Y falls 15 % (n_side 4) to 32 % (n_side 8).
   §9 measures c1 = 0.23–0.28, so the corrected prediction sits *further* from the data than
   the buggy one did. Recorded as a prediction revision; per the §9 firewall it is **not** to
   be undone by tuning.
3. **The window fix is inseparable from the channel-window fix** (the fix order's "sizing
   trap", and the one thing it got exactly right).

### The sizing trap, quantified — this is the important result

`test_charge_audit`, production kernel, leak = 1 − measured/expected with the expectation
evaluated per position (Fix 4):

| | z = 0.5 mm | 5 mm | 15 mm | 29.5 mm |
|---|---|---|---|---|
| n_side = 4, t_max = 1000 | +0.24 % | +0.92 % | +1.43 % | +1.78 % |
| n_side = 4, t_max = 3000 | **+9.68 %** | **+9.46 %** | **+9.09 %** | **+8.00 %** |
| n_side = 8, t_max = 1000 | −1.75 % | −1.03 % | −0.53 % | −0.20 % |
| **n_side = 8, t_max = 3000** | **−0.04 %** | **+0.44 %** | **+0.60 %** | **+0.51 %** |

Opening the time window at a fixed ±4 channels **loses 8–10 % of the induced charge** — the
extra time is spent spreading charge past the edge of a window that was sized for a 1 µs
kernel. n_side = 8 is the only configuration that holds every depth inside the 1 % bar.

**Production defaults are therefore now `t_max_ns = 3000` AND `n_side = 8`**
(`y_window_mm = 7.02` with it — they must grow together, since Y channels sit on the 0.78 mm
pad pitch). Shipping the window fix alone would have made the charge accounting five times
worse than the bug it fixed.

### A caveat that limits how far any of this can be pushed

The c1 numbers above integrate the whole simulated record. `event_budget`'s own docstring is
explicit that the area budget is comparable to data **only over the DAQ's 1.92 µs window**.
Re-run at `--n-samp 1920`, the honest before/after is smaller:

| at the DAQ window | c1_X | c2_X | c1_Y | c2_Y | X/Y balance |
|---|---|---|---|---|---|
| PRE (n_side 4, t 1000) | 0.2074 | 0.0432 | 0.1518 | 0.1101 | 0.514 / 0.486 |
| POST (n_side 8, t 3000) | 0.2074 | 0.0432 | 0.1344 | 0.0989 | 0.504 / 0.496 |
| change | 0.0 % | −0.0 % | **−11.5 %** | −10.2 % | → closer to 0.49/0.51 |

Which makes sense in hindsight: the kernel window only ever mattered between 1.0 and 1.92 µs,
because past 1.92 µs the DAQ frame discards the charge anyway. The X/Y balance does improve,
0.514 → 0.504 against a measured 0.49/0.51.

### The second half of A1 does not hold up

The fix order asked for the two longitudinal convolutions to be padded before truncating, on
the theory that they were losing their own tail. **They are not.** The running mean at sample
n reads only x[n−L+1..n], so every output sample inside the window is already complete, and
padding is bit-for-bit identical (verified directly); `apply_longitudinal`'s FFT is already
sized nt + len(h) − 1, so it is a linear convolution, not a wrapped one. What *is* lost is
output beyond nt, which is charge outside t_max by construction and only a longer window
recovers. A1 has one cause, not two. `test_longitudinal` gained the tail assertion that was
missing (impulse in the last transit-length: kept area 1.0000 / 0.5462 / 0.0950 against the
exact geometric fractions).

**T10 still passes** on the production kernel at the new window, and the off-grid booking
check (C2) reports 0 disagreements with the pre-fix rule differing on **1.3 %** of uniform
probes — reproducing the audit's estimate exactly.

---

## Fix 2 — production sharing block

Re-solved at production resolution (nx = 3120, ny = 512, nt = 60; 90 s on the desktop). The
regenerated table is now in plan §3. Headline: the genuinely-in-gap deposit keeps
**0.794** of its charge in the owning view at 10 µs against **0.366** for the on-strip one,
and its cross-view late feed is ~6× smaller (0.024 vs 0.143). The two rows are plainly
different objects where the mislabelled pair agreed to three digits.

⚠️ The 12 shipped products' `meta.sharing["in-gap"]` blocks predate the fix and remain
mislabelled.

---

## C6 — the ny grid does NOT pass at ny = 512

Nominal point re-solved at ny = 1024 and differenced against ny = 512
(`response/digitizer/test_ny_grid.py`):

| | max &#124;diff&#124; | vs prompt peak | bar |
|---|---|---|---|
| raw (pad edge) | 0.0162 | 1.86 % | — |
| smeared 150 µm | 0.00393 | **0.452 %** | 0.3 % → **FAIL** |
| smeared 300 µm | 0.00224 | 0.257 % | OK |
| smeared 800 µm | 0.00020 | 0.024 % | OK |

150 µm is the floor an avalanche gets (σ₀ ≈ 34 µm plus σ_T ≈ 107 µm at the shallowest
depth), and there the ny = 512 grid contributes 0.45 % — above the plan's own 0.3 % bar. By
the plan's stated rule, **production needs ny = 1024**. Cost is ~2× per point (26 → 52 s per
Y kernel) and ~2× disk.

**Acted on: the 12-point grid is running on lxplus condor**, cluster 13352202, writing to EOS
`response_sim/s1_ny1024/` — the desktop could not hold it (28 GB free on `/media/ucla`
against ~27 GB needed), EOS has 39 PB. Submission in `scripts/condor/`.

**The ny = 512 grid is kept, not retired** (decision 2026-08-07). The two directories are not
interchangeable and the difference is recorded in `scripts/condor/README.md`: `s1/` still
carries the pre-Fix-2 mislabelled `meta.sharing["in-gap"]` block and the 0.45 % transverse
grid systematic, so new work reads `s1_ny1024/`. The G arrays in `s1/` are otherwise sound —
Fix 2 was a reporting bug, never a kernel one.

---

## Fix 7 — avalanche survival: there is no bias to correct

A7 said the Polya is fitted conditional on g > 0, the per-slice `survival` is dropped in the
merge, and per-electron charge is therefore biased high by 1/P(g>0) "by an amount unknowable
from the shipped calib".

The first two are true. **The third makes no difference, because P(g>0) = 1.** Streamed from
EOS (`xrdcp`, ~264 MB per slice, parsed and deleted one at a time):

**Complete: all 56 slices, 0 download failures** (`design/report/aval_survival_2026-08-07.json`):

| V | survival | seeds | zero-gain seeds | mean gain |
|---|---|---|---|---|
| 470 | **1.0000** | 1600 | 0 | 23 107 |
| 480 | **1.0000** | 1360 | 0 | 31 387 |
| 490 | **1.0000** | 1120 | 0 | 44 472 |
| 500 | **1.0000** | 880 | 0 | 60 309 |
| 510 | **1.0000** | 640 | 0 | 81 645 |
| 520 | **1.0000** | 480 | 0 | 112 447 |
| 530 | **1.0000** | 320 | 0 | 156 731 |
| **total** | | **6400** | **0** | |

Not one seed electron in 6400 failed to multiply, at any voltage. So
`survival × E[g | g>0]` and the direct mean over all seeds are the same number — the two
routes the acceptance asks to compare agree identically, not approximately — and the shipped
calib was never biased. **A7's 1/P(g>0) correction is exactly 1, and this does not move the
P2 saturation fraction.**

Two independent consistency checks fall out of the same numbers: the gain at the 490 V bench
point is 4.45 × 10⁴, matching plan §5's "≈4.4 × 10⁴"; and gain rises ×6.78 over the 60 V
span, an e-folding every 31.3 V, against §5's "×6.8 over 60 V, e-folding every 31 V".

Implemented anyway as bookkeeping, because "it happens to be 1 here" and "it is 1 by
construction" are different statements: `collect.py` now carries `survival`, its binomial
error and the unconditional mean into an `aval_calib/3` calib and **asserts the
decomposition identity** at merge time; `digitize.py` folds it into the same single thinning
as mesh transparency and attachment, so a calib with survival = 1 is bit-identical (verified:
same 1770 survivors from the same seed) and a future non-unit value cannot be silently
ignored.

Also confirmed from the raw config blocks: `"max_avalanche": 0` — the "avalanche size cap"
of plan §5 does not exist, exactly as the audit said.

### Two by-products worth more than the fix

- **`aval_calib_v2.json` exists on EOS** and has populated `i_elec`/`i_ion`. The local
  `aval_calib.json` has those keys present but **all 2000 entries zero**, which silently
  turned the entire LUT to `nan` — `test_induce_equivalence` on HEAD reports
  `max abs difference = 0.000e+00` comparing nan to nan, which reads as a pass. Guarded now.
- **`clusters/` on EOS** holds the Stage A output (`mx17_muons_3k_spread_t0.root`) that all
  the muon numbers above depend on; it existed on neither machine.

---

## C7 — amp-gap deposits, implemented (partial gain)

Chosen over dropping them, because it makes `clusters.py`'s existing contract true and keeps
the ESL groove deposits `NEEDED_INPUTS` calls real. Verified on a synthetic z scan — mean
gain relative to a full-gap avalanche:

| z [mm] | −0.005 | 0.000 | 0.0375 | 0.075 | 0.1125 | 0.1499 |
|---|---|---|---|---|---|---|
| measured | 0.00002 | 0.00002 | 0.00034 | 0.00483 | 0.06937 | 0.98683 |
| expected | 0.00002 | 0.00002 | 0.00034 | 0.00484 | 0.06954 | 0.99292 |

with the surviving fraction exactly 1.0000 at every in-gap z (no mesh crossing) against
0.871–0.875 for drift electrons (transparency 0.873). Groove deposits (z < 0) need no special
case: the same formula gives gain < 1, i.e. no multiplication, which is right below the anode.

## C8 — Fano, implemented in Stage A, and it changes nothing observable

`SteppingAction.cc` now draws n from a truncated Normal(n̄, √(F·n̄)) with F = 0.2, falling
back to the exact floor + Bernoulli below n̄ = 5 where a Gaussian is meaningless. `MX17_FANO`
overrides it, and `MX17_FANO=0` reproduces the old behaviour exactly — which is what makes
the A/B below possible at the **same seed**.

Same seed, 3000 muons: the edep arrays are bit-identical (233 903 clusters), so only the
conversion differs. 26.3 % of clusters take the Gaussian branch.

| | conversion residual sd (Gaussian branch) | per-event electrons: mean | sd | rel |
|---|---|---|---|---|
| F = 0 (old) | 0.4095 | 327.67 | 266.70 | 81.394 % |
| F = 0.2 | **1.5263** | 327.60 | **266.75** | 81.426 % |

So the conversion statistics change exactly as specified — variance up ~14× — and the
per-event observable does not move at all. That is not a null result to be embarrassed by, it
is arithmetic: the conversion term adds √(F·n) = 8.1 electrons in quadrature to a spread of
266.7, predicting 266.82 against 266.75 measured. **The per-event electron count is dominated
by delta rays and energy straggling (81 % relative spread), which Geant4 already models.**

The plan no longer promises something the code does not do, and the honest conclusion is
recorded alongside: Fano was never a candidate explanation for any §9 tension.

C7's effect on the 400-muon table at the DAQ window is sub-percent, as expected for 0.41 % of
the electrons: c1_X −0.1 %, c1_Y −0.7 %, c2_X +2.4 %, X/Y balance +0.1 %.

## C9 / C10 / C14 — outputs, and three more real bugs

**C10.** Zero-electron events were `continue`d before the DAQ buffer, so decoded files omitted
triggers a real run records. Worse than the occupancy problem: decoded eventIds were
`arange(len(block))`, so a dropped event silently **renumbered every later one**. Events now
carry the ClusterTree's own eventID and empty ones are written as noise-only frames. Verified
by forcing every second event empty — 20 events out, 10 noise-only, ids 0–19 contiguous, join
total on both FEUs, empty frames at 231 ADC max excursion (pure common mode) against 1972 for
charged.

**C9.** `truth.py` writes a per-event sidecar keyed by eventId. Two choices worth stating: the
per-channel charge is taken **before** the shaper (with P6, ∫h = (1−β)∫h_nom, so a truth charge
read after shaping would carry β), and the format is **npz not parquet** — the response chain
runs on the system python3, which has no pyarrow or pandas, so a parquet sidecar could not be
written by the process that produces it. Free cross-check: summed per-channel true charge over
seed charge = 0.6835, where the real-pad capture (0.665) plus window and position sampling put
it.

**C14.** Two of the five were real, plus one on the Stage B side:

- **A macro without `beamOn` ran zero events.** `run_default.mac` says in its own header that
  main() launches the run after it, but main() only did that when *no* macro was given — so
  passing the default macro produced an empty file and exited 0. Now 200 events either way.
- **The impact point was randomised twice**, mislabelling provenance and (since C9) making the
  recorded truth differ from what Stage A used. Stage A now records the primary vertex and the
  spread it applied; run.py detects it and says so.
- **Charge off the 512×512 readout** was counted by the budget but dropped by the DAQ writer,
  so the two paths disagreed on out-of-area delta rays. `induce` now applies the fiducial and
  reports the count.
- Meta branch addresses hoisted off dead stack; the EventAction ×1e6 print fixed.

Note the two truth impact points are **different quantities**: Stage A's vertex is where the
particle entered, C9's is the charge-weighted arrival point at the ESL. They agree in the mean
(0.16–0.19 mm) but differ per event by up to 8.2 mm from track inclination and delta rays.
Comparing them and calling the difference an error is the obvious trap; `truth.py` says so.

## `aval_calib_v3.json` — and a provenance gap in the archived raw

v3 is on EOS, carrying v2's measured `i_elec`/`i_ion` **plus** the survival
block, with `survival = 1.0000 ± 0.0000` at all seven voltages and the
decomposition identity asserted at merge time. `collect.py` reproduces the
gains exactly (23107.4 / 31386.7 / 44471.5 / … against the streamed values) at
1.4 GB peak, so the streaming `reduce_file` design holds up on 19 GB.

⚠️ **But v3 is NOT a from-raw regeneration, and must not be made one.** Running
`collect.py` over the archived raw slices reproduces everything except the ion
templates, because **the raw slices' own `i_elec`/`i_ion` are all zeros** — every
voltage, checked directly. So v2's templates came from a different S3 pass and
are **not reproducible from what is archived**. A naive regeneration produces a
file that looks newer and is strictly worse; `digitize.py`'s template guard
(added while re-running the charge audit) refuses it, which is the guard earning
its keep on a file I had just created. `scripts/make_aval_calib_v3.py` does the
correct merge and says why in its docstring.

That archived-raw gap is worth closing separately if the ion template ever needs
re-deriving.

## The x/y sign convention — SETTLED, and the sim was wrong

Settled from the detector rather than from data (user, 2026-08-08): viewed from
the pixel side, X connectors along the bottom and Y up the right, x counting
left→right and y bottom→top, x FEU DREAM 1 at the far left and y DREAM 1 at the
bottom — channel 0 at the low end of each coordinate, exactly what
`mx17_m1_map.csv` encodes.

**`sy` was −1 and should be +1.** The header had flipped y purely to keep
(x, y, z)_local right-handed after the z flip; nothing in the response chain
takes a cross product, so handedness is free while a flipped y mislabels every
Y channel. The frame is now deliberately left-handed. The derivation needs no
agreement on which way is "up": the pads sit at larger world z than the ESL, so
+z_world points toward a viewer on the pixel side, (right, up, toward-viewer) is
right-handed for any viewer, hence (x_det, y_det, z_local) is left-handed, which
forces sx·sy = +1.

Verified on a 300 mm-spread run: physical bottom → row 74, top → row 444, far
left → col 64, far right → col 447. The connector inversion stays on the data
side (`apply_orientation`, already unconditional) and the comparison happens
after it. The bench is a third, swapped frame (bench x = detector y), which is
why the June M3 alignment needed ~90° plus a flip.

⚠️ Cluster files written before 2026-08-08 carry `sy = −1` in their Meta and stay
self-consistent read through it, but their y is mirrored against new files — do
not mix them channel-by-channel. Every observable reported here is symmetric in
y, so no number in this document changes.

## Test suite

`scripts/run_tests.sh <kernel> <calib> [ny1024]` — the whole suite in one place. On the
production kernel: `test_longitudinal`, `test_lut_vs_solver`, `test_charge_audit`,
`test_induce_equivalence`, `test_daq_wft` all **PASS**; `test_ny_grid` reports its
informational FAIL (the C6 finding).

`test_induce_equivalence` had never actually run anywhere — it hardcoded one machine's
absolute paths, and because it compared nan against nan it printed
`max abs difference = 0.000e+00`, which reads as a pass. It now takes `--kernel/--calib`, ends
with a real verdict and exit code, and passes at 3.1e-08.

## Still not done

- **The ρ_s × d_k muon scan** (12 points through Stage B) has not been re-run — only the
  nominal point. The ny=1024 products it should use now exist, so this is just compute.
- **`aval_calib_v3.json`** regenerated from the raw slices with `collect.py` (the code carries
  `survival` now; the file on EOS is still v2, and since survival = 1.0 nothing changes
  numerically).
- **The x/y sign convention** (`ActiveAreaFrame.hh`) is still the one open convention that
  could silently corrupt T14, and the beat phase vs absolute channel number is its only
  in-data check.

Everything else on the C-list is done: C1–C14 inclusive.
