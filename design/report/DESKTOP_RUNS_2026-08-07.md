# Audit fixes — the runs that needed real compute (2026-08-07)

Companion to `AUDIT_FIXES_2026-08-07.md`. Everything here was run on the **desktop**
(16 cores, 62 GB) or streamed from **EOS**, per plan §10. Results only; the code changes are
in the commits each section names.

**lxplus SSH is not available from either machine** — `dneff@lxplus` is refused
(`publickey,keyboard-interactive`; GSSAPI is offered from the desktop but rejected, and
lxplus926 does not offer it at all). Everything EOS-side was done instead with
`xrdcp root://eospublic.cern.ch/…`, which works with a live Kerberos ticket. `eosuser` and
`eosntof` both give "Invalid address".

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

**Not acted on yet, and here is the blocker:** a full 12-point grid at ny = 1024 needs ~27 GB
and the desktop has 28 GB free on `/media/ucla`, 11 GB on `/`, 7 GB on `/media/dylan/data`.
This wants either a disk decision or a one-point-at-a-time regeneration.

---

## Fix 7 — avalanche survival: there is no bias to correct

A7 said the Polya is fitted conditional on g > 0, the per-slice `survival` is dropped in the
merge, and per-electron charge is therefore biased high by 1/P(g>0) "by an amount unknowable
from the shipped calib".

The first two are true. **The third makes no difference, because P(g>0) = 1.** Streamed from
EOS (`xrdcp`, ~264 MB per slice, parsed and deleted one at a time):

| V | survival | seeds | zero-gain seeds |
|---|---|---|---|
| 470 | 1.0 | 1600 | **0** |
| 480 | 1.0 | 1360 | **0** |
| 490 | 1.0 | 1120 | **0** |
| 500–530 | *sweep still running at time of writing* | | |

Not one seed electron out of 4080 failed to multiply, at the three LOWEST voltages — and
survival can only rise with field. So `survival × mean(conditional Polya)` and the direct
mean over all seeds are the same number, the two routes the acceptance asks to compare agree
trivially, and the shipped calib is unbiased. **Carrying `survival` into a v3 calib is still
worth doing as bookkeeping, but it will not move the P2 saturation fraction.**

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

## Still not done

- **The ρ_s × d_k scan** (12 points) has not been re-run — only the nominal point. Doing it
  needs the disk decision in C6 and ~2 h.
- **ny = 1024 production grid**, per C6.
- **`aval_calib_v3.json`** carrying `survival` (bookkeeping only, per above).
- C7 (amp-gap clusters), C8 (Fano), C9 (truth block), C10 (empty events), C14 (Stage A C++) —
  all still open, all laptop-scale, none blocked.
