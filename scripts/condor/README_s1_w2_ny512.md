# `response_sim/s1_w2_ny512/` — the W2 (exposed inter-pad gap) S1 grid

Produced 2026-08-08/09 by `response/solver/w2_production.py` on lxplus condor
(66 Bloch-family jobs + 1 combine job; `run_w2_family.sh`, `run_w2_combine.sh`).
Layout:

    slabs/       66 per-family mode-coefficient slabs (intermediate; keep for
                 re-combining, not for analysis)
    slabs_mini/  ny=256 smoke-job output. NOT production. Same basenames as
                 slabs/ — never mix the two directories.
    products/    the 4 greens_comb_w2_rho{0.5,1,2,5}M_dk50um_g19um.npz

## Why this grid exists

`design/report/V6_PAD_GAPS_2026-08-08.md`: W1 grounds the 100 µm inter-pad
channel, which is not copper. That removes 24 % of the pad plane's drive and
makes **every W1 absolute amplitude ~27 % low** (prompt capture 0.67003 →
0.85194), as well as inventing a 4.5× signal swing across one pad cell that the
real board does not have (W2: 1.17×). W2 exposes the gap over the real stackup.

## THREE grids now exist and they are NOT interchangeable

| | `s1/` | `s1_ny1024/` | `s1_w2_ny512/` |
|---|---|---|---|
| boundary model | W1 (gap grounded) | W1 (gap grounded) | **W2 (gap exposed)** |
| ny | 512 (97.5 µm) | 1024 (48.8 µm) | 512 (97.5 µm) |
| absolute amplitude | **~27 % low** | **~27 % low** | **correct** |
| prompt pad-edge error after 150 µm smear | 0.452 % | reference | ~0.45 % (inherits the ny=512 term) |
| sub-pad amplitude swing | 4.51× (spurious) | 4.51× (spurious) | 1.17× |
| `meta.sharing["in-gap"]` | MISLABELLED (pre-Fix 2) | correct | correct |
| use for | reproducing pre-2026-08-07 results | shape/ratio work needing the fine grid; reproducing pre-W2 results | **everything where absolute amplitude matters** |

So the W1 products are **superseded in absolute amplitude and retained for
shape/ratio reproduction**. They are not deleted and not simply obsolete: read
this table before picking a directory. In particular `s1_w2_ny512/` reintroduces
the ~0.45 % ny=512 pad-edge shoulder term that `test_ny_grid.py` (audit C6)
rejected for W1 at a 0.3 % bar, which is exactly why `s1_ny1024/` was built.
Trading a 27 % absolute error for a 0.45 % shoulder error is a large net win,
but if a result depends on the shoulder and not on the absolute scale, the
ny=1024 W1 grid is still the finer instrument.

A W2 grid at ny = 1024 is a recorded follow-up, not an oversight: it doubles the
Bloch family to 49 920 modes (~20 GB matrices) and needs the y-parity /
mirror-family reduction first.

## Acceptance measured on the products (2026-08-09)

Identical across all four ρ_s, as expected — ρ_s rescales the relaxation rate
but not the prompt electrostatics:

| quantity | value | bar |
|---|---|---|
| `channel_capture_prompt` | **0.841977** | pre-registered 0.841977 (independent full-grid CG) |
| `channel_capture_late` | 0.845369–0.845381 | +0.40 % vs prompt — a rise, not a decay |
| `x_fraction_prompt` | **0.50000002** | 0.5 to 1e-4 |
| `x_fraction_late` | 0.5 to 1e-8 | recorded |
| vs W1 `s1_ny1024` | **+25.66 %** | V6 static +27.2 % × 0.987 ny-grid = +25.5 % predicted |
| `gd_rank` (every family) | 6200 | = N × 155/624 exactly |
| X `view_total_prompt` | 0.880258 (W1: 0.877482) | ~+0.2 % expected from V6 `pad_split` |

The prompt-capture agreement is a cross-check of the **family assembly** against
a **full-grid CG** solve — different code paths end to end — not a tautology.
Full detail, including two expectations that turned out wrong, in
`design/report/W2_NIGHT_REPORT_2026-08-09.md`.

## What these products do NOT fix

**T10 still fails on W2.** The slow-path (ion lateral shape) residual moves
8.26 % → 7.55 % at ρ_s = 2 MΩ/sq against a 2 % bar. Using a W2 kernel does not
certify the digitizer's fast path; plan §7 step 5 (build the LUT from slow-path
templates) remains required.

**The LUT caching cert reads ~1.87 % on these products, and that is a test
harness artifact, not a property of W2** — `test_lut_vs_solver` compares the
solver at 3101 ns against the LUT at 3000 ns. It reads the same on every W1
product too. Do not use that number to judge these products until the harness
is fixed.

## Reading these files

Never read a multi-GB product through the EOS fuse mount; `xrdcp` it to local
scratch. The tiny provenance/acceptance block can be read on its own without
the arrays — `scripts/condor/extract_meta.py`. Note that `/eos/experiment`
redirects to `eosexperiment.cern.ch`, which **does not resolve off-site**: the
laptop cannot xrdcp these files at all, so combine/cert work runs at CERN and
only small products travel by ssh.
