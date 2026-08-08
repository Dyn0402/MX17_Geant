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

## Acceptance measured on the products

(filled in by the combine run — see `design/report/W2_NIGHT_REPORT_2026-08-09.md`)

## Reading these files

Never read a multi-GB product through the EOS fuse mount; `xrdcp` it to local
scratch. The tiny provenance/acceptance block can be read on its own without
the arrays — `scripts/condor/extract_meta.py`. Note that `/eos/experiment`
redirects to `eosexperiment.cern.ch`, which **does not resolve off-site**: the
laptop cannot xrdcp these files at all, so combine/cert work runs at CERN and
only small products travel by ssh.
