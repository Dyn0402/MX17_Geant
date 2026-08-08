# S1 grid on lxplus condor

`s1_ny1024.sub` + `run_point.sh` re-solve the 12-point ρ_s × d_k grid at
**ny = 1024**, which `response/digitizer/test_ny_grid.py` (audit C6) showed is
required: at ny = 512 the prompt kernel's pad-edge shoulder is off by 0.452 %
after the minimum realistic 150 µm transverse smear, against a 0.3 % bar.

Output goes to EOS `response_sim/s1_ny1024/`, not AFS — a point is ~2.3 GB and
the AFS work quota is 50 GB total.

## BOTH grids are kept. Which one to use.

`response_sim/s1/` (ny = 512) is **not** retired — decision 2026-08-07. Two
directories now exist and they are not interchangeable:

| | `s1/` | `s1_ny1024/` |
|---|---|---|
| ny | 512 (97.5 µm) | 1024 (48.8 µm) |
| prompt pad-edge error after 150 µm smear | 0.452 % | reference |
| `meta.sharing["in-gap"]` | **MISLABELLED** (pre-Fix 2: selected by absolute x, so it is a second *on-strip* deposit) | correct, phase-selected |
| use for | reproducing pre-2026-08-07 results | **everything new** |

So: read the sharing block ONLY from `s1_ny1024/`, and treat any `s1/` number
that depends on the transverse grid as carrying a ~0.45 % systematic at the
shallowest depths. The G arrays in `s1/` are otherwise fine — the Fix 2 bug was
in the *reporting*, never in the kernels.

## Connecting

`ssh lxplus` works; `ssh dneff@lxplus.cern.ch` does NOT. The alias in
`~/.ssh/config` sets `GSSAPITrustDns yes`, which is what makes GSSAPI succeed
against the round-robin alias — without it the service principal does not match
the lxplusNNN node you actually land on, and every auth method is refused.

## Traps

* `MY.SendCredential = true` is required or the job cannot read AFS or `xrdcp`
  to EOS, and it fails at the very END, after a full solve.
* **LCG setup.sh reaches into the caller's shell — two measured traps
  (2026-08-08/09, both killed a job instantly):** (1) it is not `set -u`-clean
  (unbound `COMPILER` at line 18) — wrap the source in `set +u` … `set -u`;
  (2) it EXPORTS `BASE` (→ a CVMFS gcc path), clobbering any `BASE` you set
  before sourcing. Of `BASE SRC WORK OUT SLABS OUTD CALIB PROD EOSDIR EOSBASE
  RHOTAG`, `BASE` is the ONLY name LCG_105 clobbers (measured, not assumed —
  see `run_w2_cert.sh` header). Don't name a variable `BASE` in any script
  that sources LCG, and don't rename the safe ones on suspicion.
* `stream_output`/`stream_error` are no longer supported (CERN, Nov 2025) —
  submission is rejected; job stdout reaches AFS only at job END, so
  `condor_q -af MemoryUsage` is the only mid-run health signal.
* Off-site, `/eos/experiment` is unreachable by xrootd: it redirects to
  `eosexperiment.cern.ch`, which does not resolve outside CERN (the error is
  an unhelpful `[FATAL] Invalid address`). Anything that must read those files
  runs at CERN; only small products travel, by ssh/rsync or via an AFS copy.
* `OMP_NUM_THREADS` is pinned to the requested core count. The solve is
  BLAS-heavy and OpenBLAS otherwise spawns threads for cores condor did not
  give it, which is slower than single-threaded.
* Memory: `--quick` (ny = 256, nt = 12) peaks at 0.6 GB; ny = 1024 with nt = 61
  is ~20× that, hence the 24 GB request.

    ssh lxplus
    cd /afs/cern.ch/work/d/dneff/mx17_s1 && condor_submit s1_ny1024.sub
    condor_q -nobatch
