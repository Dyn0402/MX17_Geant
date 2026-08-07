# S1 grid on lxplus condor

`s1_ny1024.sub` + `run_point.sh` re-solve the 12-point ρ_s × d_k grid at
**ny = 1024**, which `response/digitizer/test_ny_grid.py` (audit C6) showed is
required: at ny = 512 the prompt kernel's pad-edge shoulder is off by 0.452 %
after the minimum realistic 150 µm transverse smear, against a 0.3 % bar.

Output goes to EOS `response_sim/s1_ny1024/`, not AFS — a point is ~2.3 GB and
the AFS work quota is 50 GB total.

## Connecting

`ssh lxplus` works; `ssh dneff@lxplus.cern.ch` does NOT. The alias in
`~/.ssh/config` sets `GSSAPITrustDns yes`, which is what makes GSSAPI succeed
against the round-robin alias — without it the service principal does not match
the lxplusNNN node you actually land on, and every auth method is refused.

## Traps

* `MY.SendCredential = true` is required or the job cannot read AFS or `xrdcp`
  to EOS, and it fails at the very END, after a full solve.
* `OMP_NUM_THREADS` is pinned to the requested core count. The solve is
  BLAS-heavy and OpenBLAS otherwise spawns threads for cores condor did not
  give it, which is slower than single-threaded.
* Memory: `--quick` (ny = 256, nt = 12) peaks at 0.6 GB; ny = 1024 with nt = 61
  is ~20× that, hence the 24 GB request.

    ssh lxplus
    cd /afs/cern.ch/work/d/dneff/mx17_s1 && condor_submit s1_ny1024.sub
    condor_q -nobatch
