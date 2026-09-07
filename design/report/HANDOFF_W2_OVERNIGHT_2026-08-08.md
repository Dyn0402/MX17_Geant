# Handoff: babysit the W2 condor fleet, then reprocess the response chain

**Written 2026-08-08 ~23:30 by the session that built and submitted the W2 re-solve.
For a fresh model running overnight. Goal: by morning, the W2 kernel products exist
and are validated, the digitizer is rebuilt on them + the T7 pooled calib, T10 has a
W2 verdict, and (stretch) T13 is done. The user is asleep; do not wake them for
anything short of "the whole fleet is dead and I cannot fix it".**

Authority for the physics chain: `design/RESPONSE_SIM_PLAN.md` (§0a for status,
next-in-order item 3 for what W2 is and why). Background you should skim first:
`design/report/V6_PAD_GAPS_2026-08-08.md` (why W1 was 27 % low),
`response/solver/wpot_w2.py` docstring (the W2 physics), `w2_production.py`
docstring (the job structure). Memory `mx17-w2-resolve-2026-08-08` has tonight in
compressed form.

## 0. Environment and access — read before touching anything

* Laptop venv: `~/PycharmProjects/nTof_x17/.venv/bin/python` from
  `~/CLionProjects/MX17_Geant` (system python silently degrades; laptop-venv
  numpy runs reference BLAS — fine for combine, which has no big matmuls).
* lxplus: `ssh lxplus` (alias ONLY — `ssh dneff@lxplus.cern.ch` fails; GSSAPI
  needs the alias's TrustDns). 2FA: auth works only through the live
  ControlMaster socket (`~/.ssh/master-*`, ControlPersist 1d). If it has died,
  GSSAPI returns "partial success" and you CANNOT re-auth without the user —
  fall back to `ssh desktop 'ssh lxplus ...'` (the desktop holds its own live
  master, verified tonight). If both are dead: the condor jobs and EOS uploads
  continue fine without you; just retry every ~30 min and note the gap.
* Kerberos ticket (laptop AND desktop) expires **08-09 11:21**, renewable to
  08-12: run `kinit -R` before big EOS pulls if you are near expiry.
* EOS: never read multi-GB files through the fuse mount — `xrdcp` via
  `root://eosuser.cern.ch/` to local scratch.
* The DESKTOP is running the T7 overnight chain (voltage ladder + HV scan,
  another session's job, `response/avalanche/overnight_chain.log`). Light use
  (ssh, file reads, an rsync) is fine; do NOT start heavy compute there and do
  NOT touch `response/avalanche/`.

## 1. What is in flight right now (submitted 2026-08-08 evening)

66 production jobs on lxplus condor = the complete W2 S1 grid at ny=512:

| cluster | jobs | what |
|---|---|---|
| 13353071 | 2 | Y-box family 0, X-box family 0 (the "canaries" — ordinary production jobs, started first ~22:20) |
| 13353072 | 64 | Y families 1–63, X family 1 (started ~23:0x) |
| 13353074 | 1 | mini smoke job (ny=256) — **DONE, PASSED end-to-end** incl. the credentialed EOS upload. Its slab lives in `slabs_mini/` and is NOT production data |

Each job: one Bloch family, ~3–3.5 h wall on 4 cores, peak ~22 GB (request
30 GB), writes `w2slab_<box><fam:03d>.npz` to EOS
`/eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/slabs/`. Every job is
independent and serves all 42 kernels × all 4 ρ_s. Expected completion
~01:30–02:30. AFS side: `/afs/cern.ch/work/d/dneff/mx17_s1/` (`run_w2_family.sh`,
`w2_canary.sub`, `w2_rest.sub`, `logs/w2_*.{out,err}`, `src/` = code copy).
stdout lands on AFS only at job END (CERN killed stream_output 11/2025) —
while a job runs, `condor_q -af MemoryUsage` is your only health signal.

Three failure classes were already found and fixed tonight (all committed):
LCG `setup.sh` vs `set -u` (`38d23d0`), drive/gap mask tie-break at
half-covered pad-edge cells, fractional constraint masks (now hard-snapped).
The code on AFS is current with laptop commit `w2: --ny smoke override…`.

## 2. Overnight monitoring loop

Poll every ~20–30 min (ScheduleWakeup or a Monitor):

    ssh lxplus 'condor_q -constraint "Owner==\"dneff\"" -af JobStatus | sort | uniq -c;
                ls /eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/slabs/ | wc -l'

* JobStatus 2 = running, 1 = idle, 5 = HELD. Held job → `condor_q -af:j
  HoldReason`; if OOM, resubmit that ONE family with more memory:
  copy `w2_canary.sub`, set the queue list to the failed `box, fam`, raise
  `request_memory` to 36000. If a job vanished without a slab, check
  `logs/w2_*.{out,err}` (they are per-canary named for 13353071; the 13353072
  jobs write `w2_y17.out` style names — actually `w2_$(box)$(fam).out`), fix if
  it is a NEW bug (read the traceback; everything found so far is fixed), and
  resubmit that family the same way. Job success ends with
  `uploaded to /eos/...` in its .out.
* Success criterion: **66 files** `w2slab_y000..y063.npz + w2slab_x000..x001.npz`
  on EOS (ignore `slabs_mini/`).

## 3. When all 66 slabs exist: combine + validate  (~1 h, laptop)

    mkdir -p /tmp/w2slabs && cd /tmp/w2slabs
    for f in $(ssh lxplus 'ls /eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/slabs/'); do
        xrdcp -f root://eosuser.cern.ch//eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/slabs/$f .
    done          # ~10 GB; do it before the 11:21 ticket expiry or kinit -R
    cd ~/CLionProjects/MX17_Geant
    ~/PycharmProjects/nTof_x17/.venv/bin/python -m response.solver.w2_production \
        combine --slabdir /tmp/w2slabs --outdir /tmp/w2products

Acceptance on each of the 4 products (`greens_comb_w2_rho{0.5,1,2,5}M_dk50um_g19um.npz`):

* `channel_capture_prompt` ≈ **0.84** (the −1.3 % ny=512 grid systematic vs
  V6's static 0.852998 — expected, recorded; if you see ~0.67 something read
  W1 masks, STOP).
* `x_fraction_prompt` ≈ 0.5 but NOT exactly — the hard x-mask gives the
  inter-column channel 110 µm vs the inter-row band's 97.5 µm, so a ~1 %
  asymmetry is geometry, not a bug. (W1's "exactly 0.5 REQUIRED" statement was
  for its tiling-symmetric fractional drives — it does not apply here.)
* Sharing report: ratios close to the W1 product's (`s1_ny1024`), absolutes
  larger. Peak kernel vs W1 should move ~+0.2 % (V6 report, pad_split).
* Capture must NOT decay to a much smaller late value than W1's does — the
  drain-free k=0 physics is unchanged.

Then upload products + a README to EOS `response_sim/s1_w2_ny512/` (xrdcp), and
put a supersession note in the repo: the W1 products (`s1/`, `s1_ny1024/`) stay
for shape/ratio reproduction but their ABSOLUTE amplitudes are 27 % low — model
the note on `scripts/condor/README.md`'s two-grid table. Update plan §0a +
T5 row + next-in-order item 3 ("W2 grid DONE, products at …").

## 4. Rebuild the fast path on W2 + recertify  (cheap)

`response/digitizer/kernel_lut.py` builds the LUT from a greens npz — read its
CLI/config to see how the product path is selected (it was pointed at
`s1_ny1024`). Rebuild from `greens_comb_w2_rho1M_dk50um_g19um.npz` (ρ_s = 1 MΩ
is the nominal production point — verify against what the digitizer config
used before; if it used a different ρ_s, match it). Re-run
`response/digitizer/test_lut_vs_solver.py` and `test_time_grid.py` — these
certify CACHING against the product they are fed, so they must pass unchanged
(~1e-4 / <0.5 %). If they read hardcoded W1 paths, repoint, don't fork.

## 5. T10 slow path, W2 verdict  (the decision fork)

Re-run `response/validation/t10_slowpath.py` against the W2 product (its
kernel-path argument; check `--help`). Note in the report that the W1 run used
ny=1024 and this uses ny=512 — a ~0.45 % pad-edge-shoulder grid effect, small
against the 8.26 % being retested. The z-lift (`zextend.py`) carries over
unchanged (gas-side physics, blind to the pad boundary). Units: t10 feeds the
S1 npz time axis in SECONDS through `t_ns*1e-9` — the guard added after the
19:12 void run will raise on any mix-up; do not bypass it.

* **PASS (<2 % worst event-normalized shaped residual):** the fast path is
  certified on W2 as-is; plan §7 step 5 (LUT from slow templates) is DEAD —
  record that in the T10 row.
* **FAIL:** record the number; §7 step 5 remains sequenced. Either way update
  the plan T10 row and drop a note in `design/report/T10_SLOWPATH_2026-08-08.md`
  (banner-style amendment; the file is owned by the ntof-x17-ed session's work,
  so amend, don't rewrite).

## 6. Stage B/C + T13  (stretch goal — only if the night is going well)

Avalanche calib decision (user, 2026-08-08, reconfirmed tonight): **the pooled
490 V Ar/iso 95/5 point is production** —
`response/avalanche/aval_calib_meshfield_pooled.json`. Do NOT wait for the HV
scan; it is for gain curves/T14 systematics later.

Re-point the digitizer chain (plan §7/§8 have the stage commands) at (a) the
W2 LUT, (b) the pooled meshfield calib, and regenerate `sim_decoded_*` ONCE
with both. Then T13 = wft reconstruction over the new sim_decoded through
`events.parquet` (plan T13 row; `wft/` package, mx_june_wft chain as the
model). Check inputs exist locally before starting (ClusterTrees; if they live
only on EOS/desktop, weigh the pull cost — do not hammer the desktop).
**T14 is NOT yours to run** — it is the blind comparison, run once, user
present. Same for anything in `response/avalanche/` (other session) and the
morning T7 collection checklist (also the other session).

## 7. Morning report

Leave `design/report/W2_NIGHT_REPORT_2026-08-09.md`: fleet statistics (jobs
resubmitted, wall times, memory peaks from `condor_history -af RemoteWallClockTime
MemoryUsage`), the 4 product acceptance numbers vs their bars, LUT cert
results, the T10 W2 verdict with its one-line consequence, what of §6 you
completed, and anything you had to fix (with commits). Update memory
`mx17-w2-resolve-2026-08-08` (it says "in flight" — flip it to the outcome)
and MEMORY.md's index line. Commit everything you touched.

## Known sharp edges (tonight's scars, do not rediscover them)

* Constraint masks are HARD by design; do not "improve" them back to
  fractional (capture collapses — see `wpot_w2.py` mask_snap comment).
* Drive snap is strictly `> 0.5`; `>=` overlaps the gap on pad-edge cells.
* `set -u` + LCG setup.sh don't mix; the job script guards it.
* lxplus system python3 has no scipy; everything condor-side sources LCG_105.
* laptop numpy `@` is reference-BLAS slow; scipy.linalg is OpenBLAS. On LCG
  both are fast. Combine does no big matmuls, so laptop is fine for it.
* Do not read the slabs through EOS fuse; xrdcp them.
* `w2slab_y000.npz` in `slabs/` (production, from 13353071) and in
  `slabs_mini/` (ny=256 smoke) share a basename — never mix the directories.
