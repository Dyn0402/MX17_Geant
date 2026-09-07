# T14 slope-hunt desktop campaign — operational log (2026-08-09 to 08-12)

This is the operational record of running the T14 HV-slope-hunt compute
campaign on the desktop: what was launched, what broke, and what is still
outstanding. It does **not** carry the physics verdict — that lives in
`design/report/DEEP_DIVE_2026-08-11.md` §21 and its follow-ups, which supersede
an earlier, incomplete read given mid-campaign (see "Superseded interpretation"
below). Treat that doc as authoritative for the physics; treat this one as
authoritative for "did the campaign actually finish, and what's left to run."

## 1. What ran, and completed cleanly

- **Dry HV-scan**, Ar/iso 95/5 (460-530V) + 90/10 (530-590V) against the real
  T6 mesh field-map ladder — 120 slices, no failures. `aval_calib_meshfield_hvscan.json`.
- **Contaminant diagnosis grid** — 48 slices (0.5/1.0/1.5% H2O + N2 co-contam
  variant, 490V; 3-point 480/490/500V mini-scan for the 1% H2O leader). No
  failures. `aval_calib_diagnosis_grid.json`.
- **Slope-hunt campaign, 144/144 slices, now fully landed** (see §2 for the
  saga to get there): iso-fraction confirmation (90/10, 24 slices), Penning rP
  A/B (0.30/0.50/0.65/0.80 x 3 voltage x 8 seed, 96 slices), field-map shape
  A/B (mesh vs uniform, 24 slices). `aval_calib_slopehunt.json`.
- **Det4 wet-quaternary Magboltz gas tables** — 3-point H2O bracket
  (1.3/1.5/1.7%) on the Ar/CF4/iC4H10 88/10/2 base, rP=0.40 manual Penning
  (same fix as the dry ternary). 6 tables built (3 points x Saclay/CERN drift
  length), no errors. **No avalanche calibration has been run against these
  yet** — gas tables only (see §4).

## 2. The rP=0.80 saga

The Penning rP A/B leg's top point (rP=0.80) pushes Ar/iso 95/5 into a
near-breakdown regime: `mean_gain` in the 8x10^5-1.8x10^6 range at 460V,
rising to 5-12x10^6 at 520V. Two failure modes followed directly from that:

1. **First attempt, full nev (230/140/60 per the standard schedule): killed
   after >11h with zero completions.** Not hung — `AvalancheMicroscopic` was
   genuinely still tracking the swarm (progress logs showed real, if glacial,
   advancement, e.g. 100/230 events after 9.2h at 460V). Killed the batch,
   rewrote the points file with nev cut ~10x (460:25, 490:15, 520:6) and moved
   to the back of the queue, behind the field-map-shape leg (which shares
   nothing with Penning and had been blocked behind the stuck batch).
2. **Second attempt, reduced nev, 16-way parallel: 20/24 completed, 4
   OOM-killed** (490V_s4, 520V_s1/s2/s5 — confirmed via `journalctl -k`,
   `Out of memory: Killed process ... python3 ... anon-rss:6-8GB`). Even at
   nev=6-15, a single event's electron swarm at this rP/voltage needs 6-8GB
   RSS, and 16 such processes landing concurrently exceeded the box's 62GB.
3. **Third attempt, same 4 slices, 2-way parallel: all 4 completed cleanly.**
   Bounding concurrency (not nev) was the actual fix — the box has plenty of
   memory for a few of these processes, just not sixteen at once.

Final merge (144/144, `python3 -m response.avalanche.collect`) succeeded
without incident.

## 3. EOS shuttle — partially blocked, not urgent

- Raw results + first merge (140/144, before the final 4-slice retry) shuttled
  to EOS successfully once a fresh Kerberos ticket was obtained (`kinit`,
  interactive — the ticket had expired past what `kinit -R` alone could
  renew).
- **The final 4 slices + the corrected 144/144 merge have NOT been shuttled.**
  lxplus now enforces a second authentication factor
  (`keyboard-interactive`) on top of Kerberos/GSSAPI, which a non-interactive
  scripted `rsync` cannot satisfy (`can't open /dev/tty`). This needs an
  interactive `ssh lxplus` session to clear, same as the Kerberos ticket did.
  **Action needed: SSH into lxplus interactively once, then the shuttle can
  be retried.** Nothing is at risk in the meantime — the complete 144/144
  results and merge are safe on the desktop
  (`/media/ucla/mx17_response_sim/avalanche/results_slopehunt/`,
  `~/CLionProjects/MX17_Geant/response/avalanche/aval_calib_slopehunt.json`).

## 4. Not yet run

- **Any avalanche/gain-scan campaign against the det4 wet-quaternary gas
  tables.** Only the Magboltz transport tables exist; no `mx17_aval_calib.py`
  points file or campaign has been built or launched for det4. No such
  request has been queued — this was gas-table prep only.
- **The final EOS shuttle** of the 4 retried rP=0.80 slices and the
  corrected 144/144 merge (§3).
- **Any rP point beyond 0.80.** Not pursued, and shouldn't be without a
  specific reason to: per `DEEP_DIVE_2026-08-11.md` §21, rP=0.80 already
  exceeds the physically-motivated gain budget, and rP≈0.91 (where the slope
  alone would match) is off the scanned range for good reason — its
  extrapolated gain is ×322 the data-implied value.
- **Uncertainty propagation on the slope-hunt slopes** beyond point estimates.
  Given the superseded-interpretation note below, this is likely moot for the
  Penning axis specifically (the DEEP_DIVE analysis already did this more
  carefully with the gain-side constraint included) but was never done for
  the iso-fraction or field-shape legs either.

## Superseded interpretation

Partway through this campaign, before the full 144/144 result was in hand, I
reported the rP=0.80 point (slope 0.426/10V vs data's 0.449±0.009/10V) as "a
genuinely promising result" for closing the T14 HV-slope discrepancy via
Penning transfer alone. **That framing was incomplete** — it only checked
whether the slope approached the data value and never checked the gain
constraint. A concurrent, more rigorous analysis in this same repo
(`design/report/OVERNIGHT_2026-08-10.md` §21, `DEEP_DIVE_2026-08-11.md`)
checked both and found gain grows exponentially in rP (×3.46 per 0.1) while
slope grows only linearly and weakly (+0.028 per 0.1) — the two constraints
cannot be satisfied by the same rP. **Penning does not close the T14
discrepancy.** See that doc for the full verdict and the current leading
suspect (Magboltz α(E) cross-sections), and its note that further closure work
on this line is deliberately deferred until after MPGD26 (Sep 3).
