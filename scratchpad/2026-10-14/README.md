# scratchpad 2026-10-14 — unbiased amplitude ledger, bG x shaper scan (FRESH_LOOK §18)

All scripts run on the laptop (yoga) with `~/PycharmProjects/nTof_x17/.venv/bin/python`, nTof_x17
imported read-only; copy them to `~` there first (they use absolute ~ paths).

- `profile_tails.py <t14 dir> <sim dir>` — transverse profile on all channels; far-channel 5σ rate.
- `hot_channels.py <sim dir> <data dir>` -> `hot_channels.log` — the sim's hot channels (§18.4).
- `peak_vs_theta.py <decoded dir> <events parquet> <tag>` -> `peak_vs_theta.log` — peak vs RECO angle.
- `t14_symwin.py` -> `t14_symwin.log` — T14 ledger with the angle window on both legs.
- `nb_ratio_vs_amp.py <t14 dir> <sim dir>` -> `nb_ratio_th00.log` — ±1/peak vs peak: no sign of
  DREAM compression below the rail (falls the same way in both legs).
- `m3_slopes_longrun.py` — adds M3 slopes and the det3 reco angle to `hv_slope/peaks_longrun.parquet`
  -> `~/x17/response_sim/hv_slope/peaks_longrun_m3.parquet`.
- `peak_vs_m3.py` -> `peak_vs_m3.log` — data peak vs TRUE angle; T14 selection factors.
- `rise_vs_m3.py` -> `rise_vs_m3.log` — data rise/FWHM by true angle via the T14 extractor.
- `ledger_m3.py [glob]` -> `ledger_m3.log/.csv` — the unbiased data/sim ledger (sim peaks cached as
  `peak_all_<feu>.npy` in each decoded dir).
- `orchestrate_bgtp.sh` — laptop driver for condor 4451431/4451432 (env CL, ANGS).

## Morning 2026-10-14 (desktop)
- `verdict_harness_table.py` (run on yoga): T14 harness table for fullphys175 / noisecm146 / verdict1. Harness amplitudes still carry the §18.1 data-selection bias.
- `kernels3d_ywings.log`: compare_kernels3d on i3d_dep_2e6_o2 (rho 2M), all Y0..Y5 channels: wings agree with the chain.
- `rho_wings.py` / `.log`: shaped point-source Y2/Y0, Y4/Y0, Y0/X0 vs rho_s (0.5/1/2/5 M) from the W2 kernels (§18.10).
