# scratchpad 2026-10-13b — moisture in the chain, and the archived two-gas ladders (FRESH_LOOK §12)

Run everything with `~/PycharmProjects/nTof_x17/.venv/bin/python`; nTof_x17 is imported read-only.

- `wf_wet_scan.py` — §10.2 toy/fit with the drift transport (v, D_T, D_L) as parameters; dry vs
  wet-1 % arms. Needs `waveforms_angled_pm3.npz` and `waveforms_pm3.npz` (laptop T14 dumps) in the
  cwd → `wf_wet_scan.log`.
- `extract_may.py RUN SUBRUN OUT.parquet` — threshold-free peak-strip amplitudes, M3-selected, for
  the May old-det3 runs in `~/x17/cosmic_bench/det_3_old/` (hv_slope `view_rows()` reused).
  The 23 parquets were written to the job tmp dir (not kept); regenerate in ~1 min with xargs -P 10.
- `ladder_may.py PARQUET_DIR` → `ladder_may.csv` (per gas/view/voltage quantiles, fiducial, rail rule).
- `fit_may.py` → `fit_may.log` (local slopes, curvature, exp/quad/free-slope feedback fits).
- `plot_may.py` → `design/report/figs/may_two_gas_ladder_2026-10-13.png`.
