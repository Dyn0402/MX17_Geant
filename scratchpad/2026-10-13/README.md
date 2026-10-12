# scratchpad 2026-10-13 — avalanche feedback test (FRESH_LOOK §9)

- `feedback_fit.py` — fits G_eff = c·G_Garfield(V)/(1 − b·G_Garfield(V)) to the SHAPE of
  the det3 27 June mesh ladder (free scale c, so no ADC calibration enters); controls:
  free-slope exponential, quadratic, Garfield-shape-only, and a 3-parameter Townsend form.
  Output in `feedback_fit.log`.
- `feedback_plot.py` — `design/report/figs/feedback_fit_2026-10-13.png`.
- `mesh_ladder_det3_2026-06-27.csv` — copy of the laptop's
  `~/x17/response_sim/hv_slope/mesh_ladder.csv` (nTof_x17 `mx17_sim_wft/hv_slope/analyse.py`
  output: threshold-free peak-strip max, M3 selection, fixed fiducial, rail flags).

Run with the nTof_x17 venv (pandas, scipy):

    ~/PycharmProjects/nTof_x17/.venv/bin/python scratchpad/2026-10-13/feedback_fit.py \
        scratchpad/2026-10-13/mesh_ladder_det3_2026-06-27.csv \
        response/avalanche/aval_calib_meshfield_hvscan.json
