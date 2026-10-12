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

Also (FRESH_LOOK §10):

- `noise_shaper.py` — the DREAM shaper's time scale from the target run's noise autocorrelation
  (det3 27 June long run, decoded_root from EOS) → ~166 ns.
- `wf_ion_T_scan166.py` — §7.5's ion-transit scan rerun at 166 ns (needs the laptop's
  `waveforms_angled_pm3.npz` in the cwd) → `wf_ion_T_scan166.log`.
- `delay_vs_rho.py` — ±1 Y-row delay from the chain's W2 kernels vs ρ_s and cloud size
  (`delay_vs_rho.log`); `--ion` varies T_ion / f_ion instead (`delay_vs_ion.log`).
  Products: `/media/ucla/x17_resp/s1_w2_ny512/` (copies of EOS `s1_w2_ny512/products`).
- `xtalk_fit.py` — vertical-track data vs chain neighbour shapes, with a crosstalk term
  (`waveforms_pm3.npz`) → `xtalk_fit.log`.
