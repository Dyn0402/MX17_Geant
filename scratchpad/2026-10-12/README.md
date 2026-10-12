# Fresh look 2026-10-12 — scripts

Written up in `design/report/FRESH_LOOK_2026-10-12.md` §6–7.

- `quality_cut_amplitude.py`, `quality_cut_rise.py` — step 0 (quality cut).
- `toy_pulse.py` — independent toy forward model (numpy only, no `response/` imports).
  `toy_scan.py`, `toy_edges.py` — its first scans.
- `laptop_dump_waveforms.py`, `laptop_dump_angled.py` — run on the laptop
  (`~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/`, nTof_x17 `.venv`);
  write `waveforms_pm3.npz` (vertical, sim TRANSPFIX + data) and
  `waveforms_angled_pm3.npz` (angled ladder th00/thx10/thx20 + matching data):
  peak strip and ±3 neighbours, no quality cut.
- `wf_*.py` — shape analyses on those npz files (run with the npz in the cwd and
  this directory on `sys.path`): `wf_shape` / `wf_sumshape` (top shape, peak vs
  ±3 sum), `wf_edge` (rise split into 10→50 / 50→90 / 90→peak), `wf_pks` /
  `wf_pop` (edge vs peak position, amplitude, angle), `wf_matched` /
  `wf_matched_ang` (edge at matched peak delay), `wf_psi_scan` (front-loaded
  ion weighting), `wf_tail20` (averaged inclined pulse shapes), `wf_peak_scan`
  (shaper peaking time × ion fraction × transit against the inclined shapes).
- `wf_ion_T_scan.py` — §7.5: ion transit T scan (= mobility scaling) at 123 and 180 ns,
  full ions, against the inclined shapes. Needs pandas (nTof_x17 `.venv`) and both npz files in the cwd.
- `dream_settings_by_run.py` — per-run DREAM settings from all archived cfgs
  (→ `design/dream_cfg/dream_settings_by_run.csv`; record in `design/DREAM_SETTINGS.md`).
