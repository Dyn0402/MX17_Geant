# aval_calib_meshfield_QUARANTINED.json — do not use

This is the raw `collect.py` output over the first T7 field-map campaign
(56 slices, run 2026-08-08), grouped by its `--voltage` label into 7
apparent points (470-530V). **Those 7 points are NOT a voltage scan.**

`ComponentGrid` loads a fixed, pre-solved field map
(`meshfield_production.txt`, solved once at V_MESH=490V in
`solve_fieldmap.py`) — `--voltage` in `mx17_aval_calib.py`'s field-map
branch only ever set `e_field` for the (bypassed) `uniform_field` path.
`mx17_aval_points.txt` was carried over unchanged from the uniform-field
campaign, where `--voltage` genuinely was the field, so every one of the
56 slices simulated the IDENTICAL 490V physics regardless of its label.
The apparent flatness of gain/theta/sigma0 across the 7 "voltages" in this
file is not a physics finding — it is confirmation that they are all the
same measurement.

The correct reduction of this campaign is `aval_calib_meshfield_pooled.json`
(all 56 slices merged as ONE 490V bench-point measurement, 6400 events).
A real voltage-dependent field-map curve requires per-voltage field maps —
see the map ladder started 2026-08-08 evening (`--v-mesh` in
solve_fieldmap.py) — followed by a re-run of the avalanche campaign against
the correct per-voltage map for each point.

`mx17_aval_calib.py` (commit 6e2aad0) now guards against this exact mistake:
it cross-checks `--voltage` against the map's own provenance JSON and warns
+ self-corrects rather than silently mislabeling.
