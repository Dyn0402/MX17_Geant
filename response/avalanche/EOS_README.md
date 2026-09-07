# S3 avalanche calibration — READ THIS BEFORE USING `raw/`

*Version-controlled at `response/avalanche/EOS_README.md`; deployed to EOS
`response_sim/avalanche/README.md`. Edit it here, not there.*

⚠️ **`raw/` in this directory is the PRE-FIX campaign and its induced-current
templates are all zeros.** Do not build a calibration from it.

Garfield's `AvalancheMicroscopic` and `AvalancheMC` both default to
`m_useWeightingPotential = true`, and `ComponentConstant` has no weighting
potential until `SetWeightingPotential()` is called. The first campaign set only
the weighting FIELD, so `WeightingPotential()` returned 0 everywhere and every
`i_elec` / `i_ion` came out identically zero — silently, with the avalanche and
the ion drift running perfectly normally and gain / sigma0 completely
unaffected, which is why it went unnoticed. Fixed in MX17_Geant commit
`42390d1` (`response/avalanche/mx17_aval_calib.py`).

| what | where | status |
|---|---|---|
| `raw/` — 56 slices, Aug 7 11:38 | here, on EOS | **PRE-FIX: zero templates, superseded** |
| `results_v2/` — 56 slices, Aug 7 14:14–14:46 | AFS `/afs/cern.ch/work/d/dneff/mx17_response/avalanche/results_v2/` | the real v2 raw |
| `aval_calib_v2.json` | here | built from `results_v2/`, correct |
| `aval_calib_v3.json` | here | v2 templates + the survival block. **Production.** |

## Verified 2026-08-08

- `python3 -m response.avalanche.collect` over `results_v2/` reproduces
  `aval_calib_v2.json` **bit-for-bit** — `max|diff| = 0.000e+00` on both
  `i_elec` and `i_ion`. So v2 is fully reproducible from archived material.
- An independent re-run with the current producer (same seeds, separate jobs,
  280 events) reproduces f_ion at 490 V to **1e-4**: 0.907809 against 0.907907.
  Its gain differs by 8.5 % on smaller statistics, which is expected — f_ion is
  a ratio and is the stable observable, which is exactly why it is the one
  checked.

### Scope of that verification — it is a snapshot, not a standing guarantee

Both checks were run against `mx17_aval_calib.py` at md5 `546b4ad2…`, the
**uniform-field** producer (`ComponentConstant`) that actually made v2. The file
has since moved under concurrent T7 work, which switches it to the **T6 field
map**. That is a deliberate physics upgrade, so a future campaign is *expected*
not to reproduce v2's numbers — the check above says the v2 templates are sound
and re-derivable, not that the current producer will reproduce them.

## Why `raw/` is kept rather than deleted or renamed

It is a real record of what was run, and renaming it would break anything
already pointing at it. `digitize.py` carries a template-content guard that
refuses an all-zero `i_elec`/`i_ion`, so a calibration accidentally built from
`raw/` fails loudly instead of turning the LUT silently to `nan`. That guard has
now caught this twice.

## T7 field-map campaigns (2026-08-08/09) — separate from the uniform-field history above

| what | where | status |
|---|---|---|
| `raw_meshfield_490V_20260808/` — 56 slices | here | first field-map run; its `--voltage` label is **fake** for every slice except 490V (see `MESHFIELD_QUARANTINE_README.md`, `response/avalanche/aval_calib_meshfield_QUARANTINED.json` in git) — `ComponentGrid` loads a fixed pre-solved map, so all 56 slices measured the identical 490V physics regardless of label |
| `aval_calib_meshfield_pooled.json` | here + git | the correct reduction of the run above: all 56 slices pooled as one 6400-event 490V point |
| `raw_meshfield_hvscan_20260808/` — 120 slices | here | real per-voltage campaign, Ar/iso 95/5 (460-530V) + Ar/iso 90/10 (530-590V), each against its own field map from `s2/meshfield_ladder/` |
| `aval_calib_meshfield_hvscan.json` | here + git | merged 15-point calib from the run above |
| `s2/meshfield_ladder/` — 41 maps | here | the gas-agnostic field-map voltage ladder (300-700V/10V, `solve_fieldmap.py --v-mesh`); pure electrostatics, so one map per voltage serves every gas — this is what the hvscan campaign's per-voltage lookup reads from |

All meshfield-campaign calib JSONs carry `voltage_V`/`gas_file`/`seed_z0_um`
per point (added 2026-08-09) rather than only in the string key or a fixed
mirrored constant downstream — see `response/avalanche/collect.py` and
`mx17_aval_calib.py`'s `config` block.

## Contaminant diagnosis grid (2026-08-09)

⚠️ **`raw_diagnosis_grid_20260809/` and `aval_calib_diagnosis_grid.json` are
DIAGNOSIS-GRID / unconstrained-contaminant-search, not a gas assay** — see
`response/avalanche/DIAGNOSIS_GRID_README.md` for the full framing. No
humidity was ever measured on the det3 bench; every water figure here is a
Magboltz fit to a slow measured drift velocity. 48 slices, 6 points: 0.5/1.5%
H2O and the June best-fit +1%N2 co-contamination at 490V, plus a 3-point
480/490/500V mini-scan for the leading 1% H2O candidate. Every point carries
`provenance.campaign_label` set to that exact string, so a consumer reading
the JSON directly sees the caveat regardless of filename.

## T14 HV-slope hunt (2026-08-09/12)

⚠️ **`raw_slopehunt_20260809/` and `aval_calib_slopehunt.json` are
DIAGNOSIS-GRID / unconstrained-slope-hunt, not a gas/field assay** — three
independent discriminators for why the sim's dry-95/5 gain-vs-voltage slope
(0.296/10V) is ~1.5x too shallow vs det3 data (0.449/10V): iso-fraction
(90/10 confirmation point, 24 slices), Penning rP A/B (0.30/0.50/0.65/0.80 x
3 voltage x 8 seed, 96 slices), field-map shape A/B (mesh vs uniform, 24
slices). **144/144 slices now landed** (see
`response/avalanche/SLOPEHUNT_OPS_LOG_2026-08-09.md` for the full run — the
rP=0.80 leg needed two reruns: a >10x nev cut after the first attempt stalled
for >11h with zero completions, then a lower-parallelism rerun of 4 slices
the kernel OOM-killed on the first rerun, all near-breakdown gain effects at
rP=0.80, not bugs).

⚠️ **Verdict: Penning does not close the T14 slope discrepancy.** An
early read of this campaign (mid-run, from the rP=0.80 slope alone) called it
"promising" — that was wrong, because it never checked the gain constraint.
See `design/report/DEEP_DIVE_2026-08-11.md` §21 for the full analysis: gain
grows exponentially in rP while the slope only grows linearly and weakly, so
no single rP satisfies both the data's gain and its slope. Field-map shape
and a gap-width scan are also ruled out as slope levers (same doc). The
leading suspect is the Magboltz α(E) cross-sections; further closure work on
this line is deliberately deferred until after MPGD26 (2026-09-03).

The final 4 rP=0.80 slices (rerun after the OOM kill) and the corrected
144/144 merge have not yet been shuttled to EOS — lxplus started enforcing a
second auth factor that a non-interactive `rsync` cannot satisfy. Results are
safe on the desktop in the meantime; see the ops log for the retry path.

## Open

Upload `results_v2/` (19 GB) here so the v2 raw lives on EOS and not only on
AFS work space.
