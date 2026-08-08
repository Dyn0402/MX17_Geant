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

## Open

Upload `results_v2/` (19 GB) here so the v2 raw lives on EOS and not only on
AFS work space.
