# T6 field map — runbook

**This is the final field-map calculation for the response chain** (plan §2
`meshfield` contract, §4 S2). It replaces S3's `uniform_field` placeholder.
Audited and rebuilt 2026-08-08; this file is the single reference for running,
gating, and shipping it. Read the *Decisions* section before touching anything
— several numbers here are deliberately different from what the first draft
(neBEM `field_map.C`, retired) assumed.

## TL;DR for the operator

```bash
# on the DESKTOP (ssh desktop) — never run production on the laptop
cd ~/CLionProjects/MX17_Geant/response/meshcell
git pull

# 1. solve + export + solver gates S1-S6  (~40 s smoke, ~15-40 min production)
~/CLionProjects/MX17_Geant/.venv-fieldmap/bin/python solve_fieldmap.py --production --jobs 14

# 2. Garfield-side gates G1-G3+G7 on the exported map (~5-10 min)
source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh
root -b -q -e 'gSystem->Load("libGarfield"); gSystem->CompileMacro("gates_check.C","k");
               gates_check("meshfield_production.txt", 2000);'

# 3. ship (see "Products & shipping") — only after BOTH steps print ALL ... PASS
```

Every gate is a hard stop: a FAIL means the map is not to be used, no
exceptions, no tolerance-widening. The gates and their meanings are listed
below.

## What is being solved

One unit cell (67 µm pitch) of the MX17 woven micromesh, approximated as a
bi-planar crossed-wire grid (x-wires at z = +9.0 µm, y-wires at −9.0 µm,
Ø19 µm, 1 µm kiss overlap at the crossing standing in for wires resting on
each other), between:

- the **ESL anode plane** at z = −168.5 µm held at +490 V (bench operating
  point, mesh = 0 V reference), and
- a **drift-slice top plane** at z = +218.5 µm, i.e. 200 µm above the mesh
  topside (plan §2 contract).

Solver: gmsh (conformal cylinder surfaces) + scikit-fem P2 tetrahedra +
pyamg, in `solve_fieldmap.py`. Side walls carry homogeneous Neumann BCs;
because the bi-planar cell is mirror-symmetric about x=0, y=0, x=±p/2,
y=±p/2, that reproduces the **infinite periodic lattice exactly** — no
periodicity machinery, no truncation error.

Output `meshfield_production.txt` is in `ComponentGrid::LoadElectricField`
"xyz" format (header `#` lines carry the mesh, so no `SetMesh` needed):

```cpp
ComponentGrid grid;
grid.LoadElectricField("meshfield_production.txt", "xyz", true, true);
grid.SetMedium(&gas);
grid.EnablePeriodicityX();
grid.EnablePeriodicityY();
```

A provenance JSON (`meshfield_production.json`) with git hash, mesh sizes,
gate values and timings is written next to it — it travels with the map,
always.

## Decisions (the audit trail)

1. **neBEM is disqualified for this geometry — by measurement.** Its
   bounding-plane (`AddPlaneZ`) and conductor-mirror interface data are read
   by `neBEMInterface.c` but never used in the solve, so plates can only be
   finite patches replicated by periodic copies. The resulting
   capacitor-fringe error on the drift side converges like 1/n²:
   +540 % at 5 copies, +12 % at 30 (`copies_scan.C`, kept in the repo as the
   evidence). ~100 copies would be needed for 1 %, with eval cost ∝ n².
   The amplification gap, enclosed by near-parallel conductors, does
   converge: **Ez(amp bulk) → 30 561 V/cm**, kept as an independent
   cross-solver check (gate S2/G2, FEM agrees within 1.5 %).
2. **The 150 µm amplification gap is pillar height** — ESL surface to mesh
   **underside** (`shared/MX17ModuleGeometry.hh`: bulk pillars stand on the
   ESL and support the mesh). The anode is therefore at −168.5 µm from the
   weave mid-plane, not −150 µm, and the mean amp field at 490 V is
   **~31.0 kV/cm, not 32.7 kV/cm**. S3 gains re-run in this map will come
   out several times lower than the uniform-field first pass — that is
   physics, not a regression; the uniform-field pass implicitly measured the
   gap from the mesh's electrical centre. (S3's absolute gain is calibrated
   against measured detector gain either way.)
3. **Mesh field penetration sets the top BC.** The 31 kV/cm amp field leaks
   through the 46 %-open mesh and offsets the drift-side potential by
   +1.95 V (effective penetration depth 0.65 µm). Over the real 30 mm drift
   that is 0.2 % of the 1000 V and irrelevant; over a 218 µm slice it is a
   28 % field error if the top plane is naively set to −E_drift·z. The
   solver therefore solves two unit problems (anode-only and top-only) and
   combines them linearly so the **bulk drift field is exactly 333 V/cm**
   (det3/T14 target). The naive single-BC approach must not be reintroduced.
4. **Region-flag semantics are erosion-compensated.** Garfield's
   `ComponentGrid` deactivates any point whose interpolation cell touches an
   inactive node — wires would look one grid step fatter (measured on the
   3 µm smoke map: transparency 0.67 instead of ~0.9). The exporter erodes
   the written flag by one grid step and fills the in-wire shell nodes with
   their radial surface projection's field, so the absorption boundary sits
   on the true wire surface and near-surface interpolation is never mixed
   with zeros. At the 1 µm production grid the residual surface bias is
   ±~0.5 µm.
5. **A local Garfield++ patch is required** (both hosts, applied 2026-08-08
   at the pin 927e5c21, `Source/ComponentGrid.cc`, marked `MX17 local
   patch`): (a) the 3D `xyz`/`ijk` loader parsed the region flag but never
   stored it — 3D maps silently lost their inactive regions (electrons
   sailed through wires, transparency exactly 1.000); (b) `m_active` was
   allocated before the header-line mesh dimensions were known, so storing
   the flag segfaulted. **Any Garfield rebuild or pin move must re-apply
   this patch until it is upstreamed** — the symptom of losing it is
   G7 = 1.000. Patch also noted in `garfield_sim/setup_garfield.sh`.
6. **The kapton+glue PCB update (commits `405e1a6`…`fb4b8ed`) does not
   enter this map.** That work refines the dielectric stack *below* the ESL
   (kapton 50 µm + lamination glue ≈ 18.8 µm in series) for the S1
   weighting-potential solver. The field map's lower boundary is the ESL
   held at the anode potential: at DC it is an equipotential and screens
   everything beneath it. If the ESL model ever stops being a uniform
   equipotential (patterned resist with exposed grooves), revisit this.
7. **Not modelled, on purpose:** bulk pillars (1.3 % of the area, Ø0.6 mm on
   4.68 mm pitch — far larger than the 67 µm cell, they belong to an
   acceptance/efficiency correction, not the unit-cell field); the true
   over/under weave (the bi-planar stack keeps the wire-to-wire focusing
   that sets transparency and funnelling; a woven-geometry refinement would
   need per-wire sinusoidal sweeps in gmsh and breaks the mirror symmetry
   — do not attempt it casually).

## The gates

Solver side (`solve_fieldmap.py`, all must PASS — the script exits nonzero
otherwise and the JSON records every value):

| gate | checks | bar |
|---|---|---|
| S1 | bulk drift field at a plane 40 µm below the top | +333 V/cm ± 1 % |
| S2 | amp bulk field vs neBEM cross-solver value 30 561 | ± 2 % |
| S3 | potential just off the wire surface | < 2 V |
| S4 | mirror symmetry of Ez (median over random points) | < 0.5 % |
| S5 | custom P2 evaluator vs scikit-fem's own interpolator | < 1e-8 |
| S6 | production only: field vs a 1.4× coarser solve, p95 | < 1 % |

Garfield side (`gates_check.C` — loads the map exactly the way S3 will):

| gate | checks | bar |
|---|---|---|
| G3 | file loads, mesh extent sane | — |
| G1 | far-field drift through the loaded grid | +333 ± 1 % |
| G2 | amp bulk through the loaded grid vs neBEM | ± 2 % (3 % smoke) |
| G7 | electron transparency at the bench point | in (0.55, 0.98) |

G7 is a catastrophe band, not a precision target: 1.000 means the region
flag was lost (see Decision 5), ~0 means a sign/BC error. **The production
G7 value is the T6 transparency deliverable and supersedes the 1D-model
0.873 (`DEFAULT_TRANSPARENCY` in Stage B) once accepted.** Smoke-map values
carry a ±0.13-scale grid bias (0.67 dilated / 0.94 eroded at 3 µm) — never
quote them.

## Products & shipping

After both gate sets pass on the desktop:

```bash
# desktop -> laptop (the .done marker discipline: rsync only after the
# producer exited; see design/report/DESKTOP_RUNS_2026-08-07.md)
rsync -a desktop:~/CLionProjects/MX17_Geant/response/meshcell/meshfield_production.{txt,json} \
      ~/x17/response_sim/meshfield/
```

Then add one manifest row (date, git hash from the JSON, gate values, where
the files are). The map is ~180 MB of text; it is gitignored — the JSON is
the record, the map is data.

## What S3 does next (the actual point of all this)

Re-run the avalanche campaign in this map instead of `uniform_field`
(plan §5, T7): load as in the snippet above, seed electrons above the mesh,
`AvalancheMicroscopic` as before. Expect: mean gain several × lower than the
uniform-field table (Decision 2), transparency and funnelling now emergent
rather than assumed. The `uniform_field` provenance nulls in
`aval_calib_v2.json` are exactly the slots this fills.

## Cheap follow-ups the map enables (good tasks for a lighter model)

- **Transparency curve ε(E_amp/E_drift)** (S2 deliverable a): loop
  `E_DRIFT` in `solve_fieldmap.py` over {50, 100, 200, 333, 500, 700, 1000},
  smoke-resolution solves are fine for the curve shape (but quote the bench
  point from the production map), G7 machinery gives ε per point. Compare
  against `mesh_transparency.csv` (1D model) and the literature bulk-MM
  curve as the plan asks.
- **Funneling map** (S2 deliverable b): drift electrons from a grid of
  (x0, y0) entry points through the loaded map with diffusion OFF
  (`AvalancheMicroscopic` seeds, or `DriftLineRKF`), record arrival (x, y)
  at the ESL → entry→seed offset map.
- **Ion endpoints** (S2 deliverable d): needs a gas table with ion
  mobility (reuse `garfield_sim` tables); `AvalancheMC` ions from the
  avalanche region, fraction terminating on mesh vs drifting up.

## Troubleshooting

- `G7 = 1.000` → the Garfield local patch is missing (Decision 5).
- Solver gate S2 drifts by >2 % after a geometry edit → you changed the
  z-datum or the overlap; the neBEM cross-check number 30 561 assumes wires
  at ±r with no overlap and pillar-height gap (re-derive with
  `copies_scan.C` if the geometry legitimately changed).
- `python3 -m venv` fails on the desktop (no ensurepip): the venv there was
  bootstrapped with `--without-pip` + get-pip.py; recreate it the same way.
- Killed mid-solve on the laptop: that is what the desktop is for. The
  original neBEM draft died exactly this way (single-threaded SVD).
