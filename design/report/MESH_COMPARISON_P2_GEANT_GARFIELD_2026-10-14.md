# Mesh and field modelling: MX17 against Alexandra's P2 simulation (2026-10-14)

Compared: this repo (`response/meshcell`, `response/comsol`) with
`github.com/akallitss/P2_Geant_Garfield` at commit `038a661` (2026-10-13, "Mesh
transparency: review, drift sweep, literature validation"), plus the COMSOL model she
sent, `~/Downloads/MeshUnitCellIon_PressureStudy.mph` (COMSOL 6.3, saved 2026-10-05).
Her repo is also Claude-assisted and written in the same style as ours (handoffs, runbooks,
audit trails), so the two are easy to compare line by line.

Only one file from her was in `~/Downloads`, the `.mph`. If a second file was meant, it did not arrive.

## 1. Headline

1. **Her main methodological finding changes our transparency number, and I tested it
   here.** Garfield's `AvalancheMicroscopic` by default reads the field once per free
   flight. In Ar or Ne with a few percent isobutane (the Ramsauer minimum) the free paths
   exceed the 18 µm wire, so electrons near the mesh are tracked with the wrong field and
   are lost on wires spuriously. `EnableNullCollisionSteps(true)` re-reads the field at
   every null collision. She found Ar/iso 95/5 going 90.6 → 96.1 %. Our `gates_check.C`
   never enabled it. Re-running our own G7 gate with it on (our maps, our gas, 333 V/cm,
   3000 electrons):

   | map | default steps (as shipped) | null-collision steps |
   |---|---|---|
   | straight wires (production) | 0.947 ± 0.004 (accepted run: 0.955; 0.960 in the weave comparison) | **0.988 ± 0.002** |
   | true weave (COMSOL) | 0.928 ± 0.005 (earlier run: 0.918 ± 0.006) | **0.988 ± 0.002** |

   Two consequences, both still preliminary (one seed each, no status breakdown):
   - The **4.2 % weave deficit in the HANDOFF and README is a stepping artefact**: with the
     better tracking both maps give the same number. The bi-planar stand-in is fine for
     transparency.
   - The **0.955 / 0.9559** that feeds `DEFAULT_TRANSPARENCY` and the avalanche
     `survival` in `aval_calib_meshfield_*.json` is probably ~3 % low. That is the
     direction the T14 amplitude tension needs (the data want more signal), but it is
     small, and the gain itself may also move because the same tracker runs in the
     avalanche campaign. I have not tested that.
   Her sweep agrees with this: Ar/iso 95/5 reaches 99.3 % at 450 V/cm, and our ratio
   E_drift/E_amp = 0.011 sits on her plateau.

2. **Her field-map route is different in kind from ours** (§3). Hers is COMSOL's
   unstructured tetrahedral mesh read by `ComponentComsol`, with the wire surface exact.
   Ours is a regular 1 µm grid with region flags, where the wire surface is represented to
   ±0.5 µm and carries an acknowledged asymmetric systematic (−0.045). Each has real
   advantages; see §5.

3. **Things she has that we lack:** validation against a measurement (Nikolopoulos et al.
   2011, 2.4 % rms), a loss-accounting by status code (wire, attachment, other), the
   electron characteristic-energy explanation of which gases lose electrons, and a
   Garfield bug fix for periodic maps.

4. **Things we have that she lacks:** a drift boundary condition by surface charge that
   gives the exact drift field in a 1-minute solve (no 2–4 mm drift gap, no penetration
   shift), the true weave verified against straight wires, a gate suite on every map
   (solver-side S1–S6, Garfield-side G1–G7), a provenance JSON shipped with each map, and
   a worked induction/resistive-layer chain on top.

## 2. Geometry and physics, side by side

| | MX17 (this repo) | P2 (her `scripts/comsol`, base `.mph`) |
|---|---|---|
| Mesh | 19 µm wire (r = 9.5), pitch 67 µm | 18 µm wire (r = 9), measured pitch 68 µm (50 µm opening) |
| Weave | wires z = ±A sin(πs/p), A = 9.0 µm, swept circle | z = meshHeight + rAmp·sin(2πs/pitch) with pitch = 2 × wire spacing, rAmp = 0.9·rWire = 8.1 µm, circle swept along a parametric curve, second wire phase-shifted by half a period |
| Wire overlap at a crossing | 1 µm (A = r − 0.5) | 1.8 µm (2·r − 2·rAmp) |
| Calendering | none (round wires, 37 µm at the crossing) | wires cut flat at ±1.5·rCut = ±14.85 µm (29.7 µm at crossings), taken from the base model and from the 30 µm of Nikolopoulos; **not measured for P2** |
| Cell | one spacing p × p, Neumann walls, mirror periodicity (translation period is 2p) | 2p × 2p periodic cell, so it needs no mirror trick |
| Amp-gap datum | 150 µm pillars: ESL to mesh **underside** (anode at −168.5 µm from the weave mid-plane). Decided at the start (runbook Decision 2) | The inherited model put the wire **centre plane** at 150 µm, so the mesh bottom was at 135.5 µm. She caught this in her 2026-10-14 review and re-solved with meshHeight 164.85 µm. Same conclusion, reached later |
| Drift side | 200 µm slice; top plane carries a surface charge so the far field is exactly E_drift (`WovenMesh.java`), or the two-solve superposition in `solve_fieldmap.py` | the real 2 mm or 4 mm drift gap with a cathode potential. Mesh penetration then shifts the drift field (625 V/cm nominal gives 632 V/cm) and each drift field needs its own solve and re-mesh |
| Dielectrics | none in the map (ESL treated as an equipotential; kapton/glue handled in the induction solver) | gas ε = 1, mesh and insulator ε = 4 (domains 4–11 in `ComponentComsol`) |
| Electrodes | mesh 0 V, anode +490 V (bench) | pads 0 V, mesh −450 V, cathode set from the drift field (SPS July 2026 polarity) |
| Weighting potentials | separate chain (`response/solver`, `Induce3D`) | three Terminals in the cell (anode, mesh, cathode); the exports are disabled in the base model |
| Gas | Ar/iC4H10 95/5, 745.83 Torr | six SPS mixtures (Ar/iso 95/5 and 90/10, Ne/iso 95/5 and 85/15, Ar/CO2/iso, Ar/CF4/iso) |

Read from the `.mph` itself (parameters, geometry feature list, analytic function in
`dmodel.xml`; nothing was executed). Not examined: the element count (the stored mesh is
138 MB and the extended mesh 494 MB, against our 419 k tets), the second study and the
`Electrostatics 2–5` interfaces behind the "Ion / Pressure study" name, and the
`xmesh1`/refinement levels 0–3 stored in the solutions.

## 3. Field-map handling

| | MX17 | P2 |
|---|---|---|
| Solver | gmsh + scikit-fem P2 tets (production), COMSOL (weave) | COMSOL, its own tets |
| Into Garfield | `ComponentGrid` text map, 1 µm grid, 1.79 M nodes, periodic (or mirror) | `ComponentComsol` from `mesh.mphtxt` + `potential.txt`, periodic X/Y, interface crossings on |
| Wire representation | region flag per node, **eroded** by one grid step, in-wire shell nodes filled by radial projection; ±0.5 µm surface bias; needs our local `ComponentGrid.cc` patch | exact tetra boundary, `IsDriftable()` tells wire from gas; needs her local `InterfacePanels` patch |
| Export route | COMSOL's sandbox blocks all file writes from model code, so we print `CSV,` lines to stdout and reassemble in Python (`woven_map.py`) | writes `mesh.mphtxt` and `potential.txt` from COMSOL export nodes straight to EOS in batch |
| Checks | gates S1–S6, G1–G3, G7, independent numpy re-derivation of the flags, provenance JSON | probe mode (prints solid z extents at chosen x, y), field-line maps, loss accounting by status code |
| Grid-dependent systematic | −0.045 on transparency, bracketed by the dilate/erode pair | none of that kind; the exact surface is in the mesh |

Notes:
- Her export route is worth trying because our memory says "a Data export with a
  coordinate file silently writes nothing". Her export nodes have **no coordinate file**:
  they write the mesh and the nodal solution under the node's own `filename`, and that
  works in batch on lxplus under 6.4. If it also works for `WovenMesh.java`, we could
  load the COMSOL tets directly, drop the 171 MB text map and the erosion/flag logic, and
  remove the grid systematic. Untested; the sandbox rule may differ for our setup.
- Her interface-crossing patch matters only with `ComponentFieldMap` classes
  (`ComponentComsol`). `ComponentGrid` has no such option, so for our current maps that
  part does not apply. Null-collision steps (§1) apply to both.
- The base model has no `ComponentComsol` dielectric table; she writes one by hand
  (`mesh_transparency.cc`, "2\n1 4\n11\n…") from the `MATERIAL` lines SolveP2MeshCell
  prints. Worth knowing if we go that way.

## 4. What she did that is new to us

1. **Null-collision tracking** (§1). Largest effect found.
2. **Status accounting.** Every electron that does not pass is classed as wire (−5),
   attachment (−7), or other (e.g. left through the top), and the top of the volume was
   moved 1 mm up so back-diffusion cannot end a track. We count "not below the mesh" as a
   loss, which is fine in Ar/iso (no attachment) but would silently mix mesh and
   attachment losses in any electronegative mixture (CF4, CO2). She shows that the entire
   Ar/CF4/iso loss is attachment (4.5 %), not the mesh.
3. **Physics explanation of which gases lose electrons.** Wire loss follows the electron
   characteristic energy ε_k = e·D_T/μ (Magboltz, 638 V/cm): Ar/CF4/iso 0.13 eV → 0 %,
   Ne/iso 85/15 0.23 → 0.1 %, Ar/iso 90/10 0.40 → 0.8 %, Ar/CO2/iso 0.46 → 2.0 %,
   Ne/iso 95/5 0.60 → 6.2 %, Ar/iso 95/5 0.76 → 5.7 %. Hot electrons diffuse off the
   converging field lines onto the top of the wires.
4. **Validation on data.** She rebuilt the Nikolopoulos et al. 2011 mesh and gas
   (Ar/CO2 85/15, 63.5 µm, 128 µm gap, 41.4 kV/cm) in the same COMSOL cell and compared
   the transparency vs E_drift/E_amp with measurements read off their figure: 2.4 % rms
   with null-collision steps, 3.4 % without. Caveat she states herself: a cold gas hardly
   tests the stepping issue, and the hot mixtures have no data validation.
5. **Garfield bug found and patched**: with interface crossings on a periodic map,
   `InterfacePanels::CrossedFolded` tested a chord across the whole cell and stopped
   about half the electrons on the cell faces (45 % instead of ~90 %).
   `scripts/garfield/patches/interfacepanels_periodic_leaf.patch`; worth reporting upstream.
6. **A drift-field sweep with a literature curve** (plateau, fall-off at
   E_drift/E_amp ≈ 0.016, comparison with a bulk Micromegas measurement she cites as
   arXiv:2602.21238; I did not check that reference).
7. **Per-gas generation of the gas JSON** with the table's own field points (`export_gas_json.C`):
   the Ar/CO2/iso Townsend coefficients reached 1e33 /cm because a regular grid sampled
   across the table's gap. Our gas tables come from different scripts, but if anyone
   exports Townsend α from a `.gas` file on a regular grid, check for this.

## 5. What we have that she could use

1. **Exact drift field without a drift gap.** The surface-charge top boundary
   (`WovenMesh.java`) gives E_drift exactly after one solve of ~1 minute (419 k tets,
   8 cores). Her sweep needs a separate 2 mm cell and a re-mesh per field because the
   field "cannot be superposed node by node". With the surface-charge condition the same
   mesh holds for any field, and the penetration offset (625 → 632 V/cm) disappears.
   Her own check that the 2 mm and 4 mm cells agree supports using the shorter slice.
2. **The mirror-symmetry reduction.** Cell p × p instead of 2p × 2p, 4× fewer nodes, with
   mirror periodicity in Garfield. We also confirmed that straight and woven wires give
   the same transparency (§1), so the weave is not what the mesh simulation is sensitive to.
3. **A gate suite and provenance.** A pass/fail list on every map (bulk drift field,
   amp-gap field against an independent solver, potential at the wire surface, mirror
   symmetry, refinement error against a 1.4× coarser solve), each run writes a JSON with
   git hash, mesh sizes and gate values.
4. **The cross-solver check.** Her amp field (28.5 kV/cm) rests on a single COMSOL
   solve. Our neBEM number for the amp gap converges independently
   (`copies_scan.C`) and agrees with the FEM to 1.5 %.
5. **Full induction chain.** COMSOL 3D strip/pad kernels, the resistive-layer
   operator, the DREAM digitizer.

## 6. Disagreements and open questions

- **Does the same null-collision issue hold for the ion/gain side of our chain?** The
  avalanche campaign (`mx17_aval_calib.py`) tracks electrons in the same maps with the
  default stepping. The `survival` value there (0.9559) is a transparency measurement.
  A rerun with `EnableNullCollisionSteps(true)` on one HV point would show whether gain
  and survival move.
- **The 0.955 −0.045/+0.005 systematic.** With better tracking the grid-surface
  systematic is probably still real, but it needs a rerun to see; the dilate/erode
  bracket was computed with the old stepping.
- **Calendering.** Our wires are round (37 µm at the crossing). Hers are cut flat to
  29.7 µm because the inherited model was. If the MX17 mesh is calendered, the crossings
  are flatter and the amp gap is slightly larger at the same pillar height. We do not
  know. Hers is not measured either.
- **Pillars.** We model none in the cell (1.3 % of the area, Ø0.6 mm on a 4.68 mm pitch);
  hers treats the pillars in Geant4 only (3.8 % / 4.9 % of the area).
- **Different operating points.** Her E_drift/E_amp = 0.020 (SPS) is in the fall-off;
  ours is 0.011, on the plateau. The two sets of numbers are consistent, but neither
  covers the other's point without a rerun.
- **The "colleague" who made the base model** is not named in her repo, only "used with
  permission"; the `.mph` says 2023 creation dates and holds a pressure-study second
  study. Ask her who owns it before we reuse it in anything we publish.

## 7. General look at the rest of her repo

- **Same layout and conventions.** Geant4 (`src/`, 12 files, 1.9 k lines) for Stage A,
  Python `stage_b/` for drift, diffusion, transparency, Polya gain, pad lookup, then a
  VMM emulator. Her README says "lift `response/digitizer/` from MX17 rather than
  writing fresh"; she took the "hard-won knowledge" instead (drift-frame trap,
  nPrimary-not-edep/W, assumptions recorded in the output).
- **Stage B is much thinner than our chain.** Her mesh transparency in Stage B is still a
  constant 1.0, the gain, induced current and pillars are placeholders (all listed in
  `stage_b/README.md` §3). Hers is a pad detector with ~1.1 pads/event, so a static Ramo
  induction is enough; ours needs the resistive-layer chain.
- **Her Geant4 mesh is an effective-density slab** like ours (hers 36 µm steel at
  ρ × 0.224, ours 38 µm at fill 0.223). The two still disagree on the mesh dimensions
  (hers default 45/18 µm, pitch 63; measured 50/18, pitch 68; ours 48/19, pitch 67).
- **Toolchain.** She uses our `setup_garfield.sh` idea (pinned Garfield `927e5c21`) and
  our observation that LCG_108/109 are too old. Her build script applies the local patch
  on checkout.
- **Documentation habits** worth borrowing: a `docs/results/mesh/README.md` that
  states what each number is, what review changed, and ships the CSVs next to the figures.
  Ours keeps the same information in HANDOFF/FRESH_LOOK and the runbook, harder to find.
- **Gas.** a Magboltz table in `stage_b/gas_tables/` (dry, 2026-10-08), gas notes in
  `docs/research/GAS_NOTES.md`. Our `response/params/wet_*` is the equivalent.

## 8. Suggested next steps (cheapest first)

1. Keep the `NULLCOLL` switch now in `gates_check.C` and rerun G7 with more statistics
   (10 k electrons), both maps, and the status breakdown (wire vs other).
2. Rerun `transparency_curve.csv` (E_drift ladder) with null-collision steps. The
   current curve (0.99 at 50–200 V/cm, 0.84 at 500, 0.55 at 1000) comes from the old
   stepping, so the fall-off will move.
3. One avalanche HV point of `mx17_aval_calib.py` with null-collision steps.
4. Try the no-coordinate-file export (`mesh` + `data` export nodes) on `WovenMesh.java`
   and load the result in `ComponentComsol`. If it works, no grid or flag systematic.
5. Send her the surface-charge top boundary and the null-collision result for her
   P2 mixtures at our ratio, and ask about the calendering.
6. After 1–3, correct `DEFAULT_TRANSPARENCY`, the HANDOFF/README weave numbers and
   FIELD_MAP_RUNBOOK's accepted values.

## Appendix: the experiment

`gates_check.C` with `NULLCOLL=1` set (adds `aval.EnableNullCollisionSteps(true)` to the G7
`AvalancheMicroscopic`), ROOT 6.40 with the pinned Garfield, 3000 electrons, gate runs
`gates_check("<map>", 3000, 333.0, <mirror>)`, all four runs in parallel:

| map | NULLCOLL | G7 | funnel rms |
|---|---|---|---|
| `meshfield_woven_production.txt` (mirror) | off | 0.928 ± 0.005 | 95.6 µm |
| `meshfield_woven_production.txt` (mirror) | on | 0.988 ± 0.002 | 94.2 µm |
| `meshfield_production.txt` | off | 0.947 ± 0.004 | 97.4 µm |
| `meshfield_production.txt` | on | 0.988 ± 0.002 | 93.3 µm |

The "on" runs print FAIL only because the old catastrophe band for G7 was (0.55, 0.98); the
band is now (0.55, 0.995) when `NULLCOLL` is set. 1.000 would still mean the region flag
was lost.
