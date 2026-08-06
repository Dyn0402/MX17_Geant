# MX17 as-built geometry vs. the Geant4 models

Derived from the mechanical CAD on 2026-08-06. This is a reference for updating
`src/DetectorConstruction.cc`; it states what the CAD says, what the two
simulations said at the time, and where they disagreed.

> **Status 2026-08-06 (later the same day):** the as-built geometry below is
> now implemented in `shared/MX17ModuleGeometry.hh` (consumed by both sims);
> the ❌ marks in the side-by-side table describe the *pre-upgrade* models and
> are kept as the historical record. See
> [`GEOMETRY_IMPLEMENTATION_NOTES.md`](GEOMETRY_IMPLEMENTATION_NOTES.md) for
> what was done and [`NEEDED_INPUTS.md`](NEEDED_INPUTS.md) for the remaining
> assumptions. Model figures: `design/figures/mx17_*.png`
> (`scripts/model/plot_mx17_model.py`).
>
> **Correction:** discrepancy 1 below originally claimed the 8 mm support
> plate is solid. That was a parser error — see the rewritten entry.

## Sources

| Source | Location | Notes |
|---|---|---|
| Full detector CAD | `~/x17/mx17_cad/MX17 2026 v1.3 05 06 2026.stp` | 14 MB, Solid Edge → STEP AP214, exported 2026-06-05. 23 products, 72 placed instances. **One** MX17 module. |
| PCB fabrication data | `design/gerbers/` (in this repo) | Readout gerbers; `DFS3498A_gasframe.gbr` independently confirms the frame aperture. |
| LS vessel | `~/x17/LS_Stuff/LS X17.step` | Shapr3D, separate file. Already modelled in `MX17_Full_Geant`. |
| Scintillator clamps | `~/x17/mx17_cad/assembly of scintillators.stp` | 4× `pince scintillateur` only; its `LS X17` children are empty placeholders. |
| Trigger paddle / PMT | `~/x17/mx17_cad/paletta_e_supportoPM.stp` | `paletta` 200×24×492, `supporto_XP2982_v2` (Photonis XP2982). Not part of the MM stack. |

## As-built stack

Along the beam axis, entrance first. Dimensions in mm, measured from the solids.

| Z range | Part | Transverse | Thickness |
|---|---|---|---|
| 307.897–307.957 | Mylar window | 440 × 440 | 0.060 |
| 302.897–307.897 | **window flange** (Al ring) | 440 outer, 400.2 aperture | 5.000 |
| 302.847–302.897 | Kapton drift cathode | 406 × 406 | 0.050 |
| 272.847–302.847 | **Alu gas flange** (the frame) | 440 outer, **410.0 aperture**, 15 mm wall | **30.000** |
| 272.847–304.447 | 4× drift field PCB (field cage) | 408.4 × 1.62, forming a 410 mm ring | 31.6 tall |
| 271.247–272.947 | PCB readout | 470 × 470 | 1.700 |
| 266.247–271.247 | rohacell | 470 × 470 | 5.000 |
| 258.247–266.247 | **Alu detector support** | 470 × 470, **402 × 402 through-aperture** | 8.000 |

Also present, outside the beam path: 4× `multi M1 PCB` (41 × 180 × 6.6) forming a
446 mm ring around the readout, the SHV/HV housing and covers, gas entry fittings,
geometer targets and fasteners.

Key derived numbers:

- **Drift gap = 30.000 mm exactly** — the Al frame *is* the spacer; its faces
  coincide with the readout PCB top and the Kapton cathode.
- **Frame aperture 410.0 mm** — cross-checked against `DFS3498A_gasframe.gbr`,
  whose routing outline runs ±205 mm. Two independent files agree.
- Window face to readout surface: **35.010 mm** (sim: 30.379 mm).
- Total module depth, Mylar front to Al support back: **49.710 mm** (sim: 35.983 mm).
- Material behind the gas, **on axis**: **6.70 mm** (1.70 PCB + 5.00 rohacell;
  the 8 mm Al plate has a 402 mm aperture clearing the active area) vs
  **5.604 mm** in both sims. Outside the aperture the plate adds its 8 mm.

## Side-by-side

`MX17_Geant` = this repo (`src/DetectorConstruction.cc`).
`MX17_Full_Geant` = `../MX17_Full_Geant/src/DetectorConstruction.cc`.
The two use **identical** MM and PCB layer constants; they differ only in
transverse face size and in that Full_Geant replicates the module across 4 arms.

| # | Feature | As-built (CAD) | MX17_Geant | MX17_Full_Geant |
|---|---|---|---|---|
| 1 | Mylar window | 60 µm, 440×440 | ❌ 40 µm | ❌ 40 µm |
| 2 | Window Al metallization | 🔶 absent from CAD | ✅ 0.1 µm | ✅ 0.1 µm |
| 3 | Window flange (Al ring, 5 mm) | present | ❌ absent | ❌ absent |
| 4 | Gas gap, window → cathode | 5.0 mm | ❌ absent | ❌ absent |
| 5 | Kapton drift cathode | 50 µm, 406×406 | ✅ 50 µm | ✅ 50 µm |
| 6 | Cu cathode cladding | 🔶 absent from CAD | ✅ 9 µm | ✅ 9 µm |
| 7 | **Drift gap** | **30.000 mm** | ✅ 30 mm | ✅ 30 mm |
| 8 | **Aluminium gas frame** | 440 out / 410 ap / 30 mm | ❌ absent | ❌ absent |
| 9 | Drift field-cage PCBs | 4× 408.4×1.62×31.6 | ❌ absent | ❌ absent |
| 10 | Micromesh | ⚠️ inside 1.7 mm board | ⚠️ 30 µm slab | ⚠️ 30 µm slab |
| 11 | Amplification gap | ⚠️ inside 1.7 mm board | ⚠️ 150 µm | ⚠️ 150 µm |
| 12 | Resistive paste | ⚠️ inside 1.7 mm board | ⚠️ 100 µm | ⚠️ 100 µm |
| 13 | Readout PCB body | 1.70 mm, 470×470 | ❌ 0.554 mm laminate | ❌ 0.554 mm laminate |
| 14 | Rohacell | 5.0 mm, 470×470 | ✅ 5 mm | ✅ 5 mm |
| 15 | **Backing plate** | **8 mm Al ring, 402 mm aperture** | ❌ 50 µm Al foil | ❌ 50 µm Al foil |
| 16 | Front-end (multi M1) PCBs | 4× 41×180×6.6 | ❌ absent | ❌ absent |
| 17 | Transverse face | 440 (window) / 470 (PCB) | ⚠️ 400 × 400 | ⚠️ 380 × 340 |
| 18 | Copper readout segmentation | pads + X/Y strips | ❌ solid Cu sheets | ❌ solid Cu sheets |
| 19 | Material behind gas | 14.70 mm | ❌ 5.604 mm | ❌ 5.604 mm |
| 20 | Aluminized mylar behind rohacell | 🔶 absent from CAD | ❌ absent | ❌ absent |
| 21 | Front-window pressure bulge | n/a (CAD is unpressurized) | ❌ flat | ❌ flat |

✅ matches · ❌ disagrees · 🔶 real, but absent from the CAD · ⚠️ not determined
by the CAD, or approximated by choice

Rows 2, 6, 20 and 21 come from Dylan's knowledge of the hardware, not from the
files. The mechanical CAD models foils as plain solids and carries no
metallization or cladding, so its silence on rows 2 and 6 is not evidence of
absence — **the simulations are right there and must not be "corrected" toward
the CAD.** See [`GEOMETRY_IMPLEMENTATION_NOTES.md`](GEOMETRY_IMPLEMENTATION_NOTES.md) §2.

## Discrepancies, worst first

**1 — The 8 mm aluminium backing plate is missing (row 15).**
`Alu detector support modifié` is a 470 × 470 × 8 mm plate directly behind the
rohacell — **with a 402 × 402 mm square through-aperture** (plus ~100
bolt/fastener holes of r ≤ 2.5 mm at the periphery).

*Correction 2026-08-06:* this entry originally claimed the plate was "verified
solid". That was wrong — the first STEP parse missed the inner boundary loops
of the two Z-normal faces. Re-analysis of the face loops shows an inner
FACE_BOUND of 402 × 402 mm on both faces (a through-cut), offset (−15, −15) in
the part frame. Cross-registration with the readout gerbers
(`DFS3498A_activearea.gbr`: 399.36 mm active area flashed at the gerber
origin; copper extents centred (+15, +15)) shows the **aperture is concentric
with the active area** — 1.32 mm clearance per side — and it is the 470 mm
plates (readout PCB, rohacell, support plate) that sit offset (+15, +15) from
the active axis. Dylan confirmed the intent: the frame screws into the plate,
but there is no aluminium behind the active area.

So the plate does *not* block on-axis Sr-90 betas; it matters for the
periphery (frame-region scattering and anything outside the 402 mm square).
Both sims previously modelled this position as 50 µm of Al foil over the full
face — wrong in both directions at once.

**2 — Material behind the gas: 14.70 mm vs 5.604 mm (row 19).**
Even setting the Al plate aside, the readout board is a single 1.70 mm solid in
CAD where the sims build a 0.554 mm laminate of 4×(26 µm Cu + 100 µm FR4) plus
50 µm Kapton.

**3 — Missing 5 mm gas gap between window and cathode (row 4).**
The window flange is 5 mm thick with a 400.2 mm aperture; the Mylar sits on its
outer face, the Kapton on its inner face. The gap between them is open gas. Both
sims stack Mylar → Al → Kapton → Cu contiguously.

**4 — Mylar window 60 µm, not 40 µm (row 1).**

**5 — All lateral structure absent (rows 3, 8, 9, 16).**
The 30 mm Al frame ring, the 5 mm window flange ring, the field cage and the
front-end boards. Irrelevant for on-axis transmission; this is exactly the
material that governs peripheral scattering and frame-induced background.

## What the CAD does *not* settle

- **Materials.** STEP carries no material assignment. "Aluminium" for the frame,
  window flange and support plate comes from part names only.
- **Mesh, amplification gap, resistive paste.** Absent as separate solids — they
  are sub-mm bulk-process features folded into the 1.70 mm `PCB readout` body.
  The sims' separate layers are not contradicted; in reality they live inside
  that 1.7 mm. Do not add their thickness on top of the CAD numbers.
- **PCB internal stackup.** Copper and dielectric thicknesses are in neither the
  CAD nor the gerbers. The 26 µm / 100 µm values in the sims have no traced
  source. The gerbers do give in-plane segmentation: 0.68 mm pads on 0.78 mm
  pitch, 512 × 512 over a 399.36 mm active area, with 3 functional copper layers
  (pads, Y strips, X strips) — not the 4 uniform sheets both sims use.

## Method

Bounding boxes and placements were extracted with a purpose-written STEP parser.
Two traps, both of which produce plausible-looking wrong answers:

1. **Bounding boxes must use `VERTEX_POINT`s only.** Taking every
   `CARTESIAN_POINT` includes arc and fillet centres that lie far outside the
   solid — the gas frame measured 931 × 927 × 598 mm instead of 440 × 440 × 30.
2. **Assembly transform direction must be resolved, not assumed.** Match
   `REPRESENTATION_RELATIONSHIP.rep_1/rep_2` against the child product's own
   shape representation to decide which `ITEM_DEFINED_TRANSFORMATION` item is
   the child frame. With the convention reversed, the four field-cage PCBs
   scatter ~2 m apart instead of closing into their 410 mm square.

Sanity checks that passed: the four field-cage PCBs close into a 410.0 mm ring;
the four front-end boards close into a 446 mm ring; the frame aperture matches
the gerber gas-frame outline; all plate thicknesses come out as round numbers.
