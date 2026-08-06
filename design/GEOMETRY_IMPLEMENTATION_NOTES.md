# MX17 geometry — implementation brief

Companion to [`GEOMETRY_FROM_CAD.md`](GEOMETRY_FROM_CAD.md), which records the
as-built dimensions. This file is the **work order**: what to change, in what
order, and which numbers are still missing.

Written 2026-08-06. Items 2–5 come from Dylan's domain knowledge and correct or
extend what the CAD alone shows.

> **Status — implemented 2026-08-06.** All six tasks are done:
>
> - §1 merge: `shared/MX17ModuleGeometry.hh` is the single module description,
>   consumed by both repos (option A, relative-path include with a CMake
>   guard). Landed as a pure refactor first — fixed-seed runs of every mode in
>   both sims were verified **bit-identical** to the pre-refactor builds.
> - §2–§6 applied on top: 60 µm window, 5 mm flange gap, copper-clad cathode,
>   440/406/410/470 faces + rings + field cage + M1 cards, woven-mesh
>   effective-density slab (19/48 µm placeholder), 1.70 mm readout body,
>   aluminized-mylar back foil (25 µm assumed), 8 mm support plate, terraced
>   front-window bulge (8 mm sag ≈ Hencky at 3 mbar, `--bulge-front`).
>   `--legacy-geometry` restores the old stack.
> - The §6 "held back" item is resolved: the plate has a **402 mm aperture in
>   the STEP file** (the earlier "verified solid" reading was a parser error —
>   see `GEOMETRY_FROM_CAD.md`, discrepancy 1) so it is modelled as a ring and
>   does not stop on-axis betas.
> - Remaining assumptions moved to [`NEEDED_INPUTS.md`](NEEDED_INPUTS.md).
> - Python mirror + figures: `scripts/model/`, output in `design/figures/`.

## Task summary

| # | Task | Blocked on input? |
|---|---|---|
| 1 | Merge the two Geant4 detector descriptions into one source of truth | no |
| 2 | Correct the CAD reading: drift copper and window aluminization | no |
| 3 | Add the aluminized mylar behind the rohacell | **yes** — thicknesses |
| 4 | Port the gas-pressure bulge on the front mylar from `p2_geant` | **yes** — overpressure |
| 5 | Port the woven-micromesh model from `p2_geant` | **yes** — mesh spec |
| 6 | Apply the CAD corrections from `GEOMETRY_FROM_CAD.md` | partly — see §6 |

---

## 1. Merge the two models

`MX17_Geant/src/DetectorConstruction.cc:206-224` and
`MX17_Full_Geant/src/DetectorConstruction.cc:193-208` currently declare the
**same MM and PCB layer constants twice**, value for value — 40 µm Mylar,
0.1 µm Al, 50 µm Kapton, 9 µm Cu, 30 mm drift, 30 µm mesh, 150 µm amp,
100 µm resistive paste, then 50 µm Kapton + 4×(26 µm Cu + 100 µm FR4) +
5 mm rohacell + 50 µm Al. Every geometry fix below therefore has to be applied
twice, and the two will drift apart the first time someone forgets.

**Goal:** one description of the MX17 module, consumed by both.

What must stay configurable per consumer (these legitimately differ):

| Property | MX17_Geant | MX17_Full_Geant |
|---|---|---|
| Transverse face | 400 × 400 mm | 380 × 340 mm (`mm_size_u_cm`, `mm_size_v_cm`) |
| Instances | 1, on axis | 4 arms, pinwheel-shifted |
| Local axes | z along beam | (u, v, w) per arm |
| Modes | 6 `SimMode` variants | full-experiment only |

Suggested shape: a header-only `MX17ModuleGeometry.hh` exposing a thickness
struct plus a `BuildModule(...)` that takes half-width, half-height and a
placement lambda, so each caller keeps its own transverse size and placement
convention. Two ways to share it across the repos:

- **A shared header referenced by relative path** in the second repo's
  `CMakeLists.txt`. Simplest; fragile if the sibling checkout moves.
- **A copy kept in sync by script**, with the header carrying a checksum
  comment. Uglier, but survives independent checkouts.

Recommend A, with a CMake `if(NOT EXISTS ...)` guard giving a clear error.

**Do the merge as a pure refactor first.** Land it with the geometry numerically
unchanged, confirm both sims produce identical output to the current build, and
only then start applying §2–§6. Mixing a refactor with physics changes makes a
regression impossible to attribute.

---

## 2. Corrections to the CAD reading

Two places where **the CAD is incomplete and the simulations are right**. Do not
"fix" the sims toward the CAD here.

| Item | CAD shows | Reality | Action |
|---|---|---|---|
| Drift cathode copper | `Kapton drift`, 50 µm, bare | Kapton is **copper-clad** | Keep `tCuCath = 9 µm`. CAD omission. |
| Window aluminization | `Mylar window`, 60 µm, bare | Mylar is **aluminized** | Keep `tAlWin = 0.1 µm`. CAD omission. |

The mechanical CAD models foils as plain solids and simply does not carry
metallization or cladding. Rows 2 and 6 of the comparison table in
`GEOMETRY_FROM_CAD.md` have been updated accordingly.

Still to confirm: which **side** each metal layer faces. The current stacking
order (Mylar → Al → Kapton → Cu → drift gas) puts both metals on their
drift-facing side, which is the physically sensible reading, but it is an
assumption. `p2_geant` hit the same question and resolved it explicitly for the
P2 cathode — see `p2_geant/docs/P2_MODEL.md`, "Al on gas side **confirmed**".

The 60 µm Mylar thickness from the CAD **does** stand — that is a mechanical
dimension, and it supersedes the 40 µm in both sims.

---

## 3. New layer: aluminized mylar behind the rohacell

Missing from the CAD **and** from both simulations. There is a second aluminized
mylar foil on the **back face of the rohacell**.

Proposed placement in the stack, after `PCB_Rohacell` and before the 8 mm
aluminium support plate:

```
... PCB_Rohacell (5 mm)
    RohacellBack_Mylar   <- thickness TBC
    RohacellBack_Al      <- thickness TBC, side TBC
    AluSupportPlate (8 mm, see GEOMETRY_FROM_CAD.md §Discrepancies)
```

**Needed before implementing:**

- mylar thickness (µm)
- aluminization thickness (µm) — 0.1 µm if it matches the window
- which side the aluminium faces
- transverse extent — assume 470 × 470 to match the rohacell unless told otherwise

Note this replaces, conceptually, what both sims currently call `PCB_AlFoil`
(50 µm Al). That 50 µm slab appears to be a stand-in for exactly this foil,
guessed rather than measured. It is **not** the 8 mm support plate. Resolve
whether the 50 µm should become this aluminized mylar, with the 8 mm plate added
separately behind it.

---

## 4. Front-window bulge from gas overpressure

Reference implementation: `p2_geant/src/DetectorConstruction.cc:870-908`
(`BuildWindow` lambda), documented in `p2_geant/docs/P2_MODEL.md` under
"Window bulge model".

### The p2_geant approach

A **terraced dome**: `N = 6` stacked gas prisms whose outline shrinks about the
opening centroid on a spherical-cap profile, each step capped by a flat mylar
ring, the last step a full cap.

```cpp
for (int k = 0; k <= N; ++k) {
    const double a = k * M_PI / (2.0 * N);
    sig[k] = std::max(std::cos(a), 0.04);   // lateral shrink factor
    h[k]   = H * std::sin(a);               // height of step k
}
```

Step `k` places a gas prism of outline `ScaleAbout(opening, centroid, sig[k])`
and height `h[k+1] - h[k]`, then a mylar terrace on top of it — an annulus
(`outline[k]` minus `outline[k+1]`) for `k < N-1`, a full cap for the last.

The property that makes this work: **every vertical path through the dome
crosses exactly one mylar thickness.** A naive nested-shell dome would
double-count material near the rim, which is precisely where a bulge changes the
path length. Preserve this property in any reimplementation.

`sign = -1` builds the dome rising upstream (front window); `+1` downstream.
MX17 needs only the front window — the readout PCB is the downstream gas
boundary, so there is no back window to bulge.

### Adapting to MX17

The p2 code extrudes an arbitrary wedge polygon via a `Prism` helper
(`p2_geant/src/DetectorConstruction.cc:851`) because the P2 opening is a wedge.
**MX17's opening is a plain 410 × 410 mm square** — the aluminium gas frame
aperture. So the polygon machinery collapses: use `G4Box` for the gas steps and
`G4Box` minus `G4Box` (`G4SubtractionSolid`) for the mylar annuli, scaling the
half-widths by `sig[k]` about the centre. Simpler and cheaper than
`G4ExtrudedSolid`.

### Sag height — must be computed, not copied

`p2_geant` uses `p2_bulge_front_mm = 10.0` from a Hencky estimate at ~3 mbar,
with the comment that 1–10 mbar spans 7–15 mm. **Do not carry that number over.**
Hencky sag scales with the membrane radius, and MX17's 410 mm square aperture is
substantially larger than the P2 wedge opening, so the same overpressure gives a
materially larger sag. Recompute for MX17's aperture, foil thickness (60 µm) and
actual operating overpressure.

Expose it as a config knob (`mx17_bulge_front_mm`) so it can be scanned, exactly
as p2 does.

**Needed:** operating overpressure, and confirmation of the Mylar elastic
modulus / foil pre-tension assumptions if a Hencky estimate is used.

---

## 5. Micromesh as a woven-mesh effective-density slab

Reference implementation: `p2_geant/src/DetectorConstruction.cc:777-791`.

Both MX17 sims currently model the mesh as a **30 µm slab of solid stainless
steel at full density** — a mesh is mostly holes, so this overstates its mass.

The p2 model replaces it with a slab of density-scaled steel whose areal mass
matches the real weave:

```cpp
const G4double dWire     = wire_diameter;
const G4double meshPitch = dWire + opening;
const G4double tMesh     = 2.0 * dWire;               // weave height
const G4double meshFill  = M_PI * dWire / (4.0 * meshPitch);
matMesh = new G4Material("MeshSteelEff", matSteel->GetDensity() * meshFill, 1);
matMesh->AddMaterial(matSteel, 1.0);
```

The fill fraction comes from the steel volume per unit area of a plain weave
with two orthogonal wire sets, `π·d²/(2p)`, divided by the slab thickness `2d`.

For P2's 48 µm opening / 19 µm wire: pitch 67 µm → 38 µm slab at fill 0.223.

**Scale of the correction:** a 30 µm solid-steel slab carries ≈ 0.024 g/cm²;
a P2-like weave carries ≈ 0.0068 g/cm². If MX17's mesh is a comparable weave,
the current model overstates mesh material by roughly a factor of 3.5.

Caveat to carry over in the comment: this reproduces the **material budget
only**. Optical transparency — `(opening/pitch)²`, ≈ 51 % for the P2 weave — is
not geometrically modelled, so the mesh does not shadow anything.

**Needed:** the MX17 mesh spec (wire diameter and opening, or the weave
designation). Also worth confirming the amplification gap while asking: both
sims use 150 µm, a bulk Micromegas is often 128 µm, and the CAD folds it
invisibly into the 1.70 mm readout board.

---

## 6. Applying the CAD corrections

From `GEOMETRY_FROM_CAD.md`. Straightforward once §1 lands:

- Mylar window 40 → 60 µm
- Add the 5 mm gas gap between window and cathode (window flange aperture 400.2 mm)
- Readout PCB laminate 0.554 → 1.70 mm single body
- Transverse faces: 440 mm at the window, 470 mm at the PCB
- Lateral structure: 30 mm Al frame ring (440 outer / 410 aperture), 5 mm Al
  window-flange ring, 4× field-cage PCB, 4× front-end board

**Held back pending confirmation:** the 8 mm aluminium backing plate. It stops
Sr-90 betas outright, so adding it changes what `kSr90Calibration` and
`kFullExperiment` predict downstream. Confirm it is in the beam path for the
runs being simulated before implementing — see `GEOMETRY_FROM_CAD.md`,
discrepancy 1.

---

## Inputs still needed

| # | Input | Blocks | Resolution (2026-08-06) |
|---|---|---|---|
| 1 | Rohacell-back mylar: mylar µm, Al µm, which side | §3 | **assumed** 25 µm + 0.1 µm Al (Dylan: use as placeholder) |
| 2 | Is the 50 µm `PCB_AlFoil` meant to be that foil? | §3 | **yes** — replaced by the aluminized mylar |
| 3 | Operating gas overpressure | §4 | **assumed** ~3 mbar → 8 mm sag (`--bulge-front`) |
| 4 | Micromesh wire Ø and opening | §5 | **placeholder** P2 weave 19/48 µm |
| 5 | Amplification gap: 150 µm or 128 µm | §5 | open — 150 µm kept |
| 6 | Which side the window Al and cathode Cu face | §2 | open — drift-facing assumed |
| 7 | Is the 8 mm Al support plate in the beam path? | §6 | **resolved** — plate has a 402 mm aperture in the STEP; modelled as a ring |
| 8 | PCB internal stackup (Cu/dielectric thicknesses) | §6 | open — historical Cu kept, FR4 filler solved to the 1.70 mm body |

Live register: [`NEEDED_INPUTS.md`](NEEDED_INPUTS.md) (modeled on
`p2_geant/docs/NEEDED_INPUTS.md`).
