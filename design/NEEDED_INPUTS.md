# Needed external inputs — the register of what the MX17 model assumes

**Written:** 2026-08-06, when the as-built geometry
(`shared/MX17ModuleGeometry.hh`) landed. One place for every number the module
description *assumes* because nobody has given us the measurement, ordered by
how much the answer would move a result. Each item states the current
assumption, why it matters, and what to ask for. Companion figure:
`design/figures/mx17_stack_status.png`.

Legend: ⬜ needed · 🟡 partially known · ✅ resolved.

---

## 1. ⬜ Micromesh weave specification

**Assume now:** the P2 weave — 19 µm wire, 48 µm opening → 38 µm
effective-density slab at fill 0.223 (areal mass of a plain weave). Chosen as
a placeholder on 2026-08-06 because the real MX17 spec was not to hand.

**Why it matters:** the mesh is directly above the amplification gap; its
areal mass sets electron backscatter/absorption at the gap boundary and the
X-ray conversion rate closest to the amp region. The previous solid-30 µm
model overstated the mass ×3.5; a wrong weave spec moves it by tens of
percent. Optical transparency (~51 % for the P2 weave) is *not* modelled —
only the material budget.

**Ask for:** the mesh designation on the MX17 build sheet (e.g. "48/19" or
wires-per-inch + wire Ø).

## 2. ⬜ Operating gas overpressure (window bulge)

**Assume now:** ~3 mbar → 8 mm sag from Hencky's membrane solution
(w₀ ≈ 0.662·a·(p·a/(E·t))^⅓, a = 200.1 mm half-span, E = 4.5 GPa, t = 60 µm).
Exposed as `--bulge-front <mm>`; the p^⅓ scaling means 1–10 mbar only spans
≈ 6–12 mm.

**Why it matters:** the bulge adds up to 8 mm of chamber gas on the beam path
ahead of the drift region and changes the window-to-mesh distance across the
face. A pressure reading, or an actual sag measurement with a ruler, pins it.

## 3. ⬜ Rohacell-back aluminized mylar thicknesses

**Assume now:** 25 µm mylar + 0.1 µm Al (Al side unknown; modelled facing
downstream), 470 mm face, replacing the old 50 µm bare-Al guess.

**Why it matters:** small on-axis effect (it sits behind the 1.70 mm board and
5 mm rohacell), but it is the last layer before the support plate and shows up
in albedo/backscatter for the Sr-90 runs. A caliper or datasheet number ends
the guess.

## 4. ✅ Amplification gap: 150 µm

**Confirmed by Dylan, 2026-08-06.** The 150 µm both sims carried is correct;
the 128 µm bulk-standard question is closed.

## 5. 🟡 Which side the window Al and cathode Cu face

**Assume now:** both metals face the drift gas (mylar → Al and kapton → Cu in
beam order). This is the physically sensible reading — P2 confirmed the
analogous question for their cathode ("Al on gas side") — but it is untraced
for MX17.

**Why it matters:** sub-µm layers; irrelevant to energy deposit, marginally
relevant to very-low-energy electron emission off the surfaces.

## 6. 🟡 Readout PCB internal stackup

**Now partly measured (2026-08-06).** The copper is modelled as the five
physical gerber layers (L3 guard ring, L4 pads, L5 Y strips, L6 X strips,
L7 fan-out; L8 carries only the outline stroke), each a 26 µm sheet
density-scaled by its measured area coverage over the 470 mm board face —
0.095 / 0.651 / 0.419 / 0.420 / 0.456 (`scripts/gerber/analyze_cu_coverage.py`,
0.05 mm/px raster; over the active area: 0.00 / 0.879 / 0.533 / 0.534 / 0.608).
FR4 filler (5 × 246.4 µm) still solves the CAD 1.70 mm body.

**Still open:** the per-layer copper *thickness* (26 µm is the historical
guess; ½-oz 18 µm is the common fab default) and the kapton/dielectric split.
A fab stackup drawing would settle it.

## 6b. 🟡 M1 front-end card internals

**Assume now (Dylan 2026-08-06):** each card is its 6.6 mm CAD envelope
(41 × 180 mm, flat on the drift-gas side of the board edge, straddling it
longways at ring 242.5 mm / tangential +90 / −93 mm), filled with the six
gerber-defined Cu layers (26 µm each, density-scaled by the panel-average
coverage 0.183/0.105/0.110/0.106/0.102/0.053) and FR4 for the remainder.

**Why it matters / caveat:** treating the whole 6.6 mm as solid laminate
means components and connectors are modelled as FR4-equivalent — that likely
*overstates* the card mass (a bare 6-layer card is ~1.6 mm). If the M1s ever
matter for a background estimate, split the envelope into laminate + connector
standoff + air. Peripheral material, outside the beam path.

## 7. 🟡 Per-arm orientation of the module asymmetries (Full_Geant)

**Assume now:** every arm has the module's +x/+y connector edges (the
(+15,+15) plate offset and the M1 card edges) oriented the same way in its
local frame.

**Why it matters:** only for peripheral scattering in the 4-arm sim; the
active regions are unaffected. The photos/survey of the actual arm mounting
would settle it. Note also two Full_Geant-only deviations, commented in its
`DetectorConstruction.cc`: window bulge off (square terraced-dome corners
would falsely hit the neighbour arm's flange) and M1 cards off (straddling
cards reach 263 mm from the active axis and interpenetrate the neighbour
arm's flange at the surveyed arm distances).

## 8. 🟡 The 0.1 mm CAD seam at the frame/readout interface

The CAD places the readout-board top 0.100 mm *inside* the gas-frame z-range
while the drift gap is nominally "exactly 30.000 mm". The model keeps
kapton → mesh = 30.000 mm (frame = spacer; the 9 µm Cu cladding eats into the
gas), which makes the model's window-face → readout-top 35.11 mm vs the CAD's
35.01 mm. Sub-0.1 mm bookkeeping — flagged so nobody chases it as a bug.
