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

## 4. 🟡 Amplification gap: 150 µm or 128 µm

**Assume now:** 150 µm — both sims' historical value, and P2 confirmed 150 for
their build. A bulk micromegas is often 128 µm, and the CAD folds the gap
invisibly into the 1.70 mm readout body, so the CAD cannot settle it.

**Why it matters:** ±15 % on amplification-gap path length; secondary for
material budget, primary for any gain modelling done downstream of Geant4.

## 5. 🟡 Which side the window Al and cathode Cu face

**Assume now:** both metals face the drift gas (mylar → Al and kapton → Cu in
beam order). This is the physically sensible reading — P2 confirmed the
analogous question for their cathode ("Al on gas side") — but it is untraced
for MX17.

**Why it matters:** sub-µm layers; irrelevant to energy deposit, marginally
relevant to very-low-energy electron emission off the surfaces.

## 6. ⬜ Readout PCB internal stackup

**Assume now:** the historical 50 µm kapton + 4 × 26 µm Cu, with the FR4
filler solved (4 × 314.5 µm) so that mesh + amp + paste + laminate = the CAD's
1.70 mm single body. In-plane, the copper is modelled as full sheets, although
the gerbers show pads + X/Y strips on 3 functional layers (0.68 mm pads on
0.78 mm pitch, 512 × 512, 399.36 mm active area).

**Why it matters:** total board mass is now pinned by the CAD (1.70 mm), but
the Cu fraction inside it is a guess. Gerber-measured area coverage (as p2 did
with `analyze_cu_coverage.py`) would turn the 4 full sheets into
density-scaled ones without new hardware input. Good candidate for the next
model iteration.

## 7. 🟡 Per-arm orientation of the module asymmetries (Full_Geant)

**Assume now:** every arm has the module's +x/+y connector edges (the
(+15,+15) plate offset and the M1 card edges) oriented the same way in its
local frame.

**Why it matters:** only for peripheral scattering in the 4-arm sim; the
active regions are unaffected. The photos/survey of the actual arm mounting
would settle it. Note also two Full_Geant-only deviations, commented in its
`DetectorConstruction.cc`: window bulge off (square terraced-dome corners
would falsely hit the neighbour arm's flange) and M1 cards off (CAD ring
radius interpenetrates the neighbour by ~1.1 mm at the surveyed arm
distances).

## 8. 🟡 The 0.1 mm CAD seam at the frame/readout interface

The CAD places the readout-board top 0.100 mm *inside* the gas-frame z-range
while the drift gap is nominally "exactly 30.000 mm". The model keeps
kapton → mesh = 30.000 mm (frame = spacer; the 9 µm Cu cladding eats into the
gas), which makes the model's window-face → readout-top 35.11 mm vs the CAD's
35.01 mm. Sub-0.1 mm bookkeeping — flagged so nobody chases it as a bug.
