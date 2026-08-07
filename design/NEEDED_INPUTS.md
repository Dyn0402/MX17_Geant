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

**Now measured (2026-08-06, rev 2).** The copper is modelled as the five
physical gerber layers (L3 guard ring, L4 pads, L5 Y strips, L6 X strips,
L7 fan-out; L8 carries only the outline stroke), each a 26 µm sheet
density-scaled by its measured area coverage, and each split into two zones —
inside and outside the 399.36 mm active window:

| layer | active | outside | board avg |
|---|---|---|---|
| L3 guard  | 0.0000 | 0.3385 | 0.0941 |
| L4 pads   | 0.7600 | 0.0568 | 0.5645 |
| L5 Y      | 0.4577 | 0.1189 | 0.3635 |
| L6 X      | 0.4568 | 0.1194 | 0.3630 |
| L7 fan-out| 0.5285 | 0.0584 | 0.3978 |

**How the FR4 filler is obtained (and the assumption buried in it).** The CAD
gives one number for the readout board: a 1.70 mm solid body, mesh front face
to laminate back face. The model does not know the real internal dielectric
split, so it treats the FR4 as the *residual*: it subtracts everything whose
thickness we do think we know — mesh, amp gap, ESL, kapton, and the five 26 µm
copper layers — and divides what is left equally over the five FR4 layers.

The consequence is that **changing any other layer silently changes the FR4**.
When the ESL went 100 → 10 µm, the residual grew by 90 µm and the FR4 filler
went 246.4 → 264.4 µm per layer, keeping the board at exactly 1.7000 mm.

That is a *choice*, and it may be the wrong one. It assumes the CAD 1.70 mm is
the reliable number and the laminate takes up the slack. The alternative — that
the dielectric stack is fixed and the real board is 1.61 mm — is equally
consistent with what we know. The two differ by 90 µm of FR4, which is not
nothing for backscatter off the board.

**Needed:** a fab stackup drawing, which would settle the per-layer dielectric
thicknesses and turn this residual into a real number. Until then, anyone
reading a board-body thickness out of this model should know it is pinned to
CAD by construction, not derived.

Two corrections went in with this revision, both of which changed the as-built
result (`--legacy-geometry` is unaffected and remains bit-identical):

1. **Rasterizer bias.** `analyze_cu_coverage.py` drew every aperture with
   PIL's endpoint-*inclusive* `rectangle`/`ellipse`, inflating each feature by
   one pixel per dimension — +16 % on the 0.68 mm pads at the old 0.05 mm/px
   default. The previous numbers (0.095/0.651/0.419/0.420/0.456) were
   therefore ~13 % high on every fine-featured layer; the L3 guard ring, being
   one large solid shape, was almost unaffected. Fixed, default raster raised
   to 0.02 mm/px, and `--selftest` now checks the rasterizer against the
   analytically exact pad grid (0.68²/0.78² = 0.76003) so it cannot regress.
2. **Radial smearing.** One board-wide coverage per layer put ~26 % too little
   copper in the active area and up to 10× too much outside it. Hence the
   two-zone build, which costs nothing.

**Still open:** the per-layer copper *thickness* (26 µm is the historical
guess; ½-oz 18 µm is the common fab default) and the kapton/dielectric split.
A fab stackup drawing would settle it — this is now the dominant uncertainty
on the board, well above the ~0.2 % residuals in the coverage measurement.

### Real copper pattern — now the DEFAULT (2026-08-06, Dylan)

The three signal layers are an exactly regular 512 × 512 grid on a 0.78 mm
pitch spanning ±199.29 mm, verified flash-by-flash against the gerbers. The
model builds their real geometry — 0.68 mm square pads (L4), Ø0.50 mm dots +
0.098 mm bus traces (L5/L6) — via nested `G4PVReplica`, at ~+3 % CPU and no
extra memory. `--homogenized-readout` falls back to density-scaled sheets.

Only the dot→bus connector stub is approximated: it exists in ~2/3 of the
cells (so it is not periodic and cannot go in a replica cell as-is), and is
therefore entered in every cell at the width that reproduces the measured
layer coverage exactly.

**L3 (guard ring) and L7 (fan-out) are NOT patterned** — their artwork is
irregular, so they stay as the two-zone density-scaled sheets above. That is
the remaining approximation in the board copper. L7 in particular is 0.53
covered inside the active area (vias/fan-out), so if board backscatter ever
turns out to matter at the percent level, L7 is the next thing to model.

## 6b. 🟡 M1 front-end cards

**Now gerber-anchored (2026-08-06, rev 2):** the M1 gerbers turn out to be
drawn in the readout-board frame at the as-mounted position — card outline
41 × 160 mm at radial 219.4..260.4, tangential centre +100 (edge-mill layer),
with the bottom pogo-pad field landing exactly on the board's connector
copper. The model places four cards (2 per connector edge, tangential ±100),
laminate flat on the board with **no standoff** (Dylan), six Cu layers at
26 µm density-scaled by card-window coverage (0.1405/0.1034/0.1078/0.1109/
0.1051/0.0362, re-measured after the rasterizer fix in §6 — the previous
0.147/0.110/0.114/0.117/0.111/0.041 were 4–12 % high) with FR4 filler.

**Still assumed:** bare-laminate thickness 1.6 mm (the CAD 6.6 mm envelope
includes the Mec8 output connectors soldered on top and the thin JZ pogo-pin
connector underneath — both omitted from the model; add simple blocks if the
M1 mass ever matters). The card inner edge is clipped 219.4 → 220.0 mm where
the plain-ring gas-frame model has its outer wall (the real frame must be
relieved there).

## 6c. 🟡 Resistive layer: ESL strips 550/250 µm — thickness still open

**Confirmed (Dylan 2026-08-06):** vertical ESL resistive strips, 550 µm wide
with 250 µm gaps (0.8 mm pitch — deliberately not the 0.78 mm pad pitch),
over 412 × 412 (the deposit boundary from `3498A_top-resist.gbr`; the strip
artwork itself is not in the gerber set — the strips are built from the
confirmed spec, not from a gerber).

**Thickness ⬜ 10 µm — STILL NEEDS CONFIRMING.** Modelled as 10 µm since
2026-08-07 (Dylan: "something like 10 µm"), replacing an earlier 100 µm guess.
It is the right *order* for a screen-printed thick-film resistive layer — the
ATLAS NSW resistive strips are ~15–20 µm, and 100 µm would be unusually thick
for a single print pass — but "right order" is not a measurement. **Wanted: the
ESL paste grade + fired thickness from the screen-printing house**, which would
also pin the sheet resistivity that `design/RESPONSE_SIM_PLAN.md` currently
scans over. A caliper measurement on a spare board would do.

Sensitivity: this sets both the ESL material budget and the depth of the
inter-strip gas grooves, and the grooves are scored as amplification gas — so
the groove contribution scales linearly with it (2.0 % of amp-gas edep at
10 µm; it was ~20 % at 100 µm).

**Now built as real strips (default).** 515 strips of 550 µm on the 0.8 mm
pitch, one `G4PVReplica` level inside a `ResistLayer` envelope. Coverage works
out to exactly 515 × 0.550 / 412 = 0.6875, i.e. the same total ESL mass as the
old density-scaled slab — but the 250 µm inter-strip grooves are now real
chamber gas rather than reduced-density paste.
`--homogenized-readout` restores the slab.

**Grooves are scored as amplification gas (Dylan, 2026-08-07).** The replica
cell is named `AmpGas`, so ionisation in the grooves counts toward the
amplification signal — right, because the grooves are only 10 µm deep at the
bottom of a 150 µm gap, so the amplification field runs through them
essentially undisturbed. The strips keep the exact name `ResistivePaste`, so
`edepResistPaste` is unchanged in meaning. Measured effect: the grooves add
**2.03 % of the amp-gas energy deposit**, matching the 2.04 % expected from
0.3125 × 10 µm / (150 + 3.125) µm.

> ⚠️ **Caveat for anything consuming `ClusterTree`.** The groove cells are
> thin, so Geant4 forces short steps there and each step writes its own row.
> Grooves are 2.0 % of amp-gas *energy* but **21.7 %** of amp-gas *rows*
> (mean 2.8 eV/row against 37.3 eV/row in the main gap). Sum `edep` or
> `nPrimary`; do **not** count rows as a proxy for ionisation. `nPrimary`
> itself is safe — `SteppingAction` carries the sub-W remainder
> probabilistically (`G4UniformRand() < frac`), so it does not lose
> ionisation to step segmentation.

**Still open:** nothing major on this layer now. Noted: the strip/pad pitch
mismatch (0.80 vs 0.78 mm) means a slowly beating moiré between strips and
pads; irrelevant for material budget, relevant for any charge-sharing model
downstream. The ESL sheet resistivity is still unknown (see
`design/RESPONSE_SIM_PLAN.md`, which scans it).

## 6d. ✅ Bulk pillars

Ø 0.6 mm pillars on an exact 85 × 85 grid at 4.68 mm pitch spanning
±196.56 mm (`3498A_bulk.gbr`), standing in the amplification gap. Modelled as
polyimide-coverlay cylinders placed as daughters of the AmpGas volume — the
scored gas is displaced (~1.3 % of the gap volume, locally dead spots).
Pillar material (coverlay ≈ kapton) is an approximation.

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
