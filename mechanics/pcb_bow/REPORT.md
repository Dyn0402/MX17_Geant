# Readout-board bow under gas overpressure: a mechanical model vs the bench dish

2026-10-28. Slide note: https://dylan-neff.web.cern.ch/notes/mx17-board-mechanics.html
(generator `make_mech_deck.py`). Every table is in `results/compare.md` (generated); figures are in `figures/`.

## Question

The bench data (nTof_x17 `rerun/model_audit/followups/fb_bow_maps/`) show the mesh centre of every chamber 1.2–2.3 mm
below a point 150 mm out:
- the report gives 1.5–2.3 ± 0.2–0.4 mm;
- the skeptic review's pivot estimator gives 1.2–1.9 mm with a ±0.4 mm systematic per chamber;
- the cathode is roughly flat.

The chambers run a few mbar over atmosphere, which visibly bulges the 60 µm window. The questions:
1. Does the same pressure bow the readout PCB by the measured amount?
2. How much does the Rohacell taped behind the board help, given doubts about how well it sticks?
3. Does mesh tension (stated as 5–10 N/m) add to it?
4. What are the true edge constraints of the board, which is larger than the active area and sealed by an O-ring on
   the Al frame?

## Answer

**Pressure gives a third of the dish at most.** A large-deflection two-layer plate model, with the mounting's own edge
stiffness, gives for a 1.6 mm PCB:

| Δp | Rohacell bonded | Rohacell loose |
|---|---|---|
| 3 mbar | 0.25 mm | 0.44 mm |
| 5 mbar | 0.38 mm | 0.59 mm |

These are bench-equivalent depths, comparable with the measured 1.2–2.3 mm.
- Across PCB 1.2–2.4 mm, every bond pattern, every edge model and mesh tension, nothing reaches 1.2 mm below about
  15 mbar (16 mbar for the softest realistic board).
- The measured 1.7 mm would need about 50 mbar (1.6 mm board, loose Rohacell). The window would then bulge about
  22 mm, against about 8.5 mm at 3 mbar.
- The measured profile shape is not a pressure shape either.

**Rohacell.** Glued full face, with the back mylar as the second skin, it cuts the bow by 35 %. Loose, it does almost
nothing (−3 %). Where it is stuck matters far more than how strongly:
- bonding the outer 60 mm captures 78 % of the full-glue gain;
- a centre that has come unstuck out to r = 100 mm costs nothing (98 % kept);
- even a soft full-face tape captures ≥ 85 %.

**Mesh tension** at 5–10 N/m adds less than 1 µm. Even at 10 N/cm it adds only 0.02–0.04 mm with the real edge.

**Edges.** The board's bending span is the 402 mm support-plate aperture. The O-ring (seal at 208.4 mm) and the frame
face sit just outside it and pinch the board onto the foam, so the edge is 96–97 % of a clamp. If the frame rides on
the O-ring rather than bottoming on the board, that drops to 89–95 % at 10 mbar. The exact edge moves the dish by only
5–10 %.

## 1. Geometry: what holds the board (from the STEP file)

`MX17 2026 v1.3 05 06 2026.stp`, global placements from gmsh's OCC importer (bounding boxes per labelled volume) plus a
direct parse of the `joint` product. Distances are from the active-area axis.

| mm | part |
|---|---|
| 0 – 201 | free span: the 8 mm Al support plate has a 402 × 402 mm through-aperture; only the Rohacell is behind the board |
| 201 – 205 | board on the Rohacell over the plate, gas above (inside the 410 mm frame aperture) |
| 205 – 206.55 | inner lip of the frame's lower face |
| 206.55 – 210.25 | O-ring groove; O-ring ("joint", 3.7 mm wide, 6 mm tall in the CAD, overlapping the uncompressed board) seals at ≈ 208.4 |
| 210.25 – 220 | frame face on the board; screws (Ø ~7.6 mm heads, 40 mm) at 214 mm on a 55 mm pitch through flange, frame and board |
| 220 – 250 | board overhang (470 mm board, Rohacell and plate offset +15, +15 mm) |

A second O-ring seals the window flange to the top of the frame (z 299–305 mm).

The screws end 1.5 mm into the Rohacell in the CAD (tips at z = 267.7, plate top at 266.2), so the drawing does not
show what they thread into. The model needs only that the frame bears on the board.

The board is larger than the active area (470 vs 399 mm), but outside 201 mm it is pressed onto the foam and plate,
and beyond the O-ring it is pinched under the frame face. So the span that bends is the plate aperture, not the board,
the frame aperture (410) or the O-ring line (417).

## 2. Models

### 2.1 Plate (`plate.py`)

Large-deflection (von Kármán) two-layer plate over the aperture.

- **Layers.** The PCB (layer 1) and the Rohacell (layer 2, optionally with the 25 µm mylar on its far face as part of
  it) share the deflection w. Each has its own in-plane displacement at its own neutral surface.
- **Energy:**
  - bending ½ΣD_i[w_xx² + w_yy² + 2νw_xx w_yy + 2(1−ν)w_xy²];
  - membrane ½Σε_iᵀA_iε_i with ε = sym∇u_i + ½∇w∇wᵀ;
  - interface ½k(x, y)(s_x² + s_y²), slip s = u₁ − u₂ + d∇w, with d the distance between the neutral surfaces;
  - loads −q w − m∇²w (gas plus weight, and the mesh-tension moment).
- **Edges.** w = 0, with rotation clamped, simply supported, an elastic spring K, or an imposed rotation. In-plane
  immovable (bolted) or free.
- **Discretisation.** Bending by finite differences (w_xy at cell centres, ghost nodes at the edge); membrane on Q4
  cells (2×2 Gauss); interface slip with one-point quadrature (the 2×2 rule locks as k → ∞). Newton iteration with
  pressure continuation; 41×41 grid in the scans.
- **Bench observable.** The surface is fitted with the bench basis [1, U, W, U²/2, W²/2, UW] (U, W in 100 mm) over
  the efficient area 399 × 360 mm, and the depth is 0.5·k̄·1.5², exactly as `fb_bow_maps/summary.py` computes it.
  For an edge-held board this is about half the true dip w(150) − w(0), because the shape is not a parabola.

### 2.2 Edge (`edge1d.py`)

The 10–20 mm clamping zone is below the 2-D grid (10 mm cells), so a 1-D cross-section in cylindrical bending at
0.25 mm resolution models it, from the centre to the overhang end:
- Rohacell foundation outside 201 mm: k_f = E/t = 1.4×10¹⁰ N/m³, compression only (the tape cannot pull the board onto
  the plate), optional precompression.
- Frame face as a rigid one-sided contact over 205–206.55 and 210.25–220 mm ("bottomed").
- O-ring sealing reaction (0 or 2000 N/m) at 208.4.
- Gas inside the seal.
- Alternative: the frame rides on the O-ring as a preloaded line spring (10⁶–10⁷ N/m²).

It returns K_eff = M/θ at the aperture edge, which the 2-D model uses as its spring edge ("real edge").

### 2.3 Window (`window_gauge.py`)

The same solver in its membrane limit: 60 µm PET, E 3.8–5.0 GPa, ν 0.38, on the 400.2 mm flange aperture, with no
pre-tension.

### 2.4 Materials and stack

| part | values |
|---|---|
| PCB | structural laminate t (FR4 + 5 Cu layers + 50 µm coverlay); nominal 1.6 mm, scanned 1.2–2.4. The CAD's 1.70 mm counts the mesh and the 150 µm gap |
| FR4 | 20 GPa (17–24) |
| Cu | 117 GPa, smeared by in-window coverage (0, 0.76, 0.458, 0.457, 0.53) |
| Rohacell 51 | E 70 MPa, 52 kg/m³ (31 IG at 36 MPa as a variant) |
| back mylar | 25 µm, 4.5 GPa |
| bench load | board weight, cathode up: ≈ 0.35 mbar equivalent |

### 2.5 Validation (`validate.py`, `results/validate.json`)

| check | model | reference |
|---|---|---|
| linear α = w₀D/(qa⁴), simply supported | 0.004062 | 0.00406 (Timoshenko, +0.05 %) |
| linear α, clamped | 0.001268 | 0.00126 (+0.6 % at N = 61, +1 % at N = 41) |
| spring edge, K → 0 / ∞ | SS / clamped | — |
| interface k = 0 | 10.29 N·m | D₁ + D₂ = 10.39 |
| interface k = 10¹² N/m³ | 18.13 N·m | composite section 18.26 |
| SS plate under edge moment m = 1 N·m/m | 1.3228 mm | exact (Poisson / Prandtl constant) 1.3231 mm |
| membrane limit (0.2 mm plate, 0.5–10 mbar) | — | Vlassak–Nix law, within 0.6–1.2 % |
| window solver | — | Vlassak–Nix, within 0.1 mm |

## 3. Results

### 3.1 Edge stiffness of the real mounting (`results/edge1d.json`, compare.md §7)

| PCB | K_eff, frame bottomed (N·m/m) | K·a/D | 1-D centre, toward pinned | frame on O-ring, 10 mbar |
|---|---|---|---|---|
| 1.3 mm | 700–770 | 70–77 | 3 % | K 370, 5 % |
| 1.6 mm | 1110–1200 | 60–64 | 3 % | K 460, 8 % |
| 2.0 mm | 1800–1880 | 50–52 | 4 % | K 590, 11 % |

- The result is almost independent of pressure and of the O-ring sealing force. A 20 µm foam precompression stiffens
  it by 40 %; a softer foam (half k_f) softens it by 15 %.
- The board lifts by less than 40 µm anywhere outside the aperture. With the frame riding on the O-ring at 10 mbar, the
  rim pivots about the seal and the overhang lifts about 0.5 mm.
- The first version of this study used K = E b³/12t ≈ 3900 N·m/m (a 15 mm foam pinch). That is 3.5× too stiff; the
  dish moves by +3–7 %.

### 3.2 Bow vs pressure, thickness and backing (real edge; `figures/pcb_and_window.png`, slide 7)

Bench depth (mm) at 5 / 10 mbar:

| PCB | Rohacell loose | Rohacell + skin bonded |
|---|---|---|
| 1.3 mm | 0.72 / 0.99 | 0.47 / 0.75 |
| 1.6 mm | 0.59 / 0.86 | 0.38 / 0.64 |
| 2.0 mm | 0.43 / 0.68 | 0.28 / 0.49 |

- **True centre deflection** is 2–2.5× larger (1.6 mm, loose: 1.55 mm at 5 mbar, 2.2 mm at 10 mbar). The bench quotes
  the fitted depth, so that is the number to compare.
- **Weight alone** gives 0.03–0.08 mm. At n_TOF (vertical) drop that term.
- **Response is sub-linear,** ∝ Δp^0.5–0.6 at 5–10 mbar, because membrane stretching takes over once w ≈ t. The
  earlier linear estimate ("0.3–0.7 mm per mbar", `physmodel.py` and `08b_bow.tex`) over-predicts by 2–3× and compares
  the true dip.
- **Other variants:**
  - FR4 17 / 24 GPa: +7 / −9 %.
  - Rohacell 31: +3 % loose, +13 % bonded.
  - Ideal clamped vs simply supported (loose, 5 mbar): 0.55 vs 0.76 mm.
  - In-plane free edges: +10–35 %.

### 3.3 Rohacell: how and where it sticks (slides 8–9; compare.md §2)

Share of the full-glue stiffening captured at 5 mbar (1.6 mm PCB), from loose (0.57 mm) to glued (0.38 mm):

| bond | captured | depth at 5 mbar |
|---|---|---|
| loose (contact only; also tape joining the foam's edge to the board outside the aperture) | 0 % | 0.57 mm |
| 25 mm band at the aperture edge only | 33 % | 0.50 |
| soft foam tape, full face (k = 10⁸ N/m³) | 40 % | 0.49 |
| central r < 100 mm patch only | 45 % | 0.48 |
| tape strips 25 mm wide every 134 mm | 68 % | 0.44 |
| 60 mm band at the aperture edge only | 78 % | 0.42 |
| bonded, centre unstuck to r = 150 mm | 81 % | 0.41 |
| thin acrylic PSA, full face (k = 10⁹) | 85 % | 0.41 |
| tape strips 50 mm wide every 100 mm | 92 % | 0.39 |
| stiff PSA, full face (k = 10¹⁰) | 97 % | 0.38 |
| bonded, centre unstuck to r = 100 mm | 98 % | 0.38 |
| glued, full face | 100 % | 0.38 |

**Why the rim matters.** The foam stiffens the board as a sandwich core between the PCB and the back mylar, carrying
shear between them. In a uniformly loaded plate that shear is largest near the supports and zero at the centre. So a
bond near the edge does the work and a debonded centre costs little. The foam is soft in-plane (70 MPa), so even a weak
adhesive can drag it along: grip barely matters, coverage does. Without the bonded mylar skin the gain is about half.

**If re-taping:** a continuous band of tape (or transfer adhesive) over at least the outer 60 mm inside the aperture,
all round, with the back mylar bonded to the foam.

The patterns were run with the 3900 N·m/m edge; the fractions do not depend on the edge.

### 3.4 Window as a pressure gauge (slide 6, `figures/gauge.png`)

Window centre sag (E = 4.5 GPa; the 3.8–5.0 GPa band is ±0.4 mm at 3 mbar):

| Δp (mbar) | 1 | 3 | 5 | 10 | 20 | 40 | 100 |
|---|---|---|---|---|---|---|---|
| sag (mm) | 5.9 | 8.5 | 10.1 | 12.7 | 16.0 | 20.2 | 27.4 |

- Sag ∝ Δp^⅓, so ±1 mm on a ~10 mm sag gives Δp to ±30 % and the board prediction to ±15 %.
- A pre-tensioned window sags less, so a measured sag means at least the tabulated Δp.
- `gauge.png` plots the expected board dish directly against the window sag.
- Overpressure for the measured dish from pressure alone (real edge):

  | scenario | 1.2 mm | 1.7 mm | 2.3 mm |
  |---|---|---|---|
  | 1.3 mm PCB, loose (softest realistic) | ≈ 16 mbar | ≈ 42 mbar | > 100 mbar |
  | 1.6 mm, loose | ≈ 21 mbar | ≈ 53 mbar | > 100 mbar |
  | 1.6 mm, bonded | ≈ 30 mbar | ≈ 63 mbar | > 100 mbar |
  | 2.0 mm, bonded | ≈ 40 mbar | ≈ 82 mbar | > 100 mbar |

### 3.5 Mesh tension (slide 10, compare.md §7)

**Bulk mesh.** The mesh is laminated on the board (bulk pillars). Its tension T is balanced by an equal compression
in the board, so the pair carries no net in-plane force and no stress stiffening or softening. Only the moment
m = T·e about the board's neutral surface is left, with e ≈ 0.8 mm. It cups the board in the measured sense.

| T | real edge, 0 / 3 / 10 mbar | simply supported edge |
|---|---|---|
| 10 N/m (the stated 5–10 N/m) | +0.3 µm | — |
| 10 N/cm (1000 N/m) | +0.04 / +0.02 / +0.01 mm | +0.23 mm at 0 mbar, +0.05 at 3 |

The moment and the pressure add almost linearly. A perfectly clamped plate suppresses a uniform moment completely.

**Floating mesh** (not bulk), for completeness. The mesh's own tension needs T·κ ≤ ε₀E²/2 ≈ 48 Pa to keep it on the
pillars. At 10 N/m that is about 0.6 Pa, so the mesh follows the board. At 10 N/cm the mesh would stop following
beyond about 0.5 mm of dish.

### 3.6 Shape test (slide 11, compare.md §5)

Each view's measured w0 profile is fitted as w0 = c + α·s_model(u). Σχ² over 9 views (63 dof):

| model shape | Σχ² |
|---|---|
| pressure, real edge | 1431 |
| pressure, clamped | 1626 |
| parabola | 1065 |
| pressure, SS | 924 |
| mesh-tension moment, SS edge | 848 |
| imposed edge rotation | 864 |

- An edge-held plate's slope turns over before the rim. The data steepen toward it (det2 Y, det3 X/Y, det4 X/Y,
  det6 Y). None of the shapes fits.
- Part of that is reconstruction near the strip ends. Read as geometry, it points to an edge load or a built-in shape,
  not a uniform pressure.
- The w0 amplitude with a flat cathode implies a dish of 2.8–3.0 mm.

## 4. What it means, and checks

1. **Pressure explains 0.25–0.6 mm of the 1.2–2.3 mm dish at 3–5 mbar.** Bond quality moves this by ±25 %, the edge
   by 5–10 %, the board thickness by about ±30 %, and mesh tension not at all.
2. **Not ruled out:**
   - The bench overpressure is really 20–60 mbar.
   - Built-in board warp: ordinary fab tolerance 0.75 % ≈ 3 mm over 400 mm. It would survive vertical mounting at
     n_TOF.
   - An edge load from the mounting (frame flatness, screw torque, cables, M1 cards); the rim-steepening shape points
     there.
   - The bench depth is over-read: the ±0.4 mm reconstruction systematic and the ≤ 0.6 mm M3 bow put the low end at
     0.8–1.2 mm.
3. **Checks:**
   - Window sag with a straightedge while flowing: about 8–10 mm means 3–5 mbar; about 20 mm means about 40 mbar.
   - Double the overpressure and re-measure the w0 route: a pressure dish grows ×1.5, a warp does not move.
   - Straightedge or dial gauge on an open, unpressurised board.
   - The fab stack-up, and how the Rohacell and its mylar are taped.

**Open hardware questions:**
- the bench overpressure (bubbler or outlet);
- any observed window sag;
- the PCB stack-up;
- the Rohacell taping (full face, strips or edge), and whether the back mylar is bonded to the foam;
- whether the frame bottoms on the board or rides on the O-ring.

## 5. Files and how to rerun

All in `MX17_Geant/mechanics/pcb_bow/`; Python from `~/PycharmProjects/nTof_x17/.venv`.

| file | what |
|---|---|
| `plate.py` | 2-D two-layer large-deflection plate solver; `fit_depth` = the bench observable |
| `validate.py` | closed-form checks → `results/validate.json` |
| `scenarios.py` | 60-scenario scan, groups A–F (thickness × backing, bond patterns, edges, materials, vertical, non-pressure shapes) → `results/scenarios.json`, `maps.npz`; per-scenario cache in `results/cache/` (git-ignored); `NPROC` workers, ≈ 45 min on 3 |
| `edge1d.py` | 1-D model of the real mounting → `results/edge1d.json` |
| `real_edge.py` | 2-D with the real-edge K and mesh tension (groups G, M) → `results/real_edge.json` |
| `high_pressure.py` | sweeps to 100 mbar (H: earlier edges, HG: real edge) → `results/high_pressure.json` |
| `window_gauge.py` | window sag vs pressure → `results/window_gauge.json` |
| `compare.py` | tables (`results/compare.md`, `compare.json`) and figures; reads nTof_x17 `fb_bow_maps/results` |
| `make_mech_deck.py` | the slide note → `deck/mx17-board-mechanics.html` |

Order: `validate.py`, `scenarios.py`, `edge1d.py`, `real_edge.py`, `high_pressure.py`, `window_gauge.py`,
`compare.py`, `make_mech_deck.py`. Publish with
`python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py deck/mx17-board-mechanics.html --slug mx17-board-mechanics --force --deploy`.

## 6. Limits

- Kirchhoff kinematics per layer: no Rohacell core shear (estimated < 1 % at this span).
- Grid ≈ 1 %.
- The board stays in contact with the foam over the aperture (true under pressure; with the weight alone and a loose
  foam they could separate, a ≤ 0.05 mm effect).
- Isotropic FR4 and a square aperture.
- No initial warp: the model gives the change from the unloaded shape.
- The window has no pre-tension.
- The 1-D edge model is linear in bending, nonlinear only in contact; its K is used as a linear spring in 2-D.
- The bond-pattern scan used the earlier edge stiffness.
