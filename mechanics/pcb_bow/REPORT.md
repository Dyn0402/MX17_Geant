# Readout-board bow under gas overpressure: a mechanical model vs the bench dish (2026-10-28)

**Question.** The bench data (nTof_x17 `rerun/model_audit/followups/fb_bow_maps/`) show the mesh centre of every
chamber 1.2–2.3 mm below a point 150 mm out (report 1.5–2.3 ± 0.2–0.4; skeptic review 1.2–1.9 ± 0.4 systematic), with
the cathode roughly flat. Is the gas overpressure (a few mbar) pushing the PCB out enough to explain this? How much
does the Rohacell help, depending on how well it sticks?

**Answer: no, not at a few mbar.** Every scenario stays below 1 mm at ≤ 10 mbar, using the bench's own definition of
the depth. Across PCB thickness 1.3–2.0 mm, every backing/bond case and every edge condition, the expected dish is
**0.2–0.75 mm at 5 mbar** and 0.5–1.1 mm at 10 mbar. To reach the measured band from pressure alone needs
**≥ 15 mbar in the softest case** (bare 1.3 mm PCB, simply supported) and **~40–100 mbar for the nominal board**.
At that pressure the 60 µm window would bulge 16–27 mm, against about 8–10 mm at 3–5 mbar. The measured profile *shape*
also disagrees with a pressure-loaded plate. Pressure contributes a fraction (roughly ¼–½ at 5 mbar) unless the
overpressure is several times "a few mbar".

## Model

`plate.py`: large-deflection (von Kármán) two-layer plate over the 402 × 402 mm support-plate aperture.
- **Layers.** The PCB and the Rohacell (+ optional 25 µm mylar skin on its far face) share the deflection w, and each
  has its own in-plane displacement. An interface shear spring k(x, y) couples them: k = 0 is loose (contact only),
  k → ∞ a rigid bond, and a map of k describes tape strips, edge bands or a debonded centre.
- **Discretisation.** Bending by finite differences; membrane and slip terms on Q4 cells; Newton iteration with
  pressure continuation.
- **Edges.** w = 0, with rotation clamped, simply supported (SS), or an elastic spring. The nominal "foam pinch" spring
  comes from the 5 mm Rohacell squeezed under the 15 mm frame wall: K = E b³/12t ≈ 3900 N·m/m (near-clamped).
  In-plane: immovable (bolted) or free.
- **PCB.** Structural laminate t (FR4 at 20 GPa with the 5 Cu layers smeared by their in-window coverage, plus the
  50 µm coverlay); nominal 1.6 mm. The CAD's 1.70 mm counts the mesh and the 150 µm gas gap.
- **Loads.** Bench weight (cathode up, about 0.35–0.4 mbar equivalent) on top of the gas load. Rohacell 51: E 70 MPa.
  The mesh-tension moment and an imposed edge rotation are run as non-pressure alternatives.
- **Matching the bench.** Each model surface is fitted with the bench's basis `[1, U, W, U²/2, W²/2, UW]` over the
  efficient area (399 × 360 mm), and the depth reported is 0.5·k̄·1.5², exactly as `summary.py` defines it.

**Validation** (`validate.py`, `results/validate.json`):
- Linear coefficients match Timoshenko: SS 0.00406 (+0.05 %), clamped 0.00126 (+0.6 % at N = 61, +1 % at N = 41, the
  grid used).
- The spring edge spans SS → clamped.
- The interface spans D₁ + D₂ (k = 0) → the composite section (k = 10¹² N/m³, within 0.7 %).
- The SS mesh-moment result matches the exact Poisson solution to 0.02 %.
- The membrane limit matches the Vlassak–Nix square-membrane law to 0.6–1.2 %.

## Results (`results/compare.md` has every table; figures in `figures/`)

**1. Pressure, thickness, backing** (foam-pinch edge; `depth_vs_pressure.png`)

| PCB | PCB only | Rohacell loose | bonded (no skin) | bonded + mylar skin |
|---|---|---|---|---|
| 1.3 mm, 5 / 10 mbar | 0.74 / 1.00 | 0.70 / 0.96 | 0.58 / 0.86 | 0.46 / 0.74 |
| 1.6 mm | 0.59 / 0.85 | 0.57 / 0.83 | 0.46 / 0.73 | 0.38 / 0.63 |
| 2.0 mm | 0.40 / 0.65 | 0.40 / 0.65 | 0.33 / 0.56 | 0.28 / 0.49 |

- Weight alone gives 0.03–0.08 mm.
- The response is sub-linear: at w ~ t membrane stretching takes over (depth ∝ p^0.5–0.6 at 5–10 mbar).
  - The existing plausibility number ("0.3–0.7 mm per mbar", `physmodel.py` and `08b_bow.tex`) is linear and
    over-predicts by 2–3× at 5 mbar.
  - It also compares the true dip, while the bench quotes the fitted depth. For a near-clamped plate the dip profile
    is far from a parabola (`cross_sections.png`), so the fitted depth is about half the true dip (5 mbar, 1.6 mm,
    loose: true dip 1.09 mm, centre 1.47 mm, fitted depth 0.57 mm).

**2. The Rohacell, and how well it sticks** (1.6 mm PCB, 5 mbar; `bonding.png`)
- **Loose Rohacell** is almost worthless structurally: D 7.4 → 8.8 N·m, depth −3 %.
- **Fully bonded** (with the back mylar as a second skin) it is a sandwich: D 16.3 N·m, depth −35 %
  (0.57 → 0.38 mm). Without the mylar skin, −20 %.
- **Adhesive strength barely matters; coverage does.** A soft full-face PSA (k = G/t ≈ 10⁹ N/m³, e.g. 0.1 MPa over
  100 µm) already captures 85 % of the rigid-bond stiffening, and 10¹⁰ captures 97 %.
- **Partial bonding**, as a fraction of the rigid-bond stiffening captured:

  | where it is bonded | captured |
  |---|---|
  | 25 mm band at the aperture edge only | 33 % |
  | 60 mm band at the aperture edge only | 78 % |
  | tape strips 50 mm wide every 100 mm | 92 % |
  | tape strips 25 mm wide every 134 mm | 68 % |
  | central r < 100 mm patch only | 45 % |
  | bonded except an unstuck centre r < 100 mm | 98 % |

  The shear that makes the sandwich work is carried near the edges, so a centre that has let go costs nothing.
  Tape that only joins the foam to the board outside the aperture is the "loose" case.
- **Bottom line.** A bad bond makes the bow at most about 1.5× the well-bonded value, never 3–5×.

**3. Edge support** (`edges.png`)
- With loose Rohacell, going from clamped to SS raises the depth by +75 % at 2 mbar and +40 % at 5 mbar (0.55 → 0.76).
  Letting the edge slide in-plane gives +10–35 % (growing with pressure). A bonded sandwich is insensitive to the edge (±5 %).
- The SS case on a 410 mm span is the same as on 402 mm.
- No edge condition reaches 1 mm at 10 mbar.

**4. How much pressure would it take** (`high_pressure.png`): for 1.2 / 1.7 / 2.3 mm,

| scenario | Δp for 1.2 mm | for 1.7 mm | for 2.3 mm |
|---|---|---|---|
| bare 1.3 mm PCB, SS (softest) | 15 mbar | ~40 mbar | > 100 mbar |
| 1.6 mm PCB, Rohacell loose, nominal edge | ~22 mbar | ~55 mbar | > 100 mbar |
| 1.6 mm PCB, Rohacell bonded | ~30 mbar | ~65 mbar | > 100 mbar |

The window sags as p^⅓ (`window_gauge.py`; solver = Vlassak–Nix to 0.1 mm): 5.9 / 8.5 / 10.1 / 12.7 mm at
1 / 3 / 5 / 10 mbar, about 20 mm at 40 mbar and about 27 mm at 100 mbar.

**5. Shape test** (`shape_test.png`): each view's measured w0(u) profile is fitted as w0 = c + α·s_model(u).

| model shape | Σχ² (9 views, 63 dof) |
|---|---|
| pressure on a near-clamped plate | 1560 |
| pressure, clamped | 1626 |
| pressure, in-plane free | 1554 |
| parabolic dish | 1065 |
| pressure, SS | 924 |
| imposed edge rotation | 864 |
| mesh-tension moment, SS edge | 848 |

- A near-clamped plate's slope turns over before the rim. The data's slope *steepens* toward the rim on most views
  (det2 Y, det3 X/Y, det4 X/Y, det6 Y). No uniform-pressure plate does that.
- None of the shapes fits (χ²/dof ≈ 13 at best). The bin errors are statistical only, so part of this is
  reconstruction systematics at the rim, which the skeptic review already flagged for the own-axis terms.
- The w0 amplitude, read through each shape with a flat cathode (r_f = ½), implies a dish of 2.8–3.0 mm (median),
  even further from the pressure numbers.

## What this means

1. **"A few mbar" explains roughly ¼–½ of the measured dish**, depending on the board and the edge: 0.4–0.75 mm of
   1.2–2.3 mm at 5 mbar. The Rohacell bond quality moves this by ±25 %, not by factors.
2. **To be all pressure, the bench overpressure must be ≥ 20–60 mbar**, with a 2 cm window bulge. That is directly
   checkable (below).
3. **The rest needs a non-pressure origin**, and the rim-steepening shape points the same way. Candidates, none
   established:
   - As-manufactured board warp. Ordinary bow/twist tolerance is 0.75 % ≈ 3 mm over 400 mm; it is not a load, would
     survive vertical mounting at n_TOF, and would differ board to board.
   - Mounting distortion at the frame/plate pinch. An imposed edge rotation would need about 50–60 mrad (3°) for
     1.7 mm, which is implausible as a pure rotation.
   - Mesh tension with a soft edge: about 0.2 mm at 10 N/cm on SS, zero when clamped.
   - Or the bench depth is over-read: the skeptic review's ±0.4 mm reconstruction systematic and the ≤ 0.6 mm M3
     common bow bring the low end toward 0.8–1.2 mm, which a thin, loose board at ~10 mbar could reach.

## Cheap checks that would settle it

1. **Measure the window sag** with a straightedge and ruler while flowing: 8–10 mm → 3–5 mbar (pressure is a minor
   part); about 20 mm → about 40 mbar (pressure can be all of it). Better still, a U-tube or bubbler height at the
   outlet.
2. **Change the overpressure and re-measure w0** (the solid route, ±0.03 µm/ns). A pressure dish scales like
   p^0.5–0.6 (doubling Δp → ×1.5); a fabrication warp doesn't move.
3. **Straightedge or dial gauge on an unpressurised board** (chamber open or spare): any dish present at 0 mbar is
   warp.
4. **The PCB stack-up from the fab** (thickness enters as t³) and how the Rohacell is taped (full face, strips or
   perimeter).

## Files

- `plate.py`: solver.
- `validate.py`: closed-form checks.
- `scenarios.py`: 60-scenario scan (`results/scenarios.json`, `results/maps.npz`; cache in `results/cache/`, not
  committed).
- `high_pressure.py`: sweep to 100 mbar.
- `window_gauge.py`: window sag vs overpressure.
- `compare.py`: tables (`results/compare.md`) and figures.
- Run: `~/PycharmProjects/nTof_x17/.venv/bin/python scenarios.py` (≈ 45 min on 3 workers), then `compare.py`.

**Limits.**
- Single-layer Kirchhoff kinematics per layer, so no core shear in the Rohacell. The estimated effect is < 1 % at a
  400 mm span.
- Grid N = 41 (≈ 1 % error).
- No contact separation. The board is pushed onto the foam, which is the case under pressure.
- Isotropic FR4 and an ideal square aperture.
- No initial warp: the model gives the change from the unloaded shape.
