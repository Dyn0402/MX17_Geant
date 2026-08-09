# Gas campaign + drift-cage simulation roadmap — recorded 2026-08-08

**What this is:** future-work spec recorded from a discussion with Dylan on the evening of
2026-08-08, immediately after the T7 voltage-label incident (see `RESPONSE_SIM_PLAN.md` §0a and
the guard in `response/avalanche/mx17_aval_calib.py`). Nothing here is launched. Fold into
`RESPONSE_SIM_PLAN.md` on its next touch — it was being edited by another session when this was
written, so this lives standalone to avoid a collision.

**Context you need first:** the T6 field map is pure electrostatics — **gas never enters the
FEM solve**, so one map per (V_mesh, E_drift) serves every gas below. The overnight 2026-08-08
job (other session) is solving the gas-agnostic mesh-voltage ladder (~300–700 V) onto the
desktop's `/media/ucla` disk; the per-gas Garfield avalanche campaigns below run *later*
against that ladder.

---

## 1. Target gases — at least these four

| # | Mixture | Where used | Amp-range `.gas` table (nTof_x17 `garfield_sim/gas_tables/`) |
|---|---------|-----------|----------------------------------------------------------------|
| 1 | Ar/iC₄H₁₀ 95/5 | cosmic bench (June) | `Ar_iC4H10_95_5_Saclay_160m.gas` ✅ exists (dry) |
| 2 | Ar/iC₄H₁₀ 90/10 | nTOF (July beam) | `Ar_iC4H10_90_10_Saclay_160m.gas` ✅ exists |
| 3 | Ar/CO₂/iC₄H₁₀ 95/3/2 | SPS — **⚠ ratios to double-check before running** | ❌ no table — needs a Magboltz job (nTof_x17 `garfield_sim` EOS/condor workflow) |
| 4 | Ar/CF₄/iC₄H₁₀ 88/10/2 | SPS | `Ar_CF4_iC4H10_88_10_2_Saclay_160m.gas` ✅ exists |

Practicalities, all four being Ar-dominant:

- **Ion mobility:** `mx17_aval_calib.py` hardcodes `IonMobility_Ar+_Ar.txt` — correct for all
  four. (It is wrong for Ne/He/CF4-dominant gases; irrelevant here but the reason the gas list
  can't be extended blindly.)
- **Penning:** `--penning auto` works for the common Ar binaries; for the two *ternary*
  mixtures Garfield's built-in table may not exist, in which case the script deliberately exits
  rather than run at rP = 0. If that happens, rerun `--penning manual --penning-rp <value>` and
  record the choice (the uniform-field campaign used auto on 95/5 only, so this is untested).
- **Voltage ranges:** per-gas, centred on each gas's actual operating point from the bench/beam
  records — this is what sized the overnight map ladder's min→max span.

## 2. Contaminants (H₂O etc.) — open question, with a proposed cheap answer

The drift-velocity work hypothesises **site-dependent water contamination** (cosmic bench ~1 %
and drying over a week — `micro-tpc` June work; nTOF ~0.2 % trace — July 90/10 study). Taken
literally that multiplies every gas above by a contamination axis, which is a lot of Magboltz
jobs and campaign slices.

**Unknown:** how much water matters for *gain* (amplification, 30–45 kV/cm) as opposed to
*drift velocity* (200–330 V/cm), where its effect is established and large. These are very
different reduced-field regimes; the drift-side sensitivity does not automatically carry over.

**Proposed strategy (not yet approved as a decision, recorded as the leading idea):**

1. Do a **sensitivity bracket first**, not the full grid: dry vs +0.5 % vs +1 % H₂O for ONE gas
   (Ar/iso 95/5) at ONE voltage (490 V), ~8 slices each. Needs wet *amp-range* Magboltz tables
   — the existing wet suites (`wet_ariso_bracket`, `wet_cf4_drift`) are drift-range only, so
   this is ~2 new Magboltz jobs before any Garfield.
2. Compare the gain shift against the calibration floor (geometry dominates the gain error
   budget; Penning is ±1–8 V equivalent — see nTof_x17 `garfield_sim` error-budget notes). If
   water moves gain by less than a few %, **drop the contaminant axis from the gain campaign
   entirely** and keep contaminants only where they are known to matter (drift transport).
3. Only if the bracket shows a real effect does the per-site contamination grid get built.

Note the field-map ladder is unaffected either way — contaminants, like gases, never enter the
FEM solve.

## 3. Drift-cage / degrador simulation — geometry recorded, study not started

A **separate, macroscopic** FEM problem, not the 67 µm mesh unit cell: the full 30 mm drift
volume with the mesh plane as a flat equipotential boundary. Same gmsh + scikit-fem toolchain
as `solve_fieldmap.py`. The two solves couple only through the electron-transparency ratio
eps(E_amp/E_drift), already measured (S2, `response/meshcell/transparency_curve.csv`).

**Geometry (Dylan, 2026-08-08 — verbal, verify against drawings before meshing):**

- 3 copper rings on a degrador PCB at the drift-gap perimeter, **evenly spaced along the 30 mm
  gap**.
- Ring widths such that the exposed copper and the PCB gaps between rings are of **roughly
  equal extent**, with **PCB (not copper) at the top and bottom** of the stack.
- Electrical: **top ring tied to the full drift voltage**, then ~1 GΩ top→middle, ~1 GΩ
  middle→bottom, ~1 GΩ bottom→ground. Equal steps ⇒ nominal ring potentials V, ⅔V, ⅓V.
- The ~1 GΩ is nominal. **Verification task (later):** fit the effective chain resistance from
  drift-HV current vs voltage readings at all three sites — nTOF, cosmic bench, SPS. One data
  point already exists: run_79 non-B drift channels drew ~0.18 µA at ~700 V ⇒ R ≈ 3.9 GΩ total,
  ~1.3 GΩ/ring (nTof_x17 memory `drift-cage-degrador-rings`).

**Planned study (when picked up):**

- Solve the cage **with and without the ring divider** (floating/removed rings ≠ graded
  solution — a floating conductor distorts the boundary field). Detector B ran the whole nTOF
  campaign with its degrador physically removed and shows a real, unexplained slowdown
  (nTof_x17 memory `detector-b-anomaly`) — this solve is the direct test of that mechanism.
- Sweep drift field. Gain is insensitive to E_drift (mesh penetration ~2 V at 490 V), so
  drift-field variations live in this solve + transparency interpolation, not in new avalanche
  campaigns; only add E_drift as a second *mesh-map* axis if transparency fidelity beyond the
  S2 curve is needed.
- Deliverables: field-uniformity map of the drift volume, edge-distortion extent with/without
  rings, and a predicted drift-time/velocity signature to compare against detector B.

**Do not start the cage solve before confirming the ring geometry against drawings/photos and
getting the chamber-wall boundary condition from Dylan.**
