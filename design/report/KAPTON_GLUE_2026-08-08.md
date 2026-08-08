# The ESL⇄pad insulator: kapton is 50 µm, and it is not alone

2026-08-08. Closes plan §1 row 4, which had carried "thickness unknown" since
2026-08-06 and a `{50, 75, 125} µm` scan.

## Verdict

The kapton is **50 µm, confirmed** — `pcbKapton_um` in
`shared/MX17ModuleGeometry.hh`, matching the `PCB_Kapton` row of the
MX17_Full_Geant stack. It is no longer a scan axis, and it is now pinned in
`assert_against_header()` so a geometry edit breaks the response chain at
import instead of quietly changing the answer.

But **50 µm is not the dielectric the pads see.** The ESL is screen-printed on
the top face of a kapton coverlay whose bottom face is laminated to the pad
copper with an adhesive, so the insulator is two dielectrics in series. Adding
it moves the prompt sum rule by 4.5 %:

| stack | series thickness | S(0)/C(0) |
|---|---|---|
| bare 50 µm kapton | 50.0 µm | 0.913043 |
| **50 µm kapton + glue (production)** | **70.5 µm** | **0.881583** |
| old 75 µm single-layer default | 75.0 µm | 0.875000 |

So the headline is not "the kapton got thinner". It is that **taking the
confirmed 50 µm at face value and dropping the glue would have been a larger
error than the 75 µm guess it replaced**, and in the opposite direction. The
old default was accidentally close to right: `c′ = 0.4985 µF/m²` and
`D = 2.006 m²/s` at 1 MΩ/sq still reproduce the figures plan §3 quotes — but
now for the right reason.

## How the glue was estimated

It is an estimate, and it is now the dominant uncertainty in the stack. In
descending order of how much each step deserves trust:

1. **25 µm as supplied.** A 50 µm (2 mil) polyimide coverlay is supplied with
   an acrylic adhesive layer, and the standard pairing for 2 mil film is 1 mil
   = 25 µm. Common alternatives are 12.5 and 35 µm — the bracket carried in
   `GLUE_THICK_SUPPLIED_UM_BRACKET`.
2. **ε_r ≈ 3.2.** Acrylic/epoxy coverlay adhesive runs 3.0–3.5, and the stack
   barely cares: across that whole range the series thickness moves ~5 %.
3. **18.8 µm over the pad, not 25.** What is supplied is not what ends up over
   the pad. Vacuum lamination drives adhesive into the 100 µm trenches between
   the 680 µm pads, which are 26 µm deep, so the over-pad layer is thinner by
   the trench volume: `t_pad = t_supplied − (1 − coverage)·t_cu`, with coverage
   = (680/780)² = 0.76003 exactly. That is the electrically relevant thickness,
   because the pad face is the electrode.

Bracket on the series thickness: **56.9 µm** (12.5 supplied) to **81.5 µm**
(35 supplied), against 70.5 µm production.

This is *not* a production scan axis. The bracket exists for sensitivity checks.

## Why the solver got a real two-layer stack instead of an effective thickness

There is a single-layer thickness with the same series capacitance,
`d_eff = d_k + d_g(ε_k/ε_g)`. It is exact at k → 0 and wrong where the kernels
live. Substituted into the one-layer form it understates S(k) by:

| k | 0.1/pitch | 0.5/pitch | 1/pitch | 2/pitch | 5/pitch |
|---|---|---|---|---|---|
| ΔS/S | −3.0e-05 | −7.6e-04 | −3.0e-03 | −1.1e-02 | −5.3e-02 |

so `stack_coeffs` is now the exact cascaded two-layer transmission-line result,
written with `tanh`/`sech` so it stays finite at every `k·d`. `d_eff` survives
only in `insulator_d_eff_m()` for reporting and `c′`.

**Layer order is fixed and it matters.** S(k) is reciprocal and does not care
which layer is which; C(k) is the admittance looking down from the ESL and does
— swapping them moves it 2.9 % at k = 2/pad-pitch. Kapton against the ESL,
glue against the pads.

Verification (`--check` equivalents all run):

- `d_glue = 0` reproduces the previous single-layer closed form to **4e-16** in
  both C(k) and S(k) over k = 1 … 1e6.
- `C(0) == c_prime()` and `S(0) ==` the series capacitance, to 1e-12.
- `prompt_sum_rule` agrees with the solver's own `S(0)/C(0)` to 1e-12, and
  `d_g=0` reproduces the legacy 0.913043 exactly.
- finite and `S ≥ 0` from k = 0 to 1e12.
- `validate.py`: **solver correctness 5/5**. (V4 still refutes A1 — a
  pre-existing result, untouched by this.)

## Consequence: every existing S1 product is legacy

All 12 products were solved on the bare-kapton stack. They are not wrong, they
are a different model, and the two are 4.5 % apart in the sum rule. Three
guards now keep them from being confused:

- the glue is **in the product filename** (`..._dk50um_g19um.npz`), so a new
  solve cannot silently overwrite a legacy one — which it would have, since
  four legacy products sit at exactly `dk50um`;
- `d_glue_m` / `glue_eps_r` / `d_insulator_eff_m` are in the product metadata,
  and `manifest.py` reports a `stack` column and **refuses to count legacy
  products toward grid completeness**;
- `KernelLUT.describe()` reports which stack a run's kernel came from.

**Re-solve in flight**: condor cluster 13352631, 4 points (ρ_s ∈ {0.5, 1, 2, 5}
MΩ/sq at the one confirmed thickness) → EOS `s1_ny1024/`. The thickness axis is
retired, so this is 4 points where the old grid was 12.

⚠️ **The ρ_s × d_k scan of 2026-08-07 is superseded for the production stack.**
It varied d_k inside a *single-layer* model; the two-layer stack differs at
finite k by up to 5 %, which is the regime c1 is built from. Its qualitative
finding (c1_X flat in ρ_s, the power is in the c1_Y/c1_X ratio) should be
re-checked against the new kernels before it is leaned on.

## What generating the README stack table caught

The README described the module as prose; MX17_Full_Geant's carries a layer
table. There is one here now, **generated** by `scripts/stack_table.py` from the
geometry header — declaration defaults with the `AsBuiltSpec()` overrides on
top, because both matter — with a `--check` mode that fails if the README has
drifted. Generating rather than typing it immediately paid for itself twice.

**`scripts/model/mx17_model.py` was stale.** It hardcoded the LEGACY 100 µm
paste declaration while `AsBuiltSpec()` has been 10 µm since 2026-08-07, so
every model figure carried 90 µm of paste that is not there — and through the
FR4 residual it took 18 µm per layer out of the laminate (246.4 vs the correct
264.4). Fixed; the board body still closes at exactly 1.7000 mm.

**A discrepancy that is NOT fixed, deliberately.** Read literally the header
stacks coverlay kapton → L3 guard ring → 264 µm FR4 → L4 pads, i.e. a quarter
of a millimetre of laminate between the ESL and the pads. The response model
instead assumes the physical adjacency for a resistive Micromegas: ESL on
kapton, kapton laminated straight onto the pads. L3 has **zero** copper inside
the active window and the FR4 is an equal-division residual, so that separation
is an artifact of an admittedly arbitrary split rather than measured geometry —
but the two models disagree about the single distance that sets the induced
signal. It is irrelevant to material budget and decisive for the response. The
missing fab stackup drawing settles it, and guessing it into the geometry would
be worse than recording it. See `design/NEEDED_INPUTS.md` §6.

## Open question closed: the ion-template provenance gap

`DESKTOP_RUNS_2026-08-07.md` recorded that "v2's templates came from a
different S3 pass and are not reproducible from what is archived". That is
true in letter and **overstated in substance**. The mechanism was already
found and fixed in this repo:

Garfield's `AvalancheMicroscopic` and `AvalancheMC` both default to
`m_useWeightingPotential = true`, and `ComponentConstant` has no weighting
potential until `SetWeightingPotential()` is called. The first campaign set
only the weighting *field*, so ψ was 0 everywhere and every induced signal came
out identically zero — silently, with the avalanche and ion drift running
normally and gain/σ₀ unaffected. Commit `42390d1` fixed it.

The timeline settles the rest:

| | when |
|---|---|
| archived raw slices (56) | Aug 7, 11:38–11:58 |
| `aval_calib.json` (all-zero templates) | Aug 7, 11:37 |
| **fix `42390d1`** | **Aug 7, 13:20** |
| `aval_calib_v2.json` (real currents) | Aug 7, 15:35 |

So the archived raw is the **pre-fix** campaign, and v2 came from a post-fix
re-run whose raw was never uploaded. The gap is **archival, not physical**, and
v2's templates are independently validated inside v2 itself: f_ion from the
current integrals (0.9079) and from the α_z histogram (0.9077) agree to four
digits, which is also the proof that no ion charge was lost off the end of the
400 ns window.

v2 is a genuinely different realization, not the same slices re-reduced — its
470 V gain is 22854.3 against the archived raw's 23107.4 (~1 %).

**Validation in flight**: condor cluster 13352632 re-runs two 490 V slices
(seeds 49001/49002, nev 140, byte-identical producer) with the fixed code. The
target is v2's `f_ion = 0.907907` at 490 V. f_ion is a ratio and stable at 280
events; the gain will differ at the ~1 % level by statistics and that is
expected, not a failure.

## Numbers to reuse

```
KAPTON_THICK_UM          = 50.0      confirmed, header-pinned
KAPTON_EPS_R             = 3.5
GLUE_THICK_SUPPLIED_UM   = 25.0      estimate, bracket (12.5, 25, 35)
GLUE_THICK_UM            = 18.76     over the pad, after squeeze-out
GLUE_EPS_R               = 3.2
insulator_d_eff_m()      = 70.52 µm  series equivalent, k->0 only
c_prime()                = 0.4985 µF/m^2
sheet_diffusivity(1e6)   = 2.006 m^2/s
prompt_sum_rule(50e-6)   = 0.881583  (legacy bare kapton: 0.913043)
```
