# COMSOL checks of the induction stage

COMSOL 6.4 runs in batch mode on lxplus from cvmfs (CERN licence, `ConductiveMedia`
is licensed). `run_lxplus.sh <Model.java> [np]` compiles and runs a Java model with
its prefs/tmp/recovery under `/tmp` (not AFS). The model reads its parameters from the
environment. It prints its table as `CSV,` lines, which the script writes to `$OUT`,
because COMSOL's security sandbox refuses file writes from the model.

Files:

| file | what |
|---|---|
| `Induce2D.java` | time-domain induction, 2D cut across strips and pad columns |
| `InduceAxi.java` | the same, 2D-axisymmetric: point avalanche over a uniform sheet that spreads |
| `compare_induce2d.py` | the chain's ion operator built from COMSOL's G, against the direct solve |
| `axi_batch.sh` | InduceAxi runs on lxplus, two at a time |
| `WovenMesh.java` | electrostatic field of the true plain-weave micromesh cell (+ Garfield grid) |
| `woven_map.py` | WovenMesh CSV to a ComponentGrid map, with an independent flag cross-check |

```
scp Induce2D.java run_lxplus.sh lxplus:MX17_comsol/
ssh lxplus 'cd ~/MX17_comsol && CASE=pair OUT=$HOME/MX17_comsol/pair.csv ./run_lxplus.sh Induce2D.java 4'
python3 -m response.comsol.compare_induce2d dep.csv pair.csv      # from the repo root
```

One case takes about 10 min on 4 cores (98k elements, 3 µs).

## Induce2D: direct time-domain induction, 2D cut (2026-10-12)

The model is a cut across the ESL strips and pad columns: gas / ESL strips at 2 MΩ/sq /
kapton / glue / pads / FR4. The avalanche is a prescribed convection current, and the
pad Terminals are held at 0 V. No weighting-field theorem is used. There are two cases:

- `deposit`: electrons only, from the mesh. This gives the surface kernel G.
- `pair`: e/ion pairs on the exp(−z/14.8 µm) profile, with ions drifting to the mesh
  over 340 ns.

The chain's operator (`kernel_lut.apply_ion_transit`: ion coupled through the same G,
f = ⟨1 − z/g⟩ = 0.901) reproduces the direct solve after the DREAM shaper:

| central pad | COMSOL pair | chain from COMSOL G |
|---|---|---|
| 10–90 % rise | 227 ns | 228 ns |
| peak time | 371 ns | 376 ns |
| peak (no free scale) | 1.000 | 1.000 |
| fitted f_eff | 0.878 | 0.901 |

Neighbours agree to about 2 %. **In 2D, resistive-layer screening does not suppress the
in-gas ion signal**, and the operator is right. The data want f_eff ≈ 0.17–0.31.

Limits of this check:

- In a 2D cut each ESL strip is an infinite line that floats. Deposited charge cannot
  spread along the strip or drain, so G is flat after landing: the cumulative charge
  on p0 holds at 0.837 from 10 ns to 3 µs. The mechanism still in question is a point
  charge spreading on the sheet during the 340 ns ion transit, and that needs 3D. The
  cheapest model that has it is 2D-axisymmetric with a uniform ESL.
- The final charge is 0.837 for deposit and 0.862 for pair, which should be equal. The
  3 % gap is most likely the 5 µm deposit sheet being under-resolved by the 8 µm gas
  mesh.
- Electrons are slowed 10× (15 µm/ns), which is invisible behind the shaper.

Raw CSVs: `/eos/experiment/ntof/data/x17/response_sim/comsol/induce2d_2026-10-12/`.

## InduceAxi: the same check with charge spreading (2026-10-12)

The model is 2D-axisymmetric, with a point avalanche (σ = 33 µm) on the axis over a
*uniform* ESL sheet that is grounded at r = 8 mm. The pads are a central disk
(r < 340 µm) plus 680 µm rings on the 780 µm pitch; p0–p2 are Terminals and the rest
are grounded. Here deposited charge really does spread off the pad. Charge integrated
on p0 for the deposit case:

| ρ_s | 10 ns | 50 ns | 100 ns | 340 ns |
|---|---|---|---|---|
| 2 MΩ/sq | 0.77 | 0.44 | 0.25 | 0.07 |
| 0.5 MΩ/sq | 0.63 | 0.14 | 0.06 | 0.01 |

These are the chain's operator built from COMSOL's G (taken from the deposit run at the
real 150 µm/ns electron velocity) against the direct pair solve, after the DREAM shaper,
on the central pad:

| ρ_s | rise pair / chain | peak t pair / chain | peak pair/chain (absolute) | fitted f_eff |
|---|---|---|---|---|
| 2 MΩ/sq | 146 / 151 ns | 250 / 258 ns | 1.11 | 0.89 |
| 0.5 MΩ/sq | 129 / 129 ns | 220 / 216 ns | 1.20 | 0.93 |

- **Shape:** the operator gets it right even while the sheet drains the pad during the
  ion transit.
- **Amplitude:** the direct solve is 11–20 % *above* the operator, and the gap grows as
  the sheet relaxes faster. So an ion moving in the gas does not couple exactly like a
  landed charge convolved with its transit. The chain *under*-predicts the ion signal,
  which is the opposite of what the data need (f_eff ≈ 0.17–0.31).
- **Neighbours:** the operator over-predicts p1/p2 relative to p0 by 2–10 %.
- **Not the cause of the 10× slower electrons:** the default deposit run used them,
  and its G gave 1.11 / 1.17. Rerunning with the real velocity (`VE_UM_NS=150`) gave
  1.11 / 1.20, so the excess is not an artefact of the slower electrons.

**Conclusion:** sheet screening of the in-gas ion cannot produce the missing ion signal,
in 2D or with radial spreading. The resistive layer, if anything, makes the ion signal
slightly larger than the chain assumes.

The remaining geometric caveat is the real ESL, which is strips (550 µm on an 800 µm
pitch) and not a uniform sheet. Charge then spreads along one strip (1D), more slowly
than radially. That sits between Induce2D (no spreading) and InduceAxi (2D spreading),
and both of those agree on the shape.

Raw CSVs: `/eos/experiment/ntof/data/x17/response_sim/comsol/induceaxi_2026-10-12/`.

## WovenMesh: the true plain-weave micromesh (2026-10-12)

The production field map (`meshcell/solve_fieldmap.py`) models the mesh as two flat
layers of straight crossed wires. `WovenMesh.java` builds the real weave instead:

- **Wires:** sinusoids z = ±A sin(πs/p), with A = 9.0 µm (the same 1 µm kiss overlap)
  and the 9.5 µm radius swept along each.
- **Cell:** 67 µm, with half-wires on all four faces. The planes through wire centres
  are exact mirror planes of a plain weave, so Neumann walls make the cell exact.
  Garfield must load it with **mirror** periodicity. The weave's translation period is
  2p, so plain periodicity would be wrong.
- **Boundaries:** the same as production, except the top carries a surface charge that
  sets the drift field to exactly 333 V/cm.
- **Mesh:** 419k tets (0.8 µm on wires, 8 µm far), solved in about 1 min.

The map uses the same grid and the same flag erosion / shell projection as production;
`woven_map.py` re-derives all 1.79M flags in numpy and they agree exactly. 199 shell
nodes at the four wire crossings project into the crossing wire and are flagged dead,
as the production exporter also does.

Solver gates:

| gate | woven | bar |
|---|---|---|
| S1 drift field | 333.00 V/cm | ± 1 % |
| S3 V 0.4 µm off a wire | 0.55 V | < 2 V |
| S6 \|E\| vs a 1.4× coarser solve, p95 | 0.65 % | < 1 % |

The ±30 µm band around the mesh is noisier on S6 (p95 1.6 %).

Garfield gates (`meshcell/gates_check.C`, now with a `mirror` argument), both maps
through the same macro, ROOT 6.40, 2000 electrons:

| | straight wires (production) | true weave |
|---|---|---|
| G2 amp field, z = −90 µm | 31 014 V/cm | 30 743 V/cm (−0.9 %) |
| G7 transparency | 0.960 ± 0.004 | **0.918 ± 0.006** |
| funnel rms | 94.6 µm | 94.3 µm |

The weave costs **4.2 ± 0.7 % of electron transparency**, about 6σ. That is a real but
small correction to the primary-electron count. It does not touch the ion/induction
question. The 0.9 % lower amp field is absorbed by the gain calibration.

To rebuild the map:

```
ssh lxplus 'cd ~/MX17_comsol && MAP=1 OUT=/tmp/$USER/map.csv LC_WIRE=0.8 LC_FAR=8 ./run_lxplus.sh WovenMesh.java 8'
python3 response/comsol/woven_map.py assemble map.csv response/meshcell/meshfield_woven_production.txt
root ... gates_check("meshfield_woven_production.txt", 2000, 333.0, true)
```

The ROOT load line is in the header of `meshcell/gates_check.C`. On the desktop,
Garfield is now built against ROOT 6.40, so source `~/Software/root_6_40/bin/thisroot.sh`
after `setup_garfield.sh`, which still selects 6.30.

The CSVs (`map_prod.csv.gz`, `map_coarse.csv.gz`) and the COMSOL logs are in
`/eos/experiment/ntof/data/x17/response_sim/comsol/woven_2026-10-12/`.

Not done yet: S3 avalanche gain and ion backflow (`meshcell/funnel_ion_endpoints.py`)
in the woven map.

COMSOL API notes:

- Model code cannot read or write files at all; no `security.external.filepermission`
  value tried allowed it. Data export with a coordinate file then silently writes
  nothing. Print to stdout instead.
- `Interp` with a `coord` matrix returns `[point][expr]`.
- Work planes and curves must be `run()` before a `Sweep` can select them.

Solver gotcha: `(s<=gap)*exp(-(gap-s)/Lb)` gives 0·inf = NaN once the exponent passes
709. COMSOL reports this as "Repeated error test failures. May have reached a
singularity" at a fixed time. Keep the exponent bounded (`abs()`).
