# COMSOL checks of the induction stage

COMSOL 6.4 runs in batch mode on lxplus from cvmfs (CERN licence, `ConductiveMedia`
is licensed). `run_lxplus.sh <Model.java> [np]` compiles and runs a Java model with
its prefs/tmp/recovery under `/tmp` (not AFS). The model reads its parameters from the
environment. It prints its table as `CSV,` lines, which the script writes to `$OUT`,
because COMSOL's security sandbox refuses file writes from the model.

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

Solver gotcha: `(s<=gap)*exp(-(gap-s)/Lb)` gives 0·inf = NaN once the exponent passes
709. COMSOL reports this as "Repeated error test failures. May have reached a
singularity" at a fixed time. Keep the exponent bounded (`abs()`).
