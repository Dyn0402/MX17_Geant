# Micromegas GEANT4 Sensitivity Simulation
## n_TOF X17 group — Dylan Neff

Simulates primary ionization and energy deposition in a 40×40 cm Micromegas
detector for the X17 sensitivity study. Fires gammas, electrons, and neutrons
across a wide energy range through five candidate gas mixtures.

---

## Detector geometry (front to back, beam along +z)

| Layer              | Material        | Thickness  |
|--------------------|-----------------|------------|
| Gas window foil    | Mylar (PET)     | 40 µm      |
| Gas window coat    | Aluminium       | 0.1 µm     |
| Drift cathode      | Kapton          | 50 µm      |
| Drift cathode Cu   | Copper          | 9 µm       |
| **Drift gas**      | **active gas**  | **3 cm**   |
| Micromesh          | Stainless steel | 30 µm      |
| **Amp gas**        | **active gas**  | **150 µm** |
| Anode              | Kapton          | 50 µm      |
| Anode Cu           | Copper          | 9 µm       |

Mesh parameters (woven SS bulk Micromegas): wire Ø 18 µm, hole 45 µm,
pitch 63 µm, total flattened thickness ~30 µm.

---

## Gas mixtures

| Tag        | Composition              | W-value |
|------------|--------------------------|---------|
| `ArCF4`    | Ar/CF₄ 90/10             | 34 eV   |
| `HeEth`    | He/C₂H₆ 96.5/3.5        | 27 eV   |
| `ArCO2`    | Ar/CO₂ 70/30             | 27 eV   |
| `ArCF4Iso` | Ar/CF₄/iC₄H₁₀ 88/10/2   | 33 eV   |
| `NeIso`    | Ne/iC₄H₁₀ 95/5           | 27 eV   |

W-values from ICRU 31 + Penning-transfer corrections (Sauli 1977, Biagi).

---

## Physics

- **EM**: `G4EmStandardPhysics_option4` (Livermore models below 100 keV,
  accurate photoelectric + Auger + delta-ray production)
- **Hadronic**: `FTFP_BERT_HP` (full high-precision neutron data below 20 MeV)
- **Production cuts**: 10 µm for e±, 100 µm for γ/p
- **Step limit**: 100 µm max step in gas volumes

Primary ionization N = Edep / W, with Poisson sampling of the remainder.

---

## Quick start on lxplus

```bash
# 1. Clone / copy the code to your lxplus home or EOS
cd ~/mm_sim   # or wherever you put the code

# 2. Source Geant4 + ROOT from CVMFS
source scripts/setup_lxplus.sh

# 3. Build
bash scripts/build.sh

# 4. Quick test (1000 events)
build/mm_sim -g ArCF4 -p gamma -e 1.0 -n 1000 -o /tmp/test_mm -v

# 5. Submit full scan to HTCondor
python3 scripts/submit_condor.py \
    --outdir /eos/experiment/ntof/data/x17/mm_sim_results \
    --nevents 50000

# 6. Monitor
condor_q

# 7. After completion: merge + summarise
python3 scripts/collect_results.py \
    --indir /eos/experiment/ntof/data/x17/mm_sim_results \
    --clusters

# 8. Plot
python3 scripts/plot_results.py \
    --summary /eos/experiment/ntof/data/x17/mm_sim_results/summary/summary.csv \
    --rootdir /eos/experiment/ntof/data/x17/mm_sim_results/summary/merged
```

---

## Output format

Each job writes (ROOT TTree):

**`EventTree`** (one entry per event):
- `eventID`, `edepDrift` [eV], `edepAmp` [eV]
- `nPrimDrift`, `nPrimAmp` — primary ion pairs
- `nClusDrift`, `nClusAmp` — number of ionising steps
- `primInDrift`, `primInAmp` — did the primary traverse the volume?

**`ClusterTree`** (one entry per ionising step):
- `eventID`, `trackID`, `parentID`
- `x`, `y`, `z` [mm] — step midpoint position
- `edep` [eV], `nPrimary` — energy and ion pairs this step
- `ke` [MeV] — particle kinetic energy
- `volume` — `DriftGas` or `AmpGas`
- `particle` — particle name

In MT mode, each thread writes `<outfile>_t<N>.root`; `collect_results.py`
runs `hadd` to merge them.

---

## Submitting a subset (testing)

```bash
# One gas, one particle, two energies
python3 scripts/submit_condor.py \
    --gases ArCF4 \
    --particles gamma \
    --nevents 5000 \
    --flavour longlunch \
    --outdir /tmp/mm_test \
    --dry-run   # remove --dry-run to actually submit
```

---

## File structure

```
mm_sim/
├── CMakeLists.txt
├── mm_sim.cc               Main entry point
├── include/
│   ├── SimConfig.hh        Shared config struct
│   ├── DetectorConstruction.hh
│   ├── PhysicsList.hh
│   ├── ActionInitialization.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── EventData.hh        IonizationCluster + EventData structs
│   ├── EventAction.hh
│   ├── SteppingAction.hh   W-value ionization scoring
│   ├── RunAction.hh        ROOT TTree / CSV output
│   └── SensitiveDetector.hh
├── src/
│   ├── DetectorConstruction.cc  All 9 layers, 5 gas mixtures
│   ├── PhysicsList.cc
│   ├── ActionInitialization.cc
│   ├── PrimaryGeneratorAction.cc
│   ├── EventAction.cc
│   ├── SteppingAction.cc
│   ├── RunAction.cc
│   └── SensitiveDetector.cc
├── macros/
│   └── run_default.mac
└── scripts/
    ├── setup_lxplus.sh     Sources G4 11.2 + ROOT from CVMFS
    ├── build.sh            CMake build
    ├── submit_condor.py    Submits full gas×particle×energy matrix
    ├── collect_results.py  hadd + uproot summary CSV
    └── plot_results.py     matplotlib sensitivity plots
```

---

## Notes and caveats

1. **Primary ionization model**: Geant4 uses continuous energy loss (Bethe-Bloch),
   not cluster-by-cluster generation. The `N = Edep/W` conversion gives the
   *expectation value* of primary ion pairs per step, with Poisson fluctuations
   on the remainder. For cross-gas comparisons this is consistent. For absolute
   counts vs Garfield++ Heed, expect ~10-20% differences due to Fano factor
   and discrete cluster statistics not modelled here.

2. **Neutron sensitivity**: For thermal/epithermal neutrons, the dominant process
   is nuclear elastic scattering producing short-range recoil nuclei. Geant4 HP
   handles this correctly. Gas compositions with H (ethane, isobutane) will show
   additional proton-recoil sensitivity.

3. **Gamma sensitivity**: Dominated by photoelectric effect at low E and Compton
   at high E. Ar K-edge at 3.2 keV will be visible as a jump in ArCF4/ArCO2/ArCF4Iso.
   He/Eth has very low Z so photoelectric cross-section is tiny; CF4 adds F (Z=9).

4. **CF4 note**: CF4 is electronegative (attaches electrons at high concentrations)
   but at 10% it's primarily a quencher. The W-value is elevated vs pure Ar.

5. **Step size**: 100 µm limit in gas. For very low-energy electrons (< 10 keV)
   where ranges are < 100 µm, Geant4 will automatically use shorter steps.
   Reduce to 10 µm if you want finer z-resolution in the cluster tree.
