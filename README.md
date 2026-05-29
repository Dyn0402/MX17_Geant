# Micromegas GEANT4 Sensitivity Simulation
## n_TOF X17 group — Dylan Neff

Two simulation modes share the same binary:

- **Vacuum mode** (default): Micromegas detector in vacuum with optional Al
  shielding upstream. Pencil beam from z = −10 cm. Used for detector sensitivity
  and gas-mixture comparison studies.
- **Full-experiment mode** (`-m full`): Complete material stack from He-3 target
  through Micromegas to the liquid-scintillator veto wall. Particles originate at
  the centre of the He-3 gas volume (simulating capture or scatter products).

---

## Vacuum mode geometry (front to back, beam along +z)

| Layer              | Material        | Thickness  | Geant4 volume name    |
|--------------------|-----------------|------------|-----------------------|
| Gas window foil    | Mylar (PET)     | 40 µm      | `GasWindow_Mylar`     |
| Gas window coat    | Aluminium       | 0.1 µm     | `GasWindow_Al`        |
| Drift cathode      | Kapton          | 50 µm      | `DriftCathode_Kapton` |
| Drift cathode Cu   | Copper          | 9 µm       | `DriftCathode_Cu`     |
| **Drift gas**      | **active gas**  | **3 cm**   | **`DriftGas`**        |
| Micromesh          | Stainless steel | 30 µm      | `Micromesh`           |
| **Amp gas**        | **active gas**  | **150 µm** | **`AmpGas`**          |
| Resistive paste    | Carbon/acrylic  | 100 µm     | `ResistivePaste`      |

Mesh parameters: wire Ø 18 µm, hole 45 µm, pitch 63 µm, flattened ~30 µm.

The anode Kapton and Cu layers are part of the PCB stack (not duplicated here).

Optional Al shielding slab (`-a <mm>`) placed 2 cm upstream of the Mylar window.

---

## Full-experiment geometry (front to back, beam along +z)

### He-3 pressurised target

The capsule is a **5 cm-diameter, 15 cm-long cylinder** with its long axis along
**Y (perpendicular to the beam)**. The beam travels along +Z and passes through
the curved surface, crossing **5 cm of He-3 gas** (the diameter).

Along the beam direction the material layers encountered are:

| Layer             | Material                   | Thickness along beam | Geant4 volume name |
|-------------------|----------------------------|----------------------|--------------------|
| CFRP curved wall  | C-fibre/epoxy (1.55 g/cm³) | 0.9 mm (fixed)       | `He3Capsule_CFRP`  |
| Al curved wall    | Aluminium                  | 0.5 mm               | `He3Capsule_Al`    |
| **He-3 gas**      | **³He at 300 bar**         | **5 cm (diameter)**  | **`He3Gas`**       |
| Al curved wall    | Aluminium                  | 0.5 mm               | (same LV, exit)    |
| CFRP curved wall  | C-fibre/epoxy              | 0.9 mm               | (same LV, exit)    |

The capsule endcaps (along Y) are also modelled as 0.5 mm Al + 0.9 mm CFRP.
Total capsule Z-extent along beam: 2 × (25 + 0.5 + 0.9) mm = **52.8 mm**.

The CFRP capsule wall thickness is hardcoded at 0.9 mm (separate from the
configurable `-c <mm>` parameter which applies only to the LS cell walls).

He-3 density: ~37.6 mg/cm³ (ideal gas law at 300 atm, 293 K).
Gun fires along +z from the centre of the He-3 gas volume.

### Air gap 1

| Layer   | Material | Thickness |
|---------|----------|-----------|
| Air gap | Air      | 20 cm     |

### Micromegas detector (same layers as vacuum mode)

| Layer              | Material        | Thickness  | Geant4 volume name    |
|--------------------|-----------------|------------|-----------------------|
| Gas window foil    | Mylar (PET)     | 40 µm      | `GasWindow_Mylar`     |
| Gas window coat    | Aluminium       | 0.1 µm     | `GasWindow_Al`        |
| Drift cathode      | Kapton          | 50 µm      | `DriftCathode_Kapton` |
| Drift cathode Cu   | Copper          | 9 µm       | `DriftCathode_Cu`     |
| **Drift gas**      | **active gas**  | **3 cm**   | **`DriftGas`**        |
| Micromesh          | Stainless steel | 30 µm      | `Micromesh`           |
| **Amp gas**        | **active gas**  | **150 µm** | **`AmpGas`**          |
| Resistive paste    | Carbon/acrylic (Saral 700A) | 100 µm | `ResistivePaste` |

Resistive paste definition: 65% C / 8% H / 27% O by mass, 1.4 g/cm³
(carbon-loaded acrylic binder approximation).

### PCB stack

| Layer        | Material                          | Thickness  | Geant4 volume name |
|--------------|-----------------------------------|------------|--------------------|
| Kapton       | Kapton (polyimide)                | 50 µm      | `PCB_Kapton`       |
| Copper       | Copper                            | 26 µm      | `PCB_Cu_1`         |
| FR4          | Epoxy-glass (1.85 g/cm³)          | 100 µm     | `PCB_FR4_1`        |
| Copper       | Copper                            | 26 µm      | `PCB_Cu_2`         |
| FR4          | Epoxy-glass                       | 100 µm     | `PCB_FR4_2`        |
| Copper       | Copper                            | 26 µm      | `PCB_Cu_3`         |
| FR4          | Epoxy-glass                       | 100 µm     | `PCB_FR4_3`        |
| Copper       | Copper                            | 26 µm      | `PCB_Cu_4`         |
| FR4          | Epoxy-glass                       | 100 µm     | `PCB_FR4_4`        |
| Rohacell 51  | PMI foam C₄H₅NO (0.052 g/cm³)    | 5 mm       | `PCB_Rohacell`     |
| Al foil      | Aluminium                         | 50 µm      | `PCB_AlFoil`       |

FR4 definition: 60% SiO₂ glass + 40% epoxy (C₁₁H₁₂O₃) by mass.
Rohacell 51 definition: C₄H₅NO, density 0.052 g/cm³.

### Air gap 2

| Layer   | Material | Thickness |
|---------|----------|-----------|
| Air gap | Air      | 2 cm      |

### Scintillator wall

| Layer          | Material                        | Thickness | Geant4 volume name      |
|----------------|---------------------------------|-----------|-------------------------|
| Black tape     | PVC (C₂H₃Cl, 1.3 g/cm³)        | 165 µm    | `ScintWall_BlackTape1`  |
| **Scint. bar** | **PVT plastic (G4_PLASTIC_SC_VINYLTOLUENE, 1.032 g/cm³)** | **3 mm** | **`PlasticScint`** |
| Black tape     | PVC                             | 165 µm    | `ScintWall_BlackTape2`  |
| Al foil        | Aluminium                       | 50 µm     | `ScintWall_AlFoil`      |

### Air gap 3

| Layer   | Material | Thickness |
|---------|----------|-----------|
| Air gap | Air      | 2 cm      |

### Liquid scintillator stack (4 × 1.5 cm layers)

| Layer           | Material                              | Thickness | Geant4 volume name |
|-----------------|---------------------------------------|-----------|--------------------|
| CFRP front wall | C-fibre/epoxy (1.55 g/cm³)            | 1.5 mm †  | `LS_CFRP_1`        |
| **LS layer 1**  | **LAB (C₁₈H₃₀, 0.86 g/cm³)**         | **1.5 cm**| **`LiqScint_1`**   |
| CFRP separator  | C-fibre/epoxy                         | 1.5 mm †  | `LS_CFRP_2`        |
| **LS layer 2**  | **LAB**                               | **1.5 cm**| **`LiqScint_2`**   |
| CFRP separator  | C-fibre/epoxy                         | 1.5 mm †  | `LS_CFRP_3`        |
| **LS layer 3**  | **LAB**                               | **1.5 cm**| **`LiqScint_3`**   |
| CFRP separator  | C-fibre/epoxy                         | 1.5 mm †  | `LS_CFRP_4`        |
| **LS layer 4**  | **LAB**                               | **1.5 cm**| **`LiqScint_4`**   |
| CFRP back wall  | C-fibre/epoxy                         | 1.5 mm †  | `LS_CFRP_5`        |

LAB definition: C₁₈H₃₀ (JUNO recipe), ignoring PPO/Bis-MSB fluors (<0.3% by mass).
† Same CFRP thickness as the He-3 capsule walls; both set by `-c <mm>`.

---

## Gas mixtures

| Tag        | Composition              | W-value | Notes                          |
|------------|--------------------------|---------|--------------------------------|
| `ArCF4`    | Ar/CF₄ 90/10             | 34 eV   |                                |
| `ArIso`    | Ar/iC₄H₁₀ 95/5           | 26 eV   | Classic MM gas; mild Penning   |
| `HeEth`    | He/C₂H₆ 96.5/3.5        | 27 eV   | Low-Z, strong Penning          |
| `ArCO2`    | Ar/CO₂ 70/30             | 27 eV   |                                |
| `ArCF4Iso` | Ar/CF₄/iC₄H₁₀ 88/10/2   | 33 eV   |                                |
| `NeIso`    | Ne/iC₄H₁₀ 95/5           | 27 eV   | Penning active                 |
| `NeCF4`    | Ne/CF₄ 90/10             | 30 eV   | Penning active                 |
| `ArCF4CO2` | Ar/CF₄/CO₂ 45/40/15      | 34 eV   |                                |

W-values from ICRU 31 + Penning-transfer corrections (Sauli 1977, Biagi).

---

## Physics

- **EM**: `G4EmStandardPhysics_option4` (Livermore models below 100 keV,
  accurate photoelectric + Auger + delta-ray production)
- **Hadronic**: `FTFP_BERT_HP` (full high-precision neutron data below 20 MeV,
  includes He-3(n,p)T capture cross section)
- **Production cuts**: 10 µm for e±, 100 µm for γ/p
- **Step limits**: 100 µm in Micromegas gas volumes; 1 mm in He-3 gas

Primary ionization N = Edep / W, with Poisson sampling of the remainder.

---

## Quick start

```bash
# Build on lxplus
source scripts/setup_lxplus.sh
bash scripts/build.sh

# Vacuum mode — 1 MeV gammas, ArCF4, 1000 events
build/mm_sim -g ArCF4 -p gamma -e 1.0 -n 1000 -o /tmp/test_vacuum

# Full-experiment mode — 3 MeV electrons from He-3 centre
build/mm_sim -m full -g ArCF4 -p electron -e 3.0 -n 10000 -o /tmp/test_full

# Full mode — tritons (He-3 capture product), 1 MeV
build/mm_sim -m full -g HeEth -p triton -e 1.0 -n 50000 -o /tmp/test_triton

# Vacuum mode with Al shielding
build/mm_sim -g ArCF4 -p gamma -e 1.0 -n 10000 -o /tmp/test_al1mm -a 1
```

---

## Command-line options

| Flag | Default | Description |
|------|---------|-------------|
| `-g <gas>`    | `ArCF4`   | Detector gas mixture (see table above) |
| `-p <particle>` | `gamma` | `gamma`, `neutron`, `electron`, `proton`, `muon`, `muon+`, `pion`, `alpha`, `triton` |
| `-e <MeV>`    | `1.0`     | Primary particle energy |
| `-n <N>`      | `10000`   | Number of events |
| `-o <path>`   | `mm_output` | Output file base (no extension) |
| `-s <seed>`   | time-based | Random seed |
| `-t <N>`      | `1`       | MT threads |
| `-m <mode>`   | `vacuum`  | `vacuum` or `full` |
| `-a <mm>`     | `0`       | Al shielding upstream (vacuum mode only) |
| `-c <mm>`     | `1.5`     | CFRP wall thickness for LS cell walls (full mode only; capsule wall is fixed at 0.9 mm) |
| `-v`          | off       | Verbose event printout |

---

## Output format

Each job writes ROOT TTrees (CSV fallback if ROOT unavailable).

**`EventTree`** (one entry per event):

| Branch | Units | Description |
|--------|-------|-------------|
| `eventID` | — | Event number |
| `edepDrift` | eV | Energy in drift gas |
| `edepAmp` | eV | Energy in amp gas |
| `nPrimDrift`, `nPrimAmp` | — | Primary ion pairs |
| `nClusDrift`, `nClusAmp` | — | Ionising steps |
| `primInDrift`, `primInAmp` | bool | Primary reached volume? |
| *(full mode only)* | | |
| `edepHe3Gas` | eV | Energy in He-3 gas |
| `edepResistPaste` | eV | Energy in resistive paste |
| `edepPCB` | eV | Total energy in PCB stack |
| `edepScintWall` | eV | Energy in plastic scintillator bar |
| `edepLS1`–`edepLS4` | eV | Energy per liquid-scintillator layer |
| `primInHe3Gas` | bool | Primary reached He-3? |
| `primInPCB` | bool | Primary reached PCB? |
| `primInScintWall` | bool | Primary reached plastic scint.? |
| `primInLS1`–`primInLS4` | bool | Primary reached each LS layer? |

**`ClusterTree`** (one entry per ionising step in gas, both modes):
`eventID`, `trackID`, `parentID`, `x/y/z` [mm], `edep` [eV], `nPrimary`,
`ke` [MeV], `volume` (`DriftGas`/`AmpGas`), `particle`

In MT mode each thread writes `<outfile>_t<N>.root`; `collect_results.py`
merges them with `hadd`.

---

## Aluminium shielding (vacuum mode)

Pass `-a <mm>` to place an Al slab 2 cm upstream of the Mylar window.
The world volume expands automatically to accommodate it.

```bash
build/mm_sim -g ArCF4 -p gamma -e 1.0 -n 10000 -o /tmp/al1mm -a 1
build/mm_sim -g ArCF4 -p gamma -e 1.0 -n 10000 -o /tmp/al4mm -a 4
```

HTCondor shielding scan via `submit_condor.py --shielding` (see script header).

---

## HTCondor submission

```bash
# Full sensitivity scan (vacuum mode)
python3 scripts/submit_condor.py \
    --outdir /eos/experiment/ntof/data/x17/mm_sim_results \
    --nevents 50000

# After jobs finish: merge + summarise
python3 scripts/collect_results.py \
    --indir /eos/experiment/ntof/data/x17/mm_sim_results

# Plot
python3 scripts/plot_results.py \
    --summary /eos/experiment/ntof/data/x17/mm_sim_results/summary/summary.csv
```

---

## Running the full-experiment stack on lxplus

### 1. Build (same as vacuum mode)

```bash
ssh lxplus.cern.ch
cd ~/mm_sim          # or wherever you checked out the code
source scripts/setup_lxplus.sh
bash scripts/build.sh
```

### 2. Single test job

The primary particles in full-experiment mode represent products generated
inside the He-3 gas — typically a **triton** (T, 2.73 MeV) and a **proton**
(p, 0.57 MeV) from the ³He(n,p)T capture reaction, or electrons/gammas from
downstream backgrounds.

```bash
OUTDIR=/eos/experiment/ntof/data/x17/mm_sim_results/full

# Triton at 2.73 MeV (He-3 capture product), ArCF4 gas, 10000 events
build/mm_sim -m full -g ArCF4 -p triton -e 2.73 -n 10000 \
    -o ${OUTDIR}/full_ArCF4_triton_2p73MeV

# Proton at 0.57 MeV
build/mm_sim -m full -g ArCF4 -p proton -e 0.57 -n 10000 \
    -o ${OUTDIR}/full_ArCF4_proton_0p57MeV

# Electron energy scan — check transmission to liquid scintillator
build/mm_sim -m full -g ArCF4 -p electron -e 5.0 -n 10000 \
    -o ${OUTDIR}/full_ArCF4_electron_5MeV
```

### 3. HTCondor energy scan

Use `scripts/submit_condor_full.py`, which mirrors the structure of the
vacuum-mode `submit_condor.py`.

**Default scan:**
- **Electrons**: 200 log-spaced points, 0.1–12 MeV (~2.4 % step per point)
- **Triton**: 9 points, 0.191–5 MeV (He-3 capture product: KE ≈ 0.191 MeV for thermal n)
- **Proton**: 7 points, 0.573–5 MeV (He-3 capture product: KE ≈ 0.573 MeV for thermal n)

```bash
# Dry run first — prints job list without submitting
python3 scripts/submit_condor_full.py \
    --outdir /eos/user/d/dneff/mx17_geant_sim_results/full \
    --jobdir /afs/cern.ch/user/d/dneff/condor/mx17_geant_full \
    --nevents 50000 \
    --dry-run

# Submit for real
python3 scripts/submit_condor_full.py \
    --outdir /eos/user/d/dneff/mx17_geant_sim_results/full \
    --jobdir /afs/cern.ch/user/d/dneff/condor/mx17_geant_full \
    --nevents 50000

condor_q
```

To restrict to specific gases or particles:

```bash
python3 scripts/submit_condor_full.py \
    --gases ArCF4 HeEth \
    --particles electron \
    --nevents 50000 \
    --outdir /eos/user/d/dneff/mx17_geant_sim_results/full
```

Key options (see `--help` for full list):

| Flag | Default | Description |
|------|---------|-------------|
| `--outdir` | EOS path | Where ROOT files are written |
| `--jobdir` | AFS path | Condor submit file, wrapper, and logs (must be AFS) |
| `--nevents` | 50000 | Events per job |
| `--gases` | `ArCF4` | Space-separated gas list |
| `--particles` | all | `electron`, `triton`, `proton` |
| `--cfrp` | 1.5 | LS cell CFRP wall thickness [mm] |
| `--flavour` | `workday` | Condor job flavour (~8 h; use `longlunch` for tests) |
| `--dry-run` | off | Print jobs without submitting |

### 4. Run the analysis script

```bash
python3 scripts/analyze_full_experiment.py \
    --indir  /eos/user/d/dneff/mx17_geant_sim_results/full \
    --gas    ArIso \
    --outfile /afs/cern.ch/user/d/dneff/x17/mm_sim_results/full_stack/full_experiment_analysis.pdf
```

This produces a multi-page PDF and a summary CSV with the following plots:

| Page | Content |
|------|---------|
| 1 | Transmission fraction per layer vs energy |
| 2 | Mean energy deposition per layer vs energy |
| 3 | LS calorimeter response, linearity, and energy resolution |
| 4 | Containment analysis (enters LS1, linearity ratio) |
| 5 | Layer-by-layer edep bar chart at 0.5, 1, 3, 5, 10 MeV |
| 6 | Angular distributions in MM drift gas at selected energies |
| 7 | Angular resolution (RMS) vs energy + Highland formula reference |

The angular analysis reads the `ClusterTree` and reconstructs a straight-line
fit to the primary electron track in the drift gas. The RMS polar angle gives
the angular uncertainty introduced by multiple scattering in the upstream
material (He-3 capsule walls, air, MM dead layers).

Skip the angular analysis if the ClusterTree files are large or not merged:

```bash
python3 scripts/analyze_full_experiment.py \
    --indir /eos/.../full --outfile /afs/cern.ch/user/d/dneff/x17/mm_sim_results/full_stack/results.pdf --no-angular
```

Key options:

| Flag | Default | Description |
|------|---------|-------------|
| `--indir` | required | Directory containing ROOT files |
| `--gas` | `ArIso` | Gas tag in filenames |
| `--particle` | `electron` | Particle tag in filenames |
| `--outfile` | `full_experiment_analysis.pdf` | Output PDF |
| `--angular-step` | `10` | Process every Nth energy for RMS plot (lower = slower) |
| `--no-angular` | off | Skip angular analysis entirely |

### 5. Merge and inspect results

Per-thread ROOT files are merged the same way as vacuum mode:

```bash
# Merge threads (if -t > 1 was used)
hadd full_ArCF4_electron_5MeV.root \
    ${OUTDIR}/full_ArCF4_electron_5MeV_t*.root

# Or use the existing collection script (reads EventTree automatically)
python3 scripts/collect_results.py --indir ${OUTDIR}
```

### 5. Reading the full-experiment output in ROOT/Python

The `EventTree` gains extra branches in full mode. Key questions and how to
answer them:

**Did electrons of a given energy reach the liquid scintillator?**

```python
import uproot, numpy as np

f = uproot.open("full_ArCF4_electron_5MeV.root")
t = f["EventTree"]
df = t.arrays(["primInLS1", "primInLS2", "primInLS3", "primInLS4",
               "edepLS1",   "edepLS2",   "edepLS3",   "edepLS4"],
              library="pd")

print("Transmission to LS1:", df["primInLS1"].mean())
print("Transmission to LS4:", df["primInLS4"].mean())
print("Mean edep in LS1 (eV):", df.loc[df.primInLS1, "edepLS1"].mean())
```

**Energy budget across layers:**

```python
cols = ["edepDrift", "edepAmp", "edepPCB", "edepScintWall",
        "edepLS1", "edepLS2", "edepLS3", "edepLS4"]
means = df[cols].mean()
print(means / 1e6)   # convert eV → MeV
```

**Transmission fraction vs energy (after collecting multiple jobs):**

```python
import glob, uproot, pandas as pd

rows = []
for fname in glob.glob("full_ArCF4_electron_*MeV.root"):
    energy = float(fname.split("electron_")[1].replace("MeV.root",""))
    t = uproot.open(fname)["EventTree"]
    arr = t.arrays(["primInLS1","primInLS4"], library="np")
    rows.append({"energy_MeV": energy,
                 "trans_LS1": arr["primInLS1"].mean(),
                 "trans_LS4": arr["primInLS4"].mean()})

df = pd.DataFrame(rows).sort_values("energy_MeV")
print(df)
```

---

## Sr-90/Y-90 calibration simulation (`-m sr90`)

A dedicated simulation mode for evaluating an Sr-90/Y-90 beta source as a
detector calibration tool.  Identical to the full-experiment mode **except
the He-3 pressurised target and its Al/CFRP capsule are absent**.  The source
sits in air at the former He-3 gas-centre position; the air path from source
to the MM entrance is 226.5 mm (the same physical distance as in the full
experiment, but previously occupied by He-3 gas + capsule walls + 200 mm air).

| What changes vs `-m full` | Detail |
|---------------------------|--------|
| He-3 gas removed          | 25 mm at 37.6 mg/cm³ |
| Al capsule wall removed   | 0.5 mm at 2.70 g/cm³ |
| CFRP capsule wall removed | 0.9 mm at 1.55 g/cm³ |
| Upstream material         | 226.5 mm of air only |
| Downstream stack          | Unchanged: MM → PCB → scint. wall → LS |

### 1. Build

Same binary; the new mode is selected at run time:

```bash
source scripts/setup_lxplus.sh && bash scripts/build.sh

# Quick test
build/mm_sim -m sr90 -g ArIso -p electron -e 1.0 -n 1000 -o /tmp/sr90_test
```

### 2. Submit HTCondor scan

Output files use the `sr90_` prefix and go into a separate EOS directory so
they cannot be confused with full-experiment results.

```bash
# Dry run
python3 scripts/submit_condor_sr90.py \
    --outdir /eos/user/d/dneff/mx17_geant_sim_results/sr90 \
    --jobdir /afs/cern.ch/user/d/dneff/condor/mx17_geant_sr90 \
    --nevents 50000 --dry-run

# Submit
python3 scripts/submit_condor_sr90.py \
    --outdir /eos/user/d/dneff/mx17_geant_sim_results/sr90 \
    --jobdir /afs/cern.ch/user/d/dneff/condor/mx17_geant_sr90 \
    --nevents 50000
```

Energy range: **0.1 – 3.0 MeV**, ~160 points (150 log-spaced + exact values
at 0.546 MeV Sr-90 endpoint, 2.28 MeV Y-90 endpoint, and round numbers).

### 3. Analyse

```bash
# Step A — Stack analysis on the sr90 ROOT files.
# Produces sr90_stack_analysis.pdf AND sr90_stack_analysis_electron.csv
# (one CSV per particle; the electron CSV is used in Step B).
python3 scripts/analyze_full_experiment.py \
    --indir   /eos/user/d/dneff/mx17_geant_sim_results/sr90 \
    --prefix  sr90 \
    --gas ArIso --particles electron positron \
    --outfile /afs/cern.ch/user/d/dneff/x17/mm_sim_results/sr90_calibration/sr90_stack_analysis.pdf --workers 16

# Step B — Feasibility convolution: fold the Sr-90/Y-90 beta spectrum
# against the per-energy summary table produced by Step A.
# --summary is the *_electron.csv written next to the PDF above.
python3 sr90_calibration/analyze_sr90_calibration.py \
    --spectrum sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \
    --summary  /afs/cern.ch/user/d/dneff/x17/mm_sim_results/sr90_calibration/sr90_stack_analysis_electron.csv \
    --outfile  /afs/cern.ch/user/d/dneff/x17/mm_sim_results/sr90_calibration/sr90_calibration.pdf
```

---

## Sr-90 without Micromegas (`-m sr90nomm`)

Identical to `-m sr90` but the Micromegas detector and PCB are removed entirely.
The source sees **only air** between itself and the plastic scintillator.  The
total source-to-scintillator distance is preserved (282.5 mm of air) so the
geometry is directly comparable to the `-m sr90` results.

**Purpose:** isolate the Micromegas + PCB material budget contribution by
comparing transmission and energy deposition in `sr90` vs `sr90nomm`.

### 1. Build

Same binary, mode selected at runtime:

```bash
build/mm_sim -m sr90nomm -g ArIso -p electron -e 1.0 -n 1000 -o /tmp/test_nomm
```

### 2. Submit HTCondor scan

Use the same `submit_condor_sr90.py` script with `--mode sr90nomm`.  Output
goes into a separate `_nomm` subdirectory automatically.

```bash
# Dry run
python3 scripts/submit_condor_sr90.py \
    --mode sr90nomm \
    --outdir /path/to/sr90 \
    --nevents 50000 --dry-run

# Submit
python3 scripts/submit_condor_sr90.py \
    --mode sr90nomm \
    --outdir /path/to/sr90 \
    --nevents 50000
```

### 3. Analyse

Use `analyze_full_experiment.py` with `--prefix sr90nomm`:

```bash
python3 scripts/analyze_full_experiment.py \
    --indir   /path/to/sr90_nomm \
    --prefix  sr90nomm \
    --gas ArIso --particles electron positron \
    --outfile sr90nomm_analysis.pdf --workers 16
```

### 4. Compare sr90 vs sr90nomm

The key comparison is transmission and edep with/without the MM stack.
The difference in `trans_primInScintWall` between the two modes directly
gives the fraction of electrons stopped by the Micromegas + PCB material.

---

## Sr-90/Y-90 calibration feasibility

`sr90_calibration/analyse_sr90_calibration.py` evaluates using an Sr-90/Y-90
source in place of the He-3 target to calibrate the detector system.

**Physics:** Sr-90 (Emax = 0.546 MeV) and Y-90 (Emax = 2.28 MeV) in secular
equilibrium emit two betas per decay. The script convolves the measured beta
spectrum with the simulated detector response to give spectrum-weighted
observables.

**Inputs required:**
- `sr90_calibration/Sr90_Y90_Beta_Spectrum.csv` — beta spectrum (shipped)
- `full_experiment_analysis_electron.csv` — simulation summary from
  `analyse_full_experiment.py` (run that first)

**Usage:**

```bash
python3 sr90_calibration/analyze_sr90_calibration.py \
    --spectrum sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \
    --summary  /afs/cern.ch/user/d/dneff/x17/mm_sim_results/full_stack/full_experiment_analysis_electron.csv \
    --outfile  /afs/cern.ch/user/d/dneff/x17/mm_sim_results/full_stack/sr90_calibration.pdf
```

**Output PDF pages:**

| Page | Content |
|------|---------|
| 1 | Text summary: efficiencies, mean signals, rate estimates |
| 2 | Beta spectrum (Sr-90 / Y-90 / total) with trigger acceptance highlighted |
| 3 | Electron fate vs energy: miss / contained / punch-through (pie + curve) |
| 4 | Energy deposition in plastic scint. and LS1 for triggered events |
| 5 | Trigger rate vs source activity (log-log) + per-isotope efficiency bars |
| 6 | MM drift gas signal for triggered events; Sr-90/Y-90 range vs experiment range |

**Key outputs reported:**
- Geometric efficiency (~14% — detector is large relative to source distance)
- Spectral trigger efficiency per isotope (Sr-90 ≈ 0%, Y-90 small %)
- Combined trigger rate at typical source activities (μCi → mCi range)
- Mean energy deposition in plastic scint. and LS1 for triggered events

---

## File structure

```
mm_sim/
├── CMakeLists.txt
├── mm_sim.cc                   Main entry point (-m vacuum | full)
├── include/
│   ├── SimConfig.hh            Config struct + SimMode enum
│   ├── DetectorConstruction.hh
│   ├── PhysicsList.hh
│   ├── ActionInitialization.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── EventData.hh            IonizationCluster + EventData structs
│   ├── EventAction.hh
│   ├── SteppingAction.hh       W-value ionization + per-layer edep
│   ├── RunAction.hh            ROOT TTree / CSV output
│   └── SensitiveDetector.hh
├── src/
│   ├── DetectorConstruction.cc Both mode geometries + all materials
│   ├── PhysicsList.cc
│   ├── ActionInitialization.cc
│   ├── PrimaryGeneratorAction.cc
│   ├── EventAction.cc
│   ├── SteppingAction.cc
│   ├── RunAction.cc
│   └── SensitiveDetector.cc
├── macros/
│   └── run_default.mac
├── scripts/
│   ├── setup_lxplus.sh
│   ├── build.sh
│   ├── submit_condor.py           Vacuum-mode sensitivity scan
│   ├── submit_condor_full.py      Full-experiment stack scan (e±, 0.1–18 MeV, 217 pts)
│   ├── submit_condor_sr90.py      Sr-90 scan: --mode sr90 (with MM) or sr90nomm (no MM)
│   ├── collect_results.py         Merge vacuum-mode results + summary CSV
│   ├── plot_results.py            Vacuum-mode sensitivity plots
│   └── analyze_full_experiment.py Full-experiment transmission + calorimetry + angles
└── sr90_calibration/
    ├── Sr90_Y90_Beta_Spectrum.csv  Measured beta spectrum (1000 points, 0–2.3 MeV)
    └── analyse_sr90_calibration.py Sr-90/Y-90 calibration feasibility analysis
```

---

## Notes and caveats

1. **Primary ionization model**: Geant4 uses continuous energy loss, not
   cluster-by-cluster generation. `N = Edep/W` with Poisson fluctuations is
   the standard approach; expect ~10–20% difference vs Garfield++/Heed absolute
   counts due to Fano factor and discrete cluster statistics.

2. **He-3 isotope**: The He-3 target is defined using `G4Isotope` (A=3, Z=2),
   not natural helium. Geant4 HP physics includes the He-3(n,p)T thermal capture
   cross section via the high-precision neutron data library.

3. **Gun position in full mode**: Particles originate at the geometric centre of
   the He-3 gas volume. This simulates products of a neutron capture or scatter
   event (e.g. triton + proton from ³He(n,p)T). Neutron transport to the target
   is not simulated; fire the relevant secondary particles directly.

4. **Neutron sensitivity**: For thermal/epithermal neutrons, nuclear elastic
   scattering produces short-range recoil nuclei. Gas mixtures with H (ethane,
   isobutane) show additional proton-recoil sensitivity.

5. **Step size**: 100 µm limit in Micromegas gas; 1 mm in He-3 gas. Reduce the
   He-3 limit if you need spatial detail of tracks within the target.

6. **Liquid scintillator optical photons**: LAB is defined for dE/dx transport
   only. To simulate scintillation photon yield and propagation, optical material
   properties (refractive index, absorption/emission spectra) must be added and
   `G4OpticalPhysics` registered in `PhysicsList`.
