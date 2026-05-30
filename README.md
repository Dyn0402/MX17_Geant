# MX17 Micromegas Simulation

Geant4 simulation of the MX17/X17 n-TOF Micromegas detector. Multiple geometry
modes support detector development, material budget studies, and calibration.

## Simulation modes

| Mode | `-m` flag | Description |
|------|-----------|-------------|
| Vacuum | `vacuum` | Micromegas in vacuum, optional Al upstream shielding |
| Full experiment | `full` | He-3 target → air → MM → PCB → scint wall → 2× LS |
| Sr-90 calibration | `sr90` | Same stack but without He-3 target |
| Sr-90 no-MM | `sr90nomm` | Source in air directly to scint wall → LS (no MM/PCB) |
| LS calibration | `lscalib` | Bare source gun → air → **1 LS layer only** |
| Back scint calibration | `backscintcalib` | Bare source gun → air → **1 back scint bar only** |

The two calibration modes (`lscalib`, `backscintcalib`) study each detector in
isolation with a bare source and are designed to be run with the Sr-90/Y-90
beta spectrum sampled at the gun (`--spectrum`).

## Detector geometry (from Full_Geant reference)

- **He-3 target**: r = 1.5 cm, L = 8 cm, 300 bar; Al (0.5 mm) + CFRP (0.9 mm) capsule
- **Micromegas**: 40 µm mylar + drift cathode + 3 cm drift gas + micromesh + 150 µm amp + resistive paste
- **PCB stack**: Kapton + 4×(Cu + FR4) + Rohacell + Al foil
- **Trigger scint wall**: 200 µm BlackMylar tape + 3 mm PVT + 50 µm Al foil
- **LS stack**: 2 mm CFRP wall + 600 µm inner CFRP liner + 40 µm Al liner + **2 cm LAB** × 2 layers
- **Back scint bar**: 30×20×2 cm PVT, wrapped in 20 µm Al foil + 200 µm BlackMylar tape

## Command-line options

```
-m <mode>          vacuum | full | sr90 | sr90nomm | lscalib | backscintcalib
-p <particle>      electron | gamma | positron | proton | muon | ...
-e <energy>        Particle energy [MeV]  (ignored when --spectrum is set)
-n <nevents>       Number of events
-g <gas>           ArIso | ArCF4 | ArCO2 | HeEth | NeCF4 | PureAr | ...
-o <output>        Output file base name
-s <seed>          Random seed
-t <nthreads>      MT threads (default: 1)
-c <mm>            CFRP wall thickness [mm] (default: 2.0)
--spectrum <csv>   Sample electron energies from Sr-90/Y-90 spectrum CSV
                   (two whitespace-separated columns: energy_MeV  probability)
--src-dist <mm>    Source-to-detector air gap [mm] (lscalib/backscintcalib, default: 100)
-v                 Verbose
-h                 Help
```

## Output columns

**All modes** (`_events.csv`): `eventID, edepDrift_eV, edepAmp_eV, ...`

**Full / Sr90 modes** (extra columns):
```
edepHe3Gas, edepMylar, edepCathode, edepMicromesh,
edepPCB, edepPCBKapton, edepPCBCu, edepPCBFR4, edepPCBRohacell, edepPCBAlFoil,
edepScintWall, edepScintTape, edepScintAlFoil,
edepLS1_eV, edepLS2_eV, edepLSCFRP_eV,
primInHe3Gas, primInPCB, primInScintWall, primInLS1, primInLS2, primInLSCFRP5
```

**kLSCalib / kBackScintCalib** (extra columns):
```
primaryKE_MeV,          ← sampled beta energy
edepLS1_eV,             ← energy in LAB layer  (lscalib)
edepLSCFRP_eV,          ← energy in LS structural walls (lscalib)
primInLS1,
edepBackScint_eV,       ← energy in PVT bar  (backscintcalib)
primInBackScint
```

## Building on lxplus

```bash
ssh lxplus.cern.ch
source scripts/setup_lxplus.sh
bash scripts/build.sh
```

The binary is written to `build/mm_sim`.

## Running on lxplus with HTCondor

### Sr-90 calibration (full material stack, energy sweep)

```bash
python scripts/submit_condor_sr90.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --outdir /eos/home-d/dneff/sim/sr90_out \
    --mode sr90nomm
```

### LS calibration (liquid scintillator, spectrum sampling)

```bash
python scripts/submit_condor_lscalib.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --spectrum /eos/home-d/dneff/sim/sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \
    --outdir /eos/home-d/dneff/sim/lscalib_out \
    --mode ls \
    --njobs 20 --nevents 100000
```

### Back scint calibration (plastic scintillator, spectrum sampling)

```bash
python scripts/submit_condor_lscalib.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --spectrum /eos/home-d/dneff/sim/sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \
    --outdir /eos/home-d/dneff/sim/backscint_out \
    --mode backscint \
    --njobs 20 --nevents 100000
```

### Merging job outputs

```bash
# ROOT (if built with USE_ROOT)
hadd lscalib_all.root lscalib_out/ls_calib_job*.root

# CSV
head -1 lscalib_out/ls_calib_job000_events.csv > ls_merged.csv
tail -n +2 -q lscalib_out/ls_calib_job*_events.csv >> ls_merged.csv
```

### Full experiment (X17 physics, energy sweep)

```bash
python scripts/submit_condor_full.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --outdir /eos/home-d/dneff/sim/full_out
```

## Analysis

### LS / back scint calibration

```bash
# After merging job outputs:
python ls_calibration/analyze_ls_calibration.py \
    --indir /path/to/merged --mode ls   --outpdf ls_calib.pdf

python ls_calibration/analyze_ls_calibration.py \
    --indir /path/to/merged --mode backscint --outpdf backscint_calib.pdf
```

Output PDF contains:
- Edep spectrum (all events + zoom)
- Mean edep and efficiency vs sampled beta energy
- 2D: beta energy vs detector edep
- Summary page with spectrum-integrated statistics

### Sr-90 calibration analysis (full material stack)

```bash
python sr90_calibration/analyze_sr90_calibration.py \
    --indir /path/to/sr90nomm/csvs --prefix sr90nomm
```

### Full experiment analysis

```bash
python scripts/analyze_full_experiment.py \
    --indir /path/to/full/csvs --gas ArIso --prefix full
```

## Calibration geometry details

### kLSCalib

```
Gun position (z = −totalZ/2)
↓  [100 mm air, configurable via --src-dist]
LS_CFRP_1       2 mm   CFRP structural front wall
LS_InnerCFRP_1  600 µm inner CFRP liner
LS_Al_1         40 µm  Al liner
LiqScint_1      20 mm  LAB liquid scintillator  ← scored: edepLS1
LS_InnerCFRP_2  600 µm
LS_Al_2         40 µm
LS_CFRP_2       2 mm   CFRP structural back wall
```

### kBackScintCalib

```
Gun position (z = −totalZ/2)
↓  [100 mm air, configurable via --src-dist]
BackScintWrap_Tape1   200 µm  black mylar tape
BackScintWrap_Al1      20 µm  Al foil
BackScint              20 mm  PVT plastic scintillator  ← scored: edepBackScint
BackScintWrap_Al2      20 µm
BackScintWrap_Tape2   200 µm
```
