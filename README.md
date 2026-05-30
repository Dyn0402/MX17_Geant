# MX17 Micromegas Simulation

Geant4 simulation of the MX17/X17 n-TOF Micromegas detector, with multiple geometry modes
for detector development and calibration studies.

## Modes

| Mode | `-m` flag | Description |
|------|-----------|-------------|
| Vacuum | `vacuum` | Micromegas in vacuum, optional Al upstream shielding |
| Full experiment | `full` | He-3 target → air → MM → PCB → scint wall → 2× LS |
| Sr-90 calibration | `sr90` | Same as full but without He-3 target; Sr-90 source in air |
| Sr-90 no-MM | `sr90nomm` | Source in air directly to scint wall → LS (no MM/PCB) |
| LS calibration | `lscalib` | Source capsule → air → 1× LS layer → air → back scint bar |

## Geometry (from Full_Geant reference)

All detector mode geometry is consistent with the 4-arm MX17 full simulation:

- **He-3 target**: r = 1.5 cm, L = 8 cm, 300 bar; Al (0.5 mm) + CFRP (0.9 mm) capsule
- **Micromegas**: 40 µm mylar window + drift cathode + 3 cm drift + micromesh + 150 µm amp + resistive paste
- **PCB stack**: Kapton + 4×(Cu + FR4) + Rohacell + Al foil
- **Trigger scint wall**: BlackMylar tape (200 µm) + 3 mm PVT + Al foil (50 µm)
- **LS stack**: 2 mm CFRP wall + 600 µm inner CFRP liner + 40 µm Al liner + 2 cm LAB × 2 layers
- **Back scint**: 30×20×2 cm PVT wrapped in 20 µm Al foil + 200 µm BlackMylar tape (kLSCalib)

## Building

```bash
# Local build (requires Geant4 + ROOT in environment)
source scripts/setup_lxplus.sh   # or your local Geant4 setup
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd ..
```

On lxplus the `scripts/build.sh` script sets the environment and builds automatically.

## Running

### Local (single job)

```bash
build/mm_sim -m lscalib -p electron -e 1.0 -n 100000 -g ArCF4 -o output/lscalib_e1p0
```

### Local (parallel sweep)

```bash
# Sweep the full Sr-90 energy grid with 8 workers
python scripts/run_local.py --mode sr90 --outdir /tmp/sr90 --workers 8 --nevents 50000

# Sweep the LS calibration energy grid
python scripts/run_local.py --mode lscalib --outdir /tmp/lscalib --workers 8 --nevents 100000
```

The script defaults to `cpu_count - 1` workers. Use `--max-jobs N` to test with fewer points.

### HTCondor (lxplus)

```bash
# Sr-90 calibration
python scripts/submit_condor_sr90.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --outdir /eos/home-d/dneff/sim/sr90_out \
    --mode sr90

# LS calibration
python scripts/submit_condor_lscalib.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --outdir /eos/home-d/dneff/sim/lscalib_out
```

## Command-line options

```
-m <mode>        vacuum | full | sr90 | sr90nomm | lscalib
-p <particle>    electron | positron | gamma | proton | muon | ...
-e <energy>      Particle energy [MeV]
-n <nevents>     Number of events
-g <gas>         Gas: ArIso | ArCF4 | ArCO2 | HeEth | NeCF4 | PureAr | ...
-o <output>      Output file base name (no extension)
-s <seed>        Random seed
-t <nthreads>    MT threads (default: 1)
-c <mm>          CFRP wall thickness for LS cells [mm] (default: 2.0)
-v               Verbose output
-h               Help
```

## Output format

CSV files (or ROOT with `-DUSE_ROOT`):

**All modes** (`_events.csv`): `eventID, edepDrift_eV, edepAmp_eV, nPrimDrift, nPrimAmp, ...`

**Full / Sr90 modes** (extra columns): `edepHe3Gas, edepMylar, edepCathode, edepMicromesh, edepPCB*, edepScintWall, edepLS1, edepLS2, edepLSCFRP, primIn*`

**kLSCalib mode** (extra columns): `edepLS1_eV, edepLSCFRP_eV, primInLS1, edepBackScint_eV, edepBackScintW_eV, edepSourceCap_eV, primInBackScint`

**Cluster tree** (`_clusters.csv`): per-ionisation-cluster detail in drift/amp gas volumes.

## Analysis

### Full / Sr90 analysis

```bash
python scripts/analyze_full_experiment.py \
    --indir /path/to/csvs --prefix sr90nomm --outpdf sr90_results.pdf
```

### Sr-90 calibration feasibility

```bash
python sr90_calibration/analyze_sr90_calibration.py \
    --indir /path/to/csvs --prefix sr90nomm
```

### LS calibration analysis

```bash
python ls_calibration/analyze_ls_calibration.py \
    --indir /path/to/csvs --prefix lscalib --outpdf lscalib_results.pdf
```

The analysis scripts load all CSV files in the directory, compute per-energy summaries,
apply the Sr-90/Y-90 beta spectrum weighting, and produce a multi-page PDF with:
- Mean energy deposit and detection efficiency vs. electron energy
- Spectrum-weighted detector response (expected experimental spectra)
- 2D scatter: LS vs. back scint edep
- Material budget waterfall

## Simulation modes — geometry details

### kLSCalib — LS + back scint calibration

```
Gun (z = -totalZ/2)
↓
SourceCap_Mylar  (25 µm)   — front window of sealed source capsule
SourceCap_Carbon (100 µm)  — carbon support disc
SourceCap_Al     (100 µm)  — aluminium housing
SourceCap_Tape   (100 µm)  — black mylar outer wrap
↓
Air gap (source_to_ls_mm = 100 mm, configurable)
↓
LS_CFRP_1       (2 mm)     — structural CFRP front wall
LS_InnerCFRP_1  (600 µm)   — inner CFRP liner
LS_Al_1         (40 µm)    — Al liner
LiqScint_1      (20 mm)    — LAB liquid scintillator  [SCORED: edepLS1]
LS_InnerCFRP_2  (600 µm)
LS_Al_2         (40 µm)
LS_CFRP_2       (2 mm)     — structural CFRP back wall
↓
Air gap (gap_ls_to_back_mm = 50 mm, configurable)
↓
BackScintWrap_Tape1 (200 µm)  — black mylar tape
BackScintWrap_Al1   (20 µm)   — Al foil
BackScint           (20 mm)   — PVT plastic scint  [SCORED: edepBackScint]
BackScintWrap_Al2   (20 µm)
BackScintWrap_Tape2 (200 µm)
```

Source capsule thicknesses are configurable via `-c` / `SimConfig` fields.
