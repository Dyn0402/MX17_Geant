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

## Detector geometry

The MM module (window → support plate) has **one description**,
[`shared/MX17ModuleGeometry.hh`](shared/MX17ModuleGeometry.hh), consumed by
this sim and by `../MX17_Full_Geant`. The default is the **as-built** module
derived from the mechanical CAD + hardware knowledge
([`design/GEOMETRY_FROM_CAD.md`](design/GEOMETRY_FROM_CAD.md),
[`design/GEOMETRY_IMPLEMENTATION_NOTES.md`](design/GEOMETRY_IMPLEMENTATION_NOTES.md);
remaining assumptions in [`design/NEEDED_INPUTS.md`](design/NEEDED_INPUTS.md)):

- **Window**: 60 µm aluminized mylar, 440×440, bulged by the gas overpressure
  (terraced dome, default 8 mm sag ≈ Hencky at ~3 mbar, `--bulge-front`)
- **Window flange**: 5 mm Al ring (440/400.2) = 5 mm gas gap before the cathode
- **Drift**: Cu-clad kapton cathode (406²), 30.000 mm gap set by the Al gas
  frame (440/410), 4× 1.62 mm field-cage PCBs lining the aperture
- **Mesh**: woven SS 19/48 µm → 38 µm effective-density slab (fill 0.223);
  bulk pillars (Ø0.6 mm, 85×85 grid at 4.68 mm) stand in the amplification gap
- **Resistive layer**: 515 discrete ESL strips, 550 µm wide on a 0.8 mm pitch,
  **10 µm** thick over 412²; the 250 µm inter-strip grooves are chamber gas and
  are scored as `AmpGas` (they add 2.0 % of the amp-gas energy deposit)
- **Readout board**: 1.70 mm single body (mesh + 150 µm amp + 10 µm ESL +
  laminate), 470², offset (+15,+15) from the active axis; copper = the five
  physical gerber layers (guard ring, pads, Y strips, X strips, fan-out),
  density-scaled by their measured coverage (`scripts/gerber/analyze_cu_coverage.py`).
  Each layer is built as **two zones** — inside and outside the 399.36 mm
  active window — because the copper is very unevenly distributed: the pad
  layer is 0.76 covered inside the active area and 0.057 outside, so a single
  board-wide average (0.564) would put ~26 % too little copper where the beam
  passes and ~10× too much at the edges. Zoning costs no CPU.
- **Backing**: 5 mm rohacell + 25 µm aluminized mylar + **8 mm Al support
  plate with a 402 mm aperture** concentric with the 399.36 mm active area
- **M1 front-end cards**: 2× per connector edge, flat on the drift side of
  the board edge and straddling it (6.6 mm envelope: 6 gerber Cu layers + FR4)

### Real readout pattern (default; `--homogenized-readout` to disable)

The three signal copper layers are built as the **real 512 × 512 pattern** read
off the production gerbers: 0.68 mm square pads on L4, and Ø0.50 mm dots +
0.098 mm bus traces on the L5/L6 X-Y track layers, all on the gerber-exact
0.78 mm pitch. The **ESL resistive layer** is likewise built as 515 discrete
550 µm strips on the 0.8 mm pitch, inside a gas-filled envelope — so the 250 µm
inter-strip grooves are real chamber gas, not reduced-density paste. Total ESL
mass is unchanged (515 × 0.550 / 412 = 0.6875, exactly the old scale factor).

The groove cells are named `AmpGas`, so groove ionisation is scored as
amplification signal: they are only 10 µm deep at the bottom of a 150 µm gap,
so the amplification field runs through them essentially undisturbed. They add
2.0 % of the amp-gas energy deposit.

> ⚠️ **If you consume `ClusterTree`**: the groove cells are thin, so Geant4
> forces short steps and each writes its own row. The grooves are 2.0 % of
> amp-gas *energy* but 21.7 % of amp-gas *rows*. Sum `edep` or `nPrimary` —
> never count rows as a proxy for ionisation. (`nPrimary` is safe: the sub-W
> remainder is carried probabilistically, so nothing is lost to segmentation.)

L3 (guard ring) and L7 (fan-out) have irregular artwork and stay as the
two-zone density-scaled sheets.

Both grids are perfectly regular, so they are built with `G4PVReplica` (two
nested levels for the copper, one for the strips). The whole 786432-feature
copper pattern plus 515 ESL strips costs **17 extra volumes**, not 786943
placements, and cell lookup stays O(1) index arithmetic:

| | volumes | total CPU, 100k sr90 events | max RSS |
|---|---|---|---|
| `--homogenized-readout` | 7342 | 27.0 s | 582 MB |
| patterned (default) | 7359 | 28.0 s | 579 MB |

So the full pattern costs about **+3–4 % of total CPU** and no extra memory.
Removing the fixed ~13 s of startup overlap-checking, the tracking-only cost is
roughly +9 %, most of it from the resist strips rather than the copper (the ESL
is 100 µm thick against 26 µm per Cu layer, and its envelope adds gas/solid
boundaries). Treat these as ±3 % — they were measured as user CPU time,
min-of-4, on a machine with other jobs running.

It is cheap because the patterned layers total ~180 µm of a ~45 mm stack, so
few steps land there — not because the pattern is approximated.

**When you can turn it off**: the signal copper sits *downstream of both gas
gaps*, so for anything scored in the gas the homogenized zones give
essentially the same answer. The resist strips are the part that could matter
even for gas scoring, since the ESL bounds the amplification gap. The pattern
would also get relatively more expensive for tracks at shallow incidence,
which cross many cells per layer.

`--legacy-geometry` is never patterned regardless of these flags, so it stays
bit-identical to the pre-2026-08 output (checked with `scripts/tree_digest.py`).

`--legacy-geometry` restores the pre-2026-08 uniform-slab stack (bit-identical
to the old build). Model figures live in `design/figures/` and regenerate with
`python scripts/model/plot_mx17_model.py`; the cross-section/3D/plan figures
render the **true constructed geometry** from `design/mx17_geometry.json`
(refresh it after geometry changes: `mm_sim -m sr90 --dump-geometry
design/mx17_geometry.json`), while the board/peel figures render the
production gerbers directly.

### The layer stack

Generated from the header — regenerate with `python3 scripts/stack_table.py`,
and `--check` fails if it has drifted from the geometry.

<!-- BEGIN STACK TABLE (generated by scripts/stack_table.py) -->

| Layer | Material | Thickness | G4 volume | Note |
|---|---|---|---|---|
| Gas window | Mylar (PET) | 60 µm | `GasWindow_Mylar` |  |
| Gas window coat | Aluminium | 0.1 µm | `GasWindow_Al` | drift side |
| Window→cathode gas | chamber gas | 5 mm | `WindowGapGas` | set by the 5 mm flange |
| Drift cathode | Kapton | 50 µm | `DriftCathode_Kapton` |  |
| Drift cathode Cu | Copper | 9 µm | `DriftCathode_Cu` | cladding, drift side |
| **Drift gas** ★ | **chamber gas** | 29.991 mm | **`DriftGas`** | kapton back → mesh front = 30.000 mm |
| Micromesh | Stainless steel | 38 µm | `Micromesh` | woven 19/48 µm, fill 0.223 |
| **Amp gas** ★ | **chamber gas** | 150 µm | **`AmpGas`** | bulk pillars Ø0.6 mm on a 4.68 mm grid stand in it |
| Resistive layer (ESL) | C/acrylic paste | 10 µm | `ResistLayer` + `ESLStrip` | 515 strips, 550 µm on 0.8 mm; grooves are `AmpGas` |
| Coverlay kapton | Kapton | 50 µm | `PCB_Kapton` | **the ESL⇄pad insulator** — see the note below |
| Lamination adhesive | Acrylic | 18.8 µm | — *not placed* | response-model estimate, not in the header |
| Cu L3 guard ring | Copper | 26 µm | `PCB_Cu_3` | density-scaled |
| FR4 under L3 | Epoxy-glass | 264.4 µm | `PCB_FR4_3` | **residual** of the CAD 1.70 mm board |
| Cu L4 pads | Copper | 26 µm | `PCB_Cu_4` | real 512×512 pattern |
| FR4 under L4 | Epoxy-glass | 264.4 µm | `PCB_FR4_4` | **residual** of the CAD 1.70 mm board |
| Cu L5 Y strips | Copper | 26 µm | `PCB_Cu_5` | real 512×512 pattern |
| FR4 under L5 | Epoxy-glass | 264.4 µm | `PCB_FR4_5` | **residual** of the CAD 1.70 mm board |
| Cu L6 X strips | Copper | 26 µm | `PCB_Cu_6` | real 512×512 pattern |
| FR4 under L6 | Epoxy-glass | 264.4 µm | `PCB_FR4_6` | **residual** of the CAD 1.70 mm board |
| Cu L7 fan-out | Copper | 26 µm | `PCB_Cu_7` | density-scaled |
| FR4 under L7 | Epoxy-glass | 264.4 µm | `PCB_FR4_7` | **residual** of the CAD 1.70 mm board |
| Rohacell 51 | PMI foam | 5 mm | `PCB_Rohacell` |  |
| Backing mylar | Mylar (PET) | 25 µm | `PCB_BackMylar` | aluminized |
| Backing Al | Aluminium | 0.1 µm | `PCB_AlFoil` |  |
| Support plate | Aluminium | 8 mm | `SupportPlate_Al` | 402 mm aperture |

Depths are in beam order from the window front face; the stack totals **49.835 mm** to the back of the support plate (the readout board body is exactly the CAD **1.700 mm** from mesh front to laminate back). Board **470²**, **399.36 mm** readout window. ★ = scored volume.

**Two cautions about the bottom half of this table.**

*The FR4 is a residual, not a measurement.* The CAD gives one 1.70 mm number for the whole board body; everything we think we know is subtracted and the remainder is divided equally over the five laminate layers. So changing any other layer silently changes the FR4 — when the ESL went 100 → 10 µm the filler went 246.4 → 264.4 µm. A fab stackup drawing would turn this into a real number (`design/NEEDED_INPUTS.md` §6).

*Do not read the ESL⇄pad capacitance off this ordering.* Taken literally the table puts a 264 µm FR4 layer between the coverlay and the pad copper, because the guard-ring layer above the pads is placed first — and that layer has **zero** copper inside the active window, so the separation there is an artifact of the equal-division residual above, not geometry anyone measured. The response model instead takes the physically motivated adjacency for a resistive Micromegas: the ESL sits on the coverlay kapton, which is laminated straight onto the pads. That is the stack `response/solver/wpot.py` solves, and the ~19 µm of adhesive is 29 % of the insulator it sees. The two models disagree here and the disagreement is unresolved; it does not matter for material budget (order barely matters for energy deposit) and it matters completely for the induced signal. The same missing fab drawing settles both.

<!-- END STACK TABLE -->

Everything downstream of the module is unchanged:

- **He-3 target**: r = 1.5 cm, L = 8 cm, 300 bar; Al (0.5 mm) + CFRP (0.9 mm) capsule
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
--bulge-front <mm> Front-window overpressure dome sag (default: 8, 0 = flat)
--legacy-geometry  Pre-2026-08 uniform-slab MM module (old stack, bit-identical;
                   never patterned, regardless of the flags below)
--homogenized-readout  Density-scaled copper sheets + ESL slab instead of the
                   default real pattern (~3 % faster; see "Real readout pattern")
--patterned-readout Explicitly request the real pattern (this is the default)
--dump-geometry <f> Write the constructed geometry to JSON and exit (feeds
                   scripts/model/plot_mx17_model.py so figures show the true geometry)
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

> **AFS vs EOS rule**: HTCondor's submit file (`executable`, `log`, `output`, `error`)
> must reference AFS paths. Simulation output files are written by the job script and
> may live on EOS. Each submit script enforces this automatically:
> `--jobdir` (AFS, default `~/condor/...`) controls job control files;
> `--outdir` (EOS) controls where ROOT/CSV output is written.

### Sr-90 calibration (full material stack, energy sweep)

```bash
python scripts/submit_condor_sr90.py \
    --exe /eos/home-d/dneff/sim/mm_sim \
    --outdir /eos/home-d/dneff/sim/sr90_out \
    --mode sr90nomm
# jobdir defaults to ~/condor/mx17_geant_sr90 (AFS)
```

### LS calibration (liquid scintillator only, Sr-90 spectrum sampling)

```bash
python scripts/submit_condor_lscalib.py \
    --exe      /eos/home-d/dneff/sim/mm_sim \
    --spectrum /eos/home-d/dneff/sim/sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \
    --outdir   /eos/home-d/dneff/sim/lscalib_out \
    --mode ls --njobs 20 --nevents 100000
# jobdir defaults to ~/condor/mx17_lscalib (AFS)
```

### Back scint calibration (plastic scint only, Sr-90 spectrum sampling)

```bash
python scripts/submit_condor_lscalib.py \
    --exe      /eos/home-d/dneff/sim/mm_sim \
    --spectrum /eos/home-d/dneff/sim/sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \
    --outdir   /eos/home-d/dneff/sim/backscint_out \
    --mode backscint --njobs 20 --nevents 100000
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
