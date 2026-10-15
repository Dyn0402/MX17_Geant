# Primary ionization: Garfield/Heed vs Geant4 on a bare gas slab (2026-10-15)

Follow-up to Alexandra's check (akallitss/P2_Geant_Garfield `scripts/analysis/heed_vs_geant4.py`,
commit bcb4734): her full-stack Geant4 gave 30-50 % more drift-gas primaries than Heed in argon.
Both her sim and MX17 make primaries inside Geant4 as `edep/W` (SteppingAction); the Python
response chain takes `nPrimary` as given. This test asks how much of the excess is the ionization
model and how much is the material around the gas.

- `heed_mip.cc`: Alexandra's Heed driver, unchanged (bare slab, delta transport on).
- `slab.cc`: the same bare slab in Geant4: vacuum world, gas slab 10x10 cm, 100 MeV e- at 0 deg,
  EM option4 with mm_sim's cuts (e+- 10 um, gamma 0.1 mm), fluo/Auger on.
  `slab <gas> <gap_mm> <n> <out> <opt4|pai|paiphot> [seed] [E_MeV]`.
  Gas = ideal-gas density at 20 C with proper mass fractions; `<gas>_mx17` builds it exactly as
  MX17 `DetectorConstruction.cc` does (0 C densities, volume fractions passed as mass fractions,
  SteppingAction W table).
- PAI: in Geant4 11.3 `G4EmParameters::AddPAIModel` is silently a no-op here (the activator checks
  `IsIonisationProcess()` before PreparePhysicsTable sets it), so `slab.cc` adds PAI through
  `G4EmConfigurator::SetExtraEmModel` on eIoni in the gas region.
- **The gas region must have its own `G4ProductionCuts` object.** G4EmModelManager maps each
  material-cuts couple to a region by cuts *pointer* (G4EmModelManager.cc:490), so a region left
  on the default cuts gives its PAI-only model set to every couple in the world. Sub-keV electrons
  in steel, copper or air then lose almost no energy and run for >1e5 steps. In `mm_sim` this
  slowed runs ~10x; in the vacuum-world slab it is invisible. Repro: see the header of `slab.cc`.
- `compare.py <heed dir> <slab dir>`: the table below.

Build: Garfield (ROOT 6.40) for heed_mip, Geant4 11.3.2 for slab; both with a 3-line CMakeLists.

## Result (100k events at 3.964 mm; 30k at 30 mm). Ratios are Geant4 / Heed, primaries per crossing

| case | opt4 mean | opt4 median | PAI mean | PAI median | MX17 gas, opt4 mean / median |
|---|---|---|---|---|---|
| Ar/iso 95/5, 3.964 mm | 1.035 | 1.10 | 0.998 | 1.00 | 1.12 / 1.21 |
| Ar/CF4/iso 88/10/2, 3.964 mm | 1.066 | 1.10 | 1.028 | 1.03 | (W=33) 0.92 / 0.97 |
| Ne/iso 95/5, 3.964 mm | 0.976 | 0.93 | 0.970 | 0.93 | (W=27) 1.38 / 1.36 |
| Ar/iso 95/5, 30 mm | 1.033 | 1.05 | 1.006 | 1.00 | 1.13 / 1.14 (PAI: 1.09 / 1.09) |

Alexandra's full-stack numbers (Ar/iso 3.964 mm): mean 53.7 vs Heed 39.7 (1.35), median 37 vs 29.
The bare slab gives 41.3 / 32 with the same physics, so ~+30 % of the mean (~+15 % of the median)
comes from the stack (deltas/backscatter from foils, mesh, PCB), which Heed's bare slab does not have.
The thin-layer peak shift (edep MPV 612 vs 338 eV at 4 mm Ar) is the Urban fluctuation model; PAI
brings it to 362 eV.

Heed's "W_eff" (sum cluster energy / sum electrons) is 29 eV at 4 mm but 27.6 eV at 30 mm: cluster
energy includes delta-electron energy that leaves the slab, so it is not a W to copy.

## Fix in mm_sim (2026-10-15)

`include/GasProperties.hh` (one gas table: ideal-gas density, mass fractions, Garfield W/Fano),
`GasRegion` + PAI for e±/mu±/p/pi± in `PhysicsList`, with the region's own cuts. MX17 module in
vacuum, Ar/iso, 6 x 5000 events, events with drift primaries:

| beam | old mean / median | new mean / median | new / old |
|---|---|---|---|
| 100 MeV e- | 382.8 / 306 | 341.8 / 273 | 0.893 / 0.892 |
| 150 GeV mu- | 402.1 / 325 | 365.3 / 294 | 0.908 / 0.905 |

The bare-slab prediction for the combined change at 30 mm was 1/1.126 = 0.888. CPU per event is
unchanged (5000 e- events: 7.15 s new vs 7.13 s old).

`HANDOFF_ALEXANDRA.md`: the claims for the P2 repo to confirm.
