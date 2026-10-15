# Handoff: primary-ionization check, Heed vs Geant4 (from Dylan's MX17 session, 2026-10-15)

For Alexandra's Claude, working in `akallitss/P2_Geant_Garfield`. Please **reproduce the
numbers below independently** rather than trusting them, and report back where you disagree.
Everything needed is in this directory of `MX17_Geant`: `heed_mip.cc` (your file, unchanged),
`slab.cc` (Geant4 bare slab), `compare.py`, and `README.md` (first-pass table).

## 1. Starting point: your result

`scripts/analysis/heed_vs_geant4.py` (commit bcb4734, slide `heedg4`): for 100 MeV e- at 0 deg,
Ar/iso 95/5, 3.964 mm, your Geant4 `nPrimDrift` gives mean 53.7 and median 37, against Heed's
39.7 and 29 (ratio 1.35). Neon agrees at the peak. Your slide already notes that the comparison
is not like-for-like: Heed is a bare gas slab, while Geant4 is the full stack.

Both repos make primaries the same way: `SteppingAction` converts each gas step's
`edep` to pairs with a fixed W. So your test also covers MX17. No Python code makes primaries.

## 2. Claims to confirm (each with how we got it)

### C1. Heed reproduces
Build `heed_mip.cc` against Garfield (ROOT 6.40 here) and run
`heed_mip ArIso 3.964 0 100000 heed_ArIso`. We get **mean 39.88, median 29**, matching yours.
`ArCF4Iso`: 41.52 / 31. `NeIso`: 18.69 / 14. 30 mm Ar/iso: 314.0 / 268.

### C2. On a bare slab, Geant4 option4 is only slightly high
`slab.cc` puts the same slab in Geant4: a vacuum world, a 10x10 cm gas box, the same 100 MeV e-
gun, `G4EmStandardPhysics_option4`, and your cuts (e± 10 um, gamma 0.1 mm), with fluo/Auger and
ignore-cut on. The gas is the ideal gas at 20 C with proper mass fractions, as in your fixed
`GasMixtures`, and W is linear in volume fraction (26.38 eV for Ar/iso).
Run: `slab ArIso 3.964 100000 g4_ArIso_4mm_opt4 opt4`.

| 3.964 mm, 100 MeV e- | Heed mean / median | G4 opt4 | G4 PAI |
|---|---|---|---|
| Ar/iso 95/5 | 39.88 / 29 | 41.27 / 32 (**1.035 / 1.10**) | 39.78 / 29 (**0.998 / 1.00**) |
| Ar/CF4/iso 88/10/2 | 41.52 / 31 | 44.28 / 34 (1.066 / 1.10) | 42.67 / 32 (1.028 / 1.03) |
| Ne/iso 95/5 | 18.69 / 14 | 18.24 / 13 (0.976 / 0.93) | 18.14 / 13 (0.970 / 0.93) |
| Ar/iso, **30 mm** | 314.0 / 268 | 324.2 / 281 (1.033 / 1.05) | 316.0 / 269 (1.006 / 1.00) |
| Ar/iso, 150 GeV mu- | 42.62 / 31 | 43.51 / 34 (1.021 / 1.10) | 42.51 / 32 (0.997 / 1.03) |

In energy terms, the most probable deposit in 4 mm of Ar/iso is 338 eV for Heed, **612 eV for
option4** and 362 eV for PAI. The option4 problem is the thin-layer peak position, which comes
from the Urban fluctuation model. The mean is only 3.5 % high.

**Implication:** of your 1.35, only ~1.035 is the ionization model. The other **~1.30 on the
mean (~1.15 on the median)** is the material around the gas: deltas and backscatter from the
cathode foil, mesh and pad copper, which Heed's bare slab doesn't have. That part is probably
real physics, and it is in the data too. Heed can't check it. Please confirm by running your
full-stack `mm_sim` with the gas-only PAI fix (C3) and comparing it with the bare slab.

### C3. The obvious way to enable PAI is a silent no-op in Geant4 11.3
`G4EmParameters::AddPAIModel("all"|"e-", "<region>", "pai")` does nothing with option4 here.
`G4EmModelActivator::ActivatePAI` loops over the energy-loss processes and skips those with
`IsIonisationProcess() == false`. In 11.3 that flag is set only in
`G4VEnergyLossProcess::PreparePhysicsTable` (G4VEnergyLossProcess.cc:221), which runs *after*
the activator. Symptom: PAI runs are bit-identical to opt4 at the same seed, and with EM verbose
2 you see "ActivatePAI for 1 regions" but no "add <pai> model" line.
**Working route:** override `ConstructProcess()`, call the base, then
`G4LossTableManager::Instance()->EmConfigurator()->SetExtraEmModel(particle, "eIoni"|"muIoni"|"hIoni", new G4PAIModel(...), region, emin, 100*TeV, sameModelAsFluct)`.
See `slab.cc` and MX17 `src/PhysicsList.cc`. Your slide's "next: PAI in the gas region" would
hit this. **Do not use this route without C4.**

### C4. The gas region needs its own production-cuts object (or every run crawls)
After C3, `mm_sim` ran **~10x slower** in some seeds. A worker thread was stuck in Goudsmit-Saunderson msc,
with 0.3-3 keV electrons in steel/copper/kapton/mylar/air taking 1e5+ nm-sized steps and losing
~1e-4 eV per step. Cause: `G4EmModelManager::Initialise` assigns each material-cuts couple to a
region by comparing `G4ProductionCuts` pointers (G4EmModelManager.cc:490 in 11.3). Our
`GasRegion` had no cuts object of its own, so it inherited the default one, and every couple in the
world matched it and got the region's PAI-only model set. The observed energy loss outside the
gas is then ~0. We assume this is because PAI's tables are built only for the region's couples,
but we did not check that step. The bare-slab gas numbers in C2 are unaffected (rerun with the
fix: PAI 39.83 / 29). **Fix:** `region->SetProductionCuts(new G4ProductionCuts(...))` with
the same values (MX17 `PhysicsList::SetCuts`).
Minimal repro: `SLAB_SHARED_CUTS=1 SLAB_WORLD=G4_STAINLESS-STEEL slab ArIso 3.964 200 z pai 1 0.001`.
Every 1 keV electron hits the 1e5-step kill; without `SLAB_SHARED_CUTS` they stop in 2 steps.
**For P2:** your gas already sits in `P2FineCut`, which has its own 1 um cuts, so you would not hit
this trap if PAI went there. But P2FineCut also holds the mesh, the pad copper and the cathode.
PAI is valid physics in solids, but it would replace the ionisation model in the layers whose
surface-escape physics you tuned the 1 um cuts for.
Cleaner: move DriftGas/AmpGas into a new region with its *own* 1 um cuts object and PAI, and leave
the solids in P2FineCut. (PAI in MX17's kapton pillars,
which are AmpGas daughters, costs only ~5 % once the cuts are right.)

### C5. Do not take a W from Heed's output
Heed's sum(cluster energy) / sum(electrons) is 29.1 eV at 4 mm but 27.6 eV at 30 mm. Cluster
energy includes delta-electron energy that leaves the slab, so it is not a W. Heed's actual W is
`MediumGas::GetW()`: Garfield's ICRU-31 table, linear in mole fraction.

### C6. W table: isobutane
Garfield (`MediumGas::GetGasInfo`) uses Ar 26.4, Ne 35.4, iC4H10 **23.4**, CF4 34.3, CO2 33.0,
C2H6 25.0, He 41.3. Your `GasMixtures.cc` uses iC4H10 26.0 (your own caveat 2) and Ne 36.4.
That moves Ar/iso 95/5 from 26.38 to 26.25 eV (-0.5 %), Ne/iso 95/5 from 35.88 to 34.8 eV (-3 %),
and Ne/iso 85/15 by about -4 %. These are small but systematic. Matching Garfield makes Geant4
and Heed consistent by construction. Fano factors per gas are in the same Garfield function.

## 3. What changed in MX17 (commit in this repo, see `git log -- include/GasProperties.hh`)

- `include/GasProperties.hh`, `src/GasProperties.cc`: one gas table (molar mass, formula,
  Garfield W and Fano), ideal-gas density at 20 C, and mass fractions from volume fractions.
  Both `DetectorConstruction` and `SteppingAction` read it. MX17 still had the 7.3 % density
  bug and the volume-as-mass-fraction bug you fixed on 2026-08-05, plus a hand-written W table
  (NeIso 27 eV, ArCF4Iso 33 eV).
- `GasRegion` = DriftGas + AmpGas, with its own cuts object (C4).
- `PhysicsList::ConstructProcess`: PAI via EmConfigurator for e±, mu±, p, pi± in GasRegion.
  Also fluo + Auger + deexcitation-ignore-cut, as in your list.
- Fano per gas from the table (Ar/iso 0.175; it was a flat 0.2). `MX17_FANO` still overrides.

MX17 `mm_sim`, vacuum mode, Ar/iso, old vs new build (6 x 5000 events each), events with drift
primaries:

| beam | old mean / median | new mean / median | new / old |
|---|---|---|---|
| 100 MeV e- | 382.8 / 306 | 341.8 / 273 | 0.893 / 0.892 |
| 150 GeV mu- | 402.1 / 325 | 365.3 / 294 | 0.908 / 0.905 |

The bare-slab prediction was 0.888. CPU per event is unchanged.

## 4. Open questions we did NOT settle

1. **Is the +30 % stack excess right?** It depends on the delta/backscatter yield from thin foils
   and copper at your 1 um cuts. A direct check would be data. For example, compare the most
   probable cluster charge on tracks with and without the cathode foil, or run a Geant4 variant
   with the mylar/Al removed.
2. **Penning / Jesse effect on W.** Ar metastables (11.55 eV) lie above the iC4H10 ionization
   potential (10.67 eV), so W in Ar/iso can sit a few % below the linear mix. Neither Heed nor
   Geant4 includes this. It is a shared systematic, and Heed agreeing with us doesn't remove it.
3. **Neon at the median** (0.93 even with PAI) is not understood. It could be the W (Heed's Ne
   35.4 vs our 35.88 here) or the PAI photoabsorption data for Ne.

## 5. How to report back

A short note to Dylan: confirm or contest C1-C6, with your own numbers, and say whether your
full-stack ratio drops to ~1.30 with PAI on the gas. When you enable PAI, check the time per
event before and after, and check for >1e5-step tracks (C4).
