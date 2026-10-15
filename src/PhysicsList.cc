// PhysicsList.cc
// Builds the physics for the Micromegas simulation.
// Key considerations:
//   - Gammas  : photoelectric, Compton, pair production (all via EM option4)
//   - Neutrons: elastic + inelastic via FTFP_BERT (QGSP or FTFP + Bertini cascade)
//              + thermal neutrons via NeutronHP
//   - Electrons: accurate low-energy EM (Livermore model in option4 goes down to eV)
//   - Ionization: G4ionIonisation tracks energy loss in gas steps

#include "PhysicsList.hh"

// Modular physics components
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4PAIModel.hh"
#include "G4ParticleTable.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "DetectorConstruction.hh"   // kGasRegionName

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
    SetVerboseLevel(0);

    // EM physics: option4 uses Livermore models below 100 keV for e-/gamma
    // This gives accurate photoelectric + Auger + ionization in low-Z gas
    RegisterPhysics(new G4EmStandardPhysics_option4(0));

    // EM extras: synchrotron, GDR etc (harmless to include)
    RegisterPhysics(new G4EmExtraPhysics(0));

    // Hadronic elastic with HP (high-precision neutron data, <20 MeV)
    RegisterPhysics(new G4HadronElasticPhysicsHP(0));

    // Hadronic inelastic with HP neutrons
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP(0));

    // Decay
    RegisterPhysics(new G4DecayPhysics(0));

    // Radioactive decay (useful for activation studies)
    RegisterPhysics(new G4RadioactiveDecayPhysics(0));

    // Step limiter (respects G4UserLimits set in detector volumes)
    RegisterPhysics(new G4StepLimiterPhysics());

    // Kill neutrons below 1 eV after tracking to avoid infinite loops
    RegisterPhysics(new G4NeutronTrackingCut(0));

    // Atomic de-excitation. Set AFTER the EM constructor, whose ctor writes
    // the same G4EmParameters singleton. IgnoreCut matters most: the 0.1 mm
    // gamma range cut is ~5-10 keV in copper, so without it the Cu-K line
    // (8.05 keV) can sit under threshold and not be produced at all.
    auto* emp = G4EmParameters::Instance();
    emp->SetFluo(true);
    emp->SetAuger(true);
    emp->SetDeexcitationIgnoreCut(true);
}

void PhysicsList::ConstructProcess() {
    G4VModularPhysicsList::ConstructProcess();

    // PAI ionisation in the counting gas (DetectorConstruction's GasRegion).
    // option4's default energy-loss fluctuations (Urban) are built for thick
    // absorbers; in a few-mm gas layer they put the primary-electron peak
    // ~10 % high (most probable deposit 612 vs 338 eV for 100 MeV e- in 4 mm
    // Ar/iso). PAI, built from the gas photoabsorption cross section like
    // Heed, matches Heed to <1 % (response/validation/heed_vs_g4/).
    //
    // Added through G4EmConfigurator on purpose: in Geant4 11.3
    // G4EmParameters::AddPAIModel is a silent no-op with this list (the
    // activator tests IsIonisationProcess() before PreparePhysicsTable has
    // set it), so the documented route does nothing.
    // Lower limits follow G4EmModelActivator::ActivatePAI.
    auto* conf = G4LossTableManager::Instance()->EmConfigurator();
    auto* table = G4ParticleTable::GetParticleTable();
    const struct { const char* particle; const char* process; G4double emin; } pai[] = {
        {"e-", "eIoni", 0.11*keV},  {"e+", "eIoni", 0.11*keV},
        {"mu-", "muIoni", 5*keV},   {"mu+", "muIoni", 5*keV},
        {"proton", "hIoni", 50*keV}, {"pi+", "hIoni", 50*keV}, {"pi-", "hIoni", 50*keV},
    };
    for (const auto& p : pai) {
        auto* model = new G4PAIModel(table->FindParticle(p.particle), "PAI");
        conf->SetExtraEmModel(p.particle, p.process, model, kGasRegionName,
                              p.emin, 100*TeV, model);
    }
}

void PhysicsList::SetCuts() {
    // Default range cuts -- 1 mm is G4 default
    // For gas detectors at low pressure we want to be more careful:
    // - Electrons: short range cut to capture delta rays in gas
    // - Photons  : default is fine
    SetCutValue(1.0 * mm,  "proton");
    SetCutValue(10.0 * um, "e-");       // short cut for electrons -- capture delta rays
    SetCutValue(10.0 * um, "e+");
    SetCutValue(0.1  * mm, "gamma");    // 0.1 mm photon production threshold

    // The gas region (PAI, see ConstructProcess) needs its OWN cuts object,
    // even with the same values. G4EmModelManager assigns a material-cuts
    // couple to a region by comparing G4ProductionCuts POINTERS
    // (G4EmModelManager.cc:490 in 11.3), so a region left on the default cuts
    // hands its PAI-only model set to every couple in the world. Then sub-keV
    // electrons in steel, copper or air lose almost no energy and are tracked
    // for >1e5 steps; runs slowed ~10x before this.
    if (auto* reg = G4RegionStore::GetInstance()->GetRegion(kGasRegionName, false)) {
        auto* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(10.0 * um, "e-");
        cuts->SetProductionCut(10.0 * um, "e+");
        cuts->SetProductionCut(0.1  * mm, "gamma");
        cuts->SetProductionCut(1.0  * mm, "proton");
        reg->SetProductionCuts(cuts);
    }
}
