// PrimaryGeneratorAction.cc
// Vacuum mode      : pencil beam at z = -10 cm, along +z.
// Full-experiment  : pencil beam from the centre of the He-3 gas volume, along +z.

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <stdexcept>
#include <map>

PrimaryGeneratorAction::PrimaryGeneratorAction(const SimConfig& cfg,
                                               const DetectorConstruction* detCon)
    : G4VUserPrimaryGeneratorAction(), fConfig(cfg), fDetCon(detCon) {

    fGun = std::make_unique<G4ParticleGun>(1);

    static const std::map<std::string, std::string> particleMap = {
        {"gamma",    "gamma"},
        {"neutron",  "neutron"},
        {"electron", "e-"},
        {"proton",   "proton"},
        {"muon",     "mu-"},
        {"muon+",    "mu+"},
        {"pion",     "pi-"},
        {"alpha",    "alpha"},
        {"triton",   "triton"},
    };

    auto it = particleMap.find(cfg.particle);
    if (it == particleMap.end()) {
        throw std::runtime_error("Unknown particle: " + cfg.particle +
            "\nAvailable: gamma, neutron, electron, proton, muon, muon+, pion, alpha, triton");
    }

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = ptable->FindParticle(it->second);
    if (!particle) {
        throw std::runtime_error("G4 particle not found: " + it->second);
    }

    fGun->SetParticleDefinition(particle);
    fGun->SetParticleEnergy(cfg.energy);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

    // Gun position: He-3 gas centre in full mode, fixed z = -10 cm in vacuum mode.
    // DetectorConstruction::Construct() is called before Build(), so
    // GetHe3GasCenterZ() is valid here.
    G4double gunZ = -10.0 * cm;
    if (cfg.mode == SimMode::kFullExperiment && fDetCon) {
        gunZ = fDetCon->GetHe3GasCenterZ();
    }
    fGun->SetParticlePosition(G4ThreeVector(0, 0, gunZ));
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    fGun->GeneratePrimaryVertex(event);
}
