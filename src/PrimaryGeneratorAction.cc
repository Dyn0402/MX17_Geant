// PrimaryGeneratorAction.cc

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <stdexcept>
#include <map>

PrimaryGeneratorAction::PrimaryGeneratorAction(const SimConfig& cfg)
    : G4VUserPrimaryGeneratorAction(), fConfig(cfg) {

    fGun = std::make_unique<G4ParticleGun>(1);

    // Map config string to G4 particle name
    static const std::map<std::string, std::string> particleMap = {
        {"gamma",    "gamma"},
        {"neutron",  "neutron"},
        {"electron", "e-"},
        {"proton",   "proton"},
        {"muon",     "mu-"},
        {"muon+",    "mu+"},
        {"pion",     "pi-"},
        {"alpha",    "alpha"},
    };

    auto it = particleMap.find(cfg.particle);
    if (it == particleMap.end()) {
        throw std::runtime_error("Unknown particle: " + cfg.particle +
            "\nAvailable: gamma, neutron, electron, proton, muon, muon+, pion, alpha");
    }

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = ptable->FindParticle(it->second);
    if (!particle) {
        throw std::runtime_error("G4 particle not found: " + it->second);
    }

    fGun->SetParticleDefinition(particle);
    fGun->SetParticleEnergy(cfg.energy);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));  // along +z
    // Start 2 mm before the detector front face
    // World extends from -worldZ/2, detector starts at -totalZ/2 ~ -16 mm
    // Just fire from -10 cm, which is safely in the world air gap
    fGun->SetParticlePosition(G4ThreeVector(0, 0, -10.0 * cm));
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    fGun->GeneratePrimaryVertex(event);
}
