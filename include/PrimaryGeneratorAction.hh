#pragma once
// PrimaryGeneratorAction.hh
// Fires a pencil beam of configurable particle type and energy
// along +z, starting 1 mm before the gas window.

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SimConfig.hh"
#include <memory>

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const SimConfig& cfg);
    ~PrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* event) override;

private:
    std::unique_ptr<G4ParticleGun> fGun;
    const SimConfig& fConfig;
};
