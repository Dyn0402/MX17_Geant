#pragma once
// PrimaryGeneratorAction.hh
// Vacuum mode : pencil beam along +z at z = -10 cm
// Full mode   : beam from centre of He-3 gas volume along +z

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SimConfig.hh"
#include <memory>

class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const SimConfig& cfg,
                           const DetectorConstruction* detCon = nullptr);
    ~PrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* event) override;

private:
    std::unique_ptr<G4ParticleGun> fGun;
    const SimConfig&               fConfig;
    const DetectorConstruction*    fDetCon;
};
