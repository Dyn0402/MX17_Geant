#pragma once
// PrimaryGeneratorAction.hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SimConfig.hh"
#include <memory>
#include <vector>
#include <string>

class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const SimConfig& cfg,
                           const DetectorConstruction* detCon = nullptr);
    ~PrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* event) override;

private:
    void LoadSpectrum(const std::string& filepath);
    double SampleSpectrum() const;

    std::unique_ptr<G4ParticleGun> fGun;
    const SimConfig&               fConfig;
    const DetectorConstruction*    fDetCon;

    // Sr-90/Y-90 spectrum: CDF table for inverse-transform sampling
    std::vector<double> fSpecEnergies;  // energy values [MeV]
    std::vector<double> fSpecCDF;       // cumulative probabilities [0,1]
    bool fUseSpectrum = false;
};
