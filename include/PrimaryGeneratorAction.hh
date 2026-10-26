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
    void LoadPrimariesTable(const std::string& filepath);
    double SampleSpectrum() const;

    std::unique_ptr<G4ParticleGun> fGun;
    G4ThreeVector fGunPos;                 // nominal gun position
    double        fBeamSpread = 0.0;       // transverse spread [G4 length]
    const SimConfig&               fConfig;
    const DetectorConstruction*    fDetCon;

    // Sr-90/Y-90 spectrum: CDF table for inverse-transform sampling
    std::vector<double> fSpecEnergies;  // energy values [MeV]
    std::vector<double> fSpecCDF;       // cumulative probabilities [0,1]
    bool fUseSpectrum = false;

    // --primaries-table: per table event, the primaries it fires
    struct Primary { G4ParticleDefinition* def; double E_MeV, ux, uy, uz, dx_mm, dy_mm; };
    std::vector<std::vector<Primary>> fTable;
};
