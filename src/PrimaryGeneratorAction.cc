// PrimaryGeneratorAction.cc
// Vacuum mode          : pencil beam at z = -10 cm, along +z.
// Full / sr90 modes    : beam from detector stack front, along +z.
// kLSCalib / kBackScintCalib : bare electron gun at detector front, optionally
//                              sampling the Sr-90/Y-90 beta spectrum from a CSV file.

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
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
        {"positron", "e+"},
        {"proton",   "proton"},
        {"muon",     "mu-"},
        {"muon+",    "mu+"},
        {"pion",     "pi-"},
        {"alpha",    "alpha"},
        {"triton",   "triton"},
    };

    auto it = particleMap.find(cfg.particle);
    if (it == particleMap.end())
        throw std::runtime_error("Unknown particle: " + cfg.particle);

    G4ParticleDefinition* particle =
        G4ParticleTable::GetParticleTable()->FindParticle(it->second);
    if (!particle)
        throw std::runtime_error("G4 particle not found: " + it->second);

    fGun->SetParticleDefinition(particle);
    fGun->SetParticleEnergy(cfg.energy);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

    // Gun position
    G4double gunZ = -10.0 * cm;
    bool useDetZ = (cfg.mode == SimMode::kFullExperiment  ||
                    cfg.mode == SimMode::kSr90Calibration ||
                    cfg.mode == SimMode::kSr90NoMM        ||
                    cfg.mode == SimMode::kLSCalib          ||
                    cfg.mode == SimMode::kBackScintCalib);
    if (useDetZ && fDetCon)
        gunZ = fDetCon->GetHe3GasCenterZ();

    fGun->SetParticlePosition(G4ThreeVector(0, 0, gunZ));

    // Load Sr-90/Y-90 spectrum if requested
    if (!cfg.spectrum_file.empty()) {
        LoadSpectrum(cfg.spectrum_file);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::LoadSpectrum(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f.is_open())
        throw std::runtime_error("Cannot open spectrum file: " + filepath);

    std::vector<double> energies, weights;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Skip header lines (contain non-numeric first token)
        std::istringstream ss(line);
        double e, w;
        if (!(ss >> e >> w)) continue;
        if (e < 0 || w < 0) continue;
        energies.push_back(e);
        weights.push_back(w);
    }

    if (energies.size() < 2)
        throw std::runtime_error("Spectrum file too short: " + filepath);

    // Build normalised CDF
    double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    fSpecEnergies = energies;
    fSpecCDF.resize(weights.size());
    double cumsum = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
        cumsum += weights[i] / total;
        fSpecCDF[i] = cumsum;
    }
    fSpecCDF.back() = 1.0;  // ensure exact 1 at end
    fUseSpectrum = true;

    G4cout << "PrimaryGeneratorAction: Loaded spectrum from " << filepath
           << "  (" << fSpecEnergies.size() << " points, "
           << "E_max=" << fSpecEnergies.back() << " MeV)" << G4endl;
}

// ─────────────────────────────────────────────────────────────────────────────
double PrimaryGeneratorAction::SampleSpectrum() const {
    double r = G4UniformRand();
    // Inverse CDF via binary search
    auto it = std::lower_bound(fSpecCDF.begin(), fSpecCDF.end(), r);
    size_t i = std::distance(fSpecCDF.begin(), it);
    if (i >= fSpecEnergies.size()) i = fSpecEnergies.size() - 1;
    return fSpecEnergies[i];
}

// ─────────────────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    if (fUseSpectrum) {
        double E = SampleSpectrum();
        fGun->SetParticleEnergy(E * MeV);
    }
    fGun->GeneratePrimaryVertex(event);
}
