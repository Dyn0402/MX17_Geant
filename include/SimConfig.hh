#pragma once
// SimConfig.hh
// Shared simulation configuration passed through the G4 user classes

#include <string>
#include "G4SystemOfUnits.hh"

enum class SimMode { kVacuum, kFullExperiment, kSr90Calibration };

struct SimConfig {
    std::string gas;            // Gas mixture name
    std::string particle;       // Particle type: gamma, neutron, electron, proton, muon
    double      energy;         // Beam energy [Geant4 internal units = MeV]
    int         nEvents;        // Number of primary events
    std::string outFile;        // Output file base (no extension)
    long        seed;           // Random seed
    int         nThreads;       // MT threads
    bool        verbose;        // Verbose flag
    double      alThickness_mm; // Al shielding thickness [mm], 0 = no shielding (vacuum mode only)
    SimMode     mode            = SimMode::kVacuum;
    double      cfrpThickness_mm = 1.5; // CFRP wall thickness [mm] for He-3 capsule and LS cells
};
