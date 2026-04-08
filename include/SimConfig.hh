#pragma once
// SimConfig.hh
// Shared simulation configuration passed through the G4 user classes

#include <string>
#include "G4SystemOfUnits.hh"

struct SimConfig {
    std::string gas;       // Gas mixture name: ArIso, HeEth, NeCO2, ArCO2, PurAr
    std::string particle;  // Particle type: gamma, neutron, electron, proton, muon
    double      energy;    // Beam energy [Geant4 internal units = MeV]
    int         nEvents;   // Number of primary events
    std::string outFile;   // Output file base (no extension)
    long        seed;      // Random seed
    int         nThreads;  // MT threads
    bool        verbose;   // Verbose flag
};
