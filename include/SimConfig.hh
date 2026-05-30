#pragma once
// SimConfig.hh
// Shared simulation configuration passed through the G4 user classes

#include <string>
#include "G4SystemOfUnits.hh"

enum class SimMode {
    kVacuum,
    kFullExperiment,
    kSr90Calibration,
    kSr90NoMM,
    kLSCalib,          // source → air → 1 LS layer
    kBackScintCalib,   // source → air → 1 back plastic scint bar
};

struct SimConfig {
    std::string gas;            // Gas mixture name
    std::string particle;       // Particle type
    double      energy;         // Beam energy [Geant4 internal units = MeV]
    int         nEvents;
    std::string outFile;
    long        seed;
    int         nThreads;
    bool        verbose;
    double      alThickness_mm;
    SimMode     mode = SimMode::kVacuum;

    // ── Spectrum sampling (kLSCalib / kBackScintCalib) ────────────────────
    // When non-empty, PrimaryGeneratorAction samples electron energies from
    // this CSV file (two columns: energy_MeV, probability) instead of using
    // the fixed 'energy' value above.
    std::string spectrum_file;

    // ── LS cell wall parameters (from Full_Geant geometry) ───────────────
    double cfrpThickness_mm  = 2.0;    // Structural CFRP wall [mm]
    double ls_inner_cfrp_um  = 600.0;  // Inner CFRP liner [µm]
    double ls_inner_al_um    = 40.0;   // Al liner [µm]
    double ls_thick_cm       = 2.0;    // LAB layer thickness [cm]

    // ── Back plastic scintillator (kBackScintCalib) ───────────────────────
    double backscint_u_cm     = 30.0;  // Back scint face: u-width [cm]
    double backscint_v_cm     = 20.0;  // Back scint face: v-height [cm]
    double backscint_thick_cm = 2.0;   // Back scint thickness [cm]
    double backscint_tape_um  = 200.0; // Outer black mylar tape [µm]
    double backscint_al_um    = 20.0;  // Al foil on scint surface [µm]

    // ── kLSCalib / kBackScintCalib: source-to-detector air gap ───────────
    double source_to_det_mm  = 100.0;  // Air gap from gun to detector front face [mm]
};
