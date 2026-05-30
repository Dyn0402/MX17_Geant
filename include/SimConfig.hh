#pragma once
// SimConfig.hh
// Shared simulation configuration passed through the G4 user classes

#include <string>
#include "G4SystemOfUnits.hh"

enum class SimMode { kVacuum, kFullExperiment, kSr90Calibration, kSr90NoMM, kLSCalib };

struct SimConfig {
    std::string gas;            // Gas mixture name
    std::string particle;       // Particle type: gamma, neutron, electron, positron, proton, muon
    double      energy;         // Beam energy [Geant4 internal units = MeV]
    int         nEvents;        // Number of primary events
    std::string outFile;        // Output file base (no extension)
    long        seed;           // Random seed
    int         nThreads;       // MT threads
    bool        verbose;        // Verbose flag
    double      alThickness_mm; // Al shielding thickness [mm], 0 = no shielding (vacuum mode only)
    SimMode     mode            = SimMode::kVacuum;

    // ── LS cell wall parameters (from Full_Geant geometry) ────────────────
    double cfrpThickness_mm  = 2.0;    // Structural CFRP wall thickness [mm]
    double ls_inner_cfrp_um  = 600.0;  // Inner CFRP liner before each LAB layer [µm]
    double ls_inner_al_um    = 40.0;   // Al liner before each LAB layer [µm]
    double ls_thick_cm       = 2.0;    // LAB layer thickness [cm] (2 layers total)

    // ── Back plastic scintillator (kLSCalib mode) ─────────────────────────
    double backscint_u_cm     = 30.0;  // Back scint face: u-width [cm]
    double backscint_v_cm     = 20.0;  // Back scint face: v-height [cm]
    double backscint_thick_cm = 2.0;   // Back scint thickness [cm]
    double backscint_tape_um  = 200.0; // Outer black mylar tape [µm]
    double backscint_al_um    = 20.0;  // Al foil on scint surface [µm]

    // ── kLSCalib geometry parameters ─────────────────────────────────────
    double source_to_ls_mm   = 100.0;  // Air gap: source capsule → LS front face [mm]
    double gap_ls_to_back_mm = 50.0;   // Air gap: LS back face → back scint front [mm]

    // ── Source capsule thin layers (kLSCalib; thickness=0 → layer skipped) ─
    double source_cap_mylar_um  = 25.0;   // Front mylar window [µm]
    double source_cap_carbon_um = 100.0;  // Carbon disc [µm]
    double source_cap_al_um     = 100.0;  // Al backing [µm]
    double source_cap_tape_um   = 100.0;  // Black mylar tape [µm]
};
