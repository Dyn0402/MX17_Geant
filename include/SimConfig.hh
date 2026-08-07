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

    // ── MX17 module geometry (shared/MX17ModuleGeometry.hh) ──────────────
    bool   legacy_mm_geometry = false;  // true: pre-2026-08 uniform-slab module
    double mx17_bulge_front_mm = 8.0;   // front-window dome sag from gas
    // Real 512x512 readout copper (pads / X / Y tracks) + discrete ESL
    // resistive strips, instead of density-scaled sheets. Default ON: it costs
    // ~3 % CPU and no memory. --homogenized-readout falls back to the cheaper
    // density-scaled build. See shared/MX17ModuleGeometry.hh.
    bool   mx17_patterned_readout = true;
                                        // overpressure (Hencky ~3 mbar); 0 = flat

    // ── Spectrum sampling (kLSCalib / kBackScintCalib) ────────────────────
    // When non-empty, PrimaryGeneratorAction samples electron energies from
    // this CSV file (two columns: energy_MeV, probability) instead of using
    // the fixed 'energy' value above.
    std::string spectrum_file;
    // Transverse beam spread [mm]. 0 = a pencil beam at (0,0), which is the
    // historical default but is a TRAP for any position-dependent readout
    // study: 512 pads is even, so the active-area centre sits exactly on a pad
    // BOUNDARY and every primary lands where the "which channel is biggest"
    // assignment is a coin flip. For response-chain work set this to at least
    // the 31.2 mm ESL/pad superperiod so every strip phase is sampled.
    double beam_spread_mm = 0.0;

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
