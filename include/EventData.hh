#pragma once
// EventData.hh — per-event data container

#include "G4ThreeVector.hh"
#include <vector>

struct IonizationCluster {
    double x, y, z;
    double edep;
    int    nPrimary;
    int    trackID;
    int    parentID;
    std::string volumeName;
    std::string particleName;
    double kineticEnergy;
};

struct EventData {
    int eventID = -1;

    // ── Micromegas gas scoring ─────────────────────────────────────────────
    double edepDrift = 0.0;
    double edepAmp   = 0.0;

    std::vector<IonizationCluster> driftClusters;
    std::vector<IonizationCluster> ampClusters;

    int  nPrimaryDrift = 0;
    int  nPrimaryAmp   = 0;
    bool primaryInDrift = false;
    bool primaryInAmp   = false;

    // ── Full-experiment per-layer edep [eV] ───────────────────────────────
    double edepHe3Gas      = 0.0;
    double edepResistPaste = 0.0;

    // ── MM entrance / dead layers ─────────────────────────────────────────
    double edepMylar     = 0.0;
    double edepCathode   = 0.0;
    double edepMicromesh = 0.0;

    // ── PCB stack ─────────────────────────────────────────────────────────
    double edepPCB         = 0.0;
    double edepPCBKapton   = 0.0;
    double edepPCBCu       = 0.0;
    double edepPCBFR4      = 0.0;
    double edepPCBRohacell = 0.0;
    double edepPCBAlFoil   = 0.0;

    // ── Scintillator wall ─────────────────────────────────────────────────
    double edepScintWall   = 0.0;
    double edepScintTape   = 0.0;
    double edepScintAlFoil = 0.0;

    // ── Liquid scintillator stack (2 layers × 2 cm) ───────────────────────
    double edepLS1    = 0.0;
    double edepLS2    = 0.0;
    double edepLSCFRP = 0.0;  // all structural + inner CFRP/Al walls combined

    // ── Back plastic scintillator (kBackScintCalib) ────────────────────────
    double edepBackScint = 0.0;  // PVT active volume [eV]

    // ── Primary beam kinetic energy (kLSCalib / kBackScintCalib) ──────────
    // Stored so the analysis can see what spectrum was sampled.
    double primaryKE_MeV = 0.0;

    // ── Transmission flags ─────────────────────────────────────────────────
    bool primInHe3Gas    = false;
    bool primInPCB       = false;
    bool primInScintWall = false;
    bool primInLS1       = false;
    bool primInLS2       = false;
    bool primInLSCFRP5   = false;  // primary exited LS stack (reached back CFRP wall)
    bool primInBackScint = false;

    void Reset() {
        eventID = -1;
        edepDrift = edepAmp = 0.0;
        nPrimaryDrift = nPrimaryAmp = 0;
        primaryInDrift = primaryInAmp = false;
        driftClusters.clear();
        ampClusters.clear();
        edepHe3Gas = edepResistPaste = 0.0;
        edepMylar = edepCathode = edepMicromesh = 0.0;
        edepPCB = edepPCBKapton = edepPCBCu = edepPCBFR4
                = edepPCBRohacell = edepPCBAlFoil = 0.0;
        edepScintWall = edepScintTape = edepScintAlFoil = 0.0;
        edepLS1 = edepLS2 = edepLSCFRP = 0.0;
        edepBackScint = 0.0;
        primaryKE_MeV = 0.0;
        primInHe3Gas = primInPCB = primInScintWall = false;
        primInLS1 = primInLS2 = primInLSCFRP5 = primInBackScint = false;
    }
};
