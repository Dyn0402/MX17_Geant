#pragma once
// EventData.hh
// Plain data container for one simulated event.
// Passed between EventAction and SteppingAction.

#include "G4ThreeVector.hh"
#include <vector>

struct IonizationCluster {
    double x, y, z;       // position [mm]
    double edep;          // energy deposited in this step [eV]
    int    nPrimary;      // number of primary ion pairs estimated for this step
    int    trackID;
    int    parentID;
    std::string volumeName;
    std::string particleName;
    double kineticEnergy;  // particle KE at step start [MeV]
};

struct EventData {
    int eventID = -1;

    // ── Micromegas gas scoring ──────────────────────────────────────────────
    double edepDrift = 0.0;  // [eV]
    double edepAmp   = 0.0;  // [eV]

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
    double edepMylar      = 0.0;
    double edepCathode    = 0.0;  // GasWin_Al + Kapton + Cu combined
    double edepMicromesh  = 0.0;

    // ── PCB stack: aggregate + individual components ──────────────────────
    double edepPCB         = 0.0;
    double edepPCBKapton   = 0.0;
    double edepPCBCu       = 0.0;
    double edepPCBFR4      = 0.0;
    double edepPCBRohacell = 0.0;
    double edepPCBAlFoil   = 0.0;

    // ── Scintillator wall breakdown ───────────────────────────────────────
    double edepScintWall   = 0.0;  // 3mm plastic scintillator
    double edepScintTape   = 0.0;  // black mylar tape layers
    double edepScintAlFoil = 0.0;

    // ── Liquid scintillator stack (2 layers × 2cm, from Full_Geant) ───────
    double edepLS1     = 0.0;
    double edepLS2     = 0.0;
    double edepLSCFRP  = 0.0;  // all structural + inner CFRP/Al lining walls

    // ── Back plastic scintillator (kLSCalib mode) ─────────────────────────
    double edepBackScint  = 0.0;  // PVT scint bar
    double edepBackScintW = 0.0;  // back scint wrapping (tape + Al foil)

    // ── Source capsule (kLSCalib mode) ────────────────────────────────────
    double edepSourceCap  = 0.0;  // total source capsule edep

    // ── Transmission flags ─────────────────────────────────────────────────
    bool primInHe3Gas    = false;
    bool primInPCB       = false;
    bool primInScintWall = false;
    bool primInLS1       = false;
    bool primInLS2       = false;
    bool primInLSCFRP5   = false;  // primary exited LS stack through back wall
    bool primInBackScint = false;  // primary reached back scint bar

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
        edepBackScint = edepBackScintW = edepSourceCap = 0.0;
        primInHe3Gas = primInPCB = primInScintWall = false;
        primInLS1 = primInLS2 = primInLSCFRP5 = primInBackScint = false;
    }
};
