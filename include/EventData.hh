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

    // ── Micromegas gas scoring (both modes) ─────────────────────────────
    double edepDrift = 0.0;  // [eV]
    double edepAmp   = 0.0;  // [eV]

    std::vector<IonizationCluster> driftClusters;
    std::vector<IonizationCluster> ampClusters;

    int  nPrimaryDrift = 0;
    int  nPrimaryAmp   = 0;
    bool primaryInDrift = false;
    bool primaryInAmp   = false;

    // ── Full-experiment per-layer edep [eV] (kFullExperiment mode only) ─
    double edepHe3Gas      = 0.0;  // He-3 gas volume
    double edepResistPaste = 0.0;  // resistive paste (100 µm, behind amp gas)
    double edepPCB         = 0.0;  // PCB stack total
    double edepScintWall   = 0.0;  // plastic scintillator bar
    double edepLS1         = 0.0;  // liquid scintillator layer 1
    double edepLS2         = 0.0;
    double edepLS3         = 0.0;
    double edepLS4         = 0.0;
    double edepLSCFRP      = 0.0;  // total edep in all 5 LS cell CFRP walls [eV]

    // ── Transmission flags (did primary particle reach each subsystem?) ─
    bool primInHe3Gas    = false;
    bool primInPCB       = false;
    bool primInScintWall = false;
    bool primInLS1       = false;
    bool primInLS2       = false;
    bool primInLS3       = false;
    bool primInLS4       = false;
    bool primInLSCFRP5   = false;  // entered back wall → exited LS stack

    void Reset() {
        eventID = -1;
        edepDrift = edepAmp = 0.0;
        nPrimaryDrift = nPrimaryAmp = 0;
        primaryInDrift = primaryInAmp = false;
        driftClusters.clear();
        ampClusters.clear();

        edepHe3Gas = edepResistPaste = edepPCB = edepScintWall = 0.0;
        edepLS1 = edepLS2 = edepLS3 = edepLS4 = edepLSCFRP = 0.0;
        primInHe3Gas = primInPCB = primInScintWall = false;
        primInLS1 = primInLS2 = primInLS3 = primInLS4 = primInLSCFRP5 = false;
    }
};
