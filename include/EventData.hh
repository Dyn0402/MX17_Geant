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

    // Total energy deposition per volume
    double edepDrift = 0.0;  // [eV]
    double edepAmp   = 0.0;  // [eV]

    // Primary ionization cluster list (one per ionizing step in gas)
    std::vector<IonizationCluster> driftClusters;
    std::vector<IonizationCluster> ampClusters;

    // Aggregate primary ion pair counts
    int nPrimaryDrift = 0;
    int nPrimaryAmp   = 0;

    // Did primary particle enter each volume?
    bool primaryInDrift = false;
    bool primaryInAmp   = false;

    void Reset() {
        eventID = -1;
        edepDrift = edepAmp = 0.0;
        nPrimaryDrift = nPrimaryAmp = 0;
        primaryInDrift = primaryInAmp = false;
        driftClusters.clear();
        ampClusters.clear();
    }
};
