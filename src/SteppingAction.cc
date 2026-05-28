// SteppingAction.cc
// Gas volumes: compute primary ion pairs via W-value, record per-cluster detail.
// Full-experiment extra volumes: record per-layer total edep and transmission flags.

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "EventData.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "Randomize.hh"

#include <cmath>

// W-values (mean energy per ion pair, eV) — unchanged from vacuum-mode code.
const std::map<std::string, double> SteppingAction::kWValues = {
    {"ArCF4",    34.0},
    {"HeEth",    27.0},
    {"ArCO2",    27.0},
    {"ArCF4Iso", 33.0},
    {"NeIso",    27.0},
    {"NeCF4",    30.0},
    {"ArCF4CO2", 34.0},
    {"ArIso",    26.0},   // Ar/iC4H10 95/5; mild Penning (Ar* 11.55 eV > iso IE 10.6 eV)
    {"PureCF4",  34.0},
    {"PureAr",   26.4},
    {"PureHe",   41.3},
    {"PureNe",   36.4},
    {"PureEthane",26.0},
    {"PureIso",  26.0},
    {"PureCO2",  33.0},
};

SteppingAction::SteppingAction(const SimConfig& cfg, EventAction* eventAction)
    : G4UserSteppingAction(), fConfig(cfg), fEventAction(eventAction) {}

double SteppingAction::GetWValue(const std::string& gas) const {
    auto it = kWValues.find(gas);
    return (it != kWValues.end()) ? it->second : 26.4;
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.0) return;

    const G4VPhysicalVolume* pv = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!pv) return;

    const std::string volName = pv->GetName();
    const G4Track*    track   = step->GetTrack();
    int trackID  = track->GetTrackID();
    bool isPrimary = (trackID == 1);

    EventData& data = fEventAction->GetEventData();

    // ── Gas volumes: detailed cluster scoring (both modes) ───
    bool inDrift = (volName == "DriftGas");
    bool inAmp   = (volName == "AmpGas");

    if (inDrift || inAmp) {
        int parentID = track->GetParentID();
        const std::string pName = track->GetDefinition()->GetParticleName();
        double ke = step->GetPreStepPoint()->GetKineticEnergy();
        G4ThreeVector pos = 0.5*(step->GetPreStepPoint()->GetPosition() +
                                  step->GetPostStepPoint()->GetPosition());

        double W = GetWValue(fConfig.gas) * eV;
        int nPrimary = static_cast<int>(std::floor(edep / W));
        double frac  = edep / W - nPrimary;
        if (G4UniformRand() < frac) nPrimary++;

        IonizationCluster cluster;
        cluster.x            = pos.x()/mm;
        cluster.y            = pos.y()/mm;
        cluster.z            = pos.z()/mm;
        cluster.edep         = edep/eV;
        cluster.nPrimary     = nPrimary;
        cluster.trackID      = trackID;
        cluster.parentID     = parentID;
        cluster.volumeName   = volName;
        cluster.particleName = pName;
        cluster.kineticEnergy = ke/MeV;

        if (inDrift) {
            data.edepDrift     += edep/eV;
            data.nPrimaryDrift += nPrimary;
            data.driftClusters.push_back(cluster);
            if (isPrimary) data.primaryInDrift = true;
        } else {
            data.edepAmp     += edep/eV;
            data.nPrimaryAmp += nPrimary;
            data.ampClusters.push_back(cluster);
            if (isPrimary) data.primaryInAmp = true;
        }
        return;
    }

    // ── Full-experiment extra volumes: per-layer edep only ───
    if (fConfig.mode != SimMode::kFullExperiment) return;

    if (volName == "ResistivePaste") {
        data.edepResistPaste += edep/eV;
        return;
    }

    if (volName == "He3Gas") {
        data.edepHe3Gas += edep/eV;
        if (isPrimary) data.primInHe3Gas = true;
        return;
    }

    // PCB stack: aggregate all sublayers
    if (volName.size() >= 4 && volName.substr(0,4) == "PCB_") {
        data.edepPCB += edep/eV;
        if (isPrimary) data.primInPCB = true;
        return;
    }

    if (volName == "PlasticScint") {
        data.edepScintWall += edep/eV;
        if (isPrimary) data.primInScintWall = true;
        return;
    }

    if (volName == "LiqScint_1") {
        data.edepLS1 += edep/eV;
        if (isPrimary) data.primInLS1 = true;
        return;
    }
    if (volName == "LiqScint_2") {
        data.edepLS2 += edep/eV;
        if (isPrimary) data.primInLS2 = true;
        return;
    }
    if (volName == "LiqScint_3") {
        data.edepLS3 += edep/eV;
        if (isPrimary) data.primInLS3 = true;
        return;
    }
    if (volName == "LiqScint_4") {
        data.edepLS4 += edep/eV;
        if (isPrimary) data.primInLS4 = true;
        return;
    }
}
