// SteppingAction.cc
// Core physics: for every step in the gas volumes, compute
//   N_primary = Edep / W_value
// where W is the mean energy per ion pair (includes excitations).
//
// W-value references (MIP, room temperature, 1 atm):
//   Ar       : 26.4 eV  (ICRU 31)
//   Ar/CO2   : ~27 eV   (interpolated)
//   Ar/Iso   : ~26 eV   (Sauli 1977, close to pure Ar)
//   Ne       : 36.4 eV  (ICRU 31)
//   Ne/CO2   : ~34 eV   (Biagi estimates)
//   He       : 41.3 eV  (ICRU 31)
//   He/Eth   : ~29 eV   (Penning-enhanced; ethane Penning transfer)
//
// Note: Geant4 does not natively simulate individual primary clusters or
// Fano-limited fluctuations in gas -- it gives continuous energy loss.
// We use W-value to convert Edep -> N_electrons, which is the standard
// approach for Garfield++-style primary ionization estimation.

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

// W-values (mean energy per ion pair, in eV) for the five gas mixtures.
// Sources: ICRU Report 31, Sauli (1977), Blum/Riegler/Rolandi "PD with Drift Chambers",
//          and Biagi's Magboltz calculations.
//
// Ar/CF4 90/10  : ~34 eV.  CF4 raises W slightly vs pure Ar; no Penning.
// He/Ethane 96.5: ~27 eV.  Strong Penning from He*(19.8 eV) -> C2H6 (IE 11.5 eV).
//                          W much lower than pure He (41.3 eV).
// Ar/CO2 70/30  : ~27 eV.  CO2 quenches; mild Penning possible.
// Ar/CF4/Iso 88/10/2: ~33 eV. CF4 dominant quencher, isobutane Penning mild.
// Ne/Iso 95/5   : ~27 eV.  Ne*(16.6 eV) > Iso IE (10.6 eV): Penning active.
// Ar/CF4/CO2 45/40/15: ~34 eV. High CF4 fraction dominates; no Penning; CO2 adds mild quenching.
//
const std::map<std::string, double> SteppingAction::kWValues = {
    {"ArCF4",    34.0},
    {"HeEth",    27.0},
    {"ArCO2",    27.0},
    {"ArCF4Iso", 33.0},
    {"NeIso",    27.0},
    {"NeCF4",    30.0},     // Ne*(16.6 eV) > CF4 IE(10.1 eV): Penning active; ~30 eV estimate
    {"ArCF4CO2", 34.0},    // CF4-rich (40%); no Penning with Ar; CO2 mild quencher
};

SteppingAction::SteppingAction(const SimConfig& cfg, EventAction* eventAction)
    : G4UserSteppingAction(), fConfig(cfg), fEventAction(eventAction) {}

double SteppingAction::GetWValue(const std::string& gas) const {
    auto it = kWValues.find(gas);
    return (it != kWValues.end()) ? it->second : 26.4;  // default to Ar
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.0) return;

    // Get volume name
    const G4VPhysicalVolume* pv =
        step->GetPreStepPoint()->GetPhysicalVolume();
    if (!pv) return;

    const std::string volName = pv->GetName();
    bool inDrift = (volName == "DriftGas");
    bool inAmp   = (volName == "AmpGas");
    if (!inDrift && !inAmp) return;

    // Track info
    const G4Track* track = step->GetTrack();
    int trackID   = track->GetTrackID();
    int parentID  = track->GetParentID();
    const std::string pName = track->GetDefinition()->GetParticleName();
    double ke = step->GetPreStepPoint()->GetKineticEnergy();

    // Position of the step midpoint
    G4ThreeVector pos = 0.5 * (step->GetPreStepPoint()->GetPosition() +
                                step->GetPostStepPoint()->GetPosition());

    // W-value for this gas
    double W = GetWValue(fConfig.gas) * eV;

    // Number of primary ion pairs (Poisson-like expectation value)
    int nPrimary = static_cast<int>(std::floor(edep / W));
    // Remainder: add 1 with probability = fractional part
    double frac = edep / W - nPrimary;
    if (G4UniformRand() < frac) nPrimary++;

    // Build cluster record
    IonizationCluster cluster;
    cluster.x = pos.x() / mm;
    cluster.y = pos.y() / mm;
    cluster.z = pos.z() / mm;
    cluster.edep           = edep / eV;
    cluster.nPrimary       = nPrimary;
    cluster.trackID        = trackID;
    cluster.parentID       = parentID;
    cluster.volumeName     = volName;
    cluster.particleName   = pName;
    cluster.kineticEnergy  = ke / MeV;

    // Flag if primary particle (trackID==1) reached the volume
    EventData& data = fEventAction->GetEventData();

    if (inDrift) {
        data.edepDrift    += edep / eV;
        data.nPrimaryDrift += nPrimary;
        data.driftClusters.push_back(cluster);
        if (trackID == 1) data.primaryInDrift = true;
    } else {
        data.edepAmp      += edep / eV;
        data.nPrimaryAmp   += nPrimary;
        data.ampClusters.push_back(cluster);
        if (trackID == 1) data.primaryInAmp = true;
    }
}
