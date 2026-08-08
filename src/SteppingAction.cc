// SteppingAction.cc
// Gas volumes: compute primary ion pairs via W-value, record per-cluster detail.
// Full-experiment extra volumes: record per-layer total edep and transmission flags.
// kLSCalib mode: score LS layer + back scint bar only.

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "EventData.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdlib>

// Fano factor for the edep -> ion-pair conversion (plan §7 step 1, audit C8).
// 0.2 is the standard value for argon-based counting gases; it is a property
// of the gas, so if the gas ever becomes a scan axis this belongs beside
// kWValues rather than here.
// Overridable with MX17_FANO so the with/without comparison can be run at the
// SAME seed — the only way to see the effect, since the per-event electron
// count is dominated by delta rays and varies >100 % event to event. Setting
// MX17_FANO=0 reproduces the pre-2026-08-08 floor+Bernoulli behaviour exactly.
static double FanoFactor() {
    static const double f = [] {
        const char* s = std::getenv("MX17_FANO");
        return s ? std::atof(s) : 0.2;
    }();
    return f;
}
// Below this mean pair count a Gaussian is not a sensible description (it
// would return negative counts and lose the sub-W remainder), so the exact
// floor + Bernoulli treatment is kept. At F = 0.2 the Gaussian sigma reaches
// 1 pair at nbar = 5, which is where the two descriptions meet.
static const double kFanoMinNbar = 5.0;

const std::map<std::string, double> SteppingAction::kWValues = {
    {"ArCF4",    34.0},
    {"HeEth",    27.0},
    {"ArCO2",    27.0},
    {"ArCF4Iso", 33.0},
    {"NeIso",    27.0},
    {"NeCF4",    30.0},
    {"ArCF4CO2", 34.0},
    {"ArIso",    26.0},
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

    // ── Gas volumes: ionisation cluster scoring (all non-vacuum modes) ───
    bool inDrift = (volName == "DriftGas");
    bool inAmp   = (volName == "AmpGas");

    if (inDrift || inAmp) {
        int parentID = track->GetParentID();
        const std::string pName = track->GetDefinition()->GetParticleName();
        double ke = step->GetPreStepPoint()->GetKineticEnergy();
        G4ThreeVector pos = 0.5*(step->GetPreStepPoint()->GetPosition() +
                                  step->GetPostStepPoint()->GetPosition());
        double tGlobal = 0.5*(step->GetPreStepPoint()->GetGlobalTime() +
                              step->GetPostStepPoint()->GetGlobalTime());
        const G4VProcess* creator = track->GetCreatorProcess();
        const std::string creatorName =
            creator ? std::string(creator->GetProcessName()) : std::string("primary");

        // ── edep -> primary electrons, with the FANO factor (audit C8) ──────
        //
        // Plan §7 step 1 specifies F ~ 0.2 and neither stage implemented it.
        // The old code was floor(edep/W) plus a Bernoulli on the remainder,
        // which is very nearly DETERMINISTIC: its variance is frac(1-frac) <=
        // 0.25 per cluster, where Fano statistics give F * nbar. The conversion
        // step was therefore simulated far too narrow, and the plan promised
        // something the code did not do.
        //
        // Fano statistics are sub-Poisson (F < 1) because the ionisations
        // along a track are anti-correlated by energy conservation: having
        // spent energy on one ion pair leaves less for the next, so the count
        // fluctuates less than a Poisson of the same mean. A Gaussian of
        // variance F*nbar is the standard representation and is what §7 asks
        // for; it is only meaningful once nbar is not tiny, hence the guard
        // below.
        double W = GetWValue(fConfig.gas) * eV;
        double nbar = edep / W;
        int nPrimary;
        const double fano = FanoFactor();
        if (fano <= 0.0 || nbar < kFanoMinNbar) {
            // Too few pairs for a Gaussian to mean anything: keep the exact
            // floor + Bernoulli remainder, which conserves <n> = nbar and is
            // the right small-number limit. Most clusters land here (a MIP
            // cluster is ~1-2 electrons), so the total is dominated by the
            // cluster-count and energy-straggling fluctuations Geant4 already
            // models -- this term matters for the dense clusters in the tail.
            nPrimary = static_cast<int>(std::floor(nbar));
            double frac = nbar - nPrimary;
            if (G4UniformRand() < frac) nPrimary++;
        } else {
            double n = G4RandGauss::shoot(nbar, std::sqrt(fano * nbar));
            nPrimary = static_cast<int>(std::lround(n));
            if (nPrimary < 0) nPrimary = 0;   // truncate, do not reflect
        }

        IonizationCluster cluster;
        cluster.x            = pos.x()/mm;
        cluster.y            = pos.y()/mm;
        cluster.z            = pos.z()/mm;
        cluster.time         = tGlobal/ns;
        cluster.edep         = edep/eV;
        cluster.nPrimary     = nPrimary;
        cluster.trackID      = trackID;
        cluster.parentID     = parentID;
        cluster.volumeName   = volName;
        cluster.particleName = pName;
        cluster.creatorProcess = creatorName;
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

    if (fConfig.mode == SimMode::kVacuum) return;

    // ── kLSCalib: score LS layer only ───────────────────────────────────
    if (fConfig.mode == SimMode::kLSCalib) {
        if (volName == "LiqScint_1") {
            data.edepLS1 += edep/eV;
            if (isPrimary) data.primInLS1 = true;
            return;
        }
        if ((volName.size() >= 8  && volName.substr(0,8)  == "LS_CFRP_")    ||
            (volName.size() >= 11 && volName.substr(0,11) == "LS_InnerCFR") ||
            (volName.size() >= 6  && volName.substr(0,6)  == "LS_Al_")) {
            data.edepLSCFRP += edep/eV;
            return;
        }
        return;
    }

    // ── kBackScintCalib: score back scint bar only ───────────────────────
    if (fConfig.mode == SimMode::kBackScintCalib) {
        if (volName == "BackScint") {
            data.edepBackScint += edep/eV;
            if (isPrimary) data.primInBackScint = true;
            return;
        }
        return;
    }

    // ── Full / Sr90 modes: per-layer edep ───────────────────────────────
    if (volName == "ResistivePaste") {
        data.edepResistPaste += edep/eV;
        return;
    }
    if (volName == "He3Gas") {
        data.edepHe3Gas += edep/eV;
        if (isPrimary) data.primInHe3Gas = true;
        return;
    }
    if (volName == "GasWindow_Mylar") {
        data.edepMylar += edep/eV;
        return;
    }
    if (volName == "GasWindow_Al" ||
        volName == "DriftCathode_Kapton" ||
        volName == "DriftCathode_Cu") {
        data.edepCathode += edep/eV;
        return;
    }
    if (volName == "Micromesh") {
        data.edepMicromesh += edep/eV;
        return;
    }
    if (volName.size() >= 4 && volName.substr(0,4) == "PCB_") {
        data.edepPCB += edep/eV;
        if (isPrimary) data.primInPCB = true;
        if      (volName == "PCB_Kapton")   data.edepPCBKapton   += edep/eV;
        else if (volName.size()>=6 && volName.substr(0,6)=="PCB_Cu")  data.edepPCBCu       += edep/eV;
        else if (volName.size()>=7 && volName.substr(0,7)=="PCB_FR4") data.edepPCBFR4      += edep/eV;
        else if (volName == "PCB_Rohacell") data.edepPCBRohacell += edep/eV;
        else if (volName == "PCB_AlFoil")   data.edepPCBAlFoil   += edep/eV;
        return;
    }
    if (volName == "ScintWall_BlackTape1" || volName == "ScintWall_BlackTape2") {
        data.edepScintTape += edep/eV;
        return;
    }
    if (volName == "ScintWall_AlFoil") {
        data.edepScintAlFoil += edep/eV;
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
    // LS structural CFRP walls + inner CFRP liners + Al liners → all into edepLSCFRP
    if (volName.size() >= 8 && volName.substr(0,8) == "LS_CFRP_") {
        data.edepLSCFRP += edep/eV;
        // LS_CFRP_3 is the back wall (exited LS2) — reusing primInLSCFRP5 flag
        if (isPrimary && volName == "LS_CFRP_3") data.primInLSCFRP5 = true;
        return;
    }
    if (volName.size() >= 11 && volName.substr(0,11) == "LS_InnerCFR") {
        data.edepLSCFRP += edep/eV;
        return;
    }
    if (volName.size() >= 6 && volName.substr(0,6) == "LS_Al_") {
        data.edepLSCFRP += edep/eV;
        return;
    }
}
