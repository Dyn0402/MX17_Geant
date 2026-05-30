// RunAction.cc
// Writes ROOT TTrees (or CSV fallback) per worker thread.

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

struct RunAction::Impl {
#ifdef USE_ROOT
    TFile* rootFile  = nullptr;
    TTree* evtTree   = nullptr;
    TTree* clusTree  = nullptr;

    // EventTree — always present
    Int_t    br_eventID;
    Double_t br_edepDrift, br_edepAmp;
    Int_t    br_nPrimDrift, br_nPrimAmp;
    Int_t    br_nClusDrift, br_nClusAmp;
    Bool_t   br_primInDrift, br_primInAmp;

    // EventTree — full/sr90 modes
    Double_t br_edepHe3Gas, br_edepResistPaste;
    Double_t br_edepMylar, br_edepCathode, br_edepMicromesh;
    Double_t br_edepPCB;
    Double_t br_edepPCBKapton, br_edepPCBCu, br_edepPCBFR4,
             br_edepPCBRohacell, br_edepPCBAlFoil;
    Double_t br_edepScintWall, br_edepScintTape, br_edepScintAlFoil;
    Double_t br_edepLS1, br_edepLS2;
    Bool_t   br_primInHe3Gas, br_primInPCB, br_primInScintWall;
    Bool_t   br_primInLS1, br_primInLS2;
    Double_t br_edepLSCFRP;
    Bool_t   br_primInLSCFRP5;

    // EventTree — kLSCalib mode
    Double_t br_edepBackScint, br_edepBackScintW, br_edepSourceCap;
    Bool_t   br_primInBackScint;

    // ClusterTree
    Int_t    cb_eventID, cb_trackID, cb_parentID, cb_nPrimary;
    Double_t cb_x, cb_y, cb_z, cb_edep, cb_ke;
    Char_t   cb_volume[32];
    Char_t   cb_particle[32];
#else
    std::ofstream evtFile;
    std::ofstream clusFile;
#endif
};

RunAction::RunAction(const SimConfig& cfg, bool isMaster)
    : G4UserRunAction(), fConfig(cfg), fIsMaster(isMaster),
      fImpl(std::make_unique<Impl>()) {}

RunAction::~RunAction() = default;

// ============================================================
void RunAction::BeginOfRunAction(const G4Run*) {
    fTotalEvents = fSumEdepDrift = fSumNprimDrift = 0.0;
    fSumEdepAmp  = fSumNprimAmp  = 0.0;

    if (fIsMaster) return;

    std::ostringstream ss;
    ss << fConfig.outFile;
    int tid = G4Threading::G4GetThreadId();
    if (tid >= 0) ss << "_t" << tid;

    bool isFull    = (fConfig.mode == SimMode::kFullExperiment  ||
                      fConfig.mode == SimMode::kSr90Calibration ||
                      fConfig.mode == SimMode::kSr90NoMM);
    bool isLSCalib = (fConfig.mode == SimMode::kLSCalib);

#ifdef USE_ROOT
    std::string fname = ss.str() + ".root";
    fImpl->rootFile = TFile::Open(fname.c_str(), "RECREATE");
    if (!fImpl->rootFile || fImpl->rootFile->IsZombie()) {
        G4cerr << "ERROR: Cannot open ROOT file " << fname << G4endl;
        return;
    }

    std::string modeStr = (fConfig.mode == SimMode::kFullExperiment  ? "full"    :
                           fConfig.mode == SimMode::kSr90Calibration ? "sr90"    :
                           fConfig.mode == SimMode::kSr90NoMM        ? "sr90nomm":
                           fConfig.mode == SimMode::kLSCalib         ? "lscalib" : "vacuum");
    std::string tag = "mode=" + modeStr
                    + "_gas=" + fConfig.gas
                    + "_particle=" + fConfig.particle
                    + "_E=" + std::to_string(fConfig.energy/MeV) + "MeV";

    fImpl->evtTree = new TTree("EventTree", "Per-event summary");
    fImpl->evtTree->SetTitle(tag.c_str());

    fImpl->evtTree->Branch("eventID",    &fImpl->br_eventID);
    fImpl->evtTree->Branch("edepDrift",  &fImpl->br_edepDrift);
    fImpl->evtTree->Branch("edepAmp",    &fImpl->br_edepAmp);
    fImpl->evtTree->Branch("nPrimDrift", &fImpl->br_nPrimDrift);
    fImpl->evtTree->Branch("nPrimAmp",   &fImpl->br_nPrimAmp);
    fImpl->evtTree->Branch("nClusDrift", &fImpl->br_nClusDrift);
    fImpl->evtTree->Branch("nClusAmp",   &fImpl->br_nClusAmp);
    fImpl->evtTree->Branch("primInDrift",&fImpl->br_primInDrift);
    fImpl->evtTree->Branch("primInAmp",  &fImpl->br_primInAmp);

    if (isFull) {
        fImpl->evtTree->Branch("edepHe3Gas",      &fImpl->br_edepHe3Gas);
        fImpl->evtTree->Branch("edepResistPaste", &fImpl->br_edepResistPaste);
        fImpl->evtTree->Branch("edepMylar",       &fImpl->br_edepMylar);
        fImpl->evtTree->Branch("edepCathode",     &fImpl->br_edepCathode);
        fImpl->evtTree->Branch("edepMicromesh",   &fImpl->br_edepMicromesh);
        fImpl->evtTree->Branch("edepPCB",         &fImpl->br_edepPCB);
        fImpl->evtTree->Branch("edepPCBKapton",   &fImpl->br_edepPCBKapton);
        fImpl->evtTree->Branch("edepPCBCu",       &fImpl->br_edepPCBCu);
        fImpl->evtTree->Branch("edepPCBFR4",      &fImpl->br_edepPCBFR4);
        fImpl->evtTree->Branch("edepPCBRohacell", &fImpl->br_edepPCBRohacell);
        fImpl->evtTree->Branch("edepPCBAlFoil",   &fImpl->br_edepPCBAlFoil);
        fImpl->evtTree->Branch("edepScintWall",   &fImpl->br_edepScintWall);
        fImpl->evtTree->Branch("edepScintTape",   &fImpl->br_edepScintTape);
        fImpl->evtTree->Branch("edepScintAlFoil", &fImpl->br_edepScintAlFoil);
        fImpl->evtTree->Branch("edepLS1",          &fImpl->br_edepLS1);
        fImpl->evtTree->Branch("edepLS2",          &fImpl->br_edepLS2);
        fImpl->evtTree->Branch("primInHe3Gas",     &fImpl->br_primInHe3Gas);
        fImpl->evtTree->Branch("primInPCB",        &fImpl->br_primInPCB);
        fImpl->evtTree->Branch("primInScintWall",  &fImpl->br_primInScintWall);
        fImpl->evtTree->Branch("primInLS1",        &fImpl->br_primInLS1);
        fImpl->evtTree->Branch("primInLS2",        &fImpl->br_primInLS2);
        fImpl->evtTree->Branch("edepLSCFRP",       &fImpl->br_edepLSCFRP);
        fImpl->evtTree->Branch("primInLSCFRP5",    &fImpl->br_primInLSCFRP5);
    }

    if (isLSCalib) {
        fImpl->evtTree->Branch("edepLS1",          &fImpl->br_edepLS1);
        fImpl->evtTree->Branch("edepLSCFRP",       &fImpl->br_edepLSCFRP);
        fImpl->evtTree->Branch("primInLS1",        &fImpl->br_primInLS1);
        fImpl->evtTree->Branch("edepBackScint",    &fImpl->br_edepBackScint);
        fImpl->evtTree->Branch("edepBackScintW",   &fImpl->br_edepBackScintW);
        fImpl->evtTree->Branch("edepSourceCap",    &fImpl->br_edepSourceCap);
        fImpl->evtTree->Branch("primInBackScint",  &fImpl->br_primInBackScint);
    }

    fImpl->clusTree = new TTree("ClusterTree", "Per-cluster ionization detail");
    fImpl->clusTree->SetTitle(tag.c_str());
    fImpl->clusTree->Branch("eventID",  &fImpl->cb_eventID);
    fImpl->clusTree->Branch("trackID",  &fImpl->cb_trackID);
    fImpl->clusTree->Branch("parentID", &fImpl->cb_parentID);
    fImpl->clusTree->Branch("x",        &fImpl->cb_x);
    fImpl->clusTree->Branch("y",        &fImpl->cb_y);
    fImpl->clusTree->Branch("z",        &fImpl->cb_z);
    fImpl->clusTree->Branch("edep",     &fImpl->cb_edep);
    fImpl->clusTree->Branch("nPrimary", &fImpl->cb_nPrimary);
    fImpl->clusTree->Branch("ke",       &fImpl->cb_ke);
    fImpl->clusTree->Branch("volume",   fImpl->cb_volume,  "volume[32]/C");
    fImpl->clusTree->Branch("particle", fImpl->cb_particle,"particle[32]/C");

    G4cout << "RunAction: Opened " << fname << G4endl;

#else
    std::string evtFname  = ss.str() + "_events.csv";
    std::string clusFname = ss.str() + "_clusters.csv";

    fImpl->evtFile.open(evtFname);
    fImpl->evtFile << "eventID,edepDrift_eV,edepAmp_eV,nPrimDrift,nPrimAmp,"
                      "nClusDrift,nClusAmp,primInDrift,primInAmp";
    if (isFull) {
        fImpl->evtFile << ",edepHe3Gas_eV,edepResistPaste_eV"
                          ",edepMylar_eV,edepCathode_eV,edepMicromesh_eV"
                          ",edepPCB_eV,edepPCBKapton_eV,edepPCBCu_eV"
                          ",edepPCBFR4_eV,edepPCBRohacell_eV,edepPCBAlFoil_eV"
                          ",edepScintWall_eV,edepScintTape_eV,edepScintAlFoil_eV"
                          ",edepLS1_eV,edepLS2_eV"
                          ",primInHe3Gas,primInPCB,primInScintWall"
                          ",primInLS1,primInLS2"
                          ",edepLSCFRP_eV,primInLSCFRP5";
    }
    if (isLSCalib) {
        fImpl->evtFile << ",edepLS1_eV,edepLSCFRP_eV,primInLS1"
                          ",edepBackScint_eV,edepBackScintW_eV"
                          ",edepSourceCap_eV,primInBackScint";
    }
    fImpl->evtFile << "\n";

    fImpl->clusFile.open(clusFname);
    fImpl->clusFile << "eventID,trackID,parentID,x_mm,y_mm,z_mm,"
                       "edep_eV,nPrimary,ke_MeV,volume,particle\n";

    G4cout << "RunAction: Writing " << evtFname << G4endl;
#endif
}

// ============================================================
void RunAction::EndOfRunAction(const G4Run* run) {
    if (!fIsMaster) {
#ifdef USE_ROOT
        if (fImpl->rootFile) {
            fImpl->rootFile->Write();
            fImpl->rootFile->Close();
            delete fImpl->rootFile;
            fImpl->rootFile = nullptr;
        }
#else
        fImpl->evtFile.close();
        fImpl->clusFile.close();
#endif
    }

    if (fIsMaster || !G4Threading::IsMultithreadedApplication()) {
        G4int nev = run->GetNumberOfEvent();
        std::string modeStr = (fConfig.mode == SimMode::kFullExperiment  ? "full-experiment"  :
                               fConfig.mode == SimMode::kSr90Calibration ? "sr90-calibration" :
                               fConfig.mode == SimMode::kSr90NoMM        ? "sr90-no-mm"       :
                               fConfig.mode == SimMode::kLSCalib         ? "ls-calibration"   : "vacuum");
        G4cout << "\n========= Run Summary =========" << G4endl;
        G4cout << "  Mode   : " << modeStr         << G4endl;
        G4cout << "  Gas    : " << fConfig.gas      << G4endl;
        G4cout << "  Events : " << nev              << G4endl;
        if (fTotalEvents > 0) {
            G4cout << std::fixed << std::setprecision(2);
            G4cout << "  <Edep_drift>  : " << fSumEdepDrift  / fTotalEvents << " eV" << G4endl;
            G4cout << "  <Edep_amp>    : " << fSumEdepAmp    / fTotalEvents << " eV" << G4endl;
        }
        G4cout << "===============================" << G4endl;
    }
}

// ============================================================
void RunAction::RecordEvent(const EventData& data) {
    fTotalEvents++;
    fSumEdepDrift  += data.edepDrift;
    fSumNprimDrift += data.nPrimaryDrift;
    fSumEdepAmp    += data.edepAmp;
    fSumNprimAmp   += data.nPrimaryAmp;

    bool isFull    = (fConfig.mode == SimMode::kFullExperiment  ||
                      fConfig.mode == SimMode::kSr90Calibration ||
                      fConfig.mode == SimMode::kSr90NoMM);
    bool isLSCalib = (fConfig.mode == SimMode::kLSCalib);

#ifdef USE_ROOT
    if (!fImpl->evtTree) return;

    fImpl->br_eventID    = data.eventID;
    fImpl->br_edepDrift  = data.edepDrift;
    fImpl->br_edepAmp    = data.edepAmp;
    fImpl->br_nPrimDrift = data.nPrimaryDrift;
    fImpl->br_nPrimAmp   = data.nPrimaryAmp;
    fImpl->br_nClusDrift = static_cast<int>(data.driftClusters.size());
    fImpl->br_nClusAmp   = static_cast<int>(data.ampClusters.size());
    fImpl->br_primInDrift = data.primaryInDrift;
    fImpl->br_primInAmp   = data.primaryInAmp;

    if (isFull) {
        fImpl->br_edepHe3Gas      = data.edepHe3Gas;
        fImpl->br_edepResistPaste = data.edepResistPaste;
        fImpl->br_edepMylar       = data.edepMylar;
        fImpl->br_edepCathode     = data.edepCathode;
        fImpl->br_edepMicromesh   = data.edepMicromesh;
        fImpl->br_edepPCB         = data.edepPCB;
        fImpl->br_edepPCBKapton   = data.edepPCBKapton;
        fImpl->br_edepPCBCu       = data.edepPCBCu;
        fImpl->br_edepPCBFR4      = data.edepPCBFR4;
        fImpl->br_edepPCBRohacell = data.edepPCBRohacell;
        fImpl->br_edepPCBAlFoil   = data.edepPCBAlFoil;
        fImpl->br_edepScintWall   = data.edepScintWall;
        fImpl->br_edepScintTape   = data.edepScintTape;
        fImpl->br_edepScintAlFoil = data.edepScintAlFoil;
        fImpl->br_edepLS1         = data.edepLS1;
        fImpl->br_edepLS2         = data.edepLS2;
        fImpl->br_primInHe3Gas    = data.primInHe3Gas;
        fImpl->br_primInPCB       = data.primInPCB;
        fImpl->br_primInScintWall = data.primInScintWall;
        fImpl->br_primInLS1       = data.primInLS1;
        fImpl->br_primInLS2       = data.primInLS2;
        fImpl->br_edepLSCFRP      = data.edepLSCFRP;
        fImpl->br_primInLSCFRP5   = data.primInLSCFRP5;
    }

    if (isLSCalib) {
        fImpl->br_edepLS1         = data.edepLS1;
        fImpl->br_edepLSCFRP      = data.edepLSCFRP;
        fImpl->br_primInLS1       = data.primInLS1;
        fImpl->br_edepBackScint   = data.edepBackScint;
        fImpl->br_edepBackScintW  = data.edepBackScintW;
        fImpl->br_edepSourceCap   = data.edepSourceCap;
        fImpl->br_primInBackScint = data.primInBackScint;
    }

    fImpl->evtTree->Fill();

    auto fillClusters = [&](const std::vector<IonizationCluster>& clusters) {
        for (const auto& c : clusters) {
            fImpl->cb_eventID  = data.eventID;
            fImpl->cb_trackID  = c.trackID;
            fImpl->cb_parentID = c.parentID;
            fImpl->cb_x        = c.x;
            fImpl->cb_y        = c.y;
            fImpl->cb_z        = c.z;
            fImpl->cb_edep     = c.edep;
            fImpl->cb_nPrimary = c.nPrimary;
            fImpl->cb_ke       = c.kineticEnergy;
            std::strncpy(fImpl->cb_volume,   c.volumeName.c_str(),   31);
            std::strncpy(fImpl->cb_particle, c.particleName.c_str(), 31);
            fImpl->cb_volume[31] = fImpl->cb_particle[31] = '\0';
            fImpl->clusTree->Fill();
        }
    };
    fillClusters(data.driftClusters);
    fillClusters(data.ampClusters);

#else
    fImpl->evtFile << data.eventID << ","
                   << data.edepDrift     << "," << data.edepAmp     << ","
                   << data.nPrimaryDrift << "," << data.nPrimaryAmp << ","
                   << data.driftClusters.size() << "," << data.ampClusters.size() << ","
                   << data.primaryInDrift << "," << data.primaryInAmp;
    if (isFull) {
        fImpl->evtFile << "," << data.edepHe3Gas
                       << "," << data.edepResistPaste
                       << "," << data.edepMylar
                       << "," << data.edepCathode
                       << "," << data.edepMicromesh
                       << "," << data.edepPCB
                       << "," << data.edepPCBKapton
                       << "," << data.edepPCBCu
                       << "," << data.edepPCBFR4
                       << "," << data.edepPCBRohacell
                       << "," << data.edepPCBAlFoil
                       << "," << data.edepScintWall
                       << "," << data.edepScintTape
                       << "," << data.edepScintAlFoil
                       << "," << data.edepLS1
                       << "," << data.edepLS2
                       << "," << data.primInHe3Gas
                       << "," << data.primInPCB
                       << "," << data.primInScintWall
                       << "," << data.primInLS1
                       << "," << data.primInLS2
                       << "," << data.edepLSCFRP
                       << "," << data.primInLSCFRP5;
    }
    if (isLSCalib) {
        fImpl->evtFile << "," << data.edepLS1
                       << "," << data.edepLSCFRP
                       << "," << data.primInLS1
                       << "," << data.edepBackScint
                       << "," << data.edepBackScintW
                       << "," << data.edepSourceCap
                       << "," << data.primInBackScint;
    }
    fImpl->evtFile << "\n";

    auto writeClusters = [&](const std::vector<IonizationCluster>& clusters) {
        for (const auto& c : clusters) {
            fImpl->clusFile << data.eventID << "," << c.trackID << "," << c.parentID << ","
                            << c.x << "," << c.y << "," << c.z << ","
                            << c.edep << "," << c.nPrimary << ","
                            << c.kineticEnergy << ","
                            << c.volumeName << "," << c.particleName << "\n";
        }
    };
    writeClusters(data.driftClusters);
    writeClusters(data.ampClusters);
#endif
}
