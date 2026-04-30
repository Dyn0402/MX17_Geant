// RunAction.cc
// Output strategy:
//   Compiled with ROOT -> writes two TTrees per run:
//     EventTree  : one entry per event  (aggregates)
//     ClusterTree: one entry per ionization cluster (full spatial detail)
//   Without ROOT  -> two CSV files (same columns).
//
// File naming: <outFile>_<threadID>.root  (merged by collect_results.py)
// For MT mode each worker writes its own file to avoid ROOT thread contention.

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

// ---- Internal state per RunAction instance ----
struct RunAction::Impl {
#ifdef USE_ROOT
    TFile* rootFile   = nullptr;
    TTree* evtTree    = nullptr;
    TTree* clusTree   = nullptr;

    // Event branch buffers
    Int_t    br_eventID;
    Double_t br_edepDrift, br_edepAmp;
    Int_t    br_nPrimDrift, br_nPrimAmp;
    Int_t    br_nClusDrift, br_nClusAmp;
    Bool_t   br_primInDrift, br_primInAmp;

    // Cluster branch buffers
    Int_t    cb_eventID, cb_trackID, cb_parentID, cb_nPrimary;
    Double_t cb_x, cb_y, cb_z, cb_edep, cb_ke;
    Char_t   cb_volume[32];
    Char_t   cb_particle[32];
#else
    std::ofstream evtFile;
    std::ofstream clusFile;
#endif
};

// ============================================================
RunAction::RunAction(const SimConfig& cfg, bool isMaster)
    : G4UserRunAction(), fConfig(cfg), fIsMaster(isMaster),
      fImpl(std::make_unique<Impl>()) {}

RunAction::~RunAction() = default;

// ============================================================
void RunAction::BeginOfRunAction(const G4Run* run) {
    fTotalEvents  = 0;
    fSumEdepDrift = fSumNprimDrift = 0.0;
    fSumEdepAmp   = fSumNprimAmp  = 0.0;

    // Master thread in MT mode does not write data
    if (fIsMaster) return;

    // Build filename: outFile_t<threadID>.root or .csv
    std::ostringstream ss;
    ss << fConfig.outFile;
    int tid = G4Threading::G4GetThreadId();
    if (tid >= 0) ss << "_t" << tid;

#ifdef USE_ROOT
    std::string fname = ss.str() + ".root";
    fImpl->rootFile = TFile::Open(fname.c_str(), "RECREATE");
    if (!fImpl->rootFile || fImpl->rootFile->IsZombie()) {
        G4cerr << "ERROR: Cannot open ROOT file " << fname << G4endl;
        return;
    }

    // --- EventTree ---
    fImpl->evtTree = new TTree("EventTree", "Per-event summary");
    fImpl->evtTree->Branch("eventID",    &fImpl->br_eventID);
    fImpl->evtTree->Branch("edepDrift",  &fImpl->br_edepDrift);   // eV
    fImpl->evtTree->Branch("edepAmp",    &fImpl->br_edepAmp);     // eV
    fImpl->evtTree->Branch("nPrimDrift", &fImpl->br_nPrimDrift);  // ion pairs
    fImpl->evtTree->Branch("nPrimAmp",   &fImpl->br_nPrimAmp);
    fImpl->evtTree->Branch("nClusDrift", &fImpl->br_nClusDrift);  // n steps
    fImpl->evtTree->Branch("nClusAmp",   &fImpl->br_nClusAmp);
    fImpl->evtTree->Branch("primInDrift",&fImpl->br_primInDrift);
    fImpl->evtTree->Branch("primInAmp",  &fImpl->br_primInAmp);

    // Metadata branches (constant per file, useful after hadd)
    fImpl->evtTree->SetTitle(
        ("gas=" + fConfig.gas + "_particle=" + fConfig.particle +
         "_E=" + std::to_string(fConfig.energy / MeV) + "MeV" +
         "_Al=" + std::to_string(fConfig.alThickness_mm) + "mm").c_str());

    // --- ClusterTree ---
    fImpl->clusTree = new TTree("ClusterTree", "Per-cluster ionization detail");
    fImpl->clusTree->Branch("eventID",   &fImpl->cb_eventID);
    fImpl->clusTree->Branch("trackID",   &fImpl->cb_trackID);
    fImpl->clusTree->Branch("parentID",  &fImpl->cb_parentID);
    fImpl->clusTree->Branch("x",         &fImpl->cb_x);           // mm
    fImpl->clusTree->Branch("y",         &fImpl->cb_y);           // mm
    fImpl->clusTree->Branch("z",         &fImpl->cb_z);           // mm
    fImpl->clusTree->Branch("edep",      &fImpl->cb_edep);        // eV
    fImpl->clusTree->Branch("nPrimary",  &fImpl->cb_nPrimary);    // ion pairs
    fImpl->clusTree->Branch("ke",        &fImpl->cb_ke);          // MeV
    fImpl->clusTree->Branch("volume",    fImpl->cb_volume,  "volume[32]/C");
    fImpl->clusTree->Branch("particle",  fImpl->cb_particle,"particle[32]/C");
    fImpl->clusTree->SetTitle(fImpl->evtTree->GetTitle());

    G4cout << "RunAction: Opened " << fname << G4endl;

#else
    // CSV fallback
    std::string evtFname  = ss.str() + "_events.csv";
    std::string clusFname = ss.str() + "_clusters.csv";

    fImpl->evtFile.open(evtFname);
    fImpl->evtFile << "eventID,edepDrift_eV,edepAmp_eV,nPrimDrift,nPrimAmp,"
                      "nClusDrift,nClusAmp,primInDrift,primInAmp\n";

    fImpl->clusFile.open(clusFname);
    fImpl->clusFile << "eventID,trackID,parentID,x_mm,y_mm,z_mm,"
                       "edep_eV,nPrimary,ke_MeV,volume,particle\n";

    G4cout << "RunAction: Writing " << evtFname << " + " << clusFname << G4endl;
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

    // Print run summary from master (or sequential mode)
    if (fIsMaster || !G4Threading::IsMultithreadedApplication()) {
        G4int nev = run->GetNumberOfEvent();
        G4cout << "\n========= Run Summary =========" << G4endl;
        G4cout << "  Gas      : " << fConfig.gas      << G4endl;
        G4cout << "  Particle : " << fConfig.particle << G4endl;
        G4cout << "  Energy   : " << fConfig.energy / MeV << " MeV" << G4endl;
        G4cout << "  Events   : " << nev << G4endl;
        if (fTotalEvents > 0) {
            G4cout << std::fixed << std::setprecision(2);
            G4cout << "  <Edep_drift>  : " << fSumEdepDrift  / fTotalEvents << " eV" << G4endl;
            G4cout << "  <Nprim_drift> : " << fSumNprimDrift / fTotalEvents << G4endl;
            G4cout << "  <Edep_amp>    : " << fSumEdepAmp    / fTotalEvents << " eV" << G4endl;
            G4cout << "  <Nprim_amp>   : " << fSumNprimAmp   / fTotalEvents << G4endl;
        }
        G4cout << "===============================" << G4endl;
    }
}

// ============================================================
void RunAction::RecordEvent(const EventData& data) {
    // Accumulate for summary
    fTotalEvents++;
    fSumEdepDrift  += data.edepDrift;
    fSumNprimDrift += data.nPrimaryDrift;
    fSumEdepAmp    += data.edepAmp;
    fSumNprimAmp   += data.nPrimaryAmp;

#ifdef USE_ROOT
    if (!fImpl->evtTree) return;

    // Fill EventTree
    fImpl->br_eventID    = data.eventID;
    fImpl->br_edepDrift  = data.edepDrift;
    fImpl->br_edepAmp    = data.edepAmp;
    fImpl->br_nPrimDrift = data.nPrimaryDrift;
    fImpl->br_nPrimAmp   = data.nPrimaryAmp;
    fImpl->br_nClusDrift = static_cast<int>(data.driftClusters.size());
    fImpl->br_nClusAmp   = static_cast<int>(data.ampClusters.size());
    fImpl->br_primInDrift = data.primaryInDrift;
    fImpl->br_primInAmp   = data.primaryInAmp;
    fImpl->evtTree->Fill();

    // Fill ClusterTree -- drift then amp
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
    // CSV fallback
    fImpl->evtFile << data.eventID << ","
                   << data.edepDrift    << "," << data.edepAmp    << ","
                   << data.nPrimaryDrift << "," << data.nPrimaryAmp << ","
                   << data.driftClusters.size() << "," << data.ampClusters.size() << ","
                   << data.primaryInDrift << "," << data.primaryInAmp << "\n";

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
