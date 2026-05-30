// EventAction.cc

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

EventAction::EventAction(const SimConfig& cfg, RunAction* runAction)
    : G4UserEventAction(), fConfig(cfg), fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fData.Reset();
    fData.eventID = event->GetEventID();

    // Record the sampled primary KE (useful when spectrum sampling is on)
    if (event->GetNumberOfPrimaryVertex() > 0) {
        const G4PrimaryVertex* vtx = event->GetPrimaryVertex(0);
        if (vtx && vtx->GetNumberOfParticle() > 0) {
            const G4PrimaryParticle* p = vtx->GetPrimary(0);
            if (p) fData.primaryKE_MeV = p->GetKineticEnergy() / MeV;
        }
    }
}

void EventAction::EndOfEventAction(const G4Event* event) {
    // Pass event summary to RunAction for accumulation / writing
    fRunAction->RecordEvent(fData);

    if (fConfig.verbose && (fData.eventID % 1000 == 0)) {
        G4cout << "[Event " << fData.eventID << "]"
               << "  Edep_drift=" << fData.edepDrift / eV << " eV"
               << "  N_prim_drift=" << fData.nPrimaryDrift
               << "  Edep_amp=" << fData.edepAmp / eV << " eV"
               << "  N_prim_amp=" << fData.nPrimaryAmp
               << G4endl;
    }
}
