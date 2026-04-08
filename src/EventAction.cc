// EventAction.cc

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

EventAction::EventAction(const SimConfig& cfg, RunAction* runAction)
    : G4UserEventAction(), fConfig(cfg), fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fData.Reset();
    fData.eventID = event->GetEventID();
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
