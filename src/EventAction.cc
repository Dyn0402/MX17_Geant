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
    // TRUTH IMPACT POINT (audit C14). The event carries its own primary
    // vertex, so record it rather than making the response chain infer the
    // impact point from a cluster centroid. `beamSpread` travels with it so a
    // consumer can tell that Stage A ALREADY randomised the impact point and
    // must not randomise it a second time in post.
    if (event->GetNumberOfPrimaryVertex() > 0) {
        const G4PrimaryVertex* v = event->GetPrimaryVertex(0);
        fData.vertexX = v->GetX0() / mm;
        fData.vertexY = v->GetY0() / mm;
        fData.vertexZ = v->GetZ0() / mm;
    }
    fData.beamSpread = fConfig.beam_spread_mm;

    // Pass event summary to RunAction for accumulation / writing
    fRunAction->RecordEvent(fData);

    if (fConfig.verbose && (fData.eventID % 1000 == 0)) {
        G4cout << "[Event " << fData.eventID << "]"
               // SteppingAction already stores these AS eV (edep/eV), so a
               // second /eV inflated the verbose print by 1e6. Files were
               // never affected -- this line is the only consumer.
               << "  Edep_drift=" << fData.edepDrift << " eV"
               << "  N_prim_drift=" << fData.nPrimaryDrift
               << "  Edep_amp=" << fData.edepAmp << " eV"
               << "  N_prim_amp=" << fData.nPrimaryAmp
               << G4endl;
    }
}
