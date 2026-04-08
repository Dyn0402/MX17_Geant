// ActionInitialization.cc

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const SimConfig& cfg)
    : G4VUserActionInitialization(), fConfig(cfg) {}

void ActionInitialization::BuildForMaster() {
    // RunAction for master thread (handles file open/close)
    SetUserAction(new RunAction(fConfig, true));
}

void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGeneratorAction(fConfig));

    auto* runAction   = new RunAction(fConfig, false);
    auto* eventAction = new EventAction(fConfig, runAction);

    SetUserAction(runAction);
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(fConfig, eventAction));
}
