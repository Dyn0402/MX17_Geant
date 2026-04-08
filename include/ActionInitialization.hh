#pragma once
// ActionInitialization.hh
// Sets up PrimaryGeneratorAction, RunAction, EventAction, SteppingAction

#include "G4VUserActionInitialization.hh"
#include "SimConfig.hh"

class ActionInitialization : public G4VUserActionInitialization {
public:
    explicit ActionInitialization(const SimConfig& cfg);
    ~ActionInitialization() override = default;

    void BuildForMaster() override;
    void Build() const override;

private:
    const SimConfig& fConfig;
};
