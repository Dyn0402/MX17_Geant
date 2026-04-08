#pragma once
// SteppingAction.hh
// Intercepts every step in the gas volumes to:
//   1. Record energy deposition
//   2. Estimate number of primary ion pairs via W-value
//   3. Record cluster position

#include "G4UserSteppingAction.hh"
#include "SimConfig.hh"
#include <map>
#include <string>

class EventAction;
class G4Step;

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(const SimConfig& cfg, EventAction* eventAction);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    // W-value (mean energy per ion pair) for each gas [eV]
    // These are well-measured values from literature
    double GetWValue(const std::string& gas) const;

    const SimConfig& fConfig;
    EventAction*     fEventAction;

    // W-value lookup
    static const std::map<std::string, double> kWValues;
};
