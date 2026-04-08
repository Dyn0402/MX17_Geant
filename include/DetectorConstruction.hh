#pragma once
// DetectorConstruction.hh
// Builds the full Micromegas detector geometry

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "SimConfig.hh"
#include <map>
#include <string>

class G4Material;
class G4Region;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const SimConfig& cfg);
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Accessors used by analysis
    G4LogicalVolume* GetDriftGasLV()  const { return fDriftGasLV; }
    G4LogicalVolume* GetAmpGasLV()    const { return fAmpGasLV; }

private:
    void DefineMaterials();
    G4Material* GetGasMixture(const std::string& name);

    const SimConfig& fConfig;
    G4LogicalVolume* fDriftGasLV = nullptr;
    G4LogicalVolume* fAmpGasLV   = nullptr;

    // Material cache
    std::map<std::string, G4Material*> fGasMaterials;
};
