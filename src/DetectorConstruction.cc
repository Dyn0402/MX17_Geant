// DetectorConstruction.cc
// Full Micromegas detector geometry:
//   Gas window   : 40 um Mylar + 0.1 um Al
//   Drift cathode: 50 um Kapton + 9 um Cu
//   Drift volume : 3 cm gas
//   Micromesh    : 30 um stainless steel (woven: 18 um wire, 45 um hole, 63 um pitch -> ~51% open)
//   Amp volume   : 150 um gas
//   Anode PCB    : 50 um Kapton + 9 um Cu (readout)
// Transverse area: 40 cm x 40 cm (your actual detector size)

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SDManager.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"

#include "SensitiveDetector.hh"

#include <stdexcept>

// ============================================================
DetectorConstruction::DetectorConstruction(const SimConfig& cfg)
    : G4VUserDetectorConstruction(), fConfig(cfg) {}

// ============================================================
void DetectorConstruction::DefineMaterials() {
    G4NistManager* nist = G4NistManager::Instance();

    // ----- Elements -----
    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elN  = nist->FindOrBuildElement("N");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elAr = nist->FindOrBuildElement("Ar");
    G4Element* elNe = nist->FindOrBuildElement("Ne");
    G4Element* elHe = nist->FindOrBuildElement("He");
    G4Element* elF  = nist->FindOrBuildElement("F");

    // =====================================================
    // Pure component gases (building blocks)
    // Densities at STP (273.15 K, 1 atm) from NIST/PDG
    // =====================================================

    // Isobutane C4H10: 2.67 kg/m3
    G4Material* isobutane = new G4Material("Isobutane", 2.67e-3 * g/cm3, 2,
                                            kStateGas, 293.15*kelvin, 1*atmosphere);
    isobutane->AddElement(elC, 4);
    isobutane->AddElement(elH, 10);

    // Ethane C2H6: 1.356 kg/m3
    G4Material* ethane = new G4Material("Ethane", 1.356e-3 * g/cm3, 2,
                                         kStateGas, 293.15*kelvin, 1*atmosphere);
    ethane->AddElement(elC, 2);
    ethane->AddElement(elH, 6);

    // CO2: 1.977 kg/m3
    G4Material* CO2 = new G4Material("CO2_gas", 1.977e-3 * g/cm3, 2,
                                      kStateGas, 293.15*kelvin, 1*atmosphere);
    CO2->AddElement(elC, 1);
    CO2->AddElement(elO, 2);

    // CF4 (tetrafluoromethane): 3.72 kg/m3
    G4Material* CF4 = new G4Material("CF4_gas", 3.72e-3 * g/cm3, 2,
                                      kStateGas, 293.15*kelvin, 1*atmosphere);
    CF4->AddElement(elC, 1);
    CF4->AddElement(elF, 4);

    // Pure Ar: 1.782 kg/m3
    G4Material* purAr = new G4Material("PureArgon", 1.782e-3 * g/cm3, 1,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
    purAr->AddElement(elAr, 1);

    // Pure He: 0.1786 kg/m3
    G4Material* pureHe = new G4Material("PureHe", 0.1786e-3 * g/cm3, 1,
                                         kStateGas, 293.15*kelvin, 1*atmosphere);
    pureHe->AddElement(elHe, 1);

    // Pure Ne: 0.900 kg/m3
    G4Material* pureNe = new G4Material("PureNe", 0.8999e-3 * g/cm3, 1,
                                         kStateGas, 293.15*kelvin, 1*atmosphere);
    pureNe->AddElement(elNe, 1);

    // =====================================================
    // Gas mixture 1: Ar/CF4 90/10 vol%
    // Density: weighted sum of component densities (ideal mixing)
    // =====================================================
    {
        G4double fAr=0.90, fCF4=0.10;
        G4double rho = fAr*1.782e-3 + fCF4*3.72e-3;
        G4Material* m = new G4Material("ArCF4", rho*g/cm3, 2,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr, fAr);
        m->AddMaterial(CF4,   fCF4);
        fGasMaterials["ArCF4"] = m;
    }

    // =====================================================
    // Gas mixture 2: He/Ethane 96.5/3.5 vol%
    // Low-Z, low-density -- good for minimising gamma sensitivity
    // Penning transfer from He metastables into ethane lowers W
    // =====================================================
    {
        G4double fHe=0.965, fEth=0.035;
        G4double rho = fHe*0.1786e-3 + fEth*1.356e-3;
        G4Material* m = new G4Material("HeEth", rho*g/cm3, 2,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(pureHe, fHe);
        m->AddMaterial(ethane,  fEth);
        fGasMaterials["HeEth"] = m;
    }

    // =====================================================
    // Gas mixture 3: Ar/CO2 70/30 vol%
    // Classic stable Micromegas mixture
    // =====================================================
    {
        G4double fAr=0.70, fCO2=0.30;
        G4double rho = fAr*1.782e-3 + fCO2*1.977e-3;
        G4Material* m = new G4Material("ArCO2", rho*g/cm3, 2,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr, fAr);
        m->AddMaterial(CO2,   fCO2);
        fGasMaterials["ArCO2"] = m;
    }

    // =====================================================
    // Gas mixture 4: Ar/CF4/Isobutane 88/10/2 vol%
    // Fast drift + good quenching
    // =====================================================
    {
        G4double fAr=0.88, fCF4=0.10, fIso=0.02;
        G4double rho = fAr*1.782e-3 + fCF4*3.72e-3 + fIso*2.67e-3;
        // 3-component mixture
        G4Material* m = new G4Material("ArCF4Iso", rho*g/cm3, 3,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr,     fAr);
        m->AddMaterial(CF4,       fCF4);
        m->AddMaterial(isobutane, fIso);
        fGasMaterials["ArCF4Iso"] = m;
    }

    // =====================================================
    // Gas mixture 5: Ne/Isobutane 95/5 vol%
    // Moderate Z, good Penning, potential gamma-flash middle ground
    // =====================================================
    {
        G4double fNe=0.95, fIso=0.05;
        G4double rho = fNe*0.8999e-3 + fIso*2.67e-3;
        G4Material* m = new G4Material("NeIso", rho*g/cm3, 2,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(pureNe,    fNe);
        m->AddMaterial(isobutane, fIso);
        fGasMaterials["NeIso"] = m;
    }
}

// ============================================================
G4Material* DetectorConstruction::GetGasMixture(const std::string& name) {
    auto it = fGasMaterials.find(name);
    if (it == fGasMaterials.end()) {
        throw std::runtime_error("Unknown gas mixture: " + name +
            "\nAvailable: ArIso, HeEth, NeCO2, ArCO2, PurAr");
    }
    return it->second;
}

// ============================================================
G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();

    G4NistManager* nist = G4NistManager::Instance();

    // ---- Materials ----
    G4Material* matAir     = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matMylar   = nist->FindOrBuildMaterial("G4_MYLAR");
    G4Material* matAl      = nist->FindOrBuildMaterial("G4_Al");
    G4Material* matKapton  = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material* matCu      = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* matSteel   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* matGas     = GetGasMixture(fConfig.gas);

    // ---- Geometry parameters ----
    // Transverse size: actual 40x40 cm detector
    G4double detXY = 40.0 * cm;

    // Layer thicknesses (along z = beam direction)
    G4double tMylar   = 40.0  * um;   // gas window foil
    G4double tAlWin   = 0.1   * um;   // Al coating on gas window
    G4double tKapCath = 50.0  * um;   // drift cathode kapton
    G4double tCuCath  = 9.0   * um;   // drift cathode copper
    G4double tDrift   = 3.0   * cm;   // drift gap (gas)
    // Micromesh: woven SS, 18 um wire diameter, 30 um total thickness after flattening
    // We model as a uniform thin SS sheet (effective thickness 30 um)
    // Optical transparency ~51% -> effective density reduction handled via thickness only
    G4double tMesh    = 30.0  * um;
    G4double tAmp     = 150.0 * um;   // amplification gap (gas)
    G4double tKapAno  = 50.0  * um;   // anode PCB kapton
    G4double tCuAno   = 9.0   * um;   // anode copper

    // World half-sizes
    G4double totalZ = (tMylar + tAlWin + tKapCath + tCuCath +
                       tDrift + tMesh + tAmp +
                       tKapAno + tCuAno);
    G4double worldZ = totalZ + 5.0 * cm;  // air gap on either side

    // ---- World volume ----
    G4Box* worldSolid = new G4Box("World", detXY/2 + 2*cm, detXY/2 + 2*cm, worldZ/2);
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
    G4VPhysicalVolume* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(),
                                                    worldLV, "World", nullptr, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ---- Helper: place a flat slab at z_center ----
    // All layers stacked along z.  We track the front face z as we go.
    // Convention: beam enters at -z, layers stacked in +z direction.
    G4double zFront = -totalZ / 2.0;

    auto PlaceSlab = [&](const std::string& name, G4double thickness,
                          G4Material* mat, G4VisAttributes* vis,
                          G4LogicalVolume*& outLV) -> G4double {
        G4double zCenter = zFront + thickness / 2.0;
        G4Box* solid = new G4Box(name, detXY/2, detXY/2, thickness/2);
        outLV = new G4LogicalVolume(solid, mat, name);
        if (vis) outLV->SetVisAttributes(vis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCenter),
                          outLV, name, worldLV, false, 0, true);
        zFront += thickness;
        return zCenter;
    };

    // Vis attributes
    auto visMylar  = new G4VisAttributes(G4Color(0.7, 0.9, 0.7, 0.5));
    auto visAl     = new G4VisAttributes(G4Color(0.7, 0.7, 0.7, 0.8));
    auto visKapton = new G4VisAttributes(G4Color(0.9, 0.7, 0.0, 0.7));
    auto visCu     = new G4VisAttributes(G4Color(0.8, 0.4, 0.1, 0.8));
    auto visDrift  = new G4VisAttributes(G4Color(0.2, 0.5, 1.0, 0.3));
    auto visMesh   = new G4VisAttributes(G4Color(0.5, 0.5, 0.5, 0.9));
    auto visAmp    = new G4VisAttributes(G4Color(1.0, 0.3, 0.3, 0.3));

    // Dummy LV for layers we don't care about scoring
    G4LogicalVolume* dummyLV = nullptr;

    // === Place layers ===

    // 1. Gas window: Mylar 40 um
    PlaceSlab("GasWindow_Mylar", tMylar, matMylar, visMylar, dummyLV);

    // 2. Gas window: Al 0.1 um
    PlaceSlab("GasWindow_Al", tAlWin, matAl, visAl, dummyLV);

    // 3. Drift cathode: Kapton 50 um
    PlaceSlab("DriftCathode_Kapton", tKapCath, matKapton, visKapton, dummyLV);

    // 4. Drift cathode: Cu 9 um
    PlaceSlab("DriftCathode_Cu", tCuCath, matCu, visCu, dummyLV);

    // 5. Drift gas volume: 3 cm  <-- PRIMARY SCORING VOLUME
    PlaceSlab("DriftGas", tDrift, matGas, visDrift, fDriftGasLV);

    // 6. Micromesh: SS 30 um
    PlaceSlab("Micromesh", tMesh, matSteel, visMesh, dummyLV);

    // 7. Amplification gas: 150 um  <-- SECONDARY SCORING VOLUME
    PlaceSlab("AmpGas", tAmp, matGas, visAmp, fAmpGasLV);

    // 8. Anode: Kapton 50 um
    PlaceSlab("Anode_Kapton", tKapAno, matKapton, visKapton, dummyLV);

    // 9. Anode: Cu 9 um
    PlaceSlab("Anode_Cu", tCuAno, matCu, visCu, dummyLV);

    G4cout << "\n=== Detector geometry built ===" << G4endl;
    G4cout << "  Gas mixture   : " << fConfig.gas << " (rho = "
           << matGas->GetDensity() / (mg/cm3) << " mg/cm3)" << G4endl;
    G4cout << "  Drift gap     : " << tDrift / cm << " cm" << G4endl;
    G4cout << "  Amp gap       : " << tAmp   / um << " um" << G4endl;
    G4cout << "  Total det. Z  : " << totalZ / mm << " mm" << G4endl;
    G4cout << "================================\n" << G4endl;

    return worldPV;
}

// ============================================================
void DetectorConstruction::ConstructSDandField() {
    // Sensitive detector for drift gas
    SensitiveDetector* driftSD = new SensitiveDetector("DriftGasSD", "DriftGasHits",
                                                         "DriftGas", fConfig);
    G4SDManager::GetSDMpointer()->AddNewDetector(driftSD);
    SetSensitiveDetector(fDriftGasLV, driftSD);

    // Sensitive detector for amp gas
    SensitiveDetector* ampSD = new SensitiveDetector("AmpGasSD", "AmpGasHits",
                                                       "AmpGas", fConfig);
    G4SDManager::GetSDMpointer()->AddNewDetector(ampSD);
    SetSensitiveDetector(fAmpGasLV, ampSD);

    // Set step size limit in gas volumes to improve ionization tracking
    // (allow Geant4 to create steps fine enough to resolve primary ionization clusters)
    G4UserLimits* stepLimit = new G4UserLimits(100 * um);  // max step in gas
    fDriftGasLV->SetUserLimits(stepLimit);
    fAmpGasLV->SetUserLimits(stepLimit);
}
