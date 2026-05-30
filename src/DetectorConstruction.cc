// DetectorConstruction.cc
//
// Modes:
//  kVacuum          : Micromegas in vacuum, optional Al shielding upstream.
//  kFullExperiment  : He-3 target → air → MM → PCB → air → scint wall → air → LS stack.
//  kSr90Calibration : Sr-90 source in air → MM → PCB → air → scint wall → air → LS stack.
//  kSr90NoMM        : Sr-90 source in air → scint wall → air → LS stack (no MM/PCB).
//  kLSCalib         : Sr-90 source capsule → air → 1 LS layer → air → back scint bar.
//
// Geometry updated to match Full_Geant (4-arm X17 sim):
//  - He-3 target: r=1.5 cm, L=8 cm (was r=2.5 cm, L=15 cm)
//  - LS: 2 layers × 2 cm LAB, each preceded by inner CFRP liner + Al liner
//  - Scint wall: BlackMylar tape (200 µm) instead of PVC tape (165 µm)

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"

#include "SensitiveDetector.hh"

#include <stdexcept>
#include <algorithm>

// ============================================================
DetectorConstruction::DetectorConstruction(const SimConfig& cfg)
    : G4VUserDetectorConstruction(), fConfig(cfg) {}

// ============================================================
void DetectorConstruction::DefineMaterials() {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elN  = nist->FindOrBuildElement("N");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elSi = nist->FindOrBuildElement("Si");
    G4Element* elAr = nist->FindOrBuildElement("Ar");
    G4Element* elNe = nist->FindOrBuildElement("Ne");
    G4Element* elHe = nist->FindOrBuildElement("He");
    G4Element* elF  = nist->FindOrBuildElement("F");

    G4Material* isobutane = new G4Material("Isobutane", 2.67e-3*g/cm3, 2,
                                            kStateGas, 293.15*kelvin, 1*atmosphere);
    isobutane->AddElement(elC, 4);
    isobutane->AddElement(elH, 10);

    G4Material* ethane = new G4Material("Ethane", 1.356e-3*g/cm3, 2,
                                         kStateGas, 293.15*kelvin, 1*atmosphere);
    ethane->AddElement(elC, 2);
    ethane->AddElement(elH, 6);

    G4Material* CO2 = new G4Material("CO2_gas", 1.977e-3*g/cm3, 2,
                                      kStateGas, 293.15*kelvin, 1*atmosphere);
    CO2->AddElement(elC, 1);
    CO2->AddElement(elO, 2);

    G4Material* CF4 = new G4Material("CF4_gas", 3.72e-3*g/cm3, 2,
                                      kStateGas, 293.15*kelvin, 1*atmosphere);
    CF4->AddElement(elC, 1);
    CF4->AddElement(elF, 4);

    G4Material* purAr = new G4Material("PureArgon", 1.782e-3*g/cm3, 1,
                                        kStateGas, 293.15*kelvin, 1*atmosphere);
    purAr->AddElement(elAr, 1);

    G4Material* pureHe = new G4Material("PureHe", 0.1786e-3*g/cm3, 1,
                                         kStateGas, 293.15*kelvin, 1*atmosphere);
    pureHe->AddElement(elHe, 1);

    G4Material* pureNe = new G4Material("PureNe", 0.8999e-3*g/cm3, 1,
                                         kStateGas, 293.15*kelvin, 1*atmosphere);
    pureNe->AddElement(elNe, 1);

    // ── Gas mixtures ─────────────────────────────────────────
    auto makeMix2 = [&](const char* nm, G4double rho,
                        G4Material* m1, G4double f1,
                        G4Material* m2, G4double f2) -> G4Material* {
        auto* m = new G4Material(nm, rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(m1, f1); m->AddMaterial(m2, f2); return m;
    };
    auto makeMix3 = [&](const char* nm, G4double rho,
                        G4Material* m1, G4double f1,
                        G4Material* m2, G4double f2,
                        G4Material* m3, G4double f3) -> G4Material* {
        auto* m = new G4Material(nm, rho*g/cm3, 3, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(m1, f1); m->AddMaterial(m2, f2); m->AddMaterial(m3, f3); return m;
    };

    fGasMaterials["ArCF4"]    = makeMix2("ArCF4",    0.90*1.782e-3+0.10*3.72e-3,
                                           purAr,0.90, CF4,0.10);
    fGasMaterials["HeEth"]    = makeMix2("HeEth",    0.965*0.1786e-3+0.035*1.356e-3,
                                           pureHe,0.965, ethane,0.035);
    fGasMaterials["ArCO2"]    = makeMix2("ArCO2",    0.70*1.782e-3+0.30*1.977e-3,
                                           purAr,0.70, CO2,0.30);
    fGasMaterials["ArIso"]    = makeMix2("ArIso",    0.95*1.782e-3+0.05*2.67e-3,
                                           purAr,0.95, isobutane,0.05);
    fGasMaterials["NeIso"]    = makeMix2("NeIso",    0.95*0.8999e-3+0.05*2.67e-3,
                                           pureNe,0.95, isobutane,0.05);
    fGasMaterials["NeCF4"]    = makeMix2("NeCF4",    0.90*0.8999e-3+0.10*3.72e-3,
                                           pureNe,0.90, CF4,0.10);
    fGasMaterials["ArCF4Iso"] = makeMix3("ArCF4Iso", 0.88*1.782e-3+0.10*3.72e-3+0.02*2.67e-3,
                                           purAr,0.88, CF4,0.10, isobutane,0.02);
    fGasMaterials["ArCF4CO2"] = makeMix3("ArCF4CO2", 0.45*1.782e-3+0.40*3.72e-3+0.15*1.977e-3,
                                           purAr,0.45, CF4,0.40, CO2,0.15);

    {
        auto* m = new G4Material("PureCF4", 3.72e-3*g/cm3, 2,
                                  kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddElement(elC,1); m->AddElement(elF,4);
        fGasMaterials["PureCF4"] = m;
    }
    fGasMaterials["PureAr"]     = purAr;
    fGasMaterials["PureHe"]     = pureHe;
    fGasMaterials["PureNe"]     = pureNe;
    fGasMaterials["PureEthane"] = ethane;
    fGasMaterials["PureIso"]    = isobutane;
    fGasMaterials["PureCO2"]    = CO2;

    // ── He-3 at 300 bar ──────────────────────────────────────
    {
        auto* isoHe3 = new G4Isotope("He3_iso", 2, 3, 3.0160293*g/mole);
        auto* elHe3  = new G4Element("Helium3_elem", "3He", 1);
        elHe3->AddIsotope(isoHe3, 1.0);
        auto* m = new G4Material("He3Gas_300bar", 37.6e-3*g/cm3, 1,
                                  kStateGas, 293.15*kelvin, 300*atmosphere);
        m->AddElement(elHe3, 1);
        fGasMaterials["He3Gas_300bar"] = m;
    }

    // ── Structural / detector materials ─────────────────────
    {
        auto* m = new G4Material("CFRP", 1.55*g/cm3, 3);
        m->AddElement(elC, 0.8968); m->AddElement(elH, 0.0207); m->AddElement(elO, 0.0826);
        fGasMaterials["CFRP"] = m;
    }
    {
        auto* m = new G4Material("ResistivePaste", 1.4*g/cm3, 3);
        m->AddElement(elC, 0.65); m->AddElement(elH, 0.08); m->AddElement(elO, 0.27);
        fGasMaterials["ResistivePaste"] = m;
    }
    {
        auto* m = new G4Material("FR4", 1.85*g/cm3, 4);
        m->AddElement(elSi, 0.2805); m->AddElement(elO,  0.4195);
        m->AddElement(elC,  0.2750); m->AddElement(elH,  0.0250);
        fGasMaterials["FR4"] = m;
    }
    {
        auto* m = new G4Material("Rohacell51", 0.052*g/cm3, 4);
        m->AddElement(elC, 0.5783); m->AddElement(elH, 0.0602);
        m->AddElement(elN, 0.1687); m->AddElement(elO, 0.1928);
        fGasMaterials["Rohacell51"] = m;
    }
    {
        auto* m = new G4Material("LAB_LiqScint", 0.86*g/cm3, 2);
        m->AddElement(elC, 0.8780); m->AddElement(elH, 0.1220);
        fGasMaterials["LAB_LiqScint"] = m;
    }
    // BlackMylar: used for scint wall wrapping and back scint tape (G4_MYLAR = PET).
    fGasMaterials["BlackMylar"] = nist->FindOrBuildMaterial("G4_MYLAR");
}

// ============================================================
G4Material* DetectorConstruction::GetGasMixture(const std::string& name) {
    auto it = fGasMaterials.find(name);
    if (it == fGasMaterials.end())
        throw std::runtime_error("Unknown gas/material: " + name);
    return it->second;
}

// ============================================================
G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* matAir     = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matMylar   = nist->FindOrBuildMaterial("G4_MYLAR");
    G4Material* matAl      = nist->FindOrBuildMaterial("G4_Al");
    G4Material* matKapton  = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material* matCu      = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* matSteel   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* matGas     = GetGasMixture(fConfig.gas);
    G4Material* matHe3     = fGasMaterials.at("He3Gas_300bar");
    G4Material* matCFRP    = fGasMaterials.at("CFRP");
    G4Material* matResPaste= fGasMaterials.at("ResistivePaste");
    G4Material* matFR4     = fGasMaterials.at("FR4");
    G4Material* matRohacell= fGasMaterials.at("Rohacell51");
    G4Material* matLAB     = fGasMaterials.at("LAB_LiqScint");
    G4Material* matBlkMylar= fGasMaterials.at("BlackMylar");
    G4Material* matPlScint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    G4double detXY = 40.0 * cm;

    // ── MM stack layer thicknesses ───────────────────────────
    G4double tMylar    = 40.0  * um;
    G4double tAlWin    = 0.1   * um;
    G4double tKapCath  = 50.0  * um;
    G4double tCuCath   = 9.0   * um;
    G4double tDrift    = 3.0   * cm;
    G4double tMesh     = 30.0  * um;
    G4double tAmp      = 150.0 * um;
    G4double tResPaste = 100.0 * um;
    G4double mmTotalZ  = tMylar + tAlWin + tKapCath + tCuCath
                       + tDrift + tMesh + tAmp + tResPaste;

    // ── PCB stack ─────────────────────────────────────────────
    G4double tPCB_Kap  = 50.0  * um;
    G4double tPCB_Cu   = 26.0  * um;
    G4double tPCB_FR4  = 100.0 * um;
    G4double tPCB_Roh  = 5.0   * mm;
    G4double tPCB_Al   = 50.0  * um;
    G4double pcbTotalZ = tPCB_Kap + 4*(tPCB_Cu + tPCB_FR4) + tPCB_Roh + tPCB_Al;

    // ── Scint wall (BlackMylar tape, from Full_Geant) ─────────
    G4double tBlkTape  = 200.0 * um;
    G4double tPlScint  = 3.0   * mm;
    G4double tScAl     = 50.0  * um;
    G4double scintWallZ = 2*tBlkTape + tPlScint + tScAl;

    // ── LS stack (from Full_Geant: 2 layers × 2cm with inner liners) ─────
    G4double tLSCfrp      = fConfig.cfrpThickness_mm  * mm;   // structural CFRP wall
    G4double tLSInnerCfrp = fConfig.ls_inner_cfrp_um  * um;   // inner CFRP liner
    G4double tLSInnerAl   = fConfig.ls_inner_al_um    * um;   // Al liner
    G4double tLS          = fConfig.ls_thick_cm        * cm;   // LAB layer
    // 3 CFRP walls + 2 inner CFRP liners + 2 Al liners + 2 LAB layers
    G4double lsStackZ = 3*tLSCfrp + 2*(tLSInnerCfrp + tLSInnerAl + tLS);

    // ── Vis attributes ───────────────────────────────────────
    auto visMylar     = new G4VisAttributes(G4Color(0.7, 0.9, 0.7, 0.5));
    auto visAl        = new G4VisAttributes(G4Color(0.7, 0.7, 0.7, 0.8));
    auto visKapton    = new G4VisAttributes(G4Color(0.9, 0.7, 0.0, 0.7));
    auto visCu        = new G4VisAttributes(G4Color(0.8, 0.4, 0.1, 0.8));
    auto visDrift     = new G4VisAttributes(G4Color(0.2, 0.5, 1.0, 0.3));
    auto visMesh      = new G4VisAttributes(G4Color(0.5, 0.5, 0.5, 0.9));
    auto visAmp       = new G4VisAttributes(G4Color(1.0, 0.3, 0.3, 0.3));
    auto visResPaste  = new G4VisAttributes(G4Color(0.2, 0.2, 0.2, 0.8));
    auto visHe3       = new G4VisAttributes(G4Color(0.6, 0.9, 1.0, 0.4));
    auto visCFRP      = new G4VisAttributes(G4Color(0.15, 0.15, 0.15, 0.9));
    auto visFR4       = new G4VisAttributes(G4Color(0.2, 0.6, 0.2, 0.8));
    auto visRohacell  = new G4VisAttributes(G4Color(0.9, 0.9, 0.6, 0.5));
    auto visScint     = new G4VisAttributes(G4Color(0.9, 0.9, 0.2, 0.7));
    auto visLAB       = new G4VisAttributes(G4Color(0.3, 0.8, 0.9, 0.4));
    auto visBlkMylar  = new G4VisAttributes(G4Color(0.1, 0.1, 0.1, 0.9));

    G4LogicalVolume*   worldLV = nullptr;
    G4VPhysicalVolume* worldPV = nullptr;

    // ── PlaceSlab helper ─────────────────────────────────────
    // Placing slab with its front face at *zFront, centred at zFront+t/2.
    // Increments zFront by t. hx,hy may differ from detXY/2 for larger detectors.
    auto MakeSlab = [&](const std::string& name, G4double t,
                         G4double hx, G4double hy,
                         G4Material* mat, G4VisAttributes* vis) -> G4LogicalVolume* {
        auto* box = new G4Box(name, hx, hy, t/2);
        auto* lv  = new G4LogicalVolume(box, mat, name);
        if (vis) lv->SetVisAttributes(vis);
        return lv;
    };

    // ═══════════════════════════════════════════════════════════
    // VACUUM MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kVacuum) {

        G4double alThickness = fConfig.alThickness_mm * mm;
        G4double alGap       = 2.0 * cm;

        G4double gunZ = 10.0 * cm;
        G4double minUpstream = 2.0*gunZ - mmTotalZ - 2.5*cm;
        G4double upstreamMargin = std::max(alGap + alThickness + 0.5*cm, minUpstream);
        G4double worldZ = mmTotalZ + upstreamMargin + 2.5*cm;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, worldZ/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        if (alThickness > 0) {
            G4double alZCenter = -mmTotalZ/2.0 - alGap - alThickness/2.0;
            auto* alBox = new G4Box("AlShield", detXY/2, detXY/2, alThickness/2);
            auto* alLV  = new G4LogicalVolume(alBox, matAl, "AlShield");
            auto* visAlS = new G4VisAttributes(G4Color(0.75, 0.75, 0.75, 0.9));
            visAlS->SetForceSolid(true);
            alLV->SetVisAttributes(visAlS);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,alZCenter),
                              alLV, "AlShield", worldLV, false, 0, true);
        }

        G4cout << "\n=== Vacuum mode ===" << G4endl;
        G4cout << "  Gas: " << fConfig.gas
               << "  (rho=" << matGas->GetDensity()/(mg/cm3) << " mg/cm3)" << G4endl;

        // Place MM layers
        G4double zFront = -mmTotalZ / 2.0;
        auto PlaceSlab = [&](const std::string& name, G4double t,
                              G4Material* mat, G4VisAttributes* vis,
                              G4LogicalVolume*& outLV) {
            auto* lv = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zFront+t/2), lv, name, worldLV, false, 0, true);
            zFront += t;
            outLV = lv;
        };
        G4LogicalVolume* dummy = nullptr;
        PlaceSlab("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,    dummy);
        PlaceSlab("GasWindow_Al",        tAlWin,    matAl,       visAl,       dummy);
        PlaceSlab("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,   dummy);
        PlaceSlab("DriftCathode_Cu",     tCuCath,   matCu,       visCu,       dummy);
        PlaceSlab("DriftGas",            tDrift,    matGas,      visDrift,    fDriftGasLV);
        PlaceSlab("Micromesh",           tMesh,     matSteel,    visMesh,     dummy);
        PlaceSlab("AmpGas",              tAmp,      matGas,      visAmp,      fAmpGasLV);
        PlaceSlab("ResistivePaste",      tResPaste, matResPaste, visResPaste, dummy);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // FULL EXPERIMENT MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kFullExperiment) {

        // He-3 capsule (from Full_Geant: r=1.5 cm, L=8 cm)
        G4double he3R          = 1.5  * cm;
        G4double he3HalfL      = 4.0  * cm;
        G4double alWallT       = 0.5  * mm;
        G4double cfrpWallT     = 0.9  * mm;
        G4double alR           = he3R + alWallT;
        G4double cfrpR         = alR  + cfrpWallT;
        G4double alHalfL       = he3HalfL + alWallT;
        G4double cfrpHalfL     = alHalfL  + cfrpWallT;
        G4double capsuleZExtent = 2.0 * cfrpR;

        G4double airGap1   = 200.0 * mm;
        G4double airGap2   = 20.0  * mm;
        G4double airGap3   = 20.0  * mm;

        G4double totalFullZ = capsuleZExtent + airGap1 + mmTotalZ + pcbTotalZ
                            + airGap2 + scintWallZ + airGap3 + lsStackZ;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, (totalFullZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4double capsuleZCenter = -totalFullZ/2.0 + cfrpR;
        fHe3GasCenterZ = capsuleZCenter;

        auto* capRot = new G4RotationMatrix();
        capRot->rotateX(-90.*deg);

        auto* cfrpSolid = new G4Tubs("He3Capsule_CFRP", 0, cfrpR, cfrpHalfL, 0, 360.*deg);
        auto* cfrpLV    = new G4LogicalVolume(cfrpSolid, matCFRP, "He3Capsule_CFRP");
        cfrpLV->SetVisAttributes(visCFRP);
        new G4PVPlacement(capRot, G4ThreeVector(0,0,capsuleZCenter), cfrpLV, "He3Capsule_CFRP", worldLV, false, 0, true);

        auto* alSolid = new G4Tubs("He3Capsule_Al", 0, alR, alHalfL, 0, 360.*deg);
        auto* alLV    = new G4LogicalVolume(alSolid, matAl, "He3Capsule_Al");
        alLV->SetVisAttributes(visAl);
        new G4PVPlacement(nullptr, G4ThreeVector(), alLV, "He3Capsule_Al", cfrpLV, false, 0, true);

        auto* he3Solid = new G4Tubs("He3Gas", 0, he3R, he3HalfL, 0, 360.*deg);
        fHe3GasLV = new G4LogicalVolume(he3Solid, matHe3, "He3Gas");
        fHe3GasLV->SetVisAttributes(visHe3);
        new G4PVPlacement(nullptr, G4ThreeVector(), fHe3GasLV, "He3Gas", alLV, false, 0, true);

        G4cout << "\n=== Full-experiment geometry ===" << G4endl;
        G4cout << "  He-3: r=" << he3R/cm << " cm, L=" << 2*he3HalfL/cm << " cm, 300 bar" << G4endl;
        G4cout << "  LS stack (2×" << fConfig.ls_thick_cm << " cm LAB, "
               << fConfig.cfrpThickness_mm << " mm CFRP walls)" << G4endl;
        G4cout << "  Total Z: " << totalFullZ/mm << " mm" << G4endl;

        G4double zF = -totalFullZ/2.0 + capsuleZExtent;
        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",             airGap1,   matAir,      nullptr,      dL);
        Place("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,     dL);
        Place("GasWindow_Al",        tAlWin,    matAl,       visAl,        dL);
        Place("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,    dL);
        Place("DriftCathode_Cu",     tCuCath,   matCu,       visCu,        dL);
        Place("DriftGas",            tDrift,    matGas,      visDrift,     fDriftGasLV);
        Place("Micromesh",           tMesh,     matSteel,    visMesh,      dL);
        Place("AmpGas",              tAmp,      matGas,      visAmp,       fAmpGasLV);
        Place("ResistivePaste",      tResPaste, matResPaste, visResPaste,  dL);

        Place("PCB_Kapton",   tPCB_Kap, matKapton,  visKapton,  dL);
        Place("PCB_Cu_1",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_1",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_2",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_2",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_3",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_3",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_4",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_4",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Rohacell", tPCB_Roh, matRohacell,visRohacell,dL);
        Place("PCB_AlFoil",   tPCB_Al,  matAl,      visAl,      dL);

        Place("AirGap2",               airGap2,   matAir,      nullptr,     dL);
        Place("ScintWall_BlackTape1",  tBlkTape,  matBlkMylar, visBlkMylar, dL);
        Place("PlasticScint",          tPlScint,  matPlScint,  visScint,    dL);
        Place("ScintWall_BlackTape2",  tBlkTape,  matBlkMylar, visBlkMylar, dL);
        Place("ScintWall_AlFoil",      tScAl,     matAl,       visAl,       dL);

        Place("AirGap3",          airGap3,    matAir, nullptr, dL);
        Place("LS_CFRP_1",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_1",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_1",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_1",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_2",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_2",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_2",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_2",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_3",        tLSCfrp,    matCFRP, visCFRP, dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // SR-90 CALIBRATION MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kSr90Calibration) {

        G4double airToMM  = 226.5 * mm;
        G4double airGap2  = 20.0  * mm;
        G4double airGap3  = 20.0  * mm;

        G4double totalZ = airToMM + mmTotalZ + pcbTotalZ
                        + airGap2 + scintWallZ + airGap3 + lsStackZ;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, (totalZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        fHe3GasCenterZ = -totalZ / 2.0;

        G4cout << "\n=== Sr-90 calibration geometry ===" << G4endl;
        G4cout << "  Air source-to-MM: " << airToMM/mm << " mm" << G4endl;
        G4cout << "  LS (2×" << fConfig.ls_thick_cm << " cm, CFRP " << fConfig.cfrpThickness_mm << " mm)" << G4endl;
        G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

        G4double zF = -totalZ / 2.0;
        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",             airToMM,   matAir,      nullptr,      dL);
        Place("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,     dL);
        Place("GasWindow_Al",        tAlWin,    matAl,       visAl,        dL);
        Place("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,    dL);
        Place("DriftCathode_Cu",     tCuCath,   matCu,       visCu,        dL);
        Place("DriftGas",            tDrift,    matGas,      visDrift,     fDriftGasLV);
        Place("Micromesh",           tMesh,     matSteel,    visMesh,      dL);
        Place("AmpGas",              tAmp,      matGas,      visAmp,       fAmpGasLV);
        Place("ResistivePaste",      tResPaste, matResPaste, visResPaste,  dL);

        Place("PCB_Kapton",   tPCB_Kap, matKapton,  visKapton,  dL);
        Place("PCB_Cu_1",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_1",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_2",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_2",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_3",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_3",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_4",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_4",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Rohacell", tPCB_Roh, matRohacell,visRohacell,dL);
        Place("PCB_AlFoil",   tPCB_Al,  matAl,      visAl,      dL);

        Place("AirGap2",               airGap2,  matAir,      nullptr,     dL);
        Place("ScintWall_BlackTape1",  tBlkTape, matBlkMylar, visBlkMylar, dL);
        Place("PlasticScint",          tPlScint, matPlScint,  visScint,    dL);
        Place("ScintWall_BlackTape2",  tBlkTape, matBlkMylar, visBlkMylar, dL);
        Place("ScintWall_AlFoil",      tScAl,    matAl,       visAl,       dL);

        Place("AirGap3",          airGap3,    matAir, nullptr, dL);
        Place("LS_CFRP_1",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_1",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_1",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_1",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_2",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_2",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_2",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_2",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_3",        tLSCfrp,    matCFRP, visCFRP, dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // SR-90 NO-MM MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kSr90NoMM) {

        G4double airToScint = 226.5*mm + mmTotalZ + pcbTotalZ + 20.0*mm;
        G4double airGap3    = 20.0 * mm;

        G4double totalZ = airToScint + scintWallZ + airGap3 + lsStackZ;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, (totalZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        fHe3GasCenterZ = -totalZ / 2.0;

        G4cout << "\n=== Sr-90 no-MM geometry ===" << G4endl;
        G4cout << "  Air source-to-scint: " << airToScint/mm << " mm" << G4endl;
        G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

        G4double zF = -totalZ / 2.0;
        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",              airToScint, matAir,      nullptr,     dL);
        Place("ScintWall_BlackTape1", tBlkTape,   matBlkMylar, visBlkMylar, dL);
        Place("PlasticScint",         tPlScint,   matPlScint,  visScint,    dL);
        Place("ScintWall_BlackTape2", tBlkTape,   matBlkMylar, visBlkMylar, dL);
        Place("ScintWall_AlFoil",     tScAl,      matAl,       visAl,       dL);

        Place("AirGap3",          airGap3,    matAir, nullptr, dL);
        Place("LS_CFRP_1",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_1",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_1",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_1",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_2",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_2",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_2",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_2",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_3",        tLSCfrp,    matCFRP, visCFRP, dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // LS CALIBRATION MODE
    // Bare source gun → air gap → 1 LS layer only.
    // No source capsule, no back scint.
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kLSCalib) {

        G4double airToLS  = fConfig.source_to_det_mm * mm;
        // LS: CFRP_front | InnerCFRP | Al | LAB | InnerCFRP | Al | CFRP_back
        G4double lsZ    = 2*tLSCfrp + 2*(tLSInnerCfrp + tLSInnerAl) + tLS;
        G4double totalZ = airToLS + lsZ;

        G4double lsHX = 22.5*cm;  // 45×45 cm LS face
        G4double lsHY = 22.5*cm;

        auto* worldSolid = new G4Box("World", lsHX+2*cm, lsHY+2*cm, (totalZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        fHe3GasCenterZ = -totalZ / 2.0;  // gun at front of world

        G4cout << "\n=== LS Calibration geometry ===" << G4endl;
        G4cout << "  Source-to-LS air gap: " << airToLS/mm << " mm" << G4endl;
        G4cout << "  LS: " << tLS/cm << " cm LAB,  CFRP walls " << fConfig.cfrpThickness_mm << " mm" << G4endl;
        G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

        G4double zF = -totalZ / 2.0;

        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, lsHX, lsHY, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",       airToLS,      matAir,  nullptr,  dL);
        Place("LS_CFRP_1",     tLSCfrp,      matCFRP, visCFRP,  dL);
        Place("LS_InnerCFRP_1",tLSInnerCfrp, matCFRP, visCFRP,  dL);
        Place("LS_Al_1",       tLSInnerAl,   matAl,   visAl,    dL);
        Place("LiqScint_1",    tLS,          matLAB,  visLAB,   dL);
        Place("LS_InnerCFRP_2",tLSInnerCfrp, matCFRP, visCFRP,  dL);
        Place("LS_Al_2",       tLSInnerAl,   matAl,   visAl,    dL);
        Place("LS_CFRP_2",     tLSCfrp,      matCFRP, visCFRP,  dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // BACK SCINT CALIBRATION MODE
    // Bare source gun → air gap → 1 back scint bar only.
    // No source capsule, no LS layer.
    // ═══════════════════════════════════════════════════════════
    // else kBackScintCalib

    G4double airToBSc = fConfig.source_to_det_mm * mm;
    G4double tBscTape = fConfig.backscint_tape_um  * um;
    G4double tBscAl   = fConfig.backscint_al_um    * um;
    G4double tBscPVT  = fConfig.backscint_thick_cm * cm;
    G4double bscZ     = 2*tBscTape + 2*tBscAl + tBscPVT;
    G4double totalZ   = airToBSc + bscZ;

    G4double bsHX = fConfig.backscint_u_cm/2 * cm;  // 15 cm half-width (30 cm total)
    G4double bsHY = fConfig.backscint_v_cm/2 * cm;  // 10 cm half-height (20 cm total)

    auto* worldSolid = new G4Box("World", bsHX+2*cm, bsHY+2*cm, (totalZ+2*cm)/2);
    worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
    worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    fHe3GasCenterZ = -totalZ / 2.0;

    G4cout << "\n=== Back Scint Calibration geometry ===" << G4endl;
    G4cout << "  Source-to-scint air gap: " << airToBSc/mm << " mm" << G4endl;
    G4cout << "  Back scint: " << tBscPVT/cm << " cm PVT, "
           << fConfig.backscint_u_cm << "×" << fConfig.backscint_v_cm << " cm face" << G4endl;
    G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

    G4double zF = -totalZ / 2.0;

    auto visBscPVT = new G4VisAttributes(G4Color(0.9, 0.5, 0.1, 0.8));

    auto PlaceB = [&](const std::string& name, G4double t,
                       G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
        out = MakeSlab(name, t, bsHX, bsHY, mat, vis);
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
        zF += t;
    };
    G4LogicalVolume* dL = nullptr;

    {
        // Air gap: use world-sized slab so gun position is inside air
        auto* airBox = new G4Box("AirGap1", bsHX+2*cm, bsHY+2*cm, airToBSc/2);
        auto* airLV  = new G4LogicalVolume(airBox, matAir, "AirGap1");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+airToBSc/2), airLV, "AirGap1", worldLV, false, 0, true);
        zF += airToBSc;
    }

    PlaceB("BackScintWrap_Tape1", tBscTape, matBlkMylar, visBlkMylar, dL);
    PlaceB("BackScintWrap_Al1",   tBscAl,   matAl,       visAl,       dL);
    PlaceB("BackScint",           tBscPVT,  matPlScint,  visBscPVT,   fBackScintLV);
    PlaceB("BackScintWrap_Al2",   tBscAl,   matAl,       visAl,       dL);
    PlaceB("BackScintWrap_Tape2", tBscTape, matBlkMylar, visBlkMylar, dL);

    return worldPV;
}

// ============================================================
void DetectorConstruction::ConstructSDandField() {
    if (fDriftGasLV) {
        auto* sd = new SensitiveDetector("DriftGasSD", "DriftGasHits", "DriftGas", fConfig);
        G4SDManager::GetSDMpointer()->AddNewDetector(sd);
        SetSensitiveDetector(fDriftGasLV, sd);
        fDriftGasLV->SetUserLimits(new G4UserLimits(100*um));
    }
    if (fAmpGasLV) {
        auto* sd = new SensitiveDetector("AmpGasSD", "AmpGasHits", "AmpGas", fConfig);
        G4SDManager::GetSDMpointer()->AddNewDetector(sd);
        SetSensitiveDetector(fAmpGasLV, sd);
        fAmpGasLV->SetUserLimits(new G4UserLimits(100*um));
    }
    if (fHe3GasLV) {
        fHe3GasLV->SetUserLimits(new G4UserLimits(1.0*mm));
    }
}
