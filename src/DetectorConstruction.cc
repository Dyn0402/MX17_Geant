// DetectorConstruction.cc
//
// Vacuum mode  : Micromegas detector in vacuum, optional Al shielding upstream.
// Full-exp mode: He-3 target → air → Micromegas → PCB stack → air →
//                scintillator wall → air → 4-layer liquid scintillator stack.

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

    // ── Elements ────────────────────────────────────────────
    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elN  = nist->FindOrBuildElement("N");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elSi = nist->FindOrBuildElement("Si");
    G4Element* elCl = nist->FindOrBuildElement("Cl");
    G4Element* elAr = nist->FindOrBuildElement("Ar");
    G4Element* elNe = nist->FindOrBuildElement("Ne");
    G4Element* elHe = nist->FindOrBuildElement("He");
    G4Element* elF  = nist->FindOrBuildElement("F");

    // ── Detector gas building-block materials ────────────────
    // (identical to original)

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
    {
        G4double fAr=0.90, fCF4=0.10;
        G4double rho = fAr*1.782e-3 + fCF4*3.72e-3;
        G4Material* m = new G4Material("ArCF4", rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr, fAr); m->AddMaterial(CF4, fCF4);
        fGasMaterials["ArCF4"] = m;
    }
    {
        G4double fHe=0.965, fEth=0.035;
        G4double rho = fHe*0.1786e-3 + fEth*1.356e-3;
        G4Material* m = new G4Material("HeEth", rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(pureHe, fHe); m->AddMaterial(ethane, fEth);
        fGasMaterials["HeEth"] = m;
    }
    {
        G4double fAr=0.70, fCO2=0.30;
        G4double rho = fAr*1.782e-3 + fCO2*1.977e-3;
        G4Material* m = new G4Material("ArCO2", rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr, fAr); m->AddMaterial(CO2, fCO2);
        fGasMaterials["ArCO2"] = m;
    }
    {
        G4double fAr=0.88, fCF4=0.10, fIso=0.02;
        G4double rho = fAr*1.782e-3 + fCF4*3.72e-3 + fIso*2.67e-3;
        G4Material* m = new G4Material("ArCF4Iso", rho*g/cm3, 3, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr, fAr); m->AddMaterial(CF4, fCF4); m->AddMaterial(isobutane, fIso);
        fGasMaterials["ArCF4Iso"] = m;
    }
    {
        G4double fNe=0.95, fIso=0.05;
        G4double rho = fNe*0.8999e-3 + fIso*2.67e-3;
        G4Material* m = new G4Material("NeIso", rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(pureNe, fNe); m->AddMaterial(isobutane, fIso);
        fGasMaterials["NeIso"] = m;
    }
    {
        G4double fNe=0.90, fCF4=0.10;
        G4double rho = fNe*0.8999e-3 + fCF4*3.72e-3;
        G4Material* m = new G4Material("NeCF4", rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(pureNe, fNe); m->AddMaterial(CF4, fCF4);
        fGasMaterials["NeCF4"] = m;
    }
    {
        G4double fAr=0.45, fCF4=0.40, fCO2=0.15;
        G4double rho = fAr*1.782e-3 + fCF4*3.72e-3 + fCO2*1.977e-3;
        G4Material* m = new G4Material("ArCF4CO2", rho*g/cm3, 3, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr, fAr); m->AddMaterial(CF4, fCF4); m->AddMaterial(CO2, fCO2);
        fGasMaterials["ArCF4CO2"] = m;
    }
    {
        G4Material* m = new G4Material("PureCF4", 3.72e-3*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddElement(elC, 1); m->AddElement(elF, 4);
        fGasMaterials["PureCF4"] = m;
    }

    // =====================================================
    // Gas mixture 8b: Ar/Isobutane 95/5 vol%
    // Classic Micromegas mixture: isobutane quenches UV photons and
    // provides mild Penning transfer (Ar* 11.55 eV > isobutane IE 10.6 eV).
    // W-value ~26 eV (close to pure Ar; Sauli 1977).
    // =====================================================
    {
        G4double fAr=0.95, fIso=0.05;
        G4double rho = fAr*1.782e-3 + fIso*2.67e-3;
        G4Material* m = new G4Material("ArIso", rho*g/cm3, 2, kStateGas, 293.15*kelvin, 1*atmosphere);
        m->AddMaterial(purAr,     fAr);
        m->AddMaterial(isobutane, fIso);
        fGasMaterials["ArIso"] = m;
    }

    fGasMaterials["PureAr"]     = purAr;
    fGasMaterials["PureHe"]     = pureHe;
    fGasMaterials["PureNe"]     = pureNe;
    fGasMaterials["PureEthane"] = ethane;
    fGasMaterials["PureIso"]    = isobutane;
    fGasMaterials["PureCO2"]    = CO2;

    // ── He-3 gas at 300 bar ──────────────────────────────────
    // Isotopically pure He-3 (not natural helium) — critical for n-capture cross section.
    // Density from ideal gas law: rho = P*M / (R*T)
    //   = (300 * 101325 Pa) * (3.016e-3 kg/mol) / (8.314 J/mol/K * 293.15 K)
    //   = 37.6 kg/m3 = 0.0376 g/cm3
    {
        G4Isotope* isoHe3 = new G4Isotope("He3_iso", 2, 3, 3.0160293*g/mole);
        G4Element* elHe3  = new G4Element("Helium3_elem", "3He", 1);
        elHe3->AddIsotope(isoHe3, 1.0);
        G4Material* m = new G4Material("He3Gas_300bar", 37.6e-3*g/cm3, 1,
                                        kStateGas, 293.15*kelvin, 300*atmosphere);
        m->AddElement(elHe3, 1);
        fGasMaterials["He3Gas_300bar"] = m;
    }

    // ── CFRP: carbon-fibre reinforced polymer ────────────────
    // 70% C-fibre (pure C, 1.75 g/cm3) + 30% epoxy (C2H4O, 1.2 g/cm3) by volume.
    // Stated overall density 1.55 g/cm3. Mass fractions computed from volume fractions.
    {
        G4Material* m = new G4Material("CFRP", 1.55*g/cm3, 3);
        m->AddElement(elC, 0.8968);
        m->AddElement(elH, 0.0207);
        m->AddElement(elO, 0.0826);
        fGasMaterials["CFRP"] = m;
    }

    // ── Resistive paste (Saral Carbon 700A analogue) ─────────
    // Carbon black in organic binder (approx. carbon-loaded acrylic).
    // Composition: C 65%, H 8%, O 27% by mass; density 1.4 g/cm3.
    {
        G4Material* m = new G4Material("ResistivePaste", 1.4*g/cm3, 3);
        m->AddElement(elC, 0.65);
        m->AddElement(elH, 0.08);
        m->AddElement(elO, 0.27);
        fGasMaterials["ResistivePaste"] = m;
    }

    // ── FR4 epoxy-glass PCB ──────────────────────────────────
    // 60% SiO2 glass + 40% C11H12O3 epoxy by weight; density 1.85 g/cm3.
    {
        G4Material* m = new G4Material("FR4", 1.85*g/cm3, 4);
        m->AddElement(elSi, 0.2805);
        m->AddElement(elO,  0.4195);
        m->AddElement(elC,  0.2750);
        m->AddElement(elH,  0.0250);
        fGasMaterials["FR4"] = m;
    }

    // ── Rohacell 51 PMI foam ─────────────────────────────────
    // Polymethacrylimide (C4H5NO), density 0.052 g/cm3.
    {
        G4Material* m = new G4Material("Rohacell51", 0.052*g/cm3, 4);
        m->AddElement(elC, 0.5783);
        m->AddElement(elH, 0.0602);
        m->AddElement(elN, 0.1687);
        m->AddElement(elO, 0.1928);
        fGasMaterials["Rohacell51"] = m;
    }

    // ── Linear Alkylbenzene liquid scintillator ───────────────
    // JUNO recipe: pure LAB (C18H30) with PPO/Bis-MSB fluors (<0.3% — ignored).
    // Density 0.86 g/cm3.
    {
        G4Material* m = new G4Material("LAB_LiqScint", 0.86*g/cm3, 2);
        m->AddElement(elC, 0.8780);
        m->AddElement(elH, 0.1220);
        fGasMaterials["LAB_LiqScint"] = m;
    }

    // ── PVC black tape ───────────────────────────────────────
    // (C2H3Cl)n, density 1.3 g/cm3.
    {
        G4Material* m = new G4Material("PVC_tape", 1.3*g/cm3, 3);
        m->AddElement(elC,  0.3843);
        m->AddElement(elH,  0.0480);
        m->AddElement(elCl, 0.5677);
        fGasMaterials["PVC_tape"] = m;
    }
}

// ============================================================
G4Material* DetectorConstruction::GetGasMixture(const std::string& name) {
    auto it = fGasMaterials.find(name);
    if (it == fGasMaterials.end()) {
        throw std::runtime_error("Unknown gas mixture: " + name +
            "\nAvailable: ArCF4, ArIso, HeEth, ArCO2, ArCF4Iso, NeIso, NeCF4, ArCF4CO2, "
            "PureCF4, PureAr, PureHe, PureNe, PureEthane, PureIso, PureCO2");
    }
    return it->second;
}

// ============================================================
G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();

    G4NistManager* nist = G4NistManager::Instance();

    G4Material* matAir    = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matMylar  = nist->FindOrBuildMaterial("G4_MYLAR");
    G4Material* matAl     = nist->FindOrBuildMaterial("G4_Al");
    G4Material* matKapton = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material* matCu     = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* matSteel  = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* matGas    = GetGasMixture(fConfig.gas);

    // Full-experiment structural materials
    G4Material* matHe3      = fGasMaterials.at("He3Gas_300bar");
    G4Material* matCFRP     = fGasMaterials.at("CFRP");
    G4Material* matResPaste = fGasMaterials.at("ResistivePaste");
    G4Material* matFR4      = fGasMaterials.at("FR4");
    G4Material* matRohacell = fGasMaterials.at("Rohacell51");
    G4Material* matLAB      = fGasMaterials.at("LAB_LiqScint");
    G4Material* matPVC      = fGasMaterials.at("PVC_tape");
    G4Material* matPlScint  = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // ── Transverse size ──────────────────────────────────────
    G4double detXY = 40.0 * cm;

    // ── Micromegas layer thicknesses (both modes) ────────────
    G4double tMylar    = 40.0  * um;
    G4double tAlWin    = 0.1   * um;
    G4double tKapCath  = 50.0  * um;
    G4double tCuCath   = 9.0   * um;
    G4double tDrift    = 3.0   * cm;
    G4double tMesh     = 30.0  * um;
    G4double tAmp      = 150.0 * um;
    G4double tResPaste = 100.0 * um;   // resistive paste (both modes)
    G4double mmTotalZ = tMylar + tAlWin + tKapCath + tCuCath +
                        tDrift + tMesh + tAmp + tResPaste;

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

    // ── Build world volume ───────────────────────────────────
    G4LogicalVolume*   worldLV = nullptr;
    G4VPhysicalVolume* worldPV = nullptr;

    if (fConfig.mode == SimMode::kVacuum) {

        // ======================================================
        // VACUUM MODE
        // ======================================================
        G4double alThickness = fConfig.alThickness_mm * mm;
        G4double alGap       = 2.0 * cm;

        G4double gunZ = 10.0 * cm;
        G4double minUpstream = 2.0*gunZ - mmTotalZ - 2.5*cm;
        G4double upstreamMargin = std::max(alGap + alThickness + 0.5*cm, minUpstream);
        G4double worldZ = mmTotalZ + upstreamMargin + 2.5*cm;

        G4Box* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, worldZ/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        // Al shielding upstream
        if (alThickness > 0) {
            G4double alZCenter = -mmTotalZ/2.0 - alGap - alThickness/2.0;
            G4Box* alBox = new G4Box("AlShield", detXY/2, detXY/2, alThickness/2);
            G4LogicalVolume* alLV = new G4LogicalVolume(alBox, matAl, "AlShield");
            auto visAlShield = new G4VisAttributes(G4Color(0.75, 0.75, 0.75, 0.9));
            visAlShield->SetForceSolid(true);
            alLV->SetVisAttributes(visAlShield);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,alZCenter),
                              alLV, "AlShield", worldLV, false, 0, true);
        }

        G4cout << "\n=== Vacuum mode geometry ===" << G4endl;
        G4cout << "  Gas       : " << fConfig.gas
               << " (rho = " << matGas->GetDensity()/(mg/cm3) << " mg/cm3)" << G4endl;
        if (alThickness > 0)
            G4cout << "  Al shield : " << fConfig.alThickness_mm << " mm" << G4endl;

    } else if (fConfig.mode == SimMode::kSr90Calibration) {

        // ======================================================
        // SR-90 CALIBRATION MODE
        // Identical to the full experiment but the He-3 pressurised
        // target and its Al/CFRP capsule are absent.  The source is
        // placed in air at the same z-position as the He-3 gas centre.
        // The combined air path from source to MM drift gas entry is
        // 226.5 mm (= 25 mm He-3 gas + 1.4 mm walls + 200 mm gap 1
        // + ~0.1 mm dead layers, all replaced by air).
        // ======================================================

        G4double cfrpT       = fConfig.cfrpThickness_mm * mm;
        G4double airToMM     = 226.5 * mm;   // source → MM drift gas entry (all air)

        // PCB layer thicknesses — identical to full mode
        G4double tPCB_Kap      = 50.0  * um;
        G4double tPCB_Cu       = 26.0  * um;
        G4double tPCB_FR4      = 100.0 * um;
        G4double tPCB_Rohacell = 5.0   * mm;
        G4double tPCB_AlFoil   = 50.0  * um;
        G4double pcbTotalZ = tPCB_Kap + 4*(tPCB_Cu + tPCB_FR4) + tPCB_Rohacell + tPCB_AlFoil;

        G4double airGap2 = 20.0 * mm;

        G4double tBlackTape   = 165.0 * um;
        G4double tPlScint     = 3.0   * mm;
        G4double tScintAlFoil = 50.0  * um;
        G4double scintWallZ   = 2*tBlackTape + tPlScint + tScintAlFoil;

        G4double airGap3 = 20.0 * mm;

        G4double tLS      = 15.0 * mm;
        G4double lsStackZ = 5*cfrpT + 4*tLS;

        G4double totalSr90Z = airToMM + mmTotalZ + pcbTotalZ
                            + airGap2 + scintWallZ + airGap3 + lsStackZ;

        G4double worldZ = totalSr90Z + 2.0*cm;

        G4Box* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, worldZ/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        // Source is in air at the upstream face of the world volume.
        // fHe3GasCenterZ stores the gun position (reuses the accessor name).
        fHe3GasCenterZ = -totalSr90Z / 2.0;

        G4cout << "\n=== Sr-90 calibration geometry ===" << G4endl;
        G4cout << "  No He-3 target or capsule" << G4endl;
        G4cout << "  Air source-to-MM : " << airToMM/mm << " mm" << G4endl;
        G4cout << "  Source z = " << fHe3GasCenterZ/mm << " mm" << G4endl;
        G4cout << "  Total stack Z : " << totalSr90Z/mm << " mm" << G4endl;

        G4double zFrontSr90 = -totalSr90Z / 2.0;

        auto PlaceSlabSr90 = [&](const std::string& name, G4double thickness,
                                  G4Material* mat, G4VisAttributes* vis,
                                  G4LogicalVolume*& outLV) {
            G4double zCenter = zFrontSr90 + thickness/2.0;
            G4Box* solid = new G4Box(name, detXY/2, detXY/2, thickness/2);
            outLV = new G4LogicalVolume(solid, mat, name);
            if (vis) outLV->SetVisAttributes(vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zCenter),
                              outLV, name, worldLV, false, 0, true);
            zFrontSr90 += thickness;
        };
        G4LogicalVolume* dummyLV = nullptr;

        // Air gap from source to MM entrance (replaces He-3 target + gap 1)
        PlaceSlabSr90("AirGap1", airToMM, matAir, nullptr, dummyLV);

        // ── Micromegas (identical to full mode) ───────────────────────────
        PlaceSlabSr90("GasWindow_Mylar",    tMylar,    matMylar,    visMylar,    dummyLV);
        PlaceSlabSr90("GasWindow_Al",       tAlWin,    matAl,       visAl,       dummyLV);
        PlaceSlabSr90("DriftCathode_Kapton",tKapCath,  matKapton,   visKapton,   dummyLV);
        PlaceSlabSr90("DriftCathode_Cu",    tCuCath,   matCu,       visCu,       dummyLV);
        PlaceSlabSr90("DriftGas",           tDrift,    matGas,      visDrift,    fDriftGasLV);
        PlaceSlabSr90("Micromesh",          tMesh,     matSteel,    visMesh,     dummyLV);
        PlaceSlabSr90("AmpGas",             tAmp,      matGas,      visAmp,      fAmpGasLV);
        PlaceSlabSr90("ResistivePaste",     tResPaste, matResPaste, visResPaste, dummyLV);

        // ── PCB stack ────────────────────────────────────────────────────
        PlaceSlabSr90("PCB_Kapton",   tPCB_Kap,      matKapton,  visKapton,  dummyLV);
        PlaceSlabSr90("PCB_Cu_1",     tPCB_Cu,       matCu,      visCu,      dummyLV);
        PlaceSlabSr90("PCB_FR4_1",    tPCB_FR4,      matFR4,     visFR4,     dummyLV);
        PlaceSlabSr90("PCB_Cu_2",     tPCB_Cu,       matCu,      visCu,      dummyLV);
        PlaceSlabSr90("PCB_FR4_2",    tPCB_FR4,      matFR4,     visFR4,     dummyLV);
        PlaceSlabSr90("PCB_Cu_3",     tPCB_Cu,       matCu,      visCu,      dummyLV);
        PlaceSlabSr90("PCB_FR4_3",    tPCB_FR4,      matFR4,     visFR4,     dummyLV);
        PlaceSlabSr90("PCB_Cu_4",     tPCB_Cu,       matCu,      visCu,      dummyLV);
        PlaceSlabSr90("PCB_FR4_4",    tPCB_FR4,      matFR4,     visFR4,     dummyLV);
        PlaceSlabSr90("PCB_Rohacell", tPCB_Rohacell, matRohacell,visRohacell,dummyLV);
        PlaceSlabSr90("PCB_AlFoil",   tPCB_AlFoil,   matAl,      visAl,      dummyLV);

        PlaceSlabSr90("AirGap2", airGap2, matAir, nullptr, dummyLV);

        // ── Scintillator wall ─────────────────────────────────────────────
        PlaceSlabSr90("ScintWall_BlackTape1", tBlackTape,   matPVC,     nullptr,  dummyLV);
        PlaceSlabSr90("PlasticScint",         tPlScint,     matPlScint, visScint, dummyLV);
        PlaceSlabSr90("ScintWall_BlackTape2", tBlackTape,   matPVC,     nullptr,  dummyLV);
        PlaceSlabSr90("ScintWall_AlFoil",     tScintAlFoil, matAl,      visAl,    dummyLV);

        PlaceSlabSr90("AirGap3", airGap3, matAir, nullptr, dummyLV);

        // ── Liquid scintillator stack ────────────────────────────────────
        PlaceSlabSr90("LS_CFRP_1",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabSr90("LiqScint_1", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabSr90("LS_CFRP_2",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabSr90("LiqScint_2", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabSr90("LS_CFRP_3",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabSr90("LiqScint_3", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabSr90("LS_CFRP_4",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabSr90("LiqScint_4", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabSr90("LS_CFRP_5",  cfrpT, matCFRP, visCFRP, dummyLV);

        return worldPV;

    } else {

        // ======================================================
        // FULL EXPERIMENT MODE
        // ======================================================

        // He-3 capsule dimensions.
        // The cylinder AXIS is along world Y (perpendicular to the beam).
        // The beam (+z) passes through the CURVED surface, traversing 2×he3R = 5 cm of gas.
        // cfrpT_capsule is fixed at the specified 0.9 mm; the configurable cfrpT below
        // applies only to the liquid-scintillator cell walls.
        G4double he3R          = 2.5  * cm;    // He-3 gas radius
        G4double he3HalfL      = 7.5  * cm;    // He-3 gas half-length along capsule axis (Y)
        G4double alWallT       = 0.5  * mm;    // Al capsule wall
        G4double cfrpT_capsule = 0.9  * mm;    // CFRP capsule wall (fixed)
        G4double cfrpT         = fConfig.cfrpThickness_mm * mm;  // LS cell CFRP (configurable)

        // Radii (xz-plane, along beam direction)
        G4double alR   = he3R + alWallT;
        G4double cfrpR = he3R + alWallT + cfrpT_capsule;  // = 26.4 mm

        // Half-lengths along capsule axis (Y)
        G4double alHalfL   = he3HalfL + alWallT;
        G4double cfrpHalfL = he3HalfL + alWallT + cfrpT_capsule;

        // Along beam (Z): the capsule spans ±cfrpR around its centre
        G4double capsuleZExtent = 2.0 * cfrpR;

        G4double airGap1 = 200.0 * mm;

        // PCB layer thicknesses
        G4double tPCB_Kap      = 50.0  * um;
        G4double tPCB_Cu       = 26.0  * um;  // ×4
        G4double tPCB_FR4      = 100.0 * um;  // ×4
        G4double tPCB_Rohacell = 5.0   * mm;
        G4double tPCB_AlFoil   = 50.0  * um;
        G4double pcbTotalZ = tPCB_Kap + 4*(tPCB_Cu + tPCB_FR4) + tPCB_Rohacell + tPCB_AlFoil;

        G4double airGap2 = 20.0 * mm;

        // Scintillator wall
        G4double tBlackTape   = 165.0 * um;
        G4double tPlScint     = 3.0   * mm;
        G4double tScintAlFoil = 50.0  * um;
        G4double scintWallZ   = 2*tBlackTape + tPlScint + tScintAlFoil;

        G4double airGap3 = 20.0 * mm;

        // LS stack: CFRP | LS1 | CFRP | LS2 | CFRP | LS3 | CFRP | LS4 | CFRP
        G4double tLS      = 15.0 * mm;
        G4double lsStackZ = 5*cfrpT + 4*tLS;

        G4double totalFullZ = capsuleZExtent + airGap1 + mmTotalZ + pcbTotalZ
                            + airGap2 + scintWallZ + airGap3 + lsStackZ;

        G4double worldZ = totalFullZ + 2.0*cm;  // 1 cm clearance each side

        G4Box* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, worldZ/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        // ── He-3 pressurised capsule (nested cylinders, axis along Y) ──
        // The capsule centre is placed so its upstream CFRP face sits at z = -totalFullZ/2.
        // Rotation: local Z (G4Tubs axis) → world Y via rotateX(-90°).
        // Daughters (Al, He3) inherit the rotation from their parent and need no extra rotation.
        G4double capsuleZCenter = -totalFullZ/2.0 + cfrpR;

        fHe3GasCenterZ = capsuleZCenter;  // He-3 gas is centred in the capsule

        G4RotationMatrix* capRot = new G4RotationMatrix();
        capRot->rotateX(-90.*deg);  // local z → world y

        // CFRP outer shell
        G4Tubs* cfrpSolid = new G4Tubs("He3Capsule_CFRP", 0, cfrpR, cfrpHalfL, 0, 360.*deg);
        G4LogicalVolume* cfrpLV = new G4LogicalVolume(cfrpSolid, matCFRP, "He3Capsule_CFRP");
        cfrpLV->SetVisAttributes(visCFRP);
        new G4PVPlacement(capRot, G4ThreeVector(0,0,capsuleZCenter),
                          cfrpLV, "He3Capsule_CFRP", worldLV, false, 0, true);

        // Al inner shell (daughter of CFRP — no additional rotation needed)
        G4Tubs* alSolid = new G4Tubs("He3Capsule_Al", 0, alR, alHalfL, 0, 360.*deg);
        G4LogicalVolume* alLV = new G4LogicalVolume(alSolid, matAl, "He3Capsule_Al");
        alLV->SetVisAttributes(visAl);
        new G4PVPlacement(nullptr, G4ThreeVector(), alLV, "He3Capsule_Al", cfrpLV, false, 0, true);

        // He-3 gas (daughter of Al shell)
        G4Tubs* he3Solid = new G4Tubs("He3Gas", 0, he3R, he3HalfL, 0, 360.*deg);
        fHe3GasLV = new G4LogicalVolume(he3Solid, matHe3, "He3Gas");
        fHe3GasLV->SetVisAttributes(visHe3);
        new G4PVPlacement(nullptr, G4ThreeVector(), fHe3GasLV, "He3Gas", alLV, false, 0, true);

        G4cout << "\n=== Full-experiment geometry ===" << G4endl;
        G4cout << "  He-3 capsule: r=" << he3R/cm << " cm (5 cm gas along beam)"
               << ", L=" << 2*he3HalfL/cm << " cm along Y"
               << ", 300 bar" << G4endl;
        G4cout << "  Capsule walls: Al " << alWallT/mm << " mm + CFRP "
               << cfrpT_capsule/mm << " mm" << G4endl;
        G4cout << "  He-3 gas centre z = " << fHe3GasCenterZ/cm << " cm" << G4endl;
        G4cout << "  LS cell CFRP walls: " << cfrpT/mm << " mm" << G4endl;
        G4cout << "  Total stack Z : " << totalFullZ/mm << " mm" << G4endl;

        // ── PlaceSlab for downstream layers ──────────────────
        // zFront starts right after the capsule (capsule rear face)
        G4double zFrontFull = -totalFullZ/2.0 + capsuleZExtent;

        auto PlaceSlabFull = [&](const std::string& name, G4double thickness,
                                  G4Material* mat, G4VisAttributes* vis,
                                  G4LogicalVolume*& outLV) {
            G4double zCenter = zFrontFull + thickness/2.0;
            G4Box* solid = new G4Box(name, detXY/2, detXY/2, thickness/2);
            outLV = new G4LogicalVolume(solid, mat, name);
            if (vis) outLV->SetVisAttributes(vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zCenter),
                              outLV, name, worldLV, false, 0, true);
            zFrontFull += thickness;
        };
        G4LogicalVolume* dummyLV = nullptr;

        // Air gap 1
        PlaceSlabFull("AirGap1", airGap1, matAir, nullptr, dummyLV);

        // ── Micromegas (in-order) ─────────────────────────────
        PlaceSlabFull("GasWindow_Mylar",    tMylar,    matMylar,   visMylar,  dummyLV);
        PlaceSlabFull("GasWindow_Al",       tAlWin,    matAl,      visAl,     dummyLV);
        PlaceSlabFull("DriftCathode_Kapton",tKapCath,  matKapton,  visKapton, dummyLV);
        PlaceSlabFull("DriftCathode_Cu",    tCuCath,   matCu,      visCu,     dummyLV);
        PlaceSlabFull("DriftGas",           tDrift,    matGas,     visDrift,  fDriftGasLV);
        PlaceSlabFull("Micromesh",          tMesh,     matSteel,   visMesh,   dummyLV);
        PlaceSlabFull("AmpGas",             tAmp,      matGas,     visAmp,    fAmpGasLV);
        PlaceSlabFull("ResistivePaste",     tResPaste, matResPaste,visResPaste,dummyLV);

        // ── PCB stack ─────────────────────────────────────────
        PlaceSlabFull("PCB_Kapton",   tPCB_Kap,     matKapton,  visKapton, dummyLV);
        PlaceSlabFull("PCB_Cu_1",     tPCB_Cu,      matCu,      visCu,     dummyLV);
        PlaceSlabFull("PCB_FR4_1",    tPCB_FR4,     matFR4,     visFR4,    dummyLV);
        PlaceSlabFull("PCB_Cu_2",     tPCB_Cu,      matCu,      visCu,     dummyLV);
        PlaceSlabFull("PCB_FR4_2",    tPCB_FR4,     matFR4,     visFR4,    dummyLV);
        PlaceSlabFull("PCB_Cu_3",     tPCB_Cu,      matCu,      visCu,     dummyLV);
        PlaceSlabFull("PCB_FR4_3",    tPCB_FR4,     matFR4,     visFR4,    dummyLV);
        PlaceSlabFull("PCB_Cu_4",     tPCB_Cu,      matCu,      visCu,     dummyLV);
        PlaceSlabFull("PCB_FR4_4",    tPCB_FR4,     matFR4,     visFR4,    dummyLV);
        PlaceSlabFull("PCB_Rohacell", tPCB_Rohacell,matRohacell,visRohacell,dummyLV);
        PlaceSlabFull("PCB_AlFoil",   tPCB_AlFoil,  matAl,      visAl,     dummyLV);

        // Air gap 2
        PlaceSlabFull("AirGap2", airGap2, matAir, nullptr, dummyLV);

        // ── Scintillator wall ─────────────────────────────────
        PlaceSlabFull("ScintWall_BlackTape1", tBlackTape,   matPVC,      nullptr,   dummyLV);
        PlaceSlabFull("PlasticScint",         tPlScint,     matPlScint,  visScint,  dummyLV);
        PlaceSlabFull("ScintWall_BlackTape2", tBlackTape,   matPVC,      nullptr,   dummyLV);
        PlaceSlabFull("ScintWall_AlFoil",     tScintAlFoil, matAl,       visAl,     dummyLV);

        // Air gap 3
        PlaceSlabFull("AirGap3", airGap3, matAir, nullptr, dummyLV);

        // ── Liquid scintillator stack ─────────────────────────
        PlaceSlabFull("LS_CFRP_1",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabFull("LiqScint_1", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabFull("LS_CFRP_2",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabFull("LiqScint_2", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabFull("LS_CFRP_3",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabFull("LiqScint_3", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabFull("LS_CFRP_4",  cfrpT, matCFRP, visCFRP, dummyLV);
        PlaceSlabFull("LiqScint_4", tLS,   matLAB,  visLAB,  dummyLV);
        PlaceSlabFull("LS_CFRP_5",  cfrpT, matCFRP, visCFRP, dummyLV);

        return worldPV;
    }

    // ======================================================
    // Shared: place Micromegas layers in vacuum mode
    // (Full mode already placed them above and returned.)
    // ======================================================

    G4double zFront = -mmTotalZ / 2.0;

    auto PlaceSlab = [&](const std::string& name, G4double thickness,
                          G4Material* mat, G4VisAttributes* vis,
                          G4LogicalVolume*& outLV) {
        G4double zCenter = zFront + thickness/2.0;
        G4Box* solid = new G4Box(name, detXY/2, detXY/2, thickness/2);
        outLV = new G4LogicalVolume(solid, mat, name);
        if (vis) outLV->SetVisAttributes(vis);
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,zCenter),
                          outLV, name, worldLV, false, 0, true);
        zFront += thickness;
    };
    G4LogicalVolume* dummyLV = nullptr;

    PlaceSlab("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,   dummyLV);
    PlaceSlab("GasWindow_Al",        tAlWin,    matAl,       visAl,      dummyLV);
    PlaceSlab("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,  dummyLV);
    PlaceSlab("DriftCathode_Cu",     tCuCath,   matCu,       visCu,      dummyLV);
    PlaceSlab("DriftGas",            tDrift,    matGas,      visDrift,   fDriftGasLV);
    PlaceSlab("Micromesh",           tMesh,     matSteel,    visMesh,    dummyLV);
    PlaceSlab("AmpGas",              tAmp,      matGas,      visAmp,     fAmpGasLV);
    PlaceSlab("ResistivePaste",      tResPaste, matResPaste, visResPaste,dummyLV);

    G4cout << "  Gas mixture  : " << fConfig.gas
           << " (rho = " << matGas->GetDensity()/(mg/cm3) << " mg/cm3)" << G4endl;
    G4cout << "  Drift gap    : " << tDrift/cm << " cm" << G4endl;
    G4cout << "  Amp gap      : " << tAmp/um << " um" << G4endl;
    G4cout << "=================================\n" << G4endl;

    return worldPV;
}

// ============================================================
void DetectorConstruction::ConstructSDandField() {
    SensitiveDetector* driftSD = new SensitiveDetector("DriftGasSD", "DriftGasHits",
                                                         "DriftGas", fConfig);
    G4SDManager::GetSDMpointer()->AddNewDetector(driftSD);
    SetSensitiveDetector(fDriftGasLV, driftSD);

    SensitiveDetector* ampSD = new SensitiveDetector("AmpGasSD", "AmpGasHits",
                                                       "AmpGas", fConfig);
    G4SDManager::GetSDMpointer()->AddNewDetector(ampSD);
    SetSensitiveDetector(fAmpGasLV, ampSD);

    G4UserLimits* gasStepLimit = new G4UserLimits(100*um);
    fDriftGasLV->SetUserLimits(gasStepLimit);
    fAmpGasLV->SetUserLimits(gasStepLimit);

    // He-3 gas: coarser step limit (not a cluster-scoring volume)
    if (fHe3GasLV) {
        fHe3GasLV->SetUserLimits(new G4UserLimits(1.0*mm));
    }
}

