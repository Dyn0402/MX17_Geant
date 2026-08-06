#pragma once
// MX17ModuleGeometry.hh
//
// Single source of truth for the MX17 Micromegas module layer stack, shared
// between MX17_Geant and MX17_Full_Geant (the latter includes this header via
// a relative path configured in its CMakeLists.txt — see the guard there).
//
// The module is described in a LOCAL frame: x,y transverse, z along the module
// normal with z = 0 at the upstream face of the flat window sheet and +z
// downstream into the detector. A bulged window extends upstream to z < 0
// (Module::frontExtent). Consumers place the returned pieces with their own
// rotation / translation, so each keeps its own axis convention (MX17_Geant:
// z along beam, one module on axis; MX17_Full_Geant: (u,v,w) per arm, 4 arms).
//
// Two factory specs:
//   LegacySpec(faceX, faceY) — reproduces the pre-2026-08 stack exactly
//     (uniform transverse face, 40 µm window, contiguous layers, solid mesh,
//     0.554 mm laminate, 50 µm Al back foil, no lateral structure).
//   AsBuiltSpec() — the CAD-derived module (design/GEOMETRY_FROM_CAD.md +
//     design/GEOMETRY_IMPLEMENTATION_NOTES.md): 60 µm window, 5 mm flange gas
//     gap, per-layer transverse faces (440/406/410/470), Al frame + flange
//     rings, field cage, woven-mesh effective-density slab, 1.70 mm readout
//     board, aluminized-mylar back foil, 8 mm support plate with the 402 mm
//     aperture, optional terraced-dome window bulge.
//
// Volume names are load-bearing: SteppingAction in both repos scores by name
// ("DriftGas", "AmpGas", "GasWindow_Mylar", "PCB_*", ...). All window-mylar
// pieces (flat sheet, annulus, dome terraces) share the name "GasWindow_Mylar"
// deliberately so the scoring stays correct; same for "GasWindow_Al".

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace MX17 {

// ─────────────────────────────────────────────────────────────────────────────
struct ModuleSpec {
    // Uniform transverse face (legacy). When > 0, every layer uses this
    // rectangle and all lateral structure (rings, cage, plate aperture,
    // offsets, bulge) is disabled regardless of the fields below.
    double uniformFaceX_mm = 0.0;
    double uniformFaceY_mm = 0.0;

    // ── Window (upstream face) ──
    double winMylar_um   = 40.0;   // mylar sheet
    double winAl_um      = 0.1;    // aluminization, drift-facing side (assumed)
    double winFace_mm    = 400.0;  // sheet transverse extent (as-built 440)
    double flangeGap_mm  = 0.0;    // window-flange thickness = window→cathode gas gap (as-built 5)
    double flangeAp_mm   = 0.0;    // flange aperture = window free span (as-built 400.2)
    double flangeOuter_mm= 0.0;    // flange ring outer (as-built 440); 0 → no flange/gap
    double bulgeSag_mm   = 0.0;    // overpressure dome sag over the aperture; 0 → flat
    int    bulgeSteps    = 6;      // terraces in the dome

    // ── Drift cathode ──
    double cathKapton_um = 50.0;
    double cathCu_um     = 9.0;    // copper cladding, drift-facing side (assumed)
    double cathFace_mm   = 400.0;  // as-built 406

    // ── Drift gap ──
    double drift_mm      = 30.0;   // kapton back face → mesh front (frame height)
    double driftFace_mm  = 400.0;  // drift gas transverse (as-built: frame aperture − 2×field cage)
    bool   cuInsideDrift = false;  // true: Cu cladding protrudes into the 30 mm (gas = drift − cathCu)
    double frameOuter_mm = 0.0;    // Al gas frame ring outer (as-built 440); 0 → none
    double frameAp_mm    = 0.0;    // frame aperture (as-built 410)
    double cageThick_mm  = 0.0;    // field-cage PCB thickness lining the aperture; 0 → none

    // ── Micromesh ──
    double meshWire_um   = 0.0;    // woven model: wire diameter; 0 → legacy solid slab
    double meshOpen_um   = 0.0;    // woven model: opening
    double meshSolid_um  = 30.0;   // legacy solid-steel slab thickness
    double meshFace_mm   = 400.0;  // as-built 410

    // ── Amplification gap + resistive layer ──
    double amp_um        = 150.0;
    double paste_um      = 100.0;
    double ampFace_mm    = 400.0;  // as-built 410
    double pasteFace_mm  = 0.0;    // 0 → ampFace; as-built 412 (resist-coat gerber)
    // ESL resistive strips: 550 µm wide / 250 µm gaps (confirmed 2026-08-06)
    // → the 100 µm slab is density-scaled by 550/800. 0 → uniform sheet.
    double pasteCoverage = 0.0;
    // Bulk pillars supporting the mesh, standing in the amplification gap:
    // Ø0.6 mm on a regular 85×85 grid at 4.68 mm pitch spanning ±196.56 mm
    // (3498A_bulk.gbr, exact). Placed as daughters of the AmpGas volume so
    // the scored gas is properly displaced (~1.3 % of the gap volume).
    double pillarD_mm    = 0.0;    // 0 → no pillars
    double pillarPitch_mm= 4.68;
    int    pillarN       = 85;

    // ── Readout laminate ──
    bool   includeReadout = true;  // false: stop after ResistivePaste (vacuum mode)
    double pcbKapton_um  = 50.0;
    double pcbCu_um      = 26.0;
    double pcbFR4_um     = 100.0;  // per layer; ignored when pcbTotal_mm > 0
    int    pcbNLayers    = 4;
    double pcbTotal_mm   = 0.0;    // >0: FR4 per-layer solved so mesh front → laminate back = pcbTotal
    // Cu coverage per readout layer (physical L3..L7; L8 is outline-only),
    // measured over the 470 mm board face from the production gerbers
    // (scripts/gerber/analyze_cu_coverage.py). When non-empty it replaces the
    // pcbNLayers full-density sheets with one density-scaled 26 µm sheet per
    // entry: L3 guard ring, L4 pads, L5 Y strips, L6 X strips, L7 fan-out.
    std::vector<double> pcbCuCoverage;
    double pcbFace_mm    = 400.0;  // as-built 470
    double pcbOffset_mm  = 0.0;    // +x,+y offset of the 470 mm plates from the active-area axis (as-built 15)

    // ── Front-end (multi M1) cards ──
    // Four cards lying FLAT on the drift-gas side of the readout-board edge,
    // straddling it longways, two per connector edge (+x and +y in the
    // active-area frame). The M1 gerbers are drawn in the BOARD frame at the
    // as-mounted position: card outline 41 × 160 mm at radial 219.4..260.4,
    // tangential centre +100 (edge-mill layer DFS3498AM1_ebm); the pogo-pad
    // field on the card bottom lands exactly on the board's connector copper
    // clusters (tangential ±83..117 on both edges). No standoff: the card
    // laminate sits directly on the board (the thin pogo-pin connector and
    // the Mec8 output connectors on top are not modelled —
    // design/NEEDED_INPUTS.md). The inner edge is clipped from 219.4 to
    // 220.0 mm where the plain-ring gas-frame model has its outer wall.
    double feThick_mm    = 0.0;    // card laminate thickness; 0 → none
    double feInner_mm    = 220.0;  // inner edge (gerber 219.4, clipped to frame)
    double feOuter_mm    = 260.4;  // outer edge (10.4 mm beyond the board edge)
    double feLen_mm      = 160.0;  // tangential length
    double feTang_mm     = 100.0;  // tangential centres at ±feTang on each edge
    double feCu_um       = 26.0;   // per Cu layer (same guess as the readout)
    // Cu coverage per M1 layer (top,L2..L5,bot) over the card outline, from
    // scripts/gerber/analyze_cu_coverage.py; empty → solid FR4 card.
    std::vector<double> feCuCoverage;

    // ── Backing ──
    double rohacell_mm   = 5.0;
    double backMylar_um  = 0.0;    // aluminized-mylar foil behind the rohacell (as-built, thickness TBC)
    double backAl_um     = 50.0;   // its aluminization — legacy uses 50 µm bare Al here
    double plate_mm      = 0.0;    // Al support plate; 0 → none (as-built 8)
    double plateAp_mm    = 0.0;    // its square aperture, centred on the active-area axis (as-built 402)

    bool Uniform() const { return uniformFaceX_mm > 0.0; }
};

// Legacy stack: numerically identical to the pre-refactor DetectorConstruction
// in both repos (only the transverse face differs per consumer).
inline ModuleSpec LegacySpec(double faceX_mm, double faceY_mm) {
    ModuleSpec s;
    s.uniformFaceX_mm = faceX_mm;
    s.uniformFaceY_mm = faceY_mm;
    return s;
}

// As-built module from the mechanical CAD + hardware knowledge.
// bulgeSag_mm default from a Hencky estimate at ~3 mbar overpressure for the
// 400.2 mm free span and 60 µm foil (design/GEOMETRY_IMPLEMENTATION_NOTES.md §4).
inline ModuleSpec AsBuiltSpec(double bulgeSag_mm = 8.0) {
    ModuleSpec s;
    s.winMylar_um    = 60.0;    // CAD (was 40)
    s.winFace_mm     = 440.0;
    s.flangeGap_mm   = 5.0;
    s.flangeAp_mm    = 400.2;
    s.flangeOuter_mm = 440.0;
    s.bulgeSag_mm    = bulgeSag_mm;
    s.cathFace_mm    = 406.0;
    s.cuInsideDrift  = true;    // frame is the 30 mm spacer, kapton→mesh = 30.000
    s.frameOuter_mm  = 440.0;
    s.frameAp_mm     = 410.0;
    s.cageThick_mm   = 1.62;    // 4× drift field PCB lining the aperture
    s.driftFace_mm   = 410.0 - 2*1.62;
    s.meshWire_um    = 19.0;    // P2-like weave, placeholder (design/NEEDED_INPUTS.md)
    s.meshOpen_um    = 48.0;
    s.meshFace_mm    = 410.0;
    s.ampFace_mm     = 410.0;
    s.pasteFace_mm   = 412.0;   // resistive coat boundary (3498A_top-resist)
    s.pasteCoverage  = 550.0 / 800.0;   // ESL strips 550/250 µm (confirmed)
    s.pillarD_mm     = 0.6;     // bulk pillars (3498A_bulk.gbr grid)
    s.feThick_mm     = 1.6;     // multi M1 cards: bare-laminate guess
    s.feCuCoverage   = {0.147, 0.110, 0.114, 0.117, 0.111, 0.041};
    s.pcbTotal_mm    = 1.70;    // CAD single-body readout board
    s.pcbCuCoverage  = {0.0945, 0.6512, 0.4194, 0.4200, 0.4558};
    s.pcbFace_mm     = 470.0;
    s.pcbOffset_mm   = 15.0;    // plates offset (+15,+15) from the active-area axis (STEP + gerbers)
    s.backMylar_um   = 25.0;    // assumed (design/NEEDED_INPUTS.md)
    s.backAl_um      = 0.1;
    s.plate_mm       = 8.0;
    s.plateAp_mm     = 402.0;   // through-aperture, concentric with the active area
    return s;
}

// ─────────────────────────────────────────────────────────────────────────────
struct Materials {
    G4Material* mylar    = nullptr;
    G4Material* al       = nullptr;
    G4Material* kapton   = nullptr;
    G4Material* cu       = nullptr;
    G4Material* steel    = nullptr;
    G4Material* fr4      = nullptr;
    G4Material* resPaste = nullptr;
    G4Material* rohacell = nullptr;
    G4Material* gas      = nullptr;   // chamber gas (drift/amp/gap/bulge)
};

struct Piece {
    G4LogicalVolume* lv;
    G4ThreeVector    pos;      // centre in the module local frame
    double           thick;    // z the piece advances the stack by (0 for rings)
    bool             advance;  // true: slab in the z-chain; false: side structure
};

// Recommended placement: accumulate a running z front over the advancing
// pieces (centre = run + thick/2), which reproduces the historical per-slab
// `zF += t` arithmetic bit-for-bit; non-advancing pieces (rings, field cage,
// dome terraces) are placed at zFront + pos.z().
//   G4double run = zFront;
//   for (auto& p : mod.pieces) {
//       G4double z = p.advance ? run + p.thick/2 : zFront + p.pos.z();
//       ... place at (p.pos.x(), p.pos.y(), z) ...
//       if (p.advance) run += p.thick;
//   }

struct Module {
    std::vector<Piece> pieces;
    double frontExtent = 0.0;  // upstream reach beyond z=0 (window bulge sag)
    double mmDepth     = 0.0;  // z=0 → back of ResistivePaste (legacy "mmTotalZ")
    double readoutDepth= 0.0;  // laminate + backing depth (legacy "pcbTotalZ")
    double totalDepth  = 0.0;  // z=0 → back of last placed layer
    double halfX       = 0.0;  // max transverse half-extent incl. placement offsets
    double halfY       = 0.0;
    G4LogicalVolume* driftGasLV = nullptr;
    G4LogicalVolume* ampGasLV   = nullptr;
};

// ─────────────────────────────────────────────────────────────────────────────
inline Module BuildModule(const ModuleSpec& s, const Materials& m) {
    Module out;

    // Vis attributes (same palette both consumers used)
    auto visMylar    = new G4VisAttributes(G4Color(0.7, 0.9, 0.7, 0.5));
    auto visAl       = new G4VisAttributes(G4Color(0.7, 0.7, 0.7, 0.8));
    auto visKapton   = new G4VisAttributes(G4Color(0.9, 0.7, 0.0, 0.7));
    auto visCu       = new G4VisAttributes(G4Color(0.8, 0.4, 0.1, 0.8));
    auto visDrift    = new G4VisAttributes(G4Color(0.2, 0.5, 1.0, 0.3));
    auto visMesh     = new G4VisAttributes(G4Color(0.5, 0.5, 0.5, 0.9));
    auto visAmp      = new G4VisAttributes(G4Color(1.0, 0.3, 0.3, 0.3));
    auto visResPaste = new G4VisAttributes(G4Color(0.2, 0.2, 0.2, 0.8));
    auto visFR4      = new G4VisAttributes(G4Color(0.2, 0.6, 0.2, 0.8));
    auto visRohacell = new G4VisAttributes(G4Color(0.9, 0.9, 0.6, 0.5));
    auto visGasV     = new G4VisAttributes(G4Color(0.55, 0.85, 0.95, 0.15));
    auto visFrame    = new G4VisAttributes(G4Color(0.75, 0.75, 0.78, 0.85));

    double zF = 0.0;   // running upstream face of the next layer

    auto grow = [&](double hx, double hy, double ox, double oy) {
        out.halfX = std::max(out.halfX, hx + std::abs(ox));
        out.halfY = std::max(out.halfY, hy + std::abs(oy));
    };
    auto addBox = [&](const std::string& name, double hx, double hy, double t,
                      G4Material* mat, G4VisAttributes* vis,
                      double ox = 0.0, double oy = 0.0) -> G4LogicalVolume* {
        auto* box = new G4Box(name, hx, hy, t/2);
        auto* lv  = new G4LogicalVolume(box, mat, name);
        if (vis) lv->SetVisAttributes(vis);
        out.pieces.push_back({lv, G4ThreeVector(ox, oy, zF + t/2), t, true});
        grow(hx, hy, ox, oy);
        zF += t;
        return lv;
    };
    // Square ring (outer minus aperture), centre of the OUTER square at (ox,oy),
    // aperture centred on the module axis. Advances zF only when `advances`.
    auto addRing = [&](const std::string& name, double outHX, double outHY,
                       double apHalf, double t, G4Material* mat,
                       G4VisAttributes* vis, double zFront,
                       double ox = 0.0, double oy = 0.0, bool advances = false) {
        auto* outer = new G4Box(name + "_o", outHX, outHY, t/2);
        auto* inner = new G4Box(name + "_i", apHalf, apHalf, t/2 + 1.0*mm);
        auto* ring  = new G4SubtractionSolid(name, outer, inner,
                                             nullptr, G4ThreeVector(-ox, -oy, 0));
        auto* lv = new G4LogicalVolume(ring, mat, name);
        if (vis) lv->SetVisAttributes(vis);
        out.pieces.push_back({lv, G4ThreeVector(ox, oy, zFront + t/2),
                              advances ? t : 0.0, advances});
        grow(outHX, outHY, ox, oy);
        if (advances) zF += t;
        return lv;
    };

    const double tWinMy  = s.winMylar_um * um;
    const double tWinAl  = s.winAl_um    * um;
    const double tKap    = s.cathKapton_um * um;
    const double tCu     = s.cathCu_um   * um;
    const double tAmp    = s.amp_um      * um;
    const double tPaste  = s.paste_um    * um;
    const double tPCBKap = s.pcbKapton_um* um;
    const double tPCBCu  = s.pcbCu_um    * um;
    const double tRoh    = s.rohacell_mm * mm;

    // Mesh: woven effective-density slab (areal mass of a plain weave: fill
    // fraction π·d/(4p), slab thickness 2d) or the legacy solid slab.
    // Optical transparency (~(open/pitch)²) is NOT modelled — material budget only.
    double tMesh;
    G4Material* matMesh;
    if (s.meshWire_um > 0.0) {
        const double dWire = s.meshWire_um * um;
        const double pitch = dWire + s.meshOpen_um * um;
        tMesh = 2.0 * dWire;
        const double fill = M_PI * dWire / (4.0 * pitch);
        matMesh = G4Material::GetMaterial("MX17MeshSteelEff", false);
        if (!matMesh) {
            matMesh = new G4Material("MX17MeshSteelEff",
                                     m.steel->GetDensity() * fill, 1);
            matMesh->AddMaterial(m.steel, 1.0);
        }
    } else {
        tMesh   = s.meshSolid_um * um;
        matMesh = m.steel;
    }

    // Readout laminate FR4 thickness: fixed per-layer (legacy), or solved so
    // that mesh front → laminate back equals the CAD 1.70 mm board body.
    const int nPCBCu = s.pcbCuCoverage.empty()
        ? s.pcbNLayers : static_cast<int>(s.pcbCuCoverage.size());
    double tPCBFR4 = s.pcbFR4_um * um;
    if (s.pcbTotal_mm > 0.0) {
        const double lamTotal = s.pcbTotal_mm*mm - tMesh - tAmp - tPaste
                              - tPCBKap - nPCBCu*tPCBCu;
        tPCBFR4 = lamTotal / nPCBCu;
    }

    // Density-scaled copper for partially covered layers (p2-style): a
    // full-thickness sheet whose density is the gerber-measured area fraction.
    auto EffCu = [&](const std::string& name, double frac) -> G4Material* {
        if (frac >= 0.999) return m.cu;
        G4Material* mat = G4Material::GetMaterial(name, false);
        if (!mat) {
            mat = new G4Material(name, m.cu->GetDensity() * frac, 1);
            mat->AddMaterial(m.cu, 1.0);
        }
        return mat;
    };

    // ── Uniform (legacy) transverse face ─────────────────────────────────────
    if (s.Uniform()) {
        const double hx = s.uniformFaceX_mm * mm / 2;
        const double hy = s.uniformFaceY_mm * mm / 2;
        auto slab = [&](const std::string& n, double t, G4Material* mat,
                        G4VisAttributes* v) { return addBox(n, hx, hy, t, mat, v); };

        slab("GasWindow_Mylar",     tWinMy, m.mylar,  visMylar);
        slab("GasWindow_Al",        tWinAl, m.al,     visAl);
        slab("DriftCathode_Kapton", tKap,   m.kapton, visKapton);
        slab("DriftCathode_Cu",     tCu,    m.cu,     visCu);
        out.driftGasLV = slab("DriftGas", s.drift_mm*mm, m.gas, visDrift);
        slab("Micromesh",           tMesh,  matMesh,  visMesh);
        out.ampGasLV = slab("AmpGas", tAmp, m.gas, visAmp);
        slab("ResistivePaste",      tPaste, m.resPaste, visResPaste);
        out.mmDepth = zF;

        if (s.includeReadout) {
            slab("PCB_Kapton", tPCBKap, m.kapton, visKapton);
            for (int i = 1; i <= s.pcbNLayers; ++i) {
                slab("PCB_Cu_"  + std::to_string(i), tPCBCu,  m.cu,  visCu);
                slab("PCB_FR4_" + std::to_string(i), tPCBFR4, m.fr4, visFR4);
            }
            slab("PCB_Rohacell", tRoh, m.rohacell, visRohacell);
            if (s.backMylar_um > 0.0)
                slab("PCB_BackMylar", s.backMylar_um*um, m.mylar, visMylar);
            if (s.backAl_um > 0.0)
                slab("PCB_AlFoil", s.backAl_um*um, m.al, visAl);
            // Same expression shape as the historical pcbTotalZ sum so the
            // value is bit-identical to the pre-refactor build.
            double rd = tPCBKap + s.pcbNLayers*(tPCBCu + tPCBFR4) + tRoh;
            if (s.backMylar_um > 0.0) rd += s.backMylar_um*um;
            if (s.backAl_um > 0.0)    rd += s.backAl_um*um;
            out.readoutDepth = rd;
        }
        out.totalDepth = zF;
        return out;
    }

    // ── As-built transverse faces ────────────────────────────────────────────
    const double winH   = s.winFace_mm   * mm / 2;
    const double apH    = s.flangeAp_mm  * mm / 2;   // window free span
    const double cathH  = s.cathFace_mm  * mm / 2;
    const double frameApH = s.frameAp_mm * mm / 2;
    const double driftH = s.driftFace_mm * mm / 2;
    const double meshH  = s.meshFace_mm  * mm / 2;
    const double ampH   = s.ampFace_mm   * mm / 2;
    const double pcbH   = s.pcbFace_mm   * mm / 2;
    const double off    = s.pcbOffset_mm * mm;

    // 1) Window: flat sheet (no bulge) or flange annulus + terraced dome.
    const bool bulged = (s.bulgeSag_mm > 0.0 && s.flangeOuter_mm > 0.0);
    if (!bulged) {
        addBox("GasWindow_Mylar", winH, winH, tWinMy, m.mylar, visMylar);
        addBox("GasWindow_Al",    winH, winH, tWinAl, m.al,    visAl);
    } else {
        // Flat annulus over the flange face; the aperture part becomes the dome.
        addRing("GasWindow_Mylar", winH, winH, apH, tWinMy, m.mylar, visMylar,
                zF, 0, 0, true);
        addRing("GasWindow_Al",    winH, winH, apH, tWinAl, m.al,    visAl,
                zF, 0, 0, true);

        // Terraced dome rising upstream from the flange front plane (= zF).
        // Spherical-cap profile: σ_k = cos(kπ/2N), h_k = H·sin(kπ/2N). Every
        // z-parallel path crosses exactly one mylar (+Al) thickness; oblique
        // tracks through terrace risers are the residual approximation error.
        // The 0.1 µm Al is stacked on the UPSTREAM face of each terrace to
        // avoid slicing the gas steps — immaterial for the energy budget.
        const double zBase = zF;
        const int    N     = s.bulgeSteps;
        const double H     = s.bulgeSag_mm * mm;
        std::vector<double> sig(N + 1), h(N + 1);
        for (int k = 0; k <= N; ++k) {
            const double a = k * M_PI / (2.0 * N);
            sig[k] = std::max(std::cos(a), 0.04);
            h[k]   = H * std::sin(a);
        }
        for (int k = 0; k < N; ++k) {
            const double hw = apH * sig[k];
            const double dz = h[k+1] - h[k];
            // gas step
            auto* gbox = new G4Box("WindowBulgeGas", hw, hw, dz/2);
            auto* glv  = new G4LogicalVolume(gbox, m.gas, "WindowBulgeGas");
            glv->SetVisAttributes(visGasV);
            out.pieces.push_back({glv, G4ThreeVector(0, 0, zBase - h[k] - dz/2),
                                  0.0, false});
            // mylar terrace: annulus, full cap on the last step
            auto terrace = [&](const std::string& nm, double t, G4Material* mat,
                               G4VisAttributes* vis, double zTop) {
                G4VSolid* sol;
                if (k < N - 1) {
                    const double ihw = apH * sig[k+1];
                    sol = new G4SubtractionSolid(nm,
                              new G4Box(nm + "_o", hw, hw, t/2),
                              new G4Box(nm + "_i", ihw, ihw, t/2 + 1.0*mm));
                } else {
                    sol = new G4Box(nm, hw, hw, t/2);
                }
                auto* lv = new G4LogicalVolume(sol, mat, nm);
                lv->SetVisAttributes(vis);
                out.pieces.push_back({lv, G4ThreeVector(0, 0, zTop + t/2),
                                      0.0, false});
            };
            terrace("GasWindow_Mylar", tWinMy, m.mylar, visMylar,
                    zBase - h[k+1] - tWinMy);
            terrace("GasWindow_Al", tWinAl, m.al, visAl,
                    zBase - h[k+1] - tWinMy - tWinAl);
        }
        out.frontExtent = H;   // dome apex reaches z ≈ -H
    }

    // 2) Window flange ring + gas gap inside its aperture.
    if (s.flangeOuter_mm > 0.0) {
        const double tGap = s.flangeGap_mm * mm;
        addRing("WindowFlange_Al", s.flangeOuter_mm*mm/2, s.flangeOuter_mm*mm/2,
                apH, tGap, m.al, visFrame, zF);
        auto* gbox = new G4Box("WindowGapGas", apH, apH, tGap/2);
        auto* glv  = new G4LogicalVolume(gbox, m.gas, "WindowGapGas");
        glv->SetVisAttributes(visGasV);
        out.pieces.push_back({glv, G4ThreeVector(0, 0, zF + tGap/2), tGap, true});
        zF += tGap;
    }

    // 3) Drift cathode (copper-clad kapton; CAD omits the cladding).
    addBox("DriftCathode_Kapton", cathH, cathH, tKap, m.kapton, visKapton);
    const double zFrameTop = zF;   // frame top face = kapton back face
    addBox("DriftCathode_Cu", cathH, cathH, tCu, m.cu, visCu);

    // 4) Gas frame ring (the 30 mm spacer) + field-cage PCBs lining its
    //    aperture + drift gas. The Cu cladding protrudes into the frame
    //    height, so the gas is drift_mm − cathCu (cuInsideDrift).
    const double tFrame = s.drift_mm * mm;
    const double tDriftGas = s.cuInsideDrift ? tFrame - tCu : tFrame;
    if (s.frameOuter_mm > 0.0)
        addRing("GasFrame_Al", s.frameOuter_mm*mm/2, s.frameOuter_mm*mm/2,
                frameApH, tFrame, m.al, visFrame, zFrameTop);
    if (s.cageThick_mm > 0.0) {
        // 4 boards pinwheeled inside the aperture: length = aperture − thickness.
        const double tC   = s.cageThick_mm * mm;
        const double lenH = frameApH - tC/2;
        auto* box = new G4Box("FieldCagePCB_x", tC/2, lenH, tFrame/2);
        auto* lvx = new G4LogicalVolume(box, m.fr4, "FieldCagePCB");
        lvx->SetVisAttributes(visFR4);
        auto* boy = new G4Box("FieldCagePCB_y", lenH, tC/2, tFrame/2);
        auto* lvy = new G4LogicalVolume(boy, m.fr4, "FieldCagePCB");
        lvy->SetVisAttributes(visFR4);
        const double p = frameApH - tC/2;   // board centreline
        const double sft = tC/2;            // pinwheel slide
        const double zC = zFrameTop + tFrame/2;
        out.pieces.push_back({lvx, G4ThreeVector(+p, -sft, zC), 0.0, false});
        out.pieces.push_back({lvx, G4ThreeVector(-p, +sft, zC), 0.0, false});
        out.pieces.push_back({lvy, G4ThreeVector(+sft, +p, zC), 0.0, false});
        out.pieces.push_back({lvy, G4ThreeVector(-sft, -p, zC), 0.0, false});
    }
    out.driftGasLV = addBox("DriftGas", driftH, driftH, tDriftGas, m.gas, visDrift);

    // 5) Front-end (multi M1) cards: four cards flat on the drift-gas side of
    //    the readout-board edge (no standoff — the card laminate sits directly
    //    on the board), two per connector edge, tangentially centred on the
    //    board's connector-copper clusters at ±feTang. The laminate holds the
    //    six gerber Cu layers (density-scaled) with FR4 filling the rest.
    const double zMeshFront = zF;
    if (s.feThick_mm > 0.0) {
        const double tCard = s.feThick_mm * mm;
        const double feWH  = (s.feOuter_mm - s.feInner_mm) * mm / 2;
        const double feR   = (s.feOuter_mm + s.feInner_mm) * mm / 2;
        const double feLH  = s.feLen_mm * mm / 2;
        // sub-stack: [Cu + FR4 gap] per gerber layer, top (upstream) first
        const int nFeCu = static_cast<int>(s.feCuCoverage.size());
        const double tFeCu  = s.feCu_um * um;
        const double tFeFR4 = nFeCu > 0 ? (tCard - nFeCu*tFeCu) / nFeCu : tCard;
        struct Sub { double z0, t; G4Material* mat; };
        std::vector<Sub> subs;
        double zs = zMeshFront - tCard;
        for (int i = 0; i < std::max(nFeCu, 1); ++i) {
            if (nFeCu > 0) {
                subs.push_back({zs, tFeCu,
                                EffCu("MX17M1CuEff_" + std::to_string(i),
                                      s.feCuCoverage[i])});
                zs += tFeCu;
            }
            subs.push_back({zs, tFeFR4, m.fr4});
            zs += tFeFR4;
        }
        // 2 cards on the +x edge, 2 on the +y edge
        for (const auto& sub : subs) {
            auto* lvx = new G4LogicalVolume(
                new G4Box("FrontEndPCB_x", feWH, feLH, sub.t/2),
                sub.mat, "FrontEndPCB");
            auto* lvy = new G4LogicalVolume(
                new G4Box("FrontEndPCB_y", feLH, feWH, sub.t/2),
                sub.mat, "FrontEndPCB");
            auto* vis = (sub.mat == m.fr4) ? visFR4 : visCu;
            lvx->SetVisAttributes(vis);
            lvy->SetVisAttributes(vis);
            const double zc = sub.z0 + sub.t/2;
            for (double tang : {+s.feTang_mm * mm, -s.feTang_mm * mm}) {
                out.pieces.push_back({lvx, G4ThreeVector(feR, tang, zc),
                                      0.0, false});
                out.pieces.push_back({lvy, G4ThreeVector(tang, feR, zc),
                                      0.0, false});
            }
        }
        grow(feR + feWH, feR + feWH, 0, 0);
    }

    // 6) Micromesh / amplification gap / resistive layer (top of the readout
    //    board body; mesh front face = frame bottom face).
    addBox("Micromesh", meshH, meshH, tMesh, matMesh, visMesh);
    out.ampGasLV = addBox("AmpGas", ampH, ampH, tAmp, m.gas, visAmp);
    // Bulk pillars: one LV placed on the measured grid inside the amp gas, so
    // wherever the module is placed the scored gas is displaced correctly.
    // Material: polyimide coverlay (kapton) — the photoimageable film used by
    // the bulk process. checkOverlaps off for the grid: 7225 regular daughters
    // with pitch >> diameter, and the O(n²) sibling check would dominate
    // startup for no benefit.
    if (s.pillarD_mm > 0.0 && out.ampGasLV) {
        auto* pilSolid = new G4Tubs("BulkPillar", 0, s.pillarD_mm*mm/2,
                                    tAmp/2, 0, 360*deg);
        auto* pilLV = new G4LogicalVolume(pilSolid, m.kapton, "BulkPillar");
        pilLV->SetVisAttributes(new G4VisAttributes(G4Color(0.9, 0.88, 0.8, 0.9)));
        const double pitch = s.pillarPitch_mm * mm;
        const double o = (s.pillarN - 1) / 2.0;
        int copy = 0;
        for (int i = 0; i < s.pillarN; ++i)
            for (int j = 0; j < s.pillarN; ++j)
                new G4PVPlacement(nullptr,
                    G4ThreeVector((i - o)*pitch, (j - o)*pitch, 0),
                    pilLV, "BulkPillar", out.ampGasLV, false, copy++, false);
    }

    // Resistive layer: ESL strips 550 µm / 250 µm gap → density-scaled slab.
    const double pasteH = (s.pasteFace_mm > 0 ? s.pasteFace_mm : s.ampFace_mm)
                          * mm / 2;
    G4Material* matPaste = m.resPaste;
    if (s.pasteCoverage > 0.0 && s.pasteCoverage < 0.999) {
        matPaste = G4Material::GetMaterial("MX17ResistEff", false);
        if (!matPaste) {
            matPaste = new G4Material("MX17ResistEff",
                                      m.resPaste->GetDensity() * s.pasteCoverage, 1);
            matPaste->AddMaterial(m.resPaste, 1.0);
        }
    }
    addBox("ResistivePaste", pasteH, pasteH, tPaste, matPaste, visResPaste);
    out.mmDepth = zF;

    if (s.includeReadout) {
        // 6) Readout laminate — with pcbTotal_mm the FR4 filler is solved so
        //    mesh + amp + paste + laminate = the CAD 1.70 mm board body. With
        //    pcbCuCoverage set, the copper is one density-scaled sheet per
        //    physical gerber layer (L3 guard ring, L4 pads, L5 Y strips,
        //    L6 X strips, L7 fan-out; L8 carries no copper).
        addBox("PCB_Kapton", pcbH, pcbH, tPCBKap, m.kapton, visKapton, off, off);
        for (int i = 1; i <= nPCBCu; ++i) {
            G4Material* cuMat = s.pcbCuCoverage.empty()
                ? m.cu
                : EffCu("MX17PCBCuEff_L" + std::to_string(i + 2),
                        s.pcbCuCoverage[i - 1]);
            addBox("PCB_Cu_"  + std::to_string(i), pcbH, pcbH, tPCBCu,
                   cuMat, visCu, off, off);
            addBox("PCB_FR4_" + std::to_string(i), pcbH, pcbH, tPCBFR4,
                   m.fr4, visFR4, off, off);
        }

        // 7) Rohacell + aluminized-mylar back foil + support plate ring
        //    (402 mm aperture concentric with the active area).
        addBox("PCB_Rohacell", pcbH, pcbH, tRoh, m.rohacell, visRohacell, off, off);
        if (s.backMylar_um > 0.0)
            addBox("PCB_BackMylar", pcbH, pcbH, s.backMylar_um*um,
                   m.mylar, visMylar, off, off);
        if (s.backAl_um > 0.0)
            addBox("PCB_AlFoil", pcbH, pcbH, s.backAl_um*um,
                   m.al, visAl, off, off);
        if (s.plate_mm > 0.0)
            addRing("SupportPlate_Al", pcbH, pcbH, s.plateAp_mm*mm/2,
                    s.plate_mm*mm, m.al, visAl, zF, off, off, true);
    }
    out.readoutDepth = zF - out.mmDepth;
    out.totalDepth   = zF;
    return out;
}

}  // namespace MX17
