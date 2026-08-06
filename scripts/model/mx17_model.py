#!/usr/bin/env python3
"""
Python mirror of the MX17 as-built module geometry
(shared/MX17ModuleGeometry.hh, AsBuiltSpec). Used by plot_mx17_model.py to
draw the model without needing Geant4.

KEEP IN SYNC with shared/MX17ModuleGeometry.hh — that header is the single
source of truth consumed by both MX17_Geant and MX17_Full_Geant;
design/GEOMETRY_FROM_CAD.md records where each number comes from and
design/NEEDED_INPUTS.md which are still assumptions.

Coordinates: module local frame, mm everywhere. x,y transverse with the
active-area axis at (0,0); z = 0 at the upstream face of the flat window
sheet, +z downstream. The bulged window rises to z < 0.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# ── AsBuiltSpec values (shared/MX17ModuleGeometry.hh) ───────────────────────
WIN_MYLAR   = 0.060    # CAD
WIN_AL      = 0.0001   # aluminization, drift side (assumed side)
WIN_FACE    = 440.0
FLANGE_GAP  = 5.0      # window flange thickness = window→cathode gas gap
FLANGE_AP   = 400.2    # window free span
FLANGE_OUT  = 440.0
BULGE_SAG   = 8.0      # Hencky ~3 mbar, 400.2 mm span, 60 um foil (assumed p)
BULGE_N     = 6

CATH_KAPTON = 0.050
CATH_CU     = 0.009    # cladding, drift side (assumed side; CAD omits it)
CATH_FACE   = 406.0

DRIFT       = 30.0     # kapton back → mesh front (frame is the spacer)
FRAME_OUT   = 440.0
FRAME_AP    = 410.0
CAGE_T      = 1.62     # 4× drift field PCB lining the aperture
DRIFT_FACE  = FRAME_AP - 2 * CAGE_T

MESH_WIRE   = 0.019    # P2-like weave — placeholder (NEEDED_INPUTS)
MESH_OPEN   = 0.048
MESH_T      = 2 * MESH_WIRE
MESH_FILL   = math.pi * MESH_WIRE / (4 * (MESH_WIRE + MESH_OPEN))
MESH_FACE   = 410.0

AMP         = 0.150
PASTE       = 0.100
AMP_FACE    = 410.0

PCB_TOTAL   = 1.70     # CAD single-body readout board (mesh → laminate back)
PCB_KAPTON  = 0.050
PCB_CU      = 0.026
PCB_N       = 4
PCB_FR4     = (PCB_TOTAL - MESH_T - AMP - PASTE - PCB_KAPTON - PCB_N * PCB_CU) / PCB_N
PCB_FACE    = 470.0
PCB_OFF     = 15.0     # 470 mm plates offset (+15,+15) from the active axis

ROHACELL    = 5.0
BACK_MYLAR  = 0.025    # assumed (NEEDED_INPUTS)
BACK_AL     = 0.0001
PLATE       = 8.0
PLATE_AP    = 402.0    # through-aperture, concentric with the active area

FE_T        = 6.6      # multi M1 cards: radial thickness
FE_LEN      = 180.0
FE_H        = 41.0     # rising upstream from the mesh plane
FE_RING     = 227.5
FE_TANG     = (75.0, -108.0)   # tangential centres on each connector edge

ACTIVE      = 399.36   # readout active area (gerbers), centred


@dataclass
class Layer:
    name: str
    z0: float                    # upstream face
    t: float
    material: str
    hx: float                    # half-extents of the outer box
    hy: float
    ox: float = 0.0              # centre offset of the outer box
    oy: float = 0.0
    hole: float | None = None    # half-width of a square aperture at (0,0)
    color: str = "0.5"
    alpha: float = 0.85


def terrace_profile(H: float, n: int = BULGE_N):
    sig = [max(math.cos(k * math.pi / (2 * n)), 0.04) for k in range(n + 1)]
    h = [H * math.sin(k * math.pi / (2 * n)) for k in range(n + 1)]
    return sig, h


def build_stack(bulge: float = BULGE_SAG):
    """Returns (layers, side, windows, key_z); all entries are Layer.

    layers  — the z-advancing slab chain
    side    — rings / field cage / front-end cards (fixed z, not advancing)
    windows — bulge terraces (gas steps + mylar/Al caps)
    """
    layers: list[Layer] = []
    side: list[Layer] = []
    windows: list[Layer] = []
    z = 0.0

    def add(name, t, mat, hx, hy, color, alpha=0.85, ox=0.0, oy=0.0, hole=None):
        nonlocal z
        layers.append(Layer(name, z, t, mat, hx, hy, ox, oy, hole, color, alpha))
        z += t

    apH = FLANGE_AP / 2
    bulged = bulge > 0

    # window sheet (annulus when bulged — the aperture part becomes the dome)
    add("GasWindow_Mylar", WIN_MYLAR, "mylar", WIN_FACE/2, WIN_FACE/2,
        color="#63b8d8", alpha=0.9, hole=apH if bulged else None)
    add("GasWindow_Al", WIN_AL, "Al", WIN_FACE/2, WIN_FACE/2,
        color="#b0b0b0", hole=apH if bulged else None)
    z_flange = z

    side.append(Layer("WindowFlange_Al", z_flange, FLANGE_GAP, "Al",
                      FLANGE_OUT/2, FLANGE_OUT/2, hole=apH,
                      color="#c9c9cf", alpha=0.9))
    add("WindowGapGas", FLANGE_GAP, "gas", apH, apH, color="#8dd8f0", alpha=0.15)

    add("DriftCathode_Kapton", CATH_KAPTON, "kapton", CATH_FACE/2, CATH_FACE/2,
        color="#e6b300")
    z_frame_top = z
    add("DriftCathode_Cu", CATH_CU, "Cu", CATH_FACE/2, CATH_FACE/2,
        color="#cc6619")

    side.append(Layer("GasFrame_Al", z_frame_top, DRIFT, "Al",
                      FRAME_OUT/2, FRAME_OUT/2, hole=FRAME_AP/2,
                      color="#c9c9cf", alpha=0.9))
    # field cage: 4 pinwheeled boards lining the aperture (drawn as one ring)
    side.append(Layer("FieldCagePCB", z_frame_top, DRIFT, "FR4",
                      FRAME_AP/2, FRAME_AP/2, hole=DRIFT_FACE/2,
                      color="#339933", alpha=0.9))
    add("DriftGas", DRIFT - CATH_CU, "gas", DRIFT_FACE/2, DRIFT_FACE/2,
        color="#3380ff", alpha=0.30)
    z_mesh = z

    # front-end M1 cards: two per connector edge, rising upstream from z_mesh
    for tang in FE_TANG:
        side.append(Layer("FrontEndPCB", z_mesh - FE_H, FE_H, "FR4",
                          FE_T/2, FE_LEN/2, ox=FE_RING, oy=tang,
                          color="#339933", alpha=0.9))
        side.append(Layer("FrontEndPCB", z_mesh - FE_H, FE_H, "FR4",
                          FE_LEN/2, FE_T/2, ox=tang, oy=FE_RING,
                          color="#339933", alpha=0.9))

    add("Micromesh", MESH_T, f"steel(x{MESH_FILL:.2f})", MESH_FACE/2, MESH_FACE/2,
        color="#808080")
    add("AmpGas", AMP, "gas", AMP_FACE/2, AMP_FACE/2, color="#ff4d4d", alpha=0.35)
    add("ResistivePaste", PASTE, "resistive", AMP_FACE/2, AMP_FACE/2,
        color="#333333")

    add("PCB_Kapton", PCB_KAPTON, "kapton", PCB_FACE/2, PCB_FACE/2,
        color="#e6b300", ox=PCB_OFF, oy=PCB_OFF)
    for i in range(1, PCB_N + 1):
        add(f"PCB_Cu_{i}", PCB_CU, "Cu", PCB_FACE/2, PCB_FACE/2,
            color="#cc6619", ox=PCB_OFF, oy=PCB_OFF)
        add(f"PCB_FR4_{i}", PCB_FR4, "FR4", PCB_FACE/2, PCB_FACE/2,
            color="#339933", ox=PCB_OFF, oy=PCB_OFF)
    z_pcb_end = z

    add("PCB_Rohacell", ROHACELL, "rohacell", PCB_FACE/2, PCB_FACE/2,
        color="#e8e89a", ox=PCB_OFF, oy=PCB_OFF)
    add("PCB_BackMylar", BACK_MYLAR, "mylar", PCB_FACE/2, PCB_FACE/2,
        color="#63b8d8", alpha=0.9, ox=PCB_OFF, oy=PCB_OFF)
    add("PCB_AlFoil", BACK_AL, "Al", PCB_FACE/2, PCB_FACE/2,
        color="#b0b0b0", ox=PCB_OFF, oy=PCB_OFF)
    add("SupportPlate_Al", PLATE, "Al", PCB_FACE/2, PCB_FACE/2,
        color="#9a9aa2", alpha=0.95, ox=PCB_OFF, oy=PCB_OFF, hole=PLATE_AP/2)
    z_back = z

    # terraced window dome over the aperture, rising upstream from z_flange
    if bulged:
        sig, h = terrace_profile(bulge)
        for k in range(BULGE_N):
            hw = apH * sig[k]
            dz = h[k+1] - h[k]
            windows.append(Layer(f"WindowBulgeGas{k}", z_flange - h[k] - dz, dz,
                                 "gas", hw, hw, color="#8dd8f0", alpha=0.15))
            hole = apH * sig[k+1] if k < BULGE_N - 1 else None
            windows.append(Layer(f"BulgeMylar{k}",
                                 z_flange - h[k+1] - WIN_MYLAR, WIN_MYLAR,
                                 "mylar", hw, hw, hole=hole,
                                 color="#63b8d8", alpha=0.9))
            windows.append(Layer(f"BulgeAl{k}",
                                 z_flange - h[k+1] - WIN_MYLAR - WIN_AL, WIN_AL,
                                 "Al", hw, hw, hole=hole,
                                 color="#b0b0b0", alpha=0.9))

    key_z = dict(window_plane=0.0, flange_front=z_flange,
                 mesh_front=z_mesh, pcb_end=z_pcb_end, module_back=z_back,
                 z_min=-bulge, z_max=z_back)
    return layers, side, windows, key_z


if __name__ == "__main__":
    layers, side, windows, key_z = build_stack()
    print(f"{'layer':22s} {'z0 [mm]':>10s} {'t [mm]':>9s}  face  material")
    for L in layers:
        face = f"{2*L.hx:.1f}x{2*L.hy:.1f}"
        off = f" @({L.ox:+.0f},{L.oy:+.0f})" if L.ox or L.oy else ""
        print(f"{L.name:22s} {L.z0:10.4f} {L.t:9.4f}  {face}{off}  {L.material}")
    for L in side:
        hol = f" ap {2*L.hole:.1f}" if L.hole else ""
        print(f"{L.name:22s} {L.z0:10.4f} {L.t:9.4f}  ring{hol}  {L.material}")
    print(f"\nwindows: {len(windows)} bulge terrace volumes, sag {BULGE_SAG} mm")
    for k, v in key_z.items():
        print(f"  {k:14s} z = {v:9.4f} mm")
    print(f"\nmesh fill {MESH_FILL:.3f}; PCB FR4 filler {PCB_FR4*1000:.1f} um x {PCB_N}")
