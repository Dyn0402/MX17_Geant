#!/usr/bin/env python3
"""
Plots of the as-built MX17 Geant4 model (scripts/model/mx17_model.py).

    python scripts/model/plot_mx17_model.py             # all figures
    python scripts/model/plot_mx17_model.py --only 3d   # xsec | 3d | status

Outputs to design/figures/ by default:
    mx17_stack_xsec.png    cross-section at y=0 (true scale + zooms + stack)
    mx17_3d_overview.png   assembled 3D views (z exaggerated where noted)
    mx17_3d_exploded.png   exploded 3D view
    mx17_stack_status.png  stack schematic with the source status of each number
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# On machines where an old distro matplotlib coexists with a pip user install,
# the distro's mpl_toolkits shadows the user install's namespace portion.
import mpl_toolkits
_user_mt = os.path.join(os.path.dirname(os.path.dirname(matplotlib.__file__)),
                        "mpl_toolkits")
if os.path.isdir(_user_mt) and _user_mt not in mpl_toolkits.__path__:
    mpl_toolkits.__path__.insert(0, _user_mt)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.projections import register_projection
register_projection(Axes3D)

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)

import mx17_model as M


# ─────────────────────────────────────────────────────────────────────────────
# 3D helpers (rectangular prisms and square rings)
# ─────────────────────────────────────────────────────────────────────────────

def rect_outline(hx, hy, ox=0.0, oy=0.0, k=16):
    """Closed rectangle boundary resampled to ~4k points."""
    c = [(ox-hx, oy-hy), (ox+hx, oy-hy), (ox+hx, oy+hy), (ox-hx, oy+hy)]
    pts = []
    for i in range(4):
        (x0, y0), (x1, y1) = c[i], c[(i+1) % 4]
        pts.extend((x0 + (x1-x0)*j/k, y0 + (y1-y0)*j/k) for j in range(k))
    return np.array(pts)


def prism_faces(outline, z0, z1):
    bot = np.column_stack([outline, np.full(len(outline), z0)])
    top = np.column_stack([outline, np.full(len(outline), z1)])
    faces = [bot, top]
    for i in range(len(outline)):
        j = (i + 1) % len(outline)
        faces.append(np.array([bot[i], bot[j], top[j], top[i]]))
    return faces


def ring_faces(outer, inner, z0, z1):
    n = len(outer)
    faces = []
    for i in range(n):
        j = (i + 1) % n
        for zz in (z0, z1):
            faces.append(np.array([[*outer[i], zz], [*outer[j], zz],
                                   [*inner[j], zz], [*inner[i], zz]]))
        faces.append(np.array([[*outer[i], z0], [*outer[j], z0],
                               [*outer[j], z1], [*outer[i], z1]]))
        faces.append(np.array([[*inner[i], z0], [*inner[j], z0],
                               [*inner[j], z1], [*inner[i], z1]]))
    return faces


def add_layer_3d(ax, L, z0, z1):
    outer = rect_outline(L.hx, L.hy, L.ox, L.oy)
    if L.hole is not None:
        faces = ring_faces(outer, rect_outline(L.hole, L.hole), z0, z1)
    else:
        faces = prism_faces(outer, z0, z1)
    ax.add_collection3d(Poly3DCollection(
        faces, facecolor=L.color, alpha=min(max(L.alpha, 0.35), 0.95),
        edgecolor="k", linewidths=0.15))


def draw_model_3d(ax, zexag=1.0, explode=0.0):
    layers, side, windows, key_z = M.build_stack()

    groups = {"GasWindow": -2, "WindowBulge": -2, "Bulge": -2,
              "WindowFlange": -1, "WindowGap": -1,
              "DriftCathode": 0, "GasFrame": 1, "FieldCage": 1, "DriftGas": 1,
              "FrontEndPCB": 1,
              "Micromesh": 2, "AmpGas": 2, "ResistivePaste": 2,
              "PCB_Kapton": 3, "PCB_Cu": 3, "PCB_FR4": 3,
              "PCB_Rohacell": 4, "PCB_BackMylar": 4, "PCB_AlFoil": 4,
              "SupportPlate": 5}

    def grp(name):
        for k in sorted(groups, key=len, reverse=True):
            if name.startswith(k):
                return groups[k]
        return 0

    def Z(z, g):
        return z * zexag + g * explode

    skip = {"WindowGapGas"}  # flat gas slab clutter; the flange shows the gap
    for L in layers + side + windows:
        if L.name in skip or (L.material == "gas" and "Bulge" in L.name and explode > 0):
            continue
        alpha = 0.10 if "BulgeGas" in L.name else L.alpha
        g = grp(L.name)
        Lc = L
        Lc.alpha = alpha
        add_layer_3d(ax, Lc, Z(L.z0, g), Z(L.z0 + L.t, g))

    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_xlim(-260, 260); ax.set_ylim(-260, 260)
    return key_z


def fig_3d(outdir):
    fig = plt.figure(figsize=(19, 6.8))
    views = [("front / window side", 28, -125, 4),
             ("edge-on — window bulge", 2, -90, 4),
             ("back / support-plate side", -30, 55, 4)]
    for i, (title, elev, azim, zex) in enumerate(views, 1):
        ax = fig.add_subplot(1, 3, i, projection="3d")
        draw_model_3d(ax, zexag=zex)
        ax.set_zlim(-60, 230)
        ax.set_box_aspect((520, 520, 300))
        ax.view_init(elev=elev, azim=azim)
        ax.invert_zaxis()   # beam travels +z into the page; show window on top
        ax.set_title(f"{title}  (z ×{zex})", fontsize=11)
        ax.set_zlabel(""); ax.set_zticks([])
        if "edge-on" in title:
            ax.set_yticks([]); ax.set_ylabel("")
    fig.suptitle("MX17 as-built Geant4 model — 3D (z exaggerated for visibility)",
                 fontsize=13)
    fig.tight_layout()
    out = os.path.join(outdir, "mx17_3d_overview.png")
    fig.savefig(out, dpi=160); plt.close(fig)
    print("wrote", out)

    fig = plt.figure(figsize=(10.5, 11))
    ax = fig.add_subplot(projection="3d")
    draw_model_3d(ax, zexag=4, explode=55)
    ax.set_zlim(-160, 560)
    ax.set_box_aspect((520, 520, 620))
    ax.view_init(elev=14, azim=-110)
    ax.invert_zaxis()
    ax.set_zticks([])
    for ztxt, lab in [(-125, "bulged window\n(60 µm mylar + Al, terraced dome)"),
                      (-40, "window flange ring + 5 mm gas gap"),
                      (85, "drift: Cu-clad kapton cathode, 30 mm frame,\n"
                           "field cage, M1 front-end cards"),
                      (200, "mesh (woven, eff. density) + amp gap + paste"),
                      (280, "readout laminate (1.70 mm body total)"),
                      (360, "rohacell + aluminized-mylar back foil"),
                      (440, "8 mm Al support plate (402 mm aperture)")]:
        ax.text(300, -240, ztxt, lab, fontsize=9.5, ha="left")
    ax.set_title("MX17 as-built model — exploded (z ×4, groups separated)",
                 fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(left=-0.15, right=0.78, top=0.98, bottom=0.02)
    out = os.path.join(outdir, "mx17_3d_exploded.png")
    fig.savefig(out, dpi=160); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Cross-section at y = 0
# ─────────────────────────────────────────────────────────────────────────────

def xsec_intervals(L):
    """x-intervals of layer L on the y=0 plane."""
    if abs(L.oy) >= L.hy:
        return []
    lo, hi = L.ox - L.hx, L.ox + L.hx
    if L.hole is None:
        return [(lo, hi)]
    h = L.hole
    out = []
    if lo < -h: out.append((lo, -h))
    if hi > h:  out.append((h, hi))
    return out


def fig_xsec(outdir):
    layers, side, windows, key_z = M.build_stack()
    items = layers + side + windows

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(19, 6.4), width_ratios=[2.2, 1.0, 1.1])

    for ax in (ax1, ax2):
        for L in items:
            for x0, x1 in xsec_intervals(L):
                ax.add_patch(Rectangle((x0, L.z0), x1 - x0, L.t,
                                       facecolor=L.color,
                                       alpha=max(L.alpha, 0.35),
                                       edgecolor="k", lw=0.2))
        ax.invert_yaxis()  # beam from the top of the plot

    ax1.annotate("", xy=(0, -9.5), xytext=(0, -19.5),
                 arrowprops=dict(arrowstyle="-|>", color="k", lw=1.6))
    ax1.text(6, -15.5, "beam", fontsize=10)
    ax1.set_xlim(-260, 260); ax1.set_ylim(52, -22)
    ax1.set_xlabel("x [mm]"); ax1.set_ylabel("z [mm]")
    ax1.set_title("True scale — bulged window, flange, frame, plate aperture")
    ax1.grid(alpha=0.25, lw=0.4)

    zs = {L.name: (L.z0, L.t) for L in layers}
    z_mesh = zs["Micromesh"][0]
    z_end = zs["PCB_FR4_4"][0] + zs["PCB_FR4_4"][1]
    ax2.set_xlim(-30, 30); ax2.set_ylim(z_end + 0.1, z_mesh - 0.1)
    ax2.set_xlabel("x [mm]")
    ax2.set_title("Zoom: mesh → laminate (1.70 mm board body)")
    ax2.grid(alpha=0.25, lw=0.4)
    for name, lab in [("Micromesh", "mesh (38 µm eff.)"),
                      ("AmpGas", "amp gap (150 µm)"),
                      ("ResistivePaste", "paste (100 µm)"),
                      ("PCB_Kapton", "kapton (50 µm)"),
                      ("PCB_FR4_2", "4× Cu 26 µm / FR4 314.5 µm")]:
        z0, t = zs[name]
        ax2.annotate(lab, xy=(30, z0 + t/2), xytext=(32, z0 + t/2),
                     fontsize=8, va="center", annotation_clip=False)

    ax3.set_title("Stack (not to scale)")
    y = 0
    for L in layers:
        ax3.add_patch(Rectangle((0, y), 1, 1, facecolor=L.color,
                                alpha=max(L.alpha, 0.4), edgecolor="k", lw=0.4))
        t_lab = f"{L.t*1000:.1f} µm" if L.t < 1 else f"{L.t:.0f} mm"
        ax3.text(1.05, y + 0.5, f"{L.name}   ({t_lab})", va="center", fontsize=8.5)
        y += 1
    ax3.text(0.0, y + 1.2,
             f"window: 60 µm aluminized mylar, sag {M.BULGE_SAG:.0f} mm\n"
             f"rings: flange 440/400.2×5, frame 440/410×30\n"
             f"plate aperture 402 mm; plates offset (+15,+15)",
             fontsize=8.5, va="top")
    ax3.set_xlim(0, 4.6); ax3.set_ylim(-0.5, y + 5)
    ax3.invert_yaxis(); ax3.axis("off")

    fig.suptitle("MX17 as-built Geant4 model — cross-section at y = 0 "
                 "(z = 0 at the window plane, beam along +z)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(outdir, "mx17_stack_xsec.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Stack schematic with the source status of every number
# ─────────────────────────────────────────────────────────────────────────────

C_CAD, C_HW, C_ASSUME, C_GUESS = "#1565c0", "#2e7d32", "#ef6c00", "#c62828"


def fig_status(outdir):
    fig, ax = plt.subplots(figsize=(14.5, 11))
    X0, X1 = 1.3, 4.3

    rows = [
        ("window",  "gas window — 60 µm mylar, 440², bulged",
         "60 µm CAD; bulge 8 mm from ~3 mbar Hencky (p unknown)", 0.95, "#63b8d8", C_ASSUME),
        ("winal",   "window aluminization — 0.1 µm Al, drift side",
         "hardware knowledge; side assumed", 0.4, "#b0b0b0", C_ASSUME),
        ("flange",  "window flange ring + 5 mm gas gap (ap 400.2)",
         "CAD", 0.9, "#dff2fa", C_CAD),
        ("cath",    "drift cathode — kapton + 9 µm Cu",
         "CAD; Cu cladding from hardware (absent in CAD)", 0.6, "#e6b300", C_HW),
        ("drift",   "DRIFT GAS — 30.00 mm (frame = spacer)",
         "CAD, exact", 1.3, "#b3ccff", C_CAD),
        ("cage",    "field cage — 4× 1.62 mm FR4",
         "CAD", 0.5, "#9fd49f", C_CAD),
        ("mesh",    "micromesh — woven SS 19/48 µm",
         "weave spec is a P2-like placeholder", 0.55, "#9e9e9e", C_GUESS),
        ("amp",     "amplification gap — 150 µm",
         "both sims' historical value; 128 µm bulk standard unresolved", 0.7, "#ffc4c4", C_ASSUME),
        ("paste",   "resistive paste — 100 µm",
         "historical value; not in CAD or gerbers", 0.45, "#555555", C_GUESS),
        ("pcb",     "readout board — 1.70 mm body, 470²",
         "CAD total; internal Cu/FR4 stackup unknown, FR4 filler solved", 0.9, "#cc6619", C_ASSUME),
        ("roh",     "rohacell — 5 mm, 470²", "CAD", 0.7, "#e8e89a", C_CAD),
        ("foil",    "back foil — 25 µm aluminized mylar",
         "existence: hardware; thickness assumed", 0.45, "#63b8d8", C_ASSUME),
        ("plate",   "Al support plate — 8 mm, 402² aperture",
         "CAD (aperture verified in STEP face loops)", 1.0, "#9a9aa2", C_CAD),
    ]

    y = 12.4
    for key, lab, note, h, fc, st in rows:
        y -= h + 0.07
        ax.add_patch(Rectangle((X0, y), X1 - X0, h, facecolor=fc,
                               edgecolor=st, lw=2.2))
        ax.text((X0+X1)/2, y + h/2, lab, ha="center", va="center",
                fontsize=10.5 if h >= 0.55 else 9.0,
                fontweight="bold" if key == "drift" else "normal")
        ax.text(X1 + 0.25, y + h/2, note, ha="left", va="center",
                fontsize=9.5, color=st)

    ax.annotate("", xy=((X0+X1)/2, 13.1), xytext=((X0+X1)/2, 14.0),
                arrowprops=dict(arrowstyle="-|>", lw=2.5, color="k"))
    ax.text((X0+X1)/2 + 0.12, 13.55, "beam", fontsize=11, va="center")

    for i, (c, lab) in enumerate([
            (C_CAD,    "mechanical CAD / gerbers"),
            (C_HW,     "hardware knowledge (confirmed)"),
            (C_ASSUME, "assumed — please check"),
            (C_GUESS,  "placeholder — number needed")]):
        ax.add_patch(Rectangle((0.10, 0.85 - 0.55*i), 0.30, 0.30,
                               facecolor="white", edgecolor=c, lw=2.2))
        ax.text(0.50, 1.0 - 0.55*i, lab, fontsize=10, va="center")

    ax.set_xlim(-0.3, 11.6)
    ax.set_ylim(-1.0, 14.4)
    ax.axis("off")
    ax.set_title("MX17 as-built stack (not to scale) — where every number comes from\n"
                 "(open items tracked in design/NEEDED_INPUTS.md)",
                 fontsize=13, pad=12)
    fig.tight_layout()
    out = os.path.join(outdir, "mx17_stack_status.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join(REPO, "design", "figures"))
    ap.add_argument("--only", choices=["xsec", "3d", "status"], default=None)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if args.only in (None, "xsec"):
        fig_xsec(args.out_dir)
    if args.only in (None, "3d"):
        fig_3d(args.out_dir)
    if args.only in (None, "status"):
        fig_status(args.out_dir)


if __name__ == "__main__":
    main()
