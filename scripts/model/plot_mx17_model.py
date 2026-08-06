#!/usr/bin/env python3
"""
Plots of the as-built MX17 Geant4 model.

The cross-section, 3D and plan figures render the TRUE constructed geometry
from design/mx17_geometry.json, produced by the simulation itself:

    ./mm_sim -m sr90 --dump-geometry design/mx17_geometry.json

Regenerate that dump after any geometry change. The board/peel figures render
the production gerbers directly; scripts/model/mx17_model.py remains as a
quick numeric reference but the figures no longer depend on it.

    python scripts/model/plot_mx17_model.py             # all figures
    python scripts/model/plot_mx17_model.py --only 3d   # board | xsec | 3d | status

Outputs to design/figures/ by default:
    mx17_board_copper.png  readout copper rendered from the production gerbers
    mx17_board_peel.png    close-up with layers peeled back band by band
    mx17_plan_views.png    plan views from upstream / downstream for alignment
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
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Circle, Rectangle

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
sys.path.insert(0, os.path.join(REPO, "scripts", "gerber"))

import json
from dataclasses import dataclass

import mx17_model as M
from gerber_outline import parse

GERBER_DIR = os.path.join(REPO, "design", "gerbers", "readout_pcb")
CU_COLOR = "#b87333"
DUMP = os.path.join(REPO, "design", "mx17_geometry.json")

# lv name -> (color, alpha); BulkPillar is skipped in the drawings.
STYLE = {
    "GasWindow_Mylar":     ("#63b8d8", 0.90),
    "GasWindow_Al":        ("#b0b0b0", 0.90),
    "WindowBulgeGas":      ("#8dd8f0", 0.15),
    "WindowFlange_Al":     ("#c9c9cf", 0.90),
    "WindowGapGas":        ("#8dd8f0", 0.15),
    "DriftCathode_Kapton": ("#e6b300", 0.85),
    "DriftCathode_Cu":     ("#cc6619", 0.85),
    "GasFrame_Al":         ("#c9c9cf", 0.90),
    "FieldCagePCB":        ("#339933", 0.90),
    "DriftGas":            ("#3380ff", 0.30),
    "FrontEndPCB":         ("#2e7d32", 0.90),
    "Micromesh":           ("#808080", 0.85),
    "AmpGas":              ("#ff4d4d", 0.35),
    "ResistivePaste":      ("#333333", 0.85),
    "ResistLayer":         ("#333333", 0.85),   # ESL envelope (patterned build)
    "PCB_Kapton":          ("#e6b300", 0.85),
    "PCB_Cu":              ("#cc6619", 0.85),
    "PCB_FR4":             ("#339933", 0.85),
    "PCB_Rohacell":        ("#e8e89a", 0.85),
    "PCB_BackMylar":       ("#63b8d8", 0.90),
    "PCB_AlFoil":          ("#b0b0b0", 0.85),
    "SupportPlate_Al":     ("#9a9aa2", 0.95),
}
DROP_PREFIX = ("AirGap", "ScintWall", "PlasticScint", "LS_", "LiqScint",
               "BulkPillar", "World", "He3", "ResistCell")
# Internals of the zoned / patterned readout copper layers. Each of these sits
# INSIDE its PCB_Cu_<n> parent at the same z, so drawing them too would stack a
# second copy of the layer on top of itself. The parent already carries the
# layer at full board extent; the copper micro-pattern is the peel figure's job.
DROP_SUFFIX = ("_Win", "_Col", "_Cell", "_Pad", "_Dot", "_Bus", "_Stub")


@dataclass
class Vol:
    name: str
    z0: float
    t: float
    hx: float
    hy: float
    ox: float
    oy: float
    hole: tuple | None      # (hcx, hcy, hhx, hhy) of a rectangular aperture
    color: str
    alpha: float


def style_of(name):
    for k in sorted(STYLE, key=len, reverse=True):
        if name.startswith(k):
            return STYLE[k]
    return ("0.6", 0.8)


def load_dump(path=DUMP):
    """The constructed Geant4 module, z-shifted so the window plane is 0."""
    if not os.path.exists(path):
        sys.exit(f"{path} not found — regenerate with\n"
                 f"    ./mm_sim -m sr90 --dump-geometry {path}")
    raw = json.load(open(path))
    raw = [v for v in raw if not v["lv"].startswith(DROP_PREFIX)
           and not v["lv"].endswith(DROP_SUFFIX)]
    # In the patterned build the ESL is 515 replicated 0.55 mm strips inside a
    # gas-filled "ResistLayer" envelope. The envelope is the thing to draw as a
    # stack layer; the strip prototype would otherwise show up as a 0.55 mm
    # sliver spanning the whole board thickness.
    if any(v["lv"] == "ResistLayer" for v in raw):
        raw = [v for v in raw if v["lv"] != "ResistivePaste"]
    vols = []
    for v in raw:
        s = v["solid"]
        x, y, z = v["pos"]
        if s["type"] == "box":
            vols.append((v["lv"], x, y, z, s["hx"], s["hy"], s["hz"], None))
        elif s["type"] == "sub" and s["a"]["type"] == "box":
            b = s["b"]
            hole = (x + s["bpos"][0], y + s["bpos"][1], b["hx"], b["hy"])
            vols.append((v["lv"], x, y, z,
                         s["a"]["hx"], s["a"]["hy"], s["a"]["hz"], hole))
        # tubs (pillars) and others: skipped in drawings
    # window plane: front of the widest GasWindow_Mylar piece
    win = max((v for v in vols if v[0] == "GasWindow_Mylar"),
              key=lambda v: v[4])
    zref = win[3] - win[6]
    out = []
    for name, x, y, z, hx, hy, hz, hole in vols:
        c, a = style_of(name)
        out.append(Vol(name, z - hz - zref, 2*hz, hx, hy, x, y, hole, c, a))
    return out


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
        hcx, hcy, hhx, hhy = L.hole
        faces = ring_faces(outer, rect_outline(hhx, hhy, hcx, hcy), z0, z1)
    else:
        faces = prism_faces(outer, z0, z1)
    ax.add_collection3d(Poly3DCollection(
        faces, facecolor=L.color, alpha=min(max(L.alpha, 0.35), 0.95),
        edgecolor="k", linewidths=0.15))


GROUPS = {"GasWindow": -2, "WindowBulge": -2,
          "WindowFlange": -1, "WindowGap": -1,
          "DriftCathode": 0, "GasFrame": 1, "FieldCage": 1, "DriftGas": 1,
          "Micromesh": 2, "AmpGas": 2, "ResistivePaste": 2, "FrontEndPCB": 2,
          "PCB_Kapton": 3, "PCB_Cu": 3, "PCB_FR4": 3,
          "PCB_Rohacell": 4, "PCB_BackMylar": 4, "PCB_AlFoil": 4,
          "SupportPlate": 5}


def grp(name):
    for k in sorted(GROUPS, key=len, reverse=True):
        if name.startswith(k):
            return GROUPS[k]
    return 0


def draw_model_3d(ax, zexag=1.0, explode=0.0):
    vols = load_dump()

    def Z(z, g):
        return z * zexag + g * explode

    skip = {"WindowGapGas"}  # flat gas slab clutter; the flange shows the gap
    for L in vols:
        if L.name in skip or (explode > 0 and L.name == "WindowBulgeGas"):
            continue
        alpha = 0.10 if L.name == "WindowBulgeGas" else L.alpha
        g = grp(L.name)
        L2 = Vol(L.name, L.z0, L.t, L.hx, L.hy, L.ox, L.oy, L.hole,
                 L.color, alpha)
        add_layer_3d(ax, L2, Z(L.z0, g), Z(L.z0 + L.t, g))

    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_xlim(-280, 280); ax.set_ylim(-280, 280)


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
                           "field cage"),
                      (200, "mesh + amp gap (bulk pillars) + ESL resist\n"
                            "+ M1 front-end cards on the board edges"),
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
# Board copper from the production gerbers (like p2_board_copper.png)
# ─────────────────────────────────────────────────────────────────────────────

def draw_gerber(ax, gf, lw_scale):
    """Additive copper render: regions, strokes at aperture width, flashes."""
    small = [r for r in gf.regions if 2 < len(r) <= 10000]
    zones = [r for r in gf.regions if len(r) > 10000]
    if small:
        polys = [plt.Polygon(np.array(r), closed=True) for r in small]
        ax.add_collection(PatchCollection(polys, facecolor=CU_COLOR,
                                          edgecolor="none", alpha=0.95, zorder=2))
    if zones:
        ax.add_collection(LineCollection(
            [np.array(z) for z in zones], colors=CU_COLOR, linewidths=0.1,
            alpha=0.35, zorder=1))
    by_w = {}
    for s in gf.segments:
        w = gf.apertures.get(s.aperture)
        w = w.size if w else 0.15
        by_w.setdefault(round(w, 3), []).append(s)
    for w, segs in by_w.items():
        lines = [s.sample(16) for s in segs]
        ax.add_collection(LineCollection(
            lines, colors=CU_COLOR, linewidths=max(w * lw_scale, 0.25),
            alpha=0.9, zorder=2, capstyle="round"))
    pats = []
    for f in gf.flashes:
        a = gf.apertures.get(f.aperture)
        if a is None:
            continue
        if a.template in ("R", "O") and len(a.params) >= 2:
            pats.append(Rectangle((f.x - a.params[0]/2, f.y - a.params[1]/2),
                                  a.params[0], a.params[1]))
        elif a.template == "C" and a.params:
            pats.append(Circle((f.x, f.y), max(a.params[0], 0.05) / 2))
        elif a.params:
            pats.append(Circle((f.x, f.y), max(a.params[0], 0.1) / 2))
    if pats:
        ax.add_collection(PatchCollection(pats, facecolor=CU_COLOR,
                                          edgecolor="none", alpha=0.95, zorder=3))


def board_overlays(ax, lw=1.2):
    b = plt.Rectangle((-220, -220), 470, 470, fill=False, ec="crimson",
                      lw=lw, zorder=5, label="board outline (L8)")
    ax.add_patch(b)
    a = M.ACTIVE / 2
    ax.add_patch(plt.Rectangle((-a, -a), 2*a, 2*a, fill=False, ec="royalblue",
                               ls=":", lw=lw, zorder=5,
                               label="active area (399.36)"))
    ax.add_patch(plt.Rectangle((-205, -205), 410, 410, fill=False, ec="0.35",
                               ls="--", lw=lw, zorder=5,
                               label="gas frame aperture (410)"))
    ax.add_patch(plt.Rectangle((-201, -201), 402, 402, fill=False, ec="0.6",
                               ls="-.", lw=0.9, zorder=5,
                               label="support-plate aperture (402)"))
    ax.plot(0, 0, "+", color="k", ms=10, mew=1.5, zorder=6)


def fig_board(outdir):
    layers = [("DFS3498A_L2-pads.gbr",   "L4 — readout pads (65% Cu on board)"),
              ("DFS3498A_L3-TrackY.gbr", "L5 — Y strips (42%)"),
              ("DFS3498A_L4-TrackX.gbr", "L6 — X strips (42%)")]
    fig = plt.figure(figsize=(19, 6.9))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.2, 1.2, 1.2, 1.0])
    axs = [fig.add_subplot(gs[i]) for i in range(4)]

    for ax, (fn, title) in zip(axs[:3], layers):
        gf = parse(os.path.join(GERBER_DIR, fn))
        draw_gerber(ax, gf, lw_scale=0.55)
        board_overlays(ax)
        ax.set_xlim(-240, 270); ax.set_ylim(-240, 270)
        ax.set_aspect("equal"); ax.grid(alpha=0.2, lw=0.4)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x [mm]")
    axs[0].set_ylabel("y [mm]")
    axs[0].legend(loc="upper left", fontsize=7)

    # zoom: pad field + strip fan-out toward the +x connector edge
    axz = axs[3]
    for fn, _ in (layers[0], layers[2]):
        gf = parse(os.path.join(GERBER_DIR, fn))
        draw_gerber(axz, gf, lw_scale=3.0)
    board_overlays(axz)
    axz.set_xlim(190, 258); axz.set_ylim(55, 130)
    axz.set_aspect("equal"); axz.grid(alpha=0.2, lw=0.4)
    axz.set_title("zoom: pads + X-strip fan-out\nat the +x connector edge (M1 cards)",
                  fontsize=10)
    axz.set_xlabel("x [mm]")

    fig.suptitle("MX17 readout board — production gerbers (DFS3498A) with model "
                 "overlays; Cu coverage feeds the density-scaled sheets",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(outdir, "mx17_board_copper.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Peel-back close-up: the board layers revealed band by band
# ─────────────────────────────────────────────────────────────────────────────

def _jag(x, y0, y1, amp=1.1, step=2.6, seed=0):
    """Jagged 'torn paper' vertical boundary from (x,y0) to (x,y1)."""
    rng = np.random.default_rng(seed)
    ys = np.arange(y0, y1 + step, step)
    xs = x + rng.uniform(-amp, amp, len(ys))
    xs[0] = xs[-1] = x
    return list(zip(xs, ys))


def draw_gerber_clipped(ax, gf, clip_patch, lw_scale, color=CU_COLOR, alpha=0.95):
    """draw_gerber, but every artist clipped to clip_patch."""
    arts = []
    small = [r for r in gf.regions if 2 < len(r) <= 10000]
    if small:
        polys = [plt.Polygon(np.array(r), closed=True) for r in small]
        arts.append(ax.add_collection(PatchCollection(
            polys, facecolor=color, edgecolor="none", alpha=alpha, zorder=3)))
    by_w = {}
    for s in gf.segments:
        w = gf.apertures.get(s.aperture)
        w = w.size if w else 0.15
        by_w.setdefault(round(w, 3), []).append(s)
    for w, segs in by_w.items():
        lines = [s.sample(8) for s in segs]
        arts.append(ax.add_collection(LineCollection(
            lines, colors=color, linewidths=max(w * lw_scale, 0.3),
            alpha=alpha, zorder=3, capstyle="round")))
    pats = []
    for f in gf.flashes:
        a = gf.apertures.get(f.aperture)
        if a is None or not a.params:
            continue
        if a.template in ("R", "O") and len(a.params) >= 2:
            pats.append(Rectangle((f.x - a.params[0]/2, f.y - a.params[1]/2),
                                  a.params[0], a.params[1]))
        elif a.template == "C":
            pats.append(Circle((f.x, f.y), max(a.params[0], 0.05) / 2))
    if pats:
        arts.append(ax.add_collection(PatchCollection(
            pats, facecolor=color, edgecolor="none", alpha=alpha, zorder=3)))
    for a in arts:
        a.set_clip_path(clip_patch)


def fig_peel(outdir):
    """One region of the board, layers peeled away left to right."""
    X0, X1, Y0, Y1 = -44.0, 44.0, -22.0, 22.0
    nb = 4
    bw = (X1 - X0) / nb
    bands = [(X0 + i*bw, X0 + (i+1)*bw) for i in range(nb)]

    # (title, substrate colour, copper colour, gerber file or None=resist)
    spec = [
        ("① top surface\nESL resist strips + bulk pillars", "#d9c48f", None, None),
        ("② coat removed\nL4 readout pads", "#d9c48f", "#e09a55",
         "DFS3498A_L2-pads.gbr"),
        ("③ pads removed\nL5 Y strips", "#c4ac74", "#b87333",
         "DFS3498A_L3-TrackY.gbr"),
        ("④ deepest\nL6 X strips", "#a8905c", "#8a5a28",
         "DFS3498A_L4-TrackX.gbr"),
    ]

    fig = plt.figure(figsize=(19, 8.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.4, 1.0])
    ax = fig.add_subplot(gs[0])
    axd = fig.add_subplot(gs[1])

    for i, ((bx0, bx1), (title, subc, cuc, fn)) in enumerate(zip(bands, spec)):
        left = _jag(bx0, Y0, Y1, seed=i) if i > 0 else \
            [(bx0, Y0), (bx0, Y1)]
        right = _jag(bx1, Y0, Y1, seed=i+1) if i < nb-1 else \
            [(bx1, Y0), (bx1, Y1)]
        poly = left + right[::-1]
        patch = plt.Polygon(np.array(poly), closed=True, facecolor=subc,
                            edgecolor="none", zorder=1)
        ax.add_patch(patch)

        if fn is not None:
            gf = parse(os.path.join(GERBER_DIR, fn))
            draw_gerber_clipped(ax, gf, patch, lw_scale=2.2, color=cuc)
        else:
            # ESL resistive strips: 550 um wide / 250 um gaps (confirmed
            # 2026-08-06; deliberately NOT the 0.78 mm pad pitch). The L4 pads
            # show through the gaps. Strip artwork is not in the gerber set,
            # so the strips are drawn from the confirmed spec.
            gp = parse(os.path.join(GERBER_DIR, "DFS3498A_L2-pads.gbr"))
            draw_gerber_clipped(ax, gp, patch, lw_scale=2.2,
                                color="#c8874a", alpha=0.8)
            strips = []
            xs = np.arange(np.floor(X0/0.8)*0.8, bx1 + 1, 0.8)
            for xstrip in xs:
                strips.append(Rectangle((xstrip - 0.275, Y0), 0.55, Y1 - Y0))
            sc = ax.add_collection(PatchCollection(
                strips, facecolor="#1c1c1c", edgecolor="none", alpha=0.92,
                zorder=4))
            sc.set_clip_path(patch)
            # bulk pillars (real: 3498A_bulk.gbr, Ø0.6 on ~4.7 mm pitch)
            gb = parse(os.path.join(REPO, "design", "gerbers", "readout_pcb",
                                    "3498A_bulk.gbr"))
            pil = [Circle((f.x, f.y), 0.3) for f in gb.flashes
                   if X0 - 2 < f.x < bx1 + 2 and Y0 - 2 < f.y < Y1 + 2]
            pc = ax.add_collection(PatchCollection(
                pil, facecolor="#e8e2d2", edgecolor="#9a9384",
                linewidths=0.4, zorder=5))
            pc.set_clip_path(patch)

        # torn-edge shadow on the left boundary of each deeper band
        if i > 0:
            sh = plt.Polygon(np.array(left + [(x + 1.6, y) for x, y in left[::-1]]),
                             closed=True, facecolor="k", alpha=0.28,
                             edgecolor="none", zorder=5)
            ax.add_patch(sh)
            sh.set_clip_path(patch)
        ax.plot([q[0] for q in left], [q[1] for q in left], color="k",
                lw=1.0, zorder=6)
        ax.text((bx0 + bx1)/2, Y1 + 1.2, title, ha="center", va="bottom",
                fontsize=10.5)

    ax.set_xlim(X0 - 1, X1 + 1)
    ax.set_ylim(Y0 - 2, Y1 + 7.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title("MX17 readout board close-up — layers peeled back "
                 "(gerber geometry; ESL strips drawn from the confirmed spec)",
                 fontsize=13, pad=14)

    # depth key: mini stack elevation linking bands to depth
    labels = [("ESL resistive strips — 550 µm wide,\n"
               "250 µm gaps (confirmed; pads show\n"
               "through). Geant4: 100 µm slab ×0.69.\n"
               "White dots: bulk pillars Ø0.6 @ 4.68 mm",
               "#1c1c1c"),
              ("L4 pads — 0.68 mm on 0.78 mm pitch\n(88 % Cu over active)",
               "#e09a55"),
              ("L5 Y strips (53 %)", "#b87333"),
              ("L6 X strips (53 %)", "#8a5a28")]
    yd = 0
    for i, (lab, c) in enumerate(labels):
        axd.add_patch(Rectangle((0, yd), 1.4, 0.55, facecolor=c,
                                edgecolor="k", lw=0.6))
        axd.text(1.55, yd + 0.27, f"{'①②③④'[i]}  {lab}", va="center",
                 fontsize=9)
        yd -= 1.0
    axd.annotate("", xy=(-0.35, yd + 1.0), xytext=(-0.35, 0.55),
                 arrowprops=dict(arrowstyle="-|>", lw=1.6, color="k"))
    axd.text(-0.75, (yd + 1.55)/2, "deeper into the board\n(beam direction)",
             rotation=90, ha="center", va="center", fontsize=9)
    axd.set_xlim(-1.2, 8.5); axd.set_ylim(yd + 0.3, 1.3)
    axd.axis("off")
    axd.set_title("depth order", fontsize=11)

    fig.tight_layout()
    out = os.path.join(outdir, "mx17_board_peel.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Plan views (from upstream / downstream) for alignment checking
# ─────────────────────────────────────────────────────────────────────────────

def fig_plan(outdir):
    """Plan views built from the geometry dump (unique footprints per lv)."""
    vols = load_dump()
    foot = {}   # name -> list of unique (ox, oy, hx, hy, hole)
    for L in vols:
        key = (round(L.ox, 3), round(L.oy, 3), round(L.hx, 3), round(L.hy, 3),
               L.hole and tuple(round(v, 3) for v in L.hole))
        foot.setdefault(L.name, set()).add(key)

    def draw(ax, name, *, fc, ec, lw=1.2, ls="-", alpha=0.8, zorder=2,
             label=None):
        for i, (ox, oy, hx, hy, hole) in enumerate(sorted(foot.get(name, []))):
            r = Rectangle((ox - hx, oy - hy), 2*hx, 2*hy, facecolor=fc,
                          edgecolor=ec, lw=lw, ls=ls, alpha=alpha,
                          zorder=zorder, label=label if i == 0 else None)
            ax.add_patch(r)
            if hole:
                hcx, hcy, hhx, hhy = hole
                ax.add_patch(Rectangle((hcx - hhx, hcy - hhy), 2*hhx, 2*hhy,
                                       facecolor="white", edgecolor=ec,
                                       lw=lw*0.7, alpha=alpha, zorder=zorder))

    def conn_clusters(ax, zorder):
        # board connector-copper clusters, measured from the gerbers
        for tang in (+99.8, -99.8):
            ax.add_patch(Rectangle((228.5, tang - 14.5), 20, 29,
                                   facecolor="#b87333", edgecolor="none",
                                   alpha=0.9, zorder=zorder))
            ax.add_patch(Rectangle((tang - 14.5, 228.5), 29, 20,
                                   facecolor="#b87333", edgecolor="none",
                                   alpha=0.9, zorder=zorder))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19, 9.6))

    # ── view from upstream (beam side, window first) ─────────────────────────
    ax = ax1
    draw(ax, "SupportPlate_Al", fc="#f2f2f2", ec="0.55", lw=1.0, zorder=1,
         label="470² plates (+15,+15) / plate aperture")
    draw(ax, "FrontEndPCB", fc="#c9e4c9", ec="#2e7d32", lw=1.3, zorder=2,
         alpha=0.75, label="M1 front-end cards (gerber-aligned)")
    draw(ax, "GasFrame_Al", fc="#e2e6ee", ec="0.3", lw=1.4, zorder=3,
         label="gas frame ring (440/410)")
    draw(ax, "FieldCagePCB", fc="#9fd49f", ec="#2e7d32", lw=0.8, zorder=4,
         label="field cage (4×, pinwheeled)")
    draw(ax, "WindowGapGas", fc="#dff2fa", ec="#4499bb", lw=1.3, zorder=5,
         alpha=0.8, label="window free span (400.2²), bulged")
    for L in sorted((L for L in vols if L.name == "WindowBulgeGas"),
                    key=lambda L: -L.hx):
        ax.add_patch(Rectangle((L.ox - L.hx, L.oy - L.hy), 2*L.hx, 2*L.hy,
                               facecolor="none", edgecolor="#4499bb", lw=0.5,
                               zorder=6))
    a = M.ACTIVE/2
    ax.add_patch(Rectangle((-a, -a), 2*a, 2*a, facecolor="none",
                           edgecolor="royalblue", ls=":", lw=1.4, zorder=7,
                           label="active area (399.36²)"))
    ax.annotate("M1 cards: 2 per connector edge,\nflat on the board edge,\n"
                "centred on the connector copper",
                xy=(240, 100), xytext=(150, 300), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.set_title("view from UPSTREAM (beam side) — window, frame, M1 cards",
                 fontsize=12)

    # ── view from downstream (support-plate side) ────────────────────────────
    ax = ax2
    draw(ax, "GasFrame_Al", fc="#f5f5f5", ec="0.7", lw=0.8, zorder=1,
         label="window/frame footprint (440²)")
    draw(ax, "FrontEndPCB", fc="#c9e4c9", ec="#2e7d32", lw=1.3, zorder=2,
         alpha=0.75, label="M1 front-end cards")
    draw(ax, "PCB_Rohacell", fc="#efe9d2", ec="#8a7a40", lw=1.4, zorder=3,
         alpha=0.85, label="readout board / rohacell (470², +15,+15)")
    draw(ax, "SupportPlate_Al", fc="#dcdce2", ec="0.35", lw=1.4, zorder=4,
         alpha=0.9, label="support plate (402² aperture)")
    conn_clusters(ax, 5)
    a = M.ACTIVE/2
    ax.add_patch(Rectangle((-a, -a), 2*a, 2*a, facecolor="none",
                           edgecolor="royalblue", ls=":", lw=1.4, zorder=6,
                           label="active area (399.36²)"))
    ax.annotate("board connector copper\n(measured: tangential ±83..117,\n"
                "radial 213..248)", xy=(238, -100), xytext=(120, -300),
                fontsize=9, arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.set_title("view from DOWNSTREAM (support-plate side) — board, plate, "
                 "connectors", fontsize=12)

    for ax in (ax1, ax2):
        ax.plot(0, 0, "+", color="k", ms=12, mew=1.6, zorder=10)
        ax.set_xlim(-345, 345); ax.set_ylim(-345, 345)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2, lw=0.4)
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
        ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    fig.suptitle("MX17 as-built model — plan views RENDERED FROM THE GEANT4 "
                 "GEOMETRY DUMP (active axis at the origin)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(outdir, "mx17_plan_views.png")
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
    hcx, hcy, hhx, hhy = L.hole
    if abs(hcy) >= hhy:
        return [(lo, hi)]
    out = []
    if lo < hcx - hhx: out.append((lo, hcx - hhx))
    if hi > hcx + hhx: out.append((hcx + hhx, hi))
    return out


def fig_xsec(outdir):
    items = load_dump()

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

    zs = {}
    for L in items:
        zs.setdefault(L.name, (L.z0, L.t))
    z_mesh = zs["Micromesh"][0]
    z_end = zs["PCB_FR4_5"][0] + zs["PCB_FR4_5"][1]
    ax2.set_xlim(-30, 30); ax2.set_ylim(z_end + 0.1, z_mesh - 0.1)
    ax2.set_xlabel("x [mm]")
    ax2.set_title("Zoom: mesh → laminate (1.70 mm board body)")
    ax2.grid(alpha=0.25, lw=0.4)
    # "ResistLayer" in the patterned build, "ResistivePaste" in the
    # homogenized one — annotate whichever the dump actually contains.
    resist = "ResistLayer" if "ResistLayer" in zs else "ResistivePaste"
    resist_lab = ("ESL strips 550/250 µm (100 µm)" if resist == "ResistLayer"
                  else "paste (100 µm)")
    for name, lab in [("Micromesh", "mesh (38 µm eff.)"),
                      ("AmpGas", "amp gap (150 µm)"),
                      (resist, resist_lab),
                      ("PCB_Kapton", "kapton (50 µm)"),
                      ("PCB_FR4_3", "5× Cu 26 µm (coverage-scaled)\n/ FR4 246.4 µm")]:
        z0, t = zs[name]
        ax2.annotate(lab, xy=(30, z0 + t/2), xytext=(32, z0 + t/2),
                     fontsize=8, va="center", annotation_clip=False)

    ax3.set_title("Stack (not to scale)")
    # z-advancing chain = unique on-axis volumes ordered by z (skip side items)
    side_names = {"WindowFlange_Al", "GasFrame_Al", "FieldCagePCB",
                  "FrontEndPCB", "WindowBulgeGas"}
    biggest = {}
    for L in items:
        if L.name in side_names:
            continue
        if L.name not in biggest or L.hx > biggest[L.name].hx:
            biggest[L.name] = L
    y = 0
    for L in sorted(biggest.values(), key=lambda L: L.z0):
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
         "confirmed (Dylan 2026-08-06)", 0.7, "#ffc4c4", C_HW),
        ("paste",   "ESL resist strips — 550/250 µm → slab ×0.69",
         "strip spec confirmed; 100 µm thickness still a guess", 0.45, "#555555", C_HW),
        ("pillars", "bulk pillars — Ø0.6 mm, 85×85 @ 4.68 mm",
         "exact grid from 3498A_bulk.gbr; in the amp gas", 0.4, "#e8e2d2", C_CAD),
        ("pcb",     "readout board — 1.70 mm body, 470²",
         "CAD total; Cu = 5 gerber layers × 26 µm, density-scaled\nby measured coverage (L3–L7); FR4 filler solved", 0.9, "#cc6619", C_CAD),
        ("m1",      "M1 cards — flat on the board edges",
         "gerber-anchored position; 1.6 mm laminate assumed,\nMec8/pogo connectors omitted", 0.5, "#9fd49f", C_ASSUME),
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
    ap.add_argument("--only", choices=["board", "peel", "plan", "xsec", "3d",
                                       "status"], default=None)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if args.only in (None, "board"):
        fig_board(args.out_dir)
    if args.only in (None, "peel"):
        fig_peel(args.out_dir)
    if args.only in (None, "plan"):
        fig_plan(args.out_dir)
    if args.only in (None, "xsec"):
        fig_xsec(args.out_dir)
    if args.only in (None, "3d"):
        fig_3d(args.out_dir)
    if args.only in (None, "status"):
        fig_status(args.out_dir)


if __name__ == "__main__":
    main()
