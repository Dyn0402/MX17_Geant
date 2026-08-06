#!/usr/bin/env python3
"""
Copper coverage of the MX17 readout-board and M1 front-end copper layers,
rasterized from the production gerbers (PIL, default 0.05 mm/px — exact to
better than the aperture sizes involved; no external geometry deps).

The Geant4 model spans each copper layer as a full-face sheet, so the sheet
density must be scaled by the real covered-area fraction, exactly as p2_geant
does for the P2 wedge (their analyze_cu_coverage.py uses shapely; the numbers
agree at the raster resolution).

Readout board (DFS3498A): physical copper layers are L3..L8 — L1/L2 are the
bulk-micromegas side (resistive layer + mesh) and have no copper gerber.
Coverage is reported over
  * the board outline   (x, y in −220..+250 mm, the model's PCB_Cu span), and
  * the active area     (399.36 mm square centred on the origin).
The `Bot` (L8) file contains only a 10 µm outline stroke — no real copper.

M1 front-end (DFS3498AM1): the gerbers are drawn in the READOUT-BOARD frame
at the card's as-mounted position — one 41 × 160 mm card at x 219.4..260.4,
y 19.6..179.6 (edge-mill layer), whose bottom pogo-pad field lands on the
board's connector cluster at tangential +100. Coverage is reported over that
card outline.

Usage:
    python scripts/gerber/analyze_cu_coverage.py [--res 0.05]

The resulting fractions are hard-coded in shared/MX17ModuleGeometry.hh
(kReadoutCuCoverage / kM1CuCoverage) — rerun this script if the gerbers
change and update the header.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gerber_outline import parse  # noqa: E402

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "..", ".."))
GDIR = os.path.join(REPO, "design", "gerbers")

# Readout layers: (file, physical layer, description)
READOUT_LAYERS = [
    ("readout_pcb/DFS3498A_top-cu.gbr",   "L3", "top guard/border ring"),
    ("readout_pcb/DFS3498A_L2-pads.gbr",  "L4", "readout pads"),
    ("readout_pcb/DFS3498A_L3-TrackY.gbr","L5", "Y strips"),
    ("readout_pcb/DFS3498A_L4-TrackX.gbr","L6", "X strips"),
    ("readout_pcb/DFS3498A_L5.gbr",       "L7", "via/fan-out"),
    ("readout_pcb/DFS3498A_Bot.gbr",      "L8", "outline only"),
]
M1_LAYERS = [(f"M1/DFS3498AM1_{n}.gbr", n) for n in
             ("top", "L2", "L3", "L4", "L5", "bot")]

BOARD = (-220.0, -220.0, 250.0, 250.0)      # readout board, active-area frame
ACTIVE_HALF = 399.36 / 2


class Raster:
    """Dark/clear gerber raster over a bbox at `res` mm/px."""

    def __init__(self, bbox, res):
        self.x0, self.y0, self.x1, self.y1 = bbox
        self.res = res
        self.w = int(math.ceil((self.x1 - self.x0) / res))
        self.h = int(math.ceil((self.y1 - self.y0) / res))
        self.img = Image.new("1", (self.w, self.h), 0)
        self.drw = ImageDraw.Draw(self.img)

    def px(self, x, y):
        return ((x - self.x0) / self.res, (y - self.y0) / self.res)

    def draw(self, gf):
        ap = gf.apertures
        # regions first (pours), then strokes, then flashes — order within a
        # polarity group does not matter; dark/clear order is preserved by
        # replaying the file order per primitive class, which is fine here
        # because the MX17 layers only use clear *after* dark (border cutouts).
        prims = ([("region", r, d) for r, d in zip(gf.regions, gf.region_dark)]
                 + [("seg", s, s.dark) for s in gf.segments]
                 + [("flash", f, f.dark) for f in gf.flashes])
        # draw dark first, then clear
        for want_dark in (True, False):
            fill = 1 if want_dark else 0
            for kind, obj, dark in prims:
                if dark != want_dark:
                    continue
                if kind == "region":
                    self.drw.polygon([self.px(x, y) for x, y in obj], fill=fill)
                elif kind == "seg":
                    a = ap.get(obj.aperture)
                    wpx = max(int(round((a.size if a else 0.1) / self.res)), 1)
                    pts = [self.px(x, y) for x, y in obj.sample(24)]
                    self.drw.line(pts, fill=fill, width=wpx, joint="curve")
                    r = wpx / 2
                    for p in (pts[0], pts[-1]):   # round caps
                        self.drw.ellipse([p[0]-r, p[1]-r, p[0]+r, p[1]+r],
                                         fill=fill)
                else:
                    a = ap.get(obj.aperture)
                    if a is None:
                        continue
                    cx, cy = self.px(obj.x, obj.y)
                    if a.template == "C" and a.params:
                        r = a.params[0] / 2 / self.res
                        self.drw.ellipse([cx-r, cy-r, cx+r, cy+r], fill=fill)
                    elif a.template in ("R", "O") and len(a.params) >= 2:
                        hw = a.params[0] / 2 / self.res
                        hh = a.params[1] / 2 / self.res
                        self.drw.rectangle([cx-hw, cy-hh, cx+hw, cy+hh],
                                           fill=fill)
                        if len(a.params) >= 4:    # rect with rect hole
                            hw2 = a.params[2] / 2 / self.res
                            hh2 = a.params[3] / 2 / self.res
                            self.drw.rectangle([cx-hw2, cy-hh2, cx+hw2, cy+hh2],
                                               fill=0 if want_dark else 1)
                    elif a.params:                # macro: treat as circle
                        r = max(a.params[0], 0.05) / 2 / self.res
                        self.drw.ellipse([cx-r, cy-r, cx+r, cy+r], fill=fill)

    def coverage(self, bbox=None):
        arr = np.asarray(self.img, dtype=np.uint8)
        if bbox is None:
            return arr.mean()
        x0, y0, x1, y1 = bbox
        j0 = max(int((x0 - self.x0) / self.res), 0)
        j1 = min(int((x1 - self.x0) / self.res), self.w)
        i0 = max(int((y0 - self.y0) / self.res), 0)
        i1 = min(int((y1 - self.y0) / self.res), self.h)
        return arr[i0:i1, j0:j1].mean()


def main():
    apar = argparse.ArgumentParser()
    apar.add_argument("--res", type=float, default=0.05, help="mm per pixel")
    args = apar.parse_args()

    print(f"raster {args.res} mm/px\n")
    print("Readout board (DFS3498A) — coverage over board / active area:")
    active = (-ACTIVE_HALF, -ACTIVE_HALF, ACTIVE_HALF, ACTIVE_HALF)
    for fn, lay, desc in READOUT_LAYERS:
        gf = parse(os.path.join(GDIR, fn))
        r = Raster(BOARD, args.res)
        r.draw(gf)
        cb, ca = r.coverage(), r.coverage(active)
        print(f"  {lay}  {os.path.basename(fn):26s} {desc:22s} "
              f"board {cb:6.4f}   active {ca:6.4f}")

    # M1: coverage over the card outline (from DFS3498AM1_ebm)
    card = (219.40, 19.59, 260.40, 179.59)
    print("\nM1 front-end (DFS3498AM1) — coverage over the 41 x 160 card outline:")
    for fn, lay in M1_LAYERS:
        gf = parse(os.path.join(GDIR, fn))
        r = Raster(card, args.res)
        r.draw(gf)
        print(f"  {lay:4s} {r.coverage():6.4f}")


if __name__ == "__main__":
    main()
