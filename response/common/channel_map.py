#!/usr/bin/env python3
"""
channel_map.py — T2: which readout pads belong to the X channels and which to
the Y channels, read out of the gerbers rather than assumed.

design/RESPONSE_SIM_PLAN.md §3 asks: "Confirm from the via pattern in
DFS3498A_pla1-2.gbr which pads belong to X vs Y — likely alternating
checkerboard". §12 item 7 lists it as needing confirmation.

How the board actually does it
------------------------------
The gerber TF.FileFunction headers give the physical stackup, which is not the
same as the file names:

    DFS3498A_L2-pads.gbr    -> Copper,L4,Inr    the 0.68 mm pads
    DFS3498A_L3-TrackY.gbr  -> Copper,L5,Inr    the Y bus layer
    DFS3498A_L4-TrackX.gbr  -> Copper,L6,Inr    the X bus layer
    DFS3498A_pla2-5.gbr     -> Plated,4,7,PTH   vias, physical L4 -> L7

So a single plated via drops from every pad through BOTH bus layers. The via
cannot be what selects the channel — it passes through both. What selects the
channel is which bus layer puts an annular ring on that via and which leaves a
clearance.

The first thing to check is the 0.5 mm annular rings (aperture C,0.500000).
They do NOT discriminate: inside the active area BOTH bus layers carry exactly
512 x 512 = 262144 rings, on the 0.78 mm pad grid to 0.0 um. Every via is
ringed on every layer it passes, as a plated through hole must be.

What discriminates is the 0.1 mm connector stub: a 390 um trace (exactly half
the pad pitch) running from a pad's ring out to the bus line. A pad is on a
layer's channel if and only if that layer draws a stub at it. Counting them
settles the checkerboard question outright:

    Y layer  131072 stubs = 512 x 256, ALL at (col+row) EVEN
    X layer  131072 stubs = 512 x 256, ALL at (col+row) ODD

131072 x 2 = 262144: every pad is claimed exactly once, by exactly one view.

The bus direction confirms which is which, independently: of the 0.1 mm
segments on the Y layer 132279 are horizontal and 1973 vertical, and on the X
layer 132248 are vertical and 1942 horizontal. So the Y layer busses along x
(a channel is a ROW, measuring y) and the X layer busses along y (a channel is
a COLUMN, measuring x).

CONSEQUENCE FOR THE RESPONSE MODEL, and it is not a small one: a channel is
every-OTHER pad along its line, 256 pads on a 1.56 mm pitch, not 512 pads on a
0.78 mm pitch. Its weighting potential is a comb, not a solid strip. Any model
that busses a full contiguous row will get the sharing wrong.

    python3 -m response.common.channel_map [--gerberdir design/gerbers/readout_pcb]
"""

from __future__ import annotations

import argparse
import os
import re
from collections import Counter

import numpy as np

from . import constants as C

DOT_DIAM_MM = 0.5            # via annular ring — present on BOTH layers
STUB_APERTURE_MM = 0.1       # the bus connector stub — this is what selects
STUB_LEN_MM = 0.39           # exactly half the pad pitch


def parse_flashes(path, want_diam_mm=DOT_DIAM_MM, tol=1e-6):
    """
    Return the (x, y) [mm] of every D03 flash made with a circular aperture of
    the requested diameter.

    A deliberately small gerber reader: enough for %FSLAX36Y36 absolute-mode
    files with circular apertures, which is what this board is. It handles the
    modal coordinate behaviour (a command may omit X or Y and inherit the
    previous value), which is not optional — most of the flashes in these files
    give only one coordinate.
    """
    fmt_int = fmt_dec = None
    unit_mm = True
    apertures = {}               # D-code -> diameter mm (circles only)
    want = set()
    cur = None
    x = y = 0
    out = []

    ap_re = re.compile(r"%ADD(\d+)([A-Z]),([0-9.X]+)\*%")
    fs_re = re.compile(r"%FSLAX(\d)(\d)Y(\d)(\d)\*%")
    dn_re = re.compile(r"^D(\d+)\*$")
    op_re = re.compile(r"(?:X(-?\d+))?(?:Y(-?\d+))?D0?([123])\*")

    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = fs_re.match(line)
            if m:
                fmt_int, fmt_dec = int(m.group(1)), int(m.group(2))
                continue
            if line.startswith("%MOIN"):
                unit_mm = False
                continue
            m = ap_re.match(line)
            if m:
                dcode, shape, params = int(m.group(1)), m.group(2), m.group(3)
                if shape == "C":
                    d = float(params.split("X")[0])
                    apertures[dcode] = d
                    if abs(d - want_diam_mm) < tol:
                        want.add(dcode)
                continue
            m = dn_re.match(line)
            if m:
                cur = int(m.group(1))
                continue
            m = op_re.match(line)
            if m:
                if m.group(1) is not None:
                    x = int(m.group(1))
                if m.group(2) is not None:
                    y = int(m.group(2))
                if m.group(3) == "3" and cur in want:
                    out.append((x, y))

    if fmt_dec is None:
        raise ValueError(f"{path}: no %FSLAX format spec found")
    scale = 10.0 ** -fmt_dec
    arr = np.asarray(out, dtype=float) * scale
    if not unit_mm:
        arr *= 25.4
    return arr, apertures


def to_indices(xy):
    """Map pad-grid positions [mm] onto integer (col, row) on the 0.78 mm grid."""
    # The grid is exactly regular and spans +/-199.29 mm (512 x 0.78 centred).
    origin = -(C.PAD_N - 1) / 2 * C.PAD_PITCH_UM * 1e-3
    i = np.rint((xy[:, 0] - origin) / (C.PAD_PITCH_UM * 1e-3)).astype(int)
    j = np.rint((xy[:, 1] - origin) / (C.PAD_PITCH_UM * 1e-3)).astype(int)
    resid = np.abs(xy[:, 0] - (origin + i * C.PAD_PITCH_UM * 1e-3)).max()
    return i, j, resid


def parse_draws(path, want_diam_mm=STUB_APERTURE_MM, tol=1e-6):
    """Every D01 draw segment (x0, y0, x1, y1) [mm] made with the given aperture."""
    fmt_dec = None
    want, apertures = set(), {}
    cur, x, y, out = None, 0, 0, []
    ap_re = re.compile(r"%ADD(\d+)([A-Z]),([0-9.X]+)\*%")
    fs_re = re.compile(r"%FSLAX(\d)(\d)Y(\d)(\d)\*%")
    dn_re = re.compile(r"^D(\d+)\*$")
    op_re = re.compile(r"(?:X(-?\d+))?(?:Y(-?\d+))?D0?([123])\*")
    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            m = fs_re.match(line)
            if m:
                fmt_dec = int(m.group(2)); continue
            m = ap_re.match(line)
            if m:
                if m.group(2) == "C":
                    d = float(m.group(3).split("X")[0])
                    apertures[int(m.group(1))] = d
                    if abs(d - want_diam_mm) < tol:
                        want.add(int(m.group(1)))
                continue
            m = dn_re.match(line)
            if m:
                cur = int(m.group(1)); continue
            m = op_re.match(line)
            if m:
                px, py = x, y
                if m.group(1) is not None:
                    x = int(m.group(1))
                if m.group(2) is not None:
                    y = int(m.group(2))
                if m.group(3) == "1" and cur in want:
                    out.append((px, py, x, y))
    return np.asarray(out, dtype=float) * 10.0 ** -fmt_dec


def channel_pads(path):
    """
    The set of (col, row) pads this bus layer actually connects to, plus the
    bus orientation. A pad is connected iff the layer draws a 390 um connector
    stub starting on its centre.
    """
    seg = parse_draws(path)
    L = np.hypot(seg[:, 2] - seg[:, 0], seg[:, 3] - seg[:, 1])
    horiz = int(((np.abs(seg[:, 3] - seg[:, 1]) < 1e-4) &
                 (np.abs(seg[:, 2] - seg[:, 0]) > 1e-4)).sum())
    vert = int(((np.abs(seg[:, 2] - seg[:, 0]) < 1e-4) &
                (np.abs(seg[:, 3] - seg[:, 1]) > 1e-4)).sum())

    stub = seg[np.abs(L - STUB_LEN_MM) < 1e-4]
    origin = -(C.PAD_N - 1) / 2 * C.PAD_PITCH_UM * 1e-3
    pitch = C.PAD_PITCH_UM * 1e-3
    i = np.rint((stub[:, 0] - origin) / pitch)
    j = np.rint((stub[:, 1] - origin) / pitch)
    on = ((np.abs(stub[:, 0] - (origin + i * pitch)) < 1e-4) &
          (np.abs(stub[:, 1] - (origin + j * pitch)) < 1e-4) &
          (i >= 0) & (i < C.PAD_N) & (j >= 0) & (j < C.PAD_N))
    return {
        "pads": set(zip(i[on].astype(int).tolist(), j[on].astype(int).tolist())),
        "n_stub": int(len(stub)), "n_on_grid": int(on.sum()),
        "horiz": horiz, "vert": vert,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gerberdir", default="design/gerbers/readout_pcb")
    a = ap.parse_args()
    print("T2 — pad -> X/Y channel assignment, read from the gerbers\n")

    r = {}
    for tag, fn in (("Y", "DFS3498A_L3-TrackY.gbr"),
                    ("X", "DFS3498A_L4-TrackX.gbr")):
        d = channel_pads(os.path.join(a.gerberdir, fn))
        r[tag] = d
        axis = "along x (rows)" if d["horiz"] > d["vert"] else "along y (columns)"
        print(f"  {tag} layer  {fn}")
        print(f"    bus runs {axis}: {d['horiz']} horizontal / {d['vert']} "
              f"vertical 0.1 mm segments")
        print(f"    {d['n_stub']} connector stubs of {STUB_LEN_MM*1e3:.0f} µm; "
              f"{d['n_on_grid']} start on a pad centre inside the active area")
        par = Counter(((i + j) % 2) for (i, j) in d["pads"])
        print(f"    (col+row) parity of connected pads: {dict(par)}")

    pY, pX = r["Y"]["pads"], r["X"]["pads"]
    print(f"\n  Y pads {len(pY)}   X pads {len(pX)}   "
          f"overlap {len(pY & pX)}   union {len(pY | pX)} of {C.PAD_N**2}")

    parY = set((i + j) % 2 for (i, j) in pY)
    parX = set((i + j) % 2 for (i, j) in pX)
    clean = (len(parY) == 1 and len(parX) == 1 and parY != parX
             and not (pY & pX))
    print()
    if clean:
        print("  => CHECKERBOARD CONFIRMED (plan §12 item 7, RESOLVED).")
        print(f"     Y channels (rows, measuring y) take (col+row) "
              f"parity {parY.pop()}")
        print(f"     X channels (columns, measuring x) take parity {parX.pop()}")
        print(f"     Each channel therefore has {C.PAD_N//2} pads on a "
              f"{2*C.PAD_PITCH_UM/1000:.2f} mm pitch —")
        print("     a COMB, not a solid strip. The response model must bus it "
              "that way.")
    else:
        print("  => not a clean parity split; carry an explicit per-pad table.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
