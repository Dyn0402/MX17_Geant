#!/usr/bin/env python3
"""
Minimal RS-274X reader focused on *geometry extraction*, not rendering.

Copied from p2_geant/scripts/gerber/gerber_outline.py (keep roughly in sync),
plus %LPD/%LPC polarity tracking, which the MX17 top-cu layer needs.

Written for the P2 wedge gerbers (KiCad 9, format 4.6, MM, absolute).  It is
deliberately small: it understands only the subset of Gerber that KiCad emits
for Edge_Cuts / copper layers, which is all we need to recover detector
dimensions for the Geant4 geometry.

Supported:
  %FSLAX<n><m>Y<n><m>*%   coordinate format
  %MOMM*% / %MOIN*%       units
  %ADD<n><template>*%     aperture definitions (C, R, O, P, and macro names)
  D01/D02/D03             draw / move / flash
  G01                     linear interpolation
  G02/G03 + G75           clockwise / counter-clockwise circular, multi-quadrant
  G36/G37                 region (polygon) fill

Not supported (silently ignored): aperture macros beyond their name, LP
polarity changes, step-and-repeat, mirroring/rotation blocks.  None of these
affect the outline geometry we care about.

Public API
----------
parse(path) -> GerberFile
    .segments : list of Segment  (straight or arc, in mm)
    .flashes  : list of Flash    (aperture flashed at a point, in mm)
    .regions  : list of list[(x, y)]  closed G36/G37 polygons
    .apertures: dict  D-code -> ApertureDef
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# Data model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ApertureDef:
    code: int
    template: str            # 'C', 'R', 'O', 'P', or a macro name
    params: list[float]
    function: str | None = None   # from %TA.AperFunction,...

    @property
    def size(self) -> float:
        """Nominal size in mm (diameter for C, long side for R/O)."""
        return max(self.params) if self.params else 0.0


@dataclass
class Segment:
    """A stroked path element.  kind is 'line' or 'arc'."""
    kind: str
    x0: float
    y0: float
    x1: float
    y1: float
    aperture: int
    dark: bool = True        # %LPD (additive) vs %LPC (clear)
    # arc only:
    cx: float | None = None
    cy: float | None = None
    ccw: bool = False

    @property
    def radius(self) -> float:
        if self.cx is None:
            return 0.0
        return math.hypot(self.x0 - self.cx, self.y0 - self.cy)

    def length(self) -> float:
        if self.kind == "line":
            return math.hypot(self.x1 - self.x0, self.y1 - self.y0)
        r = self.radius
        a0 = math.atan2(self.y0 - self.cy, self.x0 - self.cx)
        a1 = math.atan2(self.y1 - self.cy, self.x1 - self.cx)
        d = (a1 - a0) % (2 * math.pi) if self.ccw else (a0 - a1) % (2 * math.pi)
        if d < 1e-12:
            d = 2 * math.pi          # full circle
        return r * d

    def sample(self, n: int = 32) -> list[tuple[float, float]]:
        """Discretise into points (inclusive of both endpoints)."""
        if self.kind == "line":
            return [(self.x0, self.y0), (self.x1, self.y1)]
        r = self.radius
        a0 = math.atan2(self.y0 - self.cy, self.x0 - self.cx)
        a1 = math.atan2(self.y1 - self.cy, self.x1 - self.cx)
        if self.ccw:
            d = (a1 - a0) % (2 * math.pi)
        else:
            d = -((a0 - a1) % (2 * math.pi))
        if abs(d) < 1e-12:
            d = 2 * math.pi if self.ccw else -2 * math.pi
        return [
            (self.cx + r * math.cos(a0 + d * i / n),
             self.cy + r * math.sin(a0 + d * i / n))
            for i in range(n + 1)
        ]


@dataclass
class Flash:
    x: float
    y: float
    aperture: int
    dark: bool = True


@dataclass
class GerberFile:
    path: str
    unit: str = "mm"
    apertures: dict[int, ApertureDef] = field(default_factory=dict)
    segments: list[Segment] = field(default_factory=list)
    flashes: list[Flash] = field(default_factory=list)
    regions: list[list[tuple[float, float]]] = field(default_factory=list)
    region_dark: list[bool] = field(default_factory=list)
    attributes: dict[str, str] = field(default_factory=dict)

    def bbox(self) -> tuple[float, float, float, float]:
        xs, ys = [], []
        for s in self.segments:
            for x, y in s.sample(16):
                xs.append(x)
                ys.append(y)
        for f in self.flashes:
            xs.append(f.x)
            ys.append(f.y)
        for reg in self.regions:
            for x, y in reg:
                xs.append(x)
                ys.append(y)
        if not xs:
            return (0.0, 0.0, 0.0, 0.0)
        return (min(xs), min(ys), max(xs), max(ys))


# ─────────────────────────────────────────────────────────────────────────────
# Parser
# ─────────────────────────────────────────────────────────────────────────────

_COORD = re.compile(r"([XYIJ])(-?\d+)")
_AD = re.compile(r"ADD(\d+)([A-Za-z_$][\w.$]*),?(.*)")


def parse(path: str) -> GerberFile:
    gf = GerberFile(path=str(path))

    int_digits, dec_digits = 4, 6      # KiCad default; overridden by %FS
    scale = 10.0 ** -dec_digits
    unit_scale = 1.0                   # multiply to get mm

    x = y = 0.0
    cur_ap = 0
    dark = True                        # %LPD / %LPC state
    interp = "G01"                     # G01 | G02 | G03
    pending_function: str | None = None
    in_region = False
    region_pts: list[tuple[float, float]] = []

    with open(path, "r", errors="replace") as fh:
        text = fh.read()

    # Extended commands (%...%) can span lines; pull them out first.
    for m in re.finditer(r"%([^%]*)%", text, re.S):
        body = m.group(1).strip().rstrip("*")
        if body.startswith("FS"):
            fm = re.search(r"X(\d)(\d)", body)
            if fm:
                int_digits, dec_digits = int(fm.group(1)), int(fm.group(2))
                scale = 10.0 ** -dec_digits
        elif body.startswith("MOMM"):
            unit_scale, gf.unit = 1.0, "mm"
        elif body.startswith("MOIN"):
            unit_scale, gf.unit = 25.4, "in->mm"
        elif body.startswith("TF."):
            k, _, v = body[3:].partition(",")
            gf.attributes[k] = v

    # Now the streaming pass over data blocks.
    for raw in text.replace("%", "\n").split("*"):
        blk = raw.strip()
        if not blk or blk.startswith("G04"):
            continue

        if blk.startswith("LPD"):
            dark = True
            continue
        if blk.startswith("LPC"):
            dark = False
            continue

        if blk.startswith("TA.AperFunction,"):
            pending_function = blk.split(",", 1)[1]
            continue
        if blk.startswith("TD"):
            pending_function = None
            continue

        am = _AD.match(blk)
        if am:
            code = int(am.group(1))
            template = am.group(2)
            params = [float(p) for p in am.group(3).split("X") if p.strip()]
            gf.apertures[code] = ApertureDef(code, template,
                                             [p * unit_scale for p in params],
                                             pending_function)
            continue

        if blk.startswith("G36"):
            in_region, region_pts = True, []
            continue
        if blk.startswith("G37"):
            if len(region_pts) >= 3:
                gf.regions.append(region_pts)
                gf.region_dark.append(dark)
            in_region, region_pts = False, []
            continue

        if blk.startswith("G01"):
            interp = "G01"
            blk = blk[3:]
        elif blk.startswith("G02"):
            interp = "G02"
            blk = blk[3:]
        elif blk.startswith("G03"):
            interp = "G03"
            blk = blk[3:]
        elif blk.startswith("G74") or blk.startswith("G75"):
            blk = blk[3:]

        dm = re.search(r"D(\d+)$", blk)
        if dm and int(dm.group(1)) >= 10:
            cur_ap = int(dm.group(1))
            continue

        if not dm:
            continue
        op = int(dm.group(1))

        nx, ny, i_off, j_off = x, y, 0.0, 0.0
        for axis, val in _COORD.findall(blk):
            v = int(val) * scale * unit_scale
            if axis == "X":
                nx = v
            elif axis == "Y":
                ny = v
            elif axis == "I":
                i_off = v
            else:
                j_off = v

        if op == 2:                                    # move
            if in_region and len(region_pts) >= 3:
                gf.regions.append(region_pts)
                gf.region_dark.append(dark)
                region_pts = []
            x, y = nx, ny
            if in_region:
                region_pts = [(x, y)]

        elif op == 1:                                  # draw
            if in_region:
                if interp == "G01":
                    region_pts.append((nx, ny))
                else:
                    seg = Segment("arc", x, y, nx, ny, cur_ap, dark,
                                  x + i_off, y + j_off, ccw=(interp == "G03"))
                    region_pts.extend(seg.sample(24)[1:])
            else:
                if interp == "G01":
                    gf.segments.append(Segment("line", x, y, nx, ny, cur_ap,
                                               dark=dark))
                else:
                    gf.segments.append(Segment("arc", x, y, nx, ny, cur_ap,
                                               dark, x + i_off, y + j_off,
                                               ccw=(interp == "G03")))
            x, y = nx, ny

        elif op == 3:                                  # flash
            x, y = nx, ny
            gf.flashes.append(Flash(x, y, cur_ap, dark))

    return gf


# ─────────────────────────────────────────────────────────────────────────────
# Contour assembly
# ─────────────────────────────────────────────────────────────────────────────

def build_contours(segments: list[Segment], tol: float = 1e-3
                   ) -> list[list[tuple[float, float]]]:
    """Chain stroked segments end-to-end into polylines.

    Gerber emits outline segments in arbitrary order, so we greedily link each
    segment to whichever unused segment starts (or ends) at its current tip.
    """
    key = lambda p: (round(p[0] / tol), round(p[1] / tol))

    ends: dict[tuple[int, int], list[int]] = {}
    for idx, s in enumerate(segments):
        ends.setdefault(key((s.x0, s.y0)), []).append(idx)
        ends.setdefault(key((s.x1, s.y1)), []).append(idx)

    used = [False] * len(segments)
    contours = []

    for start in range(len(segments)):
        if used[start]:
            continue
        used[start] = True
        pts = segments[start].sample(24)
        # extend forwards, then backwards
        for _ in range(2):
            while True:
                tip = pts[-1]
                nxt = None
                for cand in ends.get(key(tip), []):
                    if not used[cand]:
                        nxt = cand
                        break
                if nxt is None:
                    break
                used[nxt] = True
                s = segments[nxt]
                sp = s.sample(24)
                if math.dist(sp[0], tip) > math.dist(sp[-1], tip):
                    sp = sp[::-1]
                pts.extend(sp[1:])
            pts.reverse()
        contours.append(pts)

    return contours


def polygon_area(pts: list[tuple[float, float]]) -> float:
    a = 0.0
    for i in range(len(pts)):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % len(pts)]
        a += x0 * y1 - x1 * y0
    return abs(a) / 2.0


def bbox_of(pts: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def fit_circle(pts: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    """Algebraic (Kasa) circle fit.  Returns (cx, cy, r, rms_residual)."""
    import numpy as np
    p = np.asarray(pts, dtype=float)
    A = np.column_stack([p[:, 0], p[:, 1], np.ones(len(p))])
    b = p[:, 0] ** 2 + p[:, 1] ** 2
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = sol[0] / 2, sol[1] / 2
    r = math.sqrt(sol[2] + cx * cx + cy * cy)
    resid = np.hypot(p[:, 0] - cx, p[:, 1] - cy) - r
    return cx, cy, r, float(np.sqrt((resid ** 2).mean()))
