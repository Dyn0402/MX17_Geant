"""
constants.py — the single Python source of truth for MX17 response geometry.

design/RESPONSE_SIM_PLAN.md §2 asks for exactly one place where the geometry
numbers live on the Python side, cross-checked against the C++ header so the
two cannot drift apart silently. `assert_against_header()` runs that check and
is called at import; set MX17_SKIP_HEADER_CHECK=1 to suppress it (only useful
off-machine, e.g. on a condor worker with no repo checkout).

Units: SI throughout (metres, seconds, farads). The _MM / _UM suffixed values
are the human-readable duplicates, and are what the assertions compare.
"""

from __future__ import annotations

import os
import re
import warnings

# ── Physical constants ───────────────────────────────────────────────────────
EPS0 = 8.8541878128e-12          # F/m

# ── Amplification region ─────────────────────────────────────────────────────
# CONFIRMED by the user 2026-08-06; plan §1 row 2 says explicitly: not 128, and
# not a scan axis.
AMP_GAP_UM = 150.0
AMP_GAP_M = AMP_GAP_UM * 1e-6

# ── ESL resistive strips ─────────────────────────────────────────────────────
# Screen-printed after fab, so NOT in the gerber set; geometry is the confirmed
# user/fab spec (plan §1 row 3, design/NEEDED_INPUTS.md §6c).
ESL_WIDTH_UM = 550.0
ESL_GAP_UM = 250.0
ESL_PITCH_UM = ESL_WIDTH_UM + ESL_GAP_UM     # 800.0
# As-built print thickness is 10 µm (user, 2026-08-07; AsBuiltSpec() sets
# s.paste_um = 10.0). It does NOT enter the response model at all — S1 treats
# the ESL as a zero-thickness sheet and absorbs the thickness into rho_s — so
# this constant exists only for the header cross-check below, which reads the
# ModuleSpec *declaration* default (still 100.0, kept because LegacySpec relies
# on it) rather than the as-built assignment.
ESL_THICK_UM = 100.0                          # header declaration default
ESL_THICK_AS_BUILT_UM = 10.0                  # AsBuiltSpec value; unused by S1
ESL_PITCH_M = ESL_PITCH_UM * 1e-6
ESL_WIDTH_M = ESL_WIDTH_UM * 1e-6

# Surface resistivity is UNKNOWN — it is a scan axis, not a constant.
# plan §1 row 3: scan {0.5, 1, 2, 5} MΩ/sq.
RHO_S_SCAN_OHM_SQ = (0.5e6, 1.0e6, 2.0e6, 5.0e6)

# ── Readout pads ─────────────────────────────────────────────────────────────
# Exact, measured flash-by-flash from DFS3498A_L2-pads.gbr.
PAD_N = 512
PAD_PITCH_UM = 780.0
PAD_SIZE_UM = 680.0
PAD_PITCH_M = PAD_PITCH_UM * 1e-6
ACTIVE_WIDTH_MM = PAD_N * PAD_PITCH_UM * 1e-3   # 399.36

# ── Insulator between pad copper and the ESL ─────────────────────────────────
# Material IS kapton (user 2026-08-06); the thickness is NOT known.
# plan §1 row 4: default 75 µm, scan {50, 75, 125}.
KAPTON_EPS_R = 3.5
KAPTON_THICK_UM_DEFAULT = 75.0
KAPTON_THICK_UM_SCAN = (50.0, 75.0, 125.0)

# ── The pitch beat ───────────────────────────────────────────────────────────
# ESL 800 µm against pads 780 µm: the registration phase advances 20 µm per
# pitch and repeats after LCM(800, 780) = 31 200 µm = 39 ESL pitches
# = 40 pad pitches. Every response kernel is periodic in x with THIS period,
# not with either single pitch — plan §1, and a prediction to look for in data.
SUPERPERIOD_UM = 31200.0
SUPERPERIOD_M = SUPERPERIOD_UM * 1e-6
N_ESL_PER_SUPER = 39
N_PAD_PER_SUPER = 40

assert abs(N_ESL_PER_SUPER * ESL_PITCH_UM - SUPERPERIOD_UM) < 1e-9
assert abs(N_PAD_PER_SUPER * PAD_PITCH_UM - SUPERPERIOD_UM) < 1e-9


def c_prime(gap_m: float = AMP_GAP_M,
            kapton_m: float = KAPTON_THICK_UM_DEFAULT * 1e-6,
            eps_r: float = KAPTON_EPS_R) -> float:
    """
    Areal capacitance c' seen by the resistive sheet [F/m^2].

        c' = eps0 * (eps_r / d_kapton + 1 / gap)

    This is the k -> 0 limit of the stack's total areal capacitance (see
    response/solver/wpot.py), and it is what sets the sheet diffusivity
    D = 1 / (rho_s c'). Plan §3 quotes c' ~ 5e-7 F/m^2 and D ~ 2.0 m^2/s at
    1 MΩ/sq, which these defaults reproduce.
    """
    return EPS0 * (eps_r / kapton_m + 1.0 / gap_m)


def sheet_diffusivity(rho_s_ohm_sq: float, **kw) -> float:
    """Sheet charge diffusivity D = 1/(rho_s c') [m^2/s]."""
    return 1.0 / (rho_s_ohm_sq * c_prime(**kw))


# ── Cross-check against the C++ header ───────────────────────────────────────
_HEADER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "shared", "MX17ModuleGeometry.hh")


def assert_against_header(path: str = _HEADER) -> dict:
    """
    Re-read the numbers this module duplicates straight out of
    shared/MX17ModuleGeometry.hh and fail loudly if they have drifted.

    Deliberately narrow: it checks only the handful of values the response
    chain actually depends on, matched by the assignment text in AsBuiltSpec().
    A header refactor that renames those fields will make this raise, which is
    the point — better a hard failure at import than a simulation quietly run
    on a stale gap.
    """
    path = os.path.normpath(path)
    if not os.path.exists(path):
        warnings.warn(f"[mx17.constants] header not found at {path}; "
                      "geometry cross-check SKIPPED")
        return {}

    src = open(path).read()

    def grab(field):
        # The ModuleSpec fields carry their values as member default
        # initialisers ("double amp_um = 150.0;"), not as assignments in
        # AsBuiltSpec(), so match the declaration.
        m = re.search(rf"\b{field}\s*=\s*([0-9.eE+-]+)\s*;", src)
        return float(m.group(1)) if m else None

    fields = ("amp_um", "paste_um", "padPitch_mm", "padSize_mm",
              "pasteStripW_mm", "pasteStripPitch_mm", "padN")
    found = {f: grab(f) for f in fields}

    checks = [
        ("amp_um", found["amp_um"], AMP_GAP_UM),
        ("paste_um", found["paste_um"], ESL_THICK_UM),
        ("padPitch_mm", found["padPitch_mm"], PAD_PITCH_UM * 1e-3),
        ("padSize_mm", found["padSize_mm"], PAD_SIZE_UM * 1e-3),
        ("pasteStripW_mm", found["pasteStripW_mm"], ESL_WIDTH_UM * 1e-3),
        ("pasteStripPitch_mm", found["pasteStripPitch_mm"], ESL_PITCH_UM * 1e-3),
        ("padN", found["padN"], PAD_N),
    ]
    bad = []
    for name, header_val, py_val in checks:
        if header_val is None:
            warnings.warn(f"[mx17.constants] '{name}' not found in the header; "
                          "cross-check for it SKIPPED (field renamed?)")
            continue
        if abs(header_val - py_val) > 1e-9:
            bad.append(f"{name}: header {header_val} != python {py_val}")
    if bad:
        raise AssertionError(
            "MX17 geometry has drifted between shared/MX17ModuleGeometry.hh "
            "and response/common/constants.py:\n  " + "\n  ".join(bad))
    return found


if not os.environ.get("MX17_SKIP_HEADER_CHECK"):
    assert_against_header()
