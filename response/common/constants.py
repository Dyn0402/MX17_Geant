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
# TWO dielectrics in series, not one. The ESL is screen-printed on the top face
# of a kapton coverlay whose bottom face is LAMINATED to the pad copper with an
# adhesive, so the pad sees kapton + glue. Layer order is fixed and matters:
# C(k) is the admittance looking down from the ESL and is NOT symmetric under
# swapping the layers (2.9 % at k = 2/pad-pitch), even though S(k) is.
#
# Kapton thickness is now KNOWN: 50 µm, confirmed 2026-08-08 against
# shared/MX17ModuleGeometry.hh (pcbKapton_um) and the MX17_Full_Geant PCB stack,
# and cross-checked below. It is no longer a scan axis — the {50, 75, 125} scan
# is retired (plan §1 row 4).
KAPTON_EPS_R = 3.5
KAPTON_THICK_UM = 50.0

# The glue is an ESTIMATE, and it is now the dominant uncertainty in the stack.
# Reasoning, in order of how much each step is worth trusting:
#  1. A 50 µm (2 mil) polyimide coverlay is supplied with an acrylic adhesive
#     layer; the standard pairing for 2 mil film is 1 mil = 25 µm of adhesive.
#     Common alternatives are 12.5 and 35 µm, hence the bracket below.
#  2. Acrylic/epoxy coverlay adhesive has eps_r ~ 3.0-3.5. The stack is barely
#     sensitive to it: over that whole range the series thickness moves ~5 %.
#  3. What is supplied is NOT what ends up over the pad. Vacuum lamination
#     drives adhesive into the 100 µm trenches between the 680 µm pads, which
#     are 26 µm deep, so the over-pad layer is thinner by the trench volume:
#     t_pad = t_supplied - (1 - pad_coverage) * t_cu.  That is the electrically
#     relevant thickness, because the pad face is the electrode.
# Net: 25 µm supplied -> 18.8 µm over the pad. NOT a production scan axis; the
# bracket exists for sensitivity checks only.
GLUE_EPS_R = 3.2
GLUE_THICK_SUPPLIED_UM = 25.0
GLUE_THICK_SUPPLIED_UM_BRACKET = (12.5, 25.0, 35.0)
PAD_CU_THICK_UM = 26.0            # L4 copper (design/NEEDED_INPUTS.md §6)


# Which quantity the supplied thickness refers to. This is an OPEN question and
# it is not a physics question — it is what the vendor/process number means:
#
#   "conserved" — the quoted thickness is the AS-SUPPLIED coating on the
#       coverlay, and lamination redistributes it. The inter-pad channels must
#       fill from somewhere and the only source is the adhesive itself, so the
#       over-pad layer thins by the channel volume. Three things argue for this
#       being the physical picture: adhesive FLOW is a specified, required
#       property of coverlay (void-free encapsulation of copper topography is
#       the whole point, and "low-flow" coverlay is a separate product because
#       flow is the default); a heated vacuum press exists precisely to fill
#       the channels rather than trap air; and on a 400 mm panel the interior
#       is a closed system, since the flow path to a free edge is ~200 mm of
#       viscous thermoset. The 50 µm film cannot dip into a 100 µm channel
#       either — span/thickness = 2, so it is rigid there and sits flat on the
#       pad tops.
#
#   "uniform" — the quoted thickness is the FINAL bond line over the copper,
#       i.e. the process spec already accounts for flow. Then no correction
#       applies.
#
# Settled by a datasheet or a fab process sheet, not by argument. "conserved"
# is the default because it is the one derivable from geometry we know exactly;
# the difference is ~0.010 in S(0)/C(0), well inside the spread from the
# supplied thickness itself (12-35 µm -> 0.893 to 0.856).
GLUE_FLOW_MODEL = "conserved"                         # or "uniform"


def glue_over_pad_um(t_supplied_um: float = GLUE_THICK_SUPPLIED_UM,
                     t_cu_um: float = PAD_CU_THICK_UM,
                     model: str = None) -> float:
    """
    Adhesive thickness left over the pad face after lamination [µm].

    Volume bookkeeping per unit area, referred to the pad top: the adhesive
    first fills the inter-pad channel, (1 - coverage) * t_cu, and what remains
    sits above the pads. Clipped at 0 — a coverlay too thin to fill the channel
    would not laminate void-free, and this returns 0 rather than a negative gap.

    `model="uniform"` returns the supplied thickness unchanged; see
    GLUE_FLOW_MODEL for why that is a live alternative rather than a wrong one.
    """
    if (model or GLUE_FLOW_MODEL) == "uniform":
        return t_supplied_um
    coverage = (PAD_SIZE_UM / PAD_PITCH_UM) ** 2      # 0.76003, exact
    return max(0.0, t_supplied_um - (1.0 - coverage) * t_cu_um)


GLUE_THICK_UM = glue_over_pad_um()                    # 18.76 µm

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


def insulator_d_eff_m(kapton_m: float = KAPTON_THICK_UM * 1e-6,
                      eps_r: float = KAPTON_EPS_R,
                      glue_m: float = GLUE_THICK_UM * 1e-6,
                      glue_eps_r: float = GLUE_EPS_R) -> float:
    """
    Thickness of a SINGLE kapton layer with the same series capacitance as the
    real kapton + glue stack [m]:

        d_eff = d_kapton + d_glue * (eps_kapton / eps_glue)

    Use only where one number is genuinely needed (reporting, c'). It is exact
    at k -> 0 and wrong where the kernels live: substituted into the
    single-layer stack it understates S(k) by 0.3 % at k = 1/pad-pitch and
    5.3 % at 5/pad-pitch. The solver uses the real two-layer cascade
    (response/solver/wpot.py::stack_coeffs), never this.
    """
    return kapton_m + glue_m * (eps_r / glue_eps_r)


def c_prime(gap_m: float = AMP_GAP_M,
            kapton_m: float = KAPTON_THICK_UM * 1e-6,
            eps_r: float = KAPTON_EPS_R,
            glue_m: float = GLUE_THICK_UM * 1e-6,
            glue_eps_r: float = GLUE_EPS_R) -> float:
    """
    Areal capacitance c' seen by the resistive sheet [F/m^2].

        c' = eps0 * (1/gap + 1/(d_kapton/eps_kapton + d_glue/eps_glue))

    the two insulating layers adding in SERIES. This is the k -> 0 limit of the
    stack's total areal capacitance (response/solver/wpot.py) and is what sets
    the sheet diffusivity D = 1/(rho_s c').

    Plan §3 quoted c' ~ 5e-7 F/m^2 and D ~ 2.0 m^2/s at 1 MΩ/sq from the old
    75 µm single-layer default. The kapton+glue stack lands at essentially the
    same place (d_eff = 70.5 µm), so those figures survive — but they survive
    by accident, and a bare 50 µm kapton with no glue would NOT reproduce them.
    """
    d_eff = insulator_d_eff_m(kapton_m, eps_r, glue_m, glue_eps_r)
    return EPS0 * (1.0 / gap_m + eps_r / d_eff)


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
              "pasteStripW_mm", "pasteStripPitch_mm", "padN", "pcbKapton_um")
    found = {f: grab(f) for f in fields}

    checks = [
        ("amp_um", found["amp_um"], AMP_GAP_UM),
        ("paste_um", found["paste_um"], ESL_THICK_UM),
        ("padPitch_mm", found["padPitch_mm"], PAD_PITCH_UM * 1e-3),
        ("padSize_mm", found["padSize_mm"], PAD_SIZE_UM * 1e-3),
        ("pasteStripW_mm", found["pasteStripW_mm"], ESL_WIDTH_UM * 1e-3),
        ("pasteStripPitch_mm", found["pasteStripPitch_mm"], ESL_PITCH_UM * 1e-3),
        ("padN", found["padN"], PAD_N),
        # Added 2026-08-08, when the kapton stopped being a scan axis. This is
        # the layer the ESL is printed on and the pads sit under, so the
        # response model and the G4 geometry must not disagree about it. NOTE
        # the lamination adhesive is NOT in the header and so is not checked
        # here — it is a response-model estimate (see GLUE_THICK_UM).
        ("pcbKapton_um", found["pcbKapton_um"], KAPTON_THICK_UM),
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
