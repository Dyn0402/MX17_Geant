#pragma once
// ActiveAreaFrame.hh — the world → active-area transform, recorded once at
// geometry construction and written into every output file's Meta tree.
//
// Stage A of design/RESPONSE_SIM_PLAN.md §6 item 2. The response chain
// (response/digitizer) needs cluster positions in the detector's own frame,
// not in the Geant4 world frame, and it must not have to hard-code the module
// placement — that placement changes with run mode and with any future
// geometry edit. So the geometry reports its own frame.
//
// DEFINITION of the active-area frame (plan §1 and §6):
//
//   origin : the centre of the 399.36 mm active area, ON the top surface of
//            the ESL resistive layer (= the downstream face of the
//            amplification gap). This is z = 0 for the whole response chain.
//   +z     : toward the micromesh, i.e. ANTI-parallel to the Geant4 world +z,
//            since the beam travels +z through window → mesh → gap → ESL.
//   +x,+y  : lateral axes of the active area.
//
// so, componentwise,
//
//   x_local = sx * (x_world - x0)
//   y_local = sy * (y_world - y0)
//   z_local = sz * (z_world - z0)      with sz = -1
//
// A pure z flip alone would make the frame left-handed, so one lateral axis is
// flipped with it; y is the one chosen (sy = -1, sx = +1).
//
// ── OPEN ITEM, read before trusting x/y ──────────────────────────────────────
// The plan asks for "x/y per strip-map convention
// (nTof_x17/common/Mx17StripMap.py)". That convention has NOT been checked
// against this frame yet — it is task T2, which also settles the pad↔X/Y
// bussing question. Until T2 lands, treat sx/sy and the x↔y assignment as a
// PLACEHOLDER: the magnitude of every lateral coordinate is right, the labels
// and signs may need a flip. Nothing upstream of the digitizer depends on it,
// and the Meta tree carries the numbers explicitly so a later correction is a
// one-line change in the digitizer, not a re-run of Geant4.
// ─────────────────────────────────────────────────────────────────────────────

namespace MX17 {

struct ActiveAreaFrame {
    // World coordinates of the active-area origin [mm].
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    // Axis signs, world → local.
    double sx = 1.0, sy = -1.0, sz = -1.0;
    // Active-area full width [mm] (512 × 0.78).
    double activeWidth_mm = 399.36;
    // True once the geometry has actually reported itself.
    bool valid = false;
};

// Single instance, written by DetectorConstruction::Construct() before any
// event is processed and read by RunAction at BeginOfRunAction. Not guarded:
// it is written once during single-threaded geometry construction and only
// read afterwards.
inline ActiveAreaFrame& TheActiveAreaFrame() {
    static ActiveAreaFrame frame;
    return frame;
}

}  // namespace MX17
