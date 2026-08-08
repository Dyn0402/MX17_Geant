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
// HANDEDNESS IS THE FREE CHOICE HERE, NOT THE AXIS DIRECTIONS. An earlier
// version flipped y (sy = -1) purely so that (x, y, z)_local came out
// right-handed after the z flip. That was the wrong trade: nothing in the
// response chain takes a cross product, so handedness costs nothing, while a
// flipped y silently mislabels every channel. The local frame IS left-handed
// and that is deliberate.
//
// ── SETTLED 2026-08-08 (was task T2) ────────────────────────────────────────
// The convention, from the detector itself (user, 2026-08-08):
//
//   Viewed from the PIXEL (readout pad) side:
//     * the X connectors are along the bottom, the Y connectors up the right;
//     * x counts left -> right, y counts bottom -> top;
//     * x FEU DREAM 1 is at the far LEFT, y FEU DREAM 1 at the BOTTOM,
//       i.e. channel 0 sits at the LOW end of each coordinate.
//
// That matches nTof_x17/common/Mx17StripMap.py and mx17_m1_map.csv exactly:
// axis 'x' channels carry x = 0.78 mm * channel_num increasing from 0, axis 'y'
// likewise in y. So the strip map's "x0 = 0, increasing" is not an assumption
// to be checked — it IS the wiring.
//
// DERIVING sx, sy from that, without needing to agree on which way is "up":
// the pads sit at LARGER world z than the ESL (the beam runs +z through
// window -> mesh -> gap -> ESL -> kapton -> pads), so an observer on the pixel
// side is at large +z_world looking back along -z_world, and +z_world points
// TOWARD them. For any viewer (right, up, toward-viewer) is right-handed, so
// (x_det, y_det, +z_world) is right-handed and therefore
// (x_det, y_det, z_local) is LEFT-handed, since z_local = -z_world. Writing
// x_local = sx*x_world and y_local = sy*y_world,
//
//     x_local ^ y_local = sx*sy * z_world = -sx*sy * z_local
//
// and left-handedness requires that be -z_local, i.e. sx*sy = +1. With the
// active area unrotated in world x, sx = +1, hence **sy = +1**.
//
// Channel 0 then lands where the hardware puts it: `digitize.induce` books
// col = rint((x_local - PAD_ORIGIN)/pitch) with PAD_ORIGIN = -199.29 mm, so
// col 0 is at the most negative x_local = the far left, and row 0 at the most
// negative y_local = the bottom. No offset or flip is applied anywhere.
//
// THE REMAINING RELABELLING IS ON THE DATA SIDE, and stays there: real runs
// need each 64-channel connector reversed, which `Mx17StripMap.apply_orientation`
// already does unconditionally. The simulation keeps the identity mapping and
// the comparison happens after that inversion (agreed 2026-08-08).
//
// ⚠️ FILES WRITTEN BEFORE 2026-08-08 CARRY sy = -1 IN THEIR Meta TREE and are
// still self-consistent when read through it (clusters.py takes sy from Meta,
// never from a constant). But their y is mirrored relative to new files, so
// old and new must not be mixed in any channel-by-channel comparison.
// Observables symmetric in y — the c1/c2 shares, which average d = +-1, and the
// X/Y balance — are unaffected.
//
// THE BENCH FRAME IS A THIRD FRAME, and a swapped one. On the cosmic bench the
// module goes in such that detector +y runs along bench +x and detector +x along
// bench +y (user, 2026-08-08) — which is why the June M3 alignment only
// converged with a ~90 deg rotation plus a flip. See kBenchFromDetector below;
// nothing in Stage A or B uses it, it exists so nobody re-derives it.
// ─────────────────────────────────────────────────────────────────────────────

namespace MX17 {

struct ActiveAreaFrame {
    // World coordinates of the active-area origin [mm].
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    // Axis signs, world → local. sy = +1 since 2026-08-08 (see above); the
    // frame is deliberately left-handed.
    double sx = 1.0, sy = 1.0, sz = -1.0;
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

// Detector -> cosmic-bench frame: the module is mounted rotated, so bench x is
// detector y and bench y is detector x. Recorded for whoever compares against
// M3; the response chain itself never leaves the detector frame.
struct BenchFromDetector {
    // (bench_x, bench_y) = (det_y, det_x)
    static constexpr bool swapXY = true;
};

}  // namespace MX17
