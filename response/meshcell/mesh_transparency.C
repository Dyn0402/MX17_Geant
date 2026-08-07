// mesh_transparency.C — T6 first pass: electron transparency of the MX17
// micromesh from its geometry, plan §4 items (a) and (d).
//
// GEOMETRY (design/RESPONSE_SIM_PLAN.md §12 item 10, resolved 2026-08-07):
// 19 um wires on a 67 um pitch, 150 um amplification gap, 30 mm drift gap.
// That weave is the standard 400 lpi / 18 um mesh scaled by 5.5 % — identical
// d/pitch, hence identical optical transparency 0.5133.
//
// WHAT THIS DOES AND DOES NOT DO.  A woven mesh is 3D: two orthogonal wire
// sets that interleave in z.  This first pass models ONE wire set as a
// periodic 2D array (ComponentAnalyticField, x-periodic), which is exact for
// the field of parallel wires and is the standard first approximation to a
// weave.  It gives the transparency's dependence on the field ratio and its
// funnelling scale, both of which are set by the wire-to-wire focusing, not by
// the over/under of the weave.  The genuine 3D woven solve is the neBEM
// follow-up; this exists to get a curve, to fix the run machinery, and to give
// that solve something to be checked against.
//
// Transparency here is ELECTRON transparency: the fraction of drifting
// electrons that reach the amplification side rather than terminating on a
// wire.  Scanned against E_amp/E_drift, which is the variable the literature
// bulk-MM curves use.
//
//   root -b -q mesh_transparency.C
// or compiled against the pinned Garfield (§5a).

#include <TApplication.h>
#include <fstream>
#include <iostream>
#include <cmath>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

// ── MX17 mesh, from shared/MX17ModuleGeometry.hh ────────────────────────────
const double kWireDiam = 19.0e-4;   // cm
const double kPitch    = 67.0e-4;   // cm  (19 wire + 48 opening)
const double kAmpGap   = 150.0e-4;  // cm  CONFIRMED, not a scan axis
const double kDriftTop = 0.30;      // cm  (a slice of the 30 mm gap is enough:
                                    //      the field is uniform until ~1 pitch
                                    //      from the mesh)

int mesh_transparency() {
  // Ar/iC4H10 95/5 at Saclay pressure, matching the S3 campaign (§5).
  MediumMagboltz gas;
  gas.SetComposition("ar", 95., "ic4h10", 5.);
  gas.SetTemperature(293.15);
  gas.SetPressure(745.83);
  gas.Initialise(true);

  std::ofstream out("mesh_transparency.csv");
  out << "E_drift_Vcm,E_amp_Vcm,ratio,n_launched,n_through,transparency,"
         "rms_offset_um\n";

  const int kNev = 2000;
  const double kVmesh = 490.0;                 // bench operating point (§5)
  const double eAmp = kVmesh / kAmpGap;        // V/cm across the 150 um gap

  for (double eDrift : {50., 100., 200., 333., 500., 700., 1000.}) {
    ComponentAnalyticField cmp;
    cmp.SetMedium(&gas);
    // One wire per period, periodic in x. The wire sits at the origin plane.
    cmp.AddWire(0., 0., kWireDiam, 0., "m", 100.);
    cmp.SetPeriodicityX(kPitch);
    // Plates: anode below (the ESL surface), drift cathode above.
    cmp.AddPlaneY(-kAmpGap, eAmp * kAmpGap, "anode");
    cmp.AddPlaneY(kDriftTop, -eDrift * kDriftTop, "cathode");

    Sensor sensor;
    sensor.AddComponent(&cmp);
    // MUST set the area explicitly. A 2D analytic field has NO z extent, so
    // the default bounding box is degenerate in z and the first nanometre of
    // z diffusion puts the electron "outside the drift area" — it dies after
    // ~16 um with status -1 and the transparency comes out identically zero.
    sensor.SetArea(-5 * kPitch, -kAmpGap, -1.0,
                    5 * kPitch,  kDriftTop,  1.0);

    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.SetTimeWindow(0., 1000.);
    // We want transparency, not gain. Capping the avalanche keeps the cost
    // flat: the seed electron's fate is decided at the mesh, and whatever it
    // does afterwards costs CPU without adding information.
    aval.EnableAvalancheSizeLimit(50);

    int through = 0;
    double sum2 = 0.;
    for (int i = 0; i < kNev; ++i) {
      // Launch one drift-gap-height above, uniformly across the period, so the
      // sampling is over the mesh cell rather than aimed at a hole.
      const double x0 = (RndmUniform() - 0.5) * kPitch;
      const double y0 = kDriftTop * 0.5;
      aval.AvalancheElectron(x0, y0, 0., 0., 0.1, 0., 0., 0.);
      if (aval.GetNumberOfElectronEndpoints() == 0) continue;
      double xs, ys, zs, ts, es, xe, ye, ze, te, ee;
      int status;
      aval.GetElectronEndpoint(0, xs, ys, zs, ts, es,
                               xe, ye, ze, te, ee, status);
      // "Through" means the seed got past the wire plane. Not a deep cut:
      // in a 32 kV/cm amplification field the seed ionises within microns of
      // the mesh, so its own endpoint is only just below it. An electron
      // caught by a wire instead ends on the wire surface, at |y| <= r_wire.
      if (ye < -kWireDiam) {
        ++through;
        sum2 += xe * xe;
      }
    }
    const double eps = double(through) / kNev;
    const double rms = through ? std::sqrt(sum2 / through) * 1e4 : 0.;
    std::cout << "E_drift=" << eDrift << " V/cm  ratio=" << eAmp / eDrift
              << "  transparency=" << eps << "  rms offset=" << rms << " um\n";
    out << eDrift << "," << eAmp << "," << eAmp / eDrift << "," << kNev << ","
        << through << "," << eps << "," << rms << "\n";
  }
  out.close();
  std::cout << "wrote mesh_transparency.csv\n";
  return 0;
}

int main(int argc, char** argv) { return mesh_transparency(); }
