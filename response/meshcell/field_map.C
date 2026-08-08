// field_map.C — T6: export a 3D E-field map of the real woven micromesh unit
// cell (amp gap + a slice of drift), replacing the uniform_field placeholder
// that S3 has been using. This is the genuine 3D solve that mesh_transparency.C
// (its 1D periodic-wire-array approximation) explicitly deferred to this file.
//
// GEOMETRY. Two orthogonal wire sets approximate the weave as a stacked
// bi-planar grid: all x-directed wires at z = +r, all y-directed wires at
// z = -r (r = wire radius), so they do not intersect. This is the standard
// simplification the plan calls "two orthogonal sinusoid-ish wires" — it
// gets the wire-to-wire focusing (what sets transparency and funnelling)
// right without modelling the physical over/under crossing of each wire.
// A finite patch of NCELL x NCELL periods stands in for the periodic weave
// (see plan §4; ComponentNeBem3d's own periodicity API exists but is not
// exercised by any shipped example, so — like MWPCbyneBEM.C, the one
// upstream example that builds a repeated-wire structure — this uses
// explicit finite copies rather than SetPeriodicityX/Y, which is the
// lower-risk choice to get right on a first pass).
//
// Two plates close the cell in z: an anode plate at the ESL surface
// (z = -ampGap, the reference "ground" per plan §1) and a cathode plate at
// z = +driftSlice standing in for the drift field over the last 200 um of
// drift (plan §2 file contract). Potentials follow mesh_transparency.C's
// convention exactly (mesh wires at 0 V, ESL at +E_amp*ampGap, drift plate
// at -E_drift*driftSlice) since that convention already reproduces the
// expected transparency trend.
//
// OUTPUT. Plain-text grid, one row per (x,y,z) sample:
//   x_cm,y_cm,z_cm,Ex_Vcm,Ey_Vcm,Ez_Vcm,V_volt,status
// covering one central unit cell in x,y and the full amp gap + drift slice
// in z, i.e. exactly the plan §2 `meshfield` file contract's domain. This
// feeds S3 (AvalancheMicroscopic can load a ComponentGrid built from it) in
// place of ComponentConstant's uniform_field.
//
// WHAT THIS DOES NOT DO YET, so it is not mistaken for a finished product:
//   * Pillars (1.3% coverage) are not in this cell.
//   * The wire z-offset is a fixed +-r stack, not a true woven crossing.
//   * NCELL and the grid resolution here are a SMOKE-TEST size (fast on a
//     laptop). Production resolution is a separate, larger run on the
//     desktop once this passes its own sanity check below.
//
// SANITY CHECK, printed at the end (do not trust the map without it): deep
// in the drift slice, away from every wire, the field must be uniform and
// equal to E_drift in -z (matching mesh_transparency.C's sign convention),
// to a few percent. This is the same kind of closed-form-adjacent check
// every other stage in this plan was required to pass before being trusted
// (S1's sum rule, T10's charge budget, etc.) — an ungrounded field solve is
// exactly the kind of silent-zero failure mode S3's own induced-current bug
// (plan §5) already caught once.
//
//   source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh
//   root -b -q -e 'gSystem->Load("libGarfield");
//                  gSystem->CompileMacro("field_map.C","k"); field_map();'

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/SolidWire.hh"

using namespace Garfield;

// ── MX17 mesh + gap, from shared/MX17ModuleGeometry.hh (cm) ────────────────
const double kWireDiam  = 19.0e-4;   // cm
const double kWireR     = kWireDiam / 2.;
const double kPitch     = 67.0e-4;   // cm (19 wire + 48 opening)
const double kAmpGap    = 150.0e-4;  // cm — CONFIRMED, not a scan axis
const double kDriftSlice = 200.0e-4; // cm — plan §2 file contract

// Smoke-test sizing (production sizing is a separate, larger desktop run).
const int kNCell = 3;           // wires per direction: NCell x NCell patch
const double kVmesh = 490.0;    // bench operating point (§5), mesh at 0 V ref
const double kEDrift = 333.0;   // V/cm, det3/T14 target field (§0a P1)

int field_map() {
  MediumMagboltz gas;
  gas.SetComposition("ar", 95., "ic4h10", 5.);
  gas.SetTemperature(293.15);
  gas.SetPressure(745.83);
  gas.Initialise(true);

  MediumConductor metal;

  GeometrySimple geo;
  geo.SetMedium(&gas);

  const double eAmp = kVmesh / kAmpGap;  // V/cm across the 150 um gap

  // Wire patch: half-extent covers kNCell periods either side of the origin.
  const double half = 0.5 * (kNCell - 1) * kPitch;
  const double wireHalfLength = half + kPitch;  // a bit longer than the patch

  std::vector<SolidWire*> wires;
  for (int i = 0; i < kNCell; ++i) {
    const double pos = -half + i * kPitch;
    // x-directed wires, spaced in y, at z = +r.
    auto* wx = new SolidWire(0., pos, kWireR, kWireR, wireHalfLength, 1, 0, 0);
    wx->SetBoundaryPotential(0.);
    geo.AddSolid(wx, &metal);
    wires.push_back(wx);
    // y-directed wires, spaced in x, at z = -r.
    auto* wy = new SolidWire(pos, 0., -kWireR, kWireR, wireHalfLength, 0, 1, 0);
    wy->SetBoundaryPotential(0.);
    geo.AddSolid(wy, &metal);
    wires.push_back(wy);
  }

  // Plates, sized generously larger than the wire patch so their own edges
  // do not distort the field over the sampled central cell.
  const double plateHalf = half + 3 * kPitch;
  SolidBox anode(0., 0., -kAmpGap, plateHalf, plateHalf, 1.0e-4);
  anode.SetBoundaryPotential(eAmp * kAmpGap);
  anode.SetLabel("esl");
  geo.AddSolid(&anode, &metal);

  SolidBox cathode(0., 0., kDriftSlice, plateHalf, plateHalf, 1.0e-4);
  cathode.SetBoundaryPotential(-kEDrift * kDriftSlice);
  cathode.SetLabel("drift");
  geo.AddSolid(&cathode, &metal);

  ComponentNeBem3d nebem;
  nebem.SetGeometry(&geo);
  nebem.SetTargetElementSize(kWireR);  // resolve the wire circumference
  nebem.SetMinMaxNumberOfElements(4, 30);
  nebem.UseSVDInversion();
  std::cout << "Initialising neBEM solve (" << wires.size()
            << " wires + 2 plates)...\n";
  if (!nebem.Initialise()) {
    std::cerr << "neBEM initialisation FAILED.\n";
    return 1;
  }
  std::cout << "neBEM solve done.\n";

  // ── Export grid: one central unit cell in x,y; full amp gap + drift slice
  //    in z. Smoke-test resolution — coarsen/refine independently of the
  //    geometry patch size above.
  const int nx = 11, ny = 11, nz = 15;
  const double xlo = -kPitch / 2, xhi = kPitch / 2;
  const double ylo = -kPitch / 2, yhi = kPitch / 2;
  const double zlo = -kAmpGap, zhi = kDriftSlice;

  std::ofstream out("field_map_smoketest.csv");
  out << "x_cm,y_cm,z_cm,Ex_Vcm,Ey_Vcm,Ez_Vcm,V_volt,status\n";

  double farfield_ez_sum = 0.;
  int farfield_n = 0;

  for (int iz = 0; iz < nz; ++iz) {
    const double z = zlo + (zhi - zlo) * iz / (nz - 1);
    for (int iy = 0; iy < ny; ++iy) {
      const double y = ylo + (yhi - ylo) * iy / (ny - 1);
      for (int ix = 0; ix < nx; ++ix) {
        const double x = xlo + (xhi - xlo) * ix / (nx - 1);
        double ex, ey, ez, v;
        Medium* medium = nullptr;
        int status = 0;
        nebem.ElectricField(x, y, z, ex, ey, ez, v, medium, status);
        out << x << "," << y << "," << z << "," << ex << "," << ey << ","
            << ez << "," << v << "," << status << "\n";
        // Deep in the drift slice (last grid plane), away from the mesh:
        // collect Ez for the sanity check below.
        if (iz == nz - 1) {
          farfield_ez_sum += ez;
          ++farfield_n;
        }
      }
    }
  }
  out.close();
  std::cout << "wrote field_map_smoketest.csv (" << nx * ny * nz
            << " points)\n";

  // ── Sanity check: far into the drift slice, Ez should be uniform and
  //    equal to -E_drift (mesh_transparency.C's sign convention: field
  //    points from mesh, at 0 V, toward the more-negative cathode).
  const double ez_mean = farfield_ez_sum / farfield_n;
  const double expect = -kEDrift;
  const double pct = 100. * (ez_mean - expect) / expect;
  std::cout << "\nSANITY CHECK (far-field limit at z = " << zhi << " cm):\n"
            << "  measured  <Ez> = " << ez_mean << " V/cm\n"
            << "  expected  -E_drift = " << expect << " V/cm\n"
            << "  deviation = " << pct << " %\n";
  if (std::abs(pct) > 5.0) {
    std::cout << "  *** FAILS the 5% far-field bar. Do not trust this map "
                  "until this is understood (undersized plates? patch too "
                  "small? wrong sign convention?). ***\n";
  } else {
    std::cout << "  within 5% — the boundary conditions are self-consistent."
                  " Proceed to widen the patch/grid for production.\n";
  }

  return 0;
}

int main(int argc, char** argv) { return field_map(); }
