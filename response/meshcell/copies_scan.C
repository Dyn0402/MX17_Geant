// copies_scan.C — convergence probe for field_map.C: how many neBEM periodic
// copies until the unit-cell solve reaches the infinite-lattice field?
// Prints Ez at the drift-slice top and amp-gap middle vs number of copies.
#include <cmath>
#include <iostream>

#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/SolidWire.hh"

using namespace Garfield;

const double kWireR      = 9.5e-4;
const double kPitch      = 67.0e-4;
const double kAmpGap     = 150.0e-4;
const double kDriftSlice = 200.0e-4;
const double kVmesh  = 490.0;
const double kEDrift = 333.0;
const double kZAnode = -(2. * kWireR + kAmpGap);
const double kZCath  = +(2. * kWireR + kDriftSlice);

int copies_scan() {
  MediumMagboltz gas;
  gas.SetComposition("ar", 95., "ic4h10", 5.);
  gas.SetTemperature(293.15);
  gas.SetPressure(745.83);
  gas.Initialise(false);
  MediumConductor metal;

  for (int nc : {2, 5, 10, 15, 20, 30}) {
    GeometrySimple geo;
    geo.SetMedium(&gas);
    SolidWire wx(0., 0., +kWireR, kWireR, kPitch / 2., 1, 0, 0);
    SolidWire wy(0., 0., -kWireR, kWireR, kPitch / 2., 0, 1, 0);
    wx.SetBoundaryPotential(0.);
    wy.SetBoundaryPotential(0.);
    wx.SetDiscretisationLevel(kWireR);
    wy.SetDiscretisationLevel(kWireR);
    geo.AddSolid(&wx, &metal);
    geo.AddSolid(&wy, &metal);
    const double plateHalfT = 1.0e-4;
    SolidBox anode(0., 0., kZAnode - plateHalfT, kPitch / 2., kPitch / 2.,
                   plateHalfT);
    anode.SetBoundaryPotential(kVmesh);
    anode.SetDiscretisationLevel(kPitch / 4.);
    geo.AddSolid(&anode, &metal);
    SolidBox cathode(0., 0., kZCath + plateHalfT, kPitch / 2., kPitch / 2.,
                     plateHalfT);
    cathode.SetBoundaryPotential(-kEDrift * kZCath);
    cathode.SetDiscretisationLevel(kPitch / 4.);
    geo.AddSolid(&cathode, &metal);

    ComponentNeBem3d nebem;
    nebem.SetGeometry(&geo);
    nebem.SetTargetElementSize(kPitch / 4.);
    nebem.SetMinMaxNumberOfElements(4, 30);
    nebem.SetPeriodicityX(kPitch);
    nebem.SetPeriodicityY(kPitch);
    nebem.SetPeriodicCopies(nc, nc, 0);
    nebem.UseLUInversion();
    nebem.SetNumberOfThreads(6);
    if (!nebem.Initialise()) {
      std::cerr << "init failed at nc=" << nc << "\n";
      continue;
    }
    double ex, ey, ez1, ez2, v;
    Medium* med = nullptr;
    int st = 0;
    nebem.ElectricField(3.0e-4, 5.0e-4, kZCath - 2.0e-4, ex, ey, ez1, v, med,
                        st);
    nebem.ElectricField(3.0e-4, 5.0e-4, -90.0e-4, ex, ey, ez2, v, med, st);
    std::cout << "COPIES " << nc << ": Ez(top) = " << ez1
              << " (want ~+333), Ez(amp) = " << ez2 << " (want ~+"
              << kVmesh / std::abs(kZAnode) << ")\n";
  }
  return 0;
}
int main() { return copies_scan(); }
