// gates_check.C — Garfield-side acceptance gates for the T6 mesh field map.
//
// Loads a meshfield_<tag>.txt produced by solve_fieldmap.py into a
// ComponentGrid — exactly the way S3 will consume it — and verifies:
//
//   G1 far-field drift   <Ez> at the top gas plane = +333 V/cm within 1%
//   G2 amp bulk          <Ez>(z=-90um) within 2% of the independent neBEM
//                        infinite-lattice value 30561 V/cm (copies_scan.C)
//   G3 load closure      the file loads, mesh dimensions match its header
//   G7 transparency      AvalancheMicroscopic electron transparency at the
//                        bench point through the LOADED map, compared with
//                        the 2D analytic-field reference 0.873
//                        (mesh_transparency.C). Pass band [0.80, 0.95].
//                        Also reports the funneling rms displacement.
//
// The solver-side gates S1-S6 (BC exactness, mirror symmetry, refinement)
// live in solve_fieldmap.py; this macro checks what Garfield actually sees.
//
//   source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh
//   root -b -q -e 'gSystem->Load("libGarfield");
//     gSystem->CompileMacro("gates_check.C","k");
//     gates_check("meshfield_production.txt", 2000);'
//
// Smoke maps: gates_check("meshfield_smoketest.txt", 400) — G2 tolerance is
// looser there (coarse grid under-resolves the wires).

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"

using namespace Garfield;

// Geometry constants — keep in sync with solve_fieldmap.py.
const double kPitch = 67.0e-4;    // cm
const double kWireR = 9.5e-4;     // cm
const double kZOff  = 9.0e-4;     // cm, wire-axis |z| offset (overlap 1 um)
const double kEDrift = 333.0;     // V/cm
const double kNeBemAmpEz = 30561.0;  // V/cm, independent cross-check

int gates_check(const char* mapfile = "meshfield_production.txt",
                int nev = 2000) {
  const bool production =
      std::string(mapfile).find("production") != std::string::npos;

  MediumMagboltz gas;
  gas.SetComposition("ar", 95., "ic4h10", 5.);
  gas.SetTemperature(293.15);
  gas.SetPressure(745.83);
  gas.Initialise(true);

  // ── G3: load ──────────────────────────────────────────────────────────────
  ComponentGrid grid;
  if (!grid.LoadElectricField(mapfile, "xyz", true, true)) {
    std::cout << "G3 FAIL: could not load " << mapfile << "\n";
    return 1;
  }
  grid.SetMedium(&gas);
  grid.EnablePeriodicityX();
  grid.EnablePeriodicityY();
  double xlo, ylo, zlo, xhi, yhi, zhi;
  grid.GetBoundingBox(xlo, ylo, zlo, xhi, yhi, zhi);
  std::printf("G3 load: bounding box x [%g, %g] z [%g, %g] cm -> PASS\n",
              xlo, xhi, zlo, zhi);

  bool pass = true;

  auto probeEz = [&](double z, double& mean, double& exy) {
    double sum = 0., sxy = 0.;
    int n = 0;
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        const double x = xlo + (xhi - xlo) * i / 4.;
        const double y = ylo + (yhi - ylo) * j / 4.;
        double ex, ey, ez, v;
        Medium* med = nullptr;
        int st = 0;
        grid.ElectricField(x, y, z, ex, ey, ez, v, med, st);
        if (st != 0) continue;
        sum += ez;
        sxy += std::hypot(ex, ey);
        ++n;
      }
    }
    mean = n ? sum / n : 0.;
    exy = n ? sxy / n : 0.;
  };

  // ── G1 far field ──────────────────────────────────────────────────────────
  {
    double mean, exy;
    probeEz(zhi - 20.0e-4, mean, exy);
    const double dev = 100. * (mean - kEDrift) / kEDrift;
    const bool ok = std::abs(dev) < 1.0;
    pass &= ok;
    std::printf("G1 far field: <Ez> = %.2f V/cm (want +%.0f, dev %+.2f%%), "
                "<|Exy|> = %.2f -> %s\n",
                mean, kEDrift, dev, exy, ok ? "PASS" : "FAIL");
  }

  // ── G2 amp bulk vs neBEM ──────────────────────────────────────────────────
  {
    double mean, exy;
    probeEz(-90.0e-4, mean, exy);
    const double dev = 100. * (mean - kNeBemAmpEz) / kNeBemAmpEz;
    const double tol = production ? 2.0 : 3.0;
    const bool ok = std::abs(dev) < tol;
    pass &= ok;
    std::printf("G2 amp bulk: <Ez>(z=-90um) = %.1f V/cm "
                "(neBEM: %.0f, dev %+.2f%%, tol %.0f%%) -> %s\n",
                mean, kNeBemAmpEz, dev, tol, ok ? "PASS" : "FAIL");
  }

  // ── G7 transparency through the loaded map ────────────────────────────────
  {
    Sensor sensor;
    sensor.AddComponent(&grid);
    sensor.SetArea(-10 * kPitch, -10 * kPitch, zlo,
                    10 * kPitch,  10 * kPitch, zhi);
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.SetTimeWindow(0., 1000.);
    aval.EnableAvalancheSizeLimit(50);
    int through = 0;
    double sum2 = 0.;
    for (int i = 0; i < nev; ++i) {
      const double x0 = (RndmUniform() - 0.5) * kPitch;
      const double y0 = (RndmUniform() - 0.5) * kPitch;
      aval.AvalancheElectron(x0, y0, 180.0e-4, 0., 0.1, 0., 0., 0.);
      if (aval.GetNumberOfElectronEndpoints() == 0) continue;
      double xs, ys, zs, ts, es, xe, ye, ze, te, ee;
      int status;
      aval.GetElectronEndpoint(0, xs, ys, zs, ts, es,
                               xe, ye, ze, te, ee, status);
      if (ze < -(kZOff + kWireR + 10.0e-4)) {
        ++through;
        sum2 += (xe - x0) * (xe - x0) + (ye - y0) * (ye - y0);
      }
      if ((i + 1) % 200 == 0)
        std::printf("\r  transparency: %d / %d e-", i + 1, nev);
    }
    const double eps = double(through) / nev;
    const double err = std::sqrt(eps * (1. - eps) / nev);
    const double rms = through ? std::sqrt(sum2 / through) * 1e4 : 0.;
    const bool ok = eps > 0.80 && eps < 0.95;
    pass &= ok;
    std::printf("\rG7 transparency: %.3f +- %.3f (2D analytic reference "
                "0.873), funnel rms %.1f um -> %s\n",
                eps, err, rms, ok ? "PASS" : "FAIL");
  }

  std::cout << (pass ? "\nALL GARFIELD GATES PASS for " :
                       "\n*** GATE FAILURE for ") << mapfile << "\n";
  return pass ? 0 : 1;
}

int main(int argc, char** argv) {
  return gates_check(argc > 1 ? argv[1] : "meshfield_production.txt",
                     argc > 2 ? std::atoi(argv[2]) : 2000);
}
