// heed_mip.cc -- primary ionization of a charged particle crossing the P2
// drift gas, with Garfield's Heed (PAI model, clusters, delta electrons), to
// compare with the Geant4 nPrimaryDrift spectra of the rnd2026 scan.
//
//   heed_mip <gas> <gap_mm> <theta_deg> <n_events> <out> [particle] [E_kin_eV] [seed]
//
// gas: mm_sim name (ArIso, ArIso9010, NeIso, NeIso8515, ArCO2Iso9352, ArCF4Iso).
// The gas is a slab 0 < z < gap (x, y within +-5 cm), 293.15 K, 760 Torr; the
// particle enters at z = 0 at angle theta to the normal. Delta electrons are
// transported inside the slab (their ionization counts, as in Geant4).
// Writes <out>.txt: one line per event: n_clusters n_electrons E_dep[eV].
// Unlike the Geant4 runs there is no window/cathode upstream: deltas knocked
// out of the cathode into the gas are absent here (compare the low tail).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cerr << "usage: heed_mip <gas> <gap_mm> <theta_deg> <n> <out> "
                 "[particle e-] [E_kin_eV 1e8] [seed 1]\n";
    return 1;
  }
  const std::string gasName = argv[1];
  const double gap = std::atof(argv[2]) * 0.1;  // cm
  const double theta = std::atof(argv[3]) * M_PI / 180.;
  const int n = std::atoi(argv[4]);
  const std::string out = argv[5];
  const std::string particle = argc > 6 ? argv[6] : "e-";
  const double ekin = argc > 7 ? std::atof(argv[7]) : 1.e8;
  const unsigned seed = argc > 8 ? std::atoi(argv[8]) : 1;
  Random::SetEngine(RandomEngineRoot(seed));

  // compositions as in src/GasMixtures.cc / make_gas_table.C
  const std::map<std::string, std::vector<std::pair<std::string, double>>> mix = {
      {"ArIso", {{"ar", 95.}, {"ic4h10", 5.}}},
      {"ArIso9010", {{"ar", 90.}, {"ic4h10", 10.}}},
      {"NeIso", {{"ne", 95.}, {"ic4h10", 5.}}},
      {"NeIso8515", {{"ne", 85.}, {"ic4h10", 15.}}},
      {"ArCO2Iso9352", {{"ar", 93.}, {"co2", 5.}, {"ic4h10", 2.}}},
      {"ArCF4Iso", {{"ar", 88.}, {"cf4", 10.}, {"ic4h10", 2.}}},
  };
  auto it = mix.find(gasName);
  if (it == mix.end()) {
    std::cerr << "unknown gas " << gasName << "\n";
    return 1;
  }
  const auto& p = it->second;
  MediumMagboltz gas;
  if (p.size() == 2)
    gas.SetComposition(p[0].first, p[0].second, p[1].first, p[1].second);
  else
    gas.SetComposition(p[0].first, p[0].second, p[1].first, p[1].second,
                       p[2].first, p[2].second);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);

  ComponentConstant cmp;
  cmp.SetElectricField(0., 0., 600.);  // irrelevant for primary ionization
  cmp.SetArea(-5., -5., 0., 5., 5., gap);
  cmp.SetMedium(&gas);
  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.SetArea(-5., -5., 0., 5., 5., gap);

  TrackHeed track(&sensor);
  track.SetParticle(particle);
  track.SetKineticEnergy(ekin);
  track.EnableDeltaElectronTransport();

  std::ofstream f(out + ".txt");
  f << "# heed_mip gas=" << gasName << " gap_mm=" << gap * 10 << " theta_deg="
    << theta * 180. / M_PI << " particle=" << particle << " Ekin_eV=" << ekin
    << "\n# n_clusters n_electrons E_dep_eV\n";
  double sumN = 0.;
  for (int i = 0; i < n; ++i) {
    track.NewTrack(0., 0., 0., 0., std::sin(theta), 0., std::cos(theta));
    int ncl = 0;
    long ne = 0;
    double edep = 0.;
    for (const auto& cl : track.GetClusters()) {
      ++ncl;
      ne += cl.electrons.size();
      edep += cl.energy;
    }
    f << ncl << " " << ne << " " << edep << "\n";
    sumN += ne;
    if ((i + 1) % 20000 == 0) std::printf("%d/%d  <n_e> = %.2f\n", i + 1, n, sumN / (i + 1));
  }
  std::printf("HEED gas=%s gap_mm=%.3f theta=%.0f n=%d mean_ne=%.3f\n", gasName.c_str(),
              gap * 10, theta * 180. / M_PI, n, sumN / n);
  return 0;
}
