// slab.cc -- Geant4 primary ionization in a bare gas slab (no upstream material),
// to compare with Garfield/Heed heed_mip.cc on identical geometry.
//   slab <gas> <gap_mm> <n> <out> <physics: opt4|pai|paiphot> [seed]
// Vacuum world, gas slab 10x10 cm x gap at 0<z<gap, 100 MeV e- entering at z=0 along +z.
// Output per event: edep_eV  n_W (edep/W per step, floor+Bernoulli as mm_sim)
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserEventAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"
#include "G4Positron.hh"
static std::string gPhys;
#include "G4DecayPhysics.hh"
#include "G4EmParameters.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ParticleGun.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "Randomize.hh"

static std::string gGas; static double gGap, gW, gE = 100; static double gEdep; static long gN;
static std::ofstream gOut;

struct Comp { const char* el; std::vector<std::pair<const char*,int>> atoms; double M, W; };
static const std::map<std::string, Comp> kComp = {
  {"Ar", {"Ar", {{"Ar",1}}, 39.948, 26.4}},
  {"Ne", {"Ne", {{"Ne",1}}, 20.1797, 36.4}},
  {"iso",{"iso",{{"C",4},{"H",10}}, 58.1222, 26.0}},
  {"CF4",{"CF4",{{"C",1},{"F",4}}, 88.0043, 34.0}},
};
static const std::map<std::string, std::vector<std::pair<std::string,double>>> kMix = {
  {"ArIso", {{"Ar",.95},{"iso",.05}}},
  {"NeIso", {{"Ne",.95},{"iso",.05}}},
  {"ArCF4Iso", {{"Ar",.88},{"CF4",.10},{"iso",.02}}},
};
static double rho(double M){ return 101325.*M*1e-3/(8.314462618*293.15)*1e-3; } // g/cm3

class Det : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    // "<gas>_mx17": the gas exactly as MX17_Geant DetectorConstruction builds it:
    // 0 C STP densities declared at 293 K, volume fractions passed as mass fractions,
    // and SteppingAction's W table.
    const bool mx = gGas.size() > 5 && gGas.substr(gGas.size()-5) == "_mx17";
    const std::string base = mx ? gGas.substr(0, gGas.size()-5) : gGas;
    static const std::map<std::string,double> stp = {{"Ar",1.782e-3},{"Ne",0.8999e-3},{"iso",2.67e-3},{"CF4",3.72e-3}};
    static const std::map<std::string,double> wmx = {{"ArIso",26.0},{"NeIso",27.0},{"ArCF4Iso",33.0}};
    auto& mix = kMix.at(base);
    double rmix = 0, wsum = 0;
    for (auto& [c,f] : mix) rmix += f*(mx ? stp.at(c) : rho(kComp.at(c).M));
    auto* gas = new G4Material("gas", rmix*g/cm3, (int)mix.size(), kStateGas, 293.15*kelvin, 1*atmosphere);
    for (auto& [c,f] : mix) {
      auto& cc = kComp.at(c);
      auto* m = new G4Material(c, rho(cc.M)*g/cm3, (int)cc.atoms.size(), kStateGas, 293.15*kelvin, 1*atmosphere);
      for (auto& [e,n] : cc.atoms) m->AddElement(nist->FindOrBuildElement(e), n);
      gas->AddMaterial(m, mx ? f : f*rho(cc.M)/rmix);   // mass fraction (mx17: volume fraction, as there)
      wsum += f*cc.W;                          // linear volume-weighted W (as P2 GasMixtures)
    }
    gW = mx ? wmx.at(base) : wsum;
    auto* vac = nist->FindOrBuildMaterial("G4_Galactic");
    auto* lw = new G4LogicalVolume(new G4Box("w", 1*m, 1*m, 1*m), vac, "World");
    auto* pw = new G4PVPlacement(nullptr, {}, lw, "World", nullptr, false, 0);
    auto* ls = new G4LogicalVolume(new G4Box("s", 5*cm, 5*cm, 0.5*gGap*mm), gas, "DriftGas");
    new G4PVPlacement(nullptr, {0,0,0.5*gGap*mm}, ls, "DriftGas", lw, false, 0);
    auto* reg = new G4Region("GasRegion"); reg->AddRootLogicalVolume(ls);
    return pw;
  }
};
class Phys : public G4VModularPhysicsList {
 public:
  Phys(){ if (gPhys.rfind("opt0",0)==0) RegisterPhysics(new G4EmStandardPhysics(0)); else RegisterPhysics(new G4EmStandardPhysics_option4(0)); RegisterPhysics(new G4DecayPhysics(0));
    auto* p = G4EmParameters::Instance(); p->SetFluo(true); p->SetAuger(true); p->SetDeexcitationIgnoreCut(true); }
  void ConstructProcess() override {
    G4VModularPhysicsList::ConstructProcess();
    if (gPhys.find("pai") == std::string::npos) return;
    auto* conf = G4LossTableManager::Instance()->EmConfigurator();
    for (G4ParticleDefinition* p : std::vector<G4ParticleDefinition*>{G4Electron::Definition(), G4Positron::Definition()}) {
      G4VEmModel* m; G4VEmFluctuationModel* f;
      if (gPhys.find("paiphot") != std::string::npos) { auto* x = new G4PAIPhotModel(p, "PAIPhot"); m = x; f = x; }
      else { auto* x = new G4PAIModel(p, "PAI"); m = x; f = x; }
      conf->SetExtraEmModel(p->GetParticleName(), "eIoni", m, "GasRegion", 0.11*keV, 100*TeV, f);
    }
  }
  void SetCuts() override {   // as mm_sim (MX17 + P2): e+- 10 um, gamma 0.1 mm
    SetCutValue(10*um,"e-"); SetCutValue(10*um,"e+"); SetCutValue(0.1*mm,"gamma"); SetCutValue(1*mm,"proton"); }
};
class Gun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun g{1};
 public:
  Gun(){ g.SetParticleDefinition(G4Electron::Definition()); g.SetParticleEnergy(gE*MeV);
         g.SetParticlePosition({0,0,-1*mm}); g.SetParticleMomentumDirection({0,0,1}); }
  void GeneratePrimaries(G4Event* e) override { g.GeneratePrimaryVertex(e); }
};
class Step : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* s) override {
    double ed = s->GetTotalEnergyDeposit(); if (ed <= 0) return;
    auto* pv = s->GetPreStepPoint()->GetPhysicalVolume();
    if (!pv || pv->GetName() != "DriftGas") return;
    gEdep += ed/eV;
    double nb = ed/(gW*eV); long n = (long)std::floor(nb); if (G4UniformRand() < nb-n) ++n; gN += n;
  }
};
class Evt : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override { gEdep = 0; gN = 0; }
  void EndOfEventAction(const G4Event*) override { gOut << gEdep << " " << gN << "\n"; }
};

int main(int argc, char** argv) {
  gGas = argv[1]; gGap = std::atof(argv[2]); int n = std::atoi(argv[3]); std::string out = argv[4];
  std::string phys = argv[5]; gPhys = phys; long seed = argc > 6 ? std::atol(argv[6]) : 1; if (argc > 7) gE = std::atof(argv[7]);
  G4Random::setTheSeed(seed);
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new Det); rm->SetUserInitialization(new Phys);
        rm->SetUserAction(new Gun); rm->SetUserAction(new Step); rm->SetUserAction(new Evt);
  gOut.open(out + ".txt");
  rm->Initialize();
  gOut << "# gas=" << gGas << " gap_mm=" << gGap << " phys=" << phys << " W=" << gW << "\n";
  rm->BeamOn(n);
  delete rm; return 0;
}
