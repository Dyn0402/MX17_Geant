// mm_sim.cc
// Micromegas + full-experiment GEANT4 simulation
// Author: Dylan Neff / n_TOF X17 group

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SimConfig.hh"

#include <iostream>
#include <string>
#include <ctime>

void PrintUsage() {
    std::cerr << "Usage: mm_sim [options] [macro_file]\n";
    std::cerr << "Options:\n";
    std::cerr << "  -g <gas>         Gas: ArCF4, ArIso, HeEth, ArCO2, ArCF4Iso, NeIso, NeCF4,\n";
    std::cerr << "                       ArCF4CO2, PureCF4, PureAr, PureHe, PureNe,\n";
    std::cerr << "                       PureEthane, PureIso, PureCO2  (default: ArCF4)\n";
    std::cerr << "  -p <particle>    gamma, neutron, electron, proton, muon, muon+,\n";
    std::cerr << "                   pion, alpha, triton  (default: gamma)\n";
    std::cerr << "  -e <energy>      Particle energy [MeV]  (default: 1.0)\n";
    std::cerr << "  -n <nevents>     Number of events  (default: 10000)\n";
    std::cerr << "  -o <output>      Output file base name  (default: mm_output)\n";
    std::cerr << "  -s <seed>        Random seed  (default: time-based)\n";
    std::cerr << "  -t <nthreads>    MT threads  (default: 1)\n";
    std::cerr << "  -m <mode>        vacuum | full  (default: vacuum)\n";
    std::cerr << "  -a <mm>          Al shielding [mm], vacuum mode only  (default: 0)\n";
    std::cerr << "  -c <mm>          CFRP wall thickness [mm] for LS cells, full mode only  (default: 1.5)\n";
    std::cerr << "  -v               Verbose output\n";
    std::cerr << "  -h               Print this help\n";
}

int main(int argc, char** argv) {
    SimConfig config;
    config.gas            = "ArCF4";
    config.particle       = "gamma";
    config.energy         = 1.0 * MeV;
    config.nEvents        = 10000;
    config.outFile        = "mm_output";
    config.seed           = static_cast<long>(std::time(nullptr));
    config.nThreads       = 1;
    config.verbose        = false;
    config.alThickness_mm = 0.0;
    config.mode           = SimMode::kVacuum;
    config.cfrpThickness_mm = 1.5;

    std::string macroFile = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "-h")              { PrintUsage(); return 0; }
        else if (arg == "-g" && i+1<argc)  config.gas              = argv[++i];
        else if (arg == "-p" && i+1<argc)  config.particle         = argv[++i];
        else if (arg == "-e" && i+1<argc)  config.energy           = std::stod(argv[++i]) * MeV;
        else if (arg == "-n" && i+1<argc)  config.nEvents          = std::stoi(argv[++i]);
        else if (arg == "-o" && i+1<argc)  config.outFile          = argv[++i];
        else if (arg == "-s" && i+1<argc)  config.seed             = std::stol(argv[++i]);
        else if (arg == "-t" && i+1<argc)  config.nThreads         = std::stoi(argv[++i]);
        else if (arg == "-a" && i+1<argc)  config.alThickness_mm   = std::stod(argv[++i]);
        else if (arg == "-c" && i+1<argc)  config.cfrpThickness_mm = std::stod(argv[++i]);
        else if (arg == "-v")              config.verbose           = true;
        else if (arg == "-m" && i+1<argc) {
            std::string mode = argv[++i];
            if      (mode == "vacuum") config.mode = SimMode::kVacuum;
            else if (mode == "full")   config.mode = SimMode::kFullExperiment;
            else { std::cerr << "Unknown mode: " << mode << " (use 'vacuum' or 'full')\n"; return 1; }
        }
        else if (arg[0] != '-') macroFile = arg;
        else { std::cerr << "Unknown option: " << arg << "\n"; PrintUsage(); return 1; }
    }

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(config.seed);

    std::cout << "=== Micromegas Simulation ===" << "\n";
    std::cout << "  Geant4 version : " << G4Version << "\n";
    std::cout << "  Mode           : " << (config.mode == SimMode::kFullExperiment ? "full-experiment" : "vacuum") << "\n";
    std::cout << "  Gas            : " << config.gas << "\n";
    std::cout << "  Particle       : " << config.particle << "\n";
    std::cout << "  Energy         : " << config.energy/MeV << " MeV" << "\n";
    std::cout << "  Events         : " << config.nEvents << "\n";
    std::cout << "  Output         : " << config.outFile << "\n";
    std::cout << "  Seed           : " << config.seed << "\n";
    std::cout << "  Threads        : " << config.nThreads << "\n";
    if (config.mode == SimMode::kVacuum)
        std::cout << "  Al shielding   : " << config.alThickness_mm << " mm\n";
    else
        std::cout << "  CFRP thickness : " << config.cfrpThickness_mm << " mm\n";
    std::cout << "=============================" << "\n";

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(config.nThreads);
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // DetectorConstruction must outlive ActionInitialization so that
    // PrimaryGeneratorAction can read GetHe3GasCenterZ() after Construct().
    auto* detCon = new DetectorConstruction(config);
    runManager->SetUserInitialization(detCon);
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(config, detCon));

    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (macroFile.empty()) {
        UI->ApplyCommand("/run/beamOn " + std::to_string(config.nEvents));
    } else {
        UI->ApplyCommand("/control/execute " + macroFile);
    }

    delete runManager;
    std::cout << "Simulation complete. Output: " << config.outFile << "\n";
    return 0;
}
