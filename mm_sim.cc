// mm_sim.cc
// Main entry point for the Micromegas GEANT4 simulation
// Simulates energy deposition and primary ionization in a Micromegas detector
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
    std::cerr << "Usage: mm_sim [options] [macro_file]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -g <gas>         Gas mixture: ArIso, HeEth, NeCO2, ArCO2, PurAr (default: ArIso)" << std::endl;
    std::cerr << "  -p <particle>    Particle: gamma, neutron, electron, proton, muon (default: gamma)" << std::endl;
    std::cerr << "  -e <energy>      Particle energy in MeV (default: 1.0)" << std::endl;
    std::cerr << "  -n <nevents>     Number of events (default: 10000)" << std::endl;
    std::cerr << "  -o <output>      Output file base name (default: mm_output)" << std::endl;
    std::cerr << "  -s <seed>        Random seed (default: time-based)" << std::endl;
    std::cerr << "  -t <nthreads>    Number of threads for MT mode (default: 1)" << std::endl;
    std::cerr << "  -a <mm>          Al shielding thickness in mm, 0 = none (default: 0)" << std::endl;
    std::cerr << "  -v               Verbose output" << std::endl;
    std::cerr << "  -h               Print this help" << std::endl;
}

int main(int argc, char** argv) {
    // ---- Parse command-line arguments ----
    SimConfig config;
    config.gas           = "ArIso";
    config.particle      = "gamma";
    config.energy        = 1.0 * MeV;
    config.nEvents       = 10000;
    config.outFile       = "mm_output";
    config.seed          = static_cast<long>(std::time(nullptr));
    config.nThreads      = 1;
    config.verbose       = false;
    config.alThickness_mm = 0.0;

    std::string macroFile = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") { PrintUsage(); return 0; }
        else if (arg == "-g" && i+1 < argc) config.gas      = argv[++i];
        else if (arg == "-p" && i+1 < argc) config.particle = argv[++i];
        else if (arg == "-e" && i+1 < argc) config.energy   = std::stod(argv[++i]) * MeV;
        else if (arg == "-n" && i+1 < argc) config.nEvents  = std::stoi(argv[++i]);
        else if (arg == "-o" && i+1 < argc) config.outFile  = argv[++i];
        else if (arg == "-s" && i+1 < argc) config.seed     = std::stol(argv[++i]);
        else if (arg == "-t" && i+1 < argc) config.nThreads      = std::stoi(argv[++i]);
        else if (arg == "-a" && i+1 < argc) config.alThickness_mm = std::stod(argv[++i]);
        else if (arg == "-v")               config.verbose        = true;
        else if (arg[0] != '-')             macroFile       = arg;
        else { std::cerr << "Unknown option: " << arg << std::endl; PrintUsage(); return 1; }
    }

    // ---- Random engine ----
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(config.seed);

    std::cout << "=== Micromegas Simulation ===" << std::endl;
    std::cout << "  Geant4 version : " << G4Version << std::endl;
    std::cout << "  Gas            : " << config.gas << std::endl;
    std::cout << "  Particle       : " << config.particle << std::endl;
    std::cout << "  Energy         : " << config.energy / MeV << " MeV" << std::endl;
    std::cout << "  Events         : " << config.nEvents << std::endl;
    std::cout << "  Output         : " << config.outFile << std::endl;
    std::cout << "  Seed           : " << config.seed << std::endl;
    std::cout << "  Threads        : " << config.nThreads << std::endl;
    std::cout << "  Al shielding   : " << config.alThickness_mm << " mm" << std::endl;
    std::cout << "=============================" << std::endl;

    // ---- Run Manager ----
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(config.nThreads);
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // ---- Mandatory initializations ----
    runManager->SetUserInitialization(new DetectorConstruction(config));
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(config));

    runManager->Initialize();

    // ---- UI manager ----
    G4UImanager* UI = G4UImanager::GetUIpointer();

    if (macroFile.empty()) {
        // Batch mode: run directly
        UI->ApplyCommand("/run/beamOn " + std::to_string(config.nEvents));
    } else {
        // Run macro file
        std::string command = "/control/execute " + macroFile;
        UI->ApplyCommand(command);
    }

    delete runManager;
    std::cout << "Simulation complete. Output: " << config.outFile << std::endl;
    return 0;
}
