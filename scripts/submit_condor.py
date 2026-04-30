#!/usr/bin/env python3
"""
submit_condor.py
Submits an HTCondor job array covering the full
  gas x particle x energy
matrix for the Micromegas sensitivity study.

Usage:
    python3 scripts/submit_condor.py [--dry-run] [--outdir /eos/...] [--nevents N]

Each job runs: mm_sim -g <gas> -p <particle> -e <E_MeV> -n <N> -o <outdir>/<tag>

Output ROOT files land in <outdir>/<tag>_t<thread>.root
After all jobs finish, run: python3 scripts/collect_results.py --indir <outdir>
"""

import argparse
import itertools
import os
import sys
import stat
import textwrap
from pathlib import Path

# ============================================================
# SCAN GRID -- edit here to add/remove points
# ============================================================

# GASES = ["ArCF4", "HeEth", "ArCO2", "ArCF4Iso", "NeIso", "NeCF4", "ArCF4CO2", "PureCF4"]
GASES = ["PureCF4"]

# Al shielding scan: thicknesses in mm (0 = no shielding baseline).
# Applied only to SHIELDING_GAS when --shielding flag is used.
SHIELDING_GAS        = "ArCF4"
AL_THICKNESSES_MM    = [1, 4]  # 0 mm (no shield) is always included

# Energies in MeV for each particle type
PARTICLE_ENERGIES = {
    # Gammas: fine log-spaced scan from 3 keV to 10 MeV (~55 points, ~5x denser)
    "gamma": [
        3.0e-3, 3.5e-3, 4.1e-3, 4.8e-3, 5.5e-3, 6.4e-3, 7.5e-3, 8.7e-3,  # 3–8.7 keV
        1.0e-2, 1.16e-2, 1.35e-2, 1.57e-2, 1.82e-2, 2.12e-2, 2.46e-2, 2.86e-2,  # 10–29 keV
        3.32e-2, 3.86e-2, 4.49e-2, 5.21e-2, 6.06e-2, 7.04e-2, 8.18e-2, 9.51e-2,  # 33–95 keV
        1.10e-1, 1.22e-1, 1.49e-1, 1.73e-1, 2.01e-1, 2.33e-1, 2.71e-1, 3.15e-1,  # 110–315 keV; 122=Co-57
        3.66e-1, 4.25e-1, 4.94e-1, 5.74e-1, 6.67e-1, 7.75e-1, 9.01e-1,  # 366–901 keV
        1.05, 1.22, 1.41, 1.64, 1.91, 2.22, 2.58, 3.00,  # 1.05–3 MeV
        3.49, 4.05, 4.71, 5.47, 6.36, 7.39, 8.59, 10.0, 20.0,  # 3.5–20 MeV
        # Extended range: 20 log-spaced points from ~28 MeV to 20 GeV
        28.3, 40.0, 56.4, 79.7, 113, 159, 224, 317, 448, 633,  # 28–633 MeV
        894, 1260, 1780, 2510, 3550, 5010, 7080, 10000, 14100, 20000,  # 894 MeV–20 GeV
    ],
    # Electrons: 1–10 MeV (relativistic MIP regime)
    "electron": [
        1.0,
        2.0,
        3.0,
        5.0,
        7.0,
        10.0,
    ],
    # Neutrons: low-energy points kept sparse; fine scan from 50 keV to 10 MeV
    "neutron": [
        1e-8,   # ~10 meV  -- cold neutron
        2.5e-8, # ~25 meV  -- thermal
        1e-6,   # ~1 eV
        1e-4,   # ~100 eV  -- epithermal
        1e-3,   # ~1 keV
        1e-2,   # ~10 keV
        # Fine scan 50 keV–10 MeV (~20 log-spaced points, ~5x denser than before)
        0.050, 0.066, 0.087, 0.115, 0.150, 0.200, 0.270, 0.350,
        0.470, 0.620, 0.810, 1.07, 1.42, 1.88, 2.48, 3.28,
        4.34, 5.74, 7.59, 10.0,
        14.0,   # D-T fusion neutrons
        20.0,
        # Extended range: 20 log-spaced points from ~28 MeV to 20 GeV
        28.3, 40.0, 56.4, 79.7, 113, 159, 224, 317, 448, 633,  # 28–633 MeV
        894, 1260, 1780, 2510, 3550, 5010, 7080, 10000, 14100, 20000,  # 894 MeV–20 GeV
    ],
}

# Per-particle event count multiplier applied on top of --nevents
NEVENTS_SCALE = {
    "gamma":    10,
    "electron": 1,
    "neutron":  200,
}

# ============================================================
# Condor configuration
# ============================================================
NEVENTS_DEFAULT = 50000
JOB_FLAVOUR    = "workday"   # ~8 h; use "longlunch" (2h) for tests
NTHREADS       = 1           # single-threaded per job; Condor provides parallelism


def parse_args():
    p = argparse.ArgumentParser(description="Submit Micromegas scan to HTCondor")
    p.add_argument("--dry-run",  action="store_true",
                   help="Print jobs without submitting")
    p.add_argument("--outdir",   default="/eos/user/d/dneff/mx17_geant_sim_results",
                   help="Output directory for ROOT files (default: EOS results dir)")
    p.add_argument("--jobdir",   default="/afs/cern.ch/user/d/dneff/condor/mx17_geant_sim",
                   help="Directory for condor submit file, wrapper, and logs (must be on AFS)")
    p.add_argument("--nevents",  type=int, default=NEVENTS_DEFAULT,
                   help=f"Events per job (default: {NEVENTS_DEFAULT})")
    p.add_argument("--exe",      default=None,
                   help="Path to mm_sim executable (default: auto-detect)")
    p.add_argument("--flavour",  default=JOB_FLAVOUR,
                   help=f"Condor job flavour (default: {JOB_FLAVOUR})")
    p.add_argument("--gases",    nargs="+", default=GASES,
                   help="Gas mixtures to include")
    p.add_argument("--particles",nargs="+", default=list(PARTICLE_ENERGIES.keys()),
                   help="Particle types to include")
    p.add_argument("--shielding", action="store_true",
                   help=f"Also submit Al shielding jobs for {SHIELDING_GAS} "
                        f"at {AL_THICKNESSES_MM} mm (plus 0 mm baseline if not already present)")
    return p.parse_args()


def find_exe():
    """Find the mm_sim executable relative to this script."""
    script_dir = Path(__file__).parent
    candidates = [
        script_dir.parent / "build" / "mm_sim",
        script_dir.parent / "install" / "bin" / "mm_sim",
        Path(os.environ.get("MM_SIM_EXE", "")),
    ]
    for c in candidates:
        if c.is_file():
            return str(c.resolve())
    return None


def make_job_tag(gas, particle, energy_mev, al_mm=0):
    """Create a short unique tag for a job."""
    e_str = f"{energy_mev:.6g}MeV".replace(".", "p")
    tag = f"{gas}_{particle}_{e_str}"
    if al_mm > 0:
        al_str = f"{al_mm:g}mm".replace(".", "p")
        tag += f"_Al{al_str}"
    return tag


def write_wrapper_script(job_dir: Path, exe: str, setup_script: str) -> Path:
    """
    Write a bash wrapper that re-sources Geant4 on the worker node
    (Condor workers may not inherit the login shell environment).
    The wrapper receives job arguments via environment variables set in the
    Condor submit file.
    """
    wrapper = job_dir / "run_job.sh"
    content = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -e
        # Re-source Geant4 environment on worker node
        source "{setup_script}"

        # Arguments passed via Condor environment / positional
        GAS="$1"
        PARTICLE="$2"
        ENERGY="$3"
        NEVENTS="$4"
        OUTFILE="$5"
        SEED="$6"
        AL_MM="${{7:-0}}"

        echo "=== mm_sim job ==="
        echo "  Node     : $(hostname)"
        echo "  Gas      : $GAS"
        echo "  Particle : $PARTICLE"
        echo "  Energy   : $ENERGY MeV"
        echo "  Events   : $NEVENTS"
        echo "  Output   : $OUTFILE"
        echo "  Seed     : $SEED"
        echo "  Al (mm)  : $AL_MM"
        echo "=================="

        "{exe}" \\
            -g "$GAS"      \\
            -p "$PARTICLE" \\
            -e "$ENERGY"   \\
            -n "$NEVENTS"  \\
            -o "$OUTFILE"  \\
            -s "$SEED"     \\
            -a "$AL_MM"

        echo "Job done: $(date)"
    """)
    wrapper.write_text(content)
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def write_condor_submit(job_dir: Path, wrapper: Path, jobs: list,
                        outdir: Path, flavour: str) -> Path:
    """
    Write a single Condor submit file with one queue entry per job.
    Uses queue ... from a list of arguments (Condor 8.8+ feature).
    """
    submit_file = job_dir / "mm_scan.sub"

    log_dir  = job_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    lines = [
        f"executable            = {wrapper}",
        f"output                = {log_dir}/$(tag).out",
        f"error                 = {log_dir}/$(tag).err",
        f"log                   = {log_dir}/condor.log",
        f'+JobFlavour           = "{flavour}"',
        "request_cpus          = 1",
        "request_memory        = 2048",   # MB -- G4 with HP data uses ~1 GB
        "request_disk          = 2048",   # MB
        "should_transfer_files = YES",
        "when_to_transfer_output = ON_EXIT",
        "",
        "# G4 HP neutron data lives on CVMFS -- needs network access on worker",
        'requirements          = (OpSysAndVer =?= "AlmaLinux9")',
        "",
        "# Arguments: gas particle energy_MeV nevents outfile_base seed al_mm",
        "arguments             = $(gas) $(particle) $(energy) $(nevents) $(outfile) $(seed) $(al_mm)",
        "",
        "queue gas,particle,energy,nevents,outfile,seed,al_mm,tag from (",
    ]

    import random
    rng = random.Random(42)

    for (gas, particle, energy_mev, nevents, al_mm) in jobs:
        tag     = make_job_tag(gas, particle, energy_mev, al_mm)
        outfile = str(outdir / tag)
        seed    = rng.randint(1, 2**31 - 1)
        lines.append(f"  {gas}, {particle}, {energy_mev}, {nevents}, {outfile}, {seed}, {al_mm}, {tag}")

    lines.append(")")

    submit_file.write_text("\n".join(lines) + "\n")
    return submit_file


def main():
    args = parse_args()

    # Resolve paths
    exe = args.exe or find_exe()
    if not exe or not Path(exe).is_file():
        print("ERROR: mm_sim executable not found.")
        print("  Build first: bash scripts/build.sh")
        print("  Or specify with --exe /path/to/mm_sim")
        sys.exit(1)
    exe = str(Path(exe).resolve())

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Job submission directory (logs, sub file, wrapper) -- must be on AFS, not EOS
    job_dir = Path(args.jobdir)
    job_dir.mkdir(parents=True, exist_ok=True)

    # Setup script (must be accessible from worker nodes -- use EOS or AFS path)
    setup_script = str(Path(__file__).parent.resolve() / "setup_lxplus.sh")

    # Build job list: (gas, particle, energy, nevents, al_mm)
    jobs = []
    for gas in args.gases:
        for particle in args.particles:
            if particle not in PARTICLE_ENERGIES:
                print(f"WARNING: Unknown particle '{particle}', skipping")
                continue
            for energy in PARTICLE_ENERGIES[particle]:
                n = args.nevents * NEVENTS_SCALE.get(particle, 1)
                jobs.append((gas, particle, energy, n, 0))

    # Al shielding jobs: ArCF4 with 1 mm and 4 mm Al (+ 0 mm baseline if not already present)
    if args.shielding:
        shielding_gases = set(j[0] for j in jobs)
        al_thicknesses = AL_THICKNESSES_MM[:]
        if SHIELDING_GAS not in shielding_gases:
            al_thicknesses = [0] + al_thicknesses  # include no-shield baseline too
        for al_mm in al_thicknesses:
            for particle in args.particles:
                if particle not in PARTICLE_ENERGIES:
                    continue
                for energy in PARTICLE_ENERGIES[particle]:
                    n = args.nevents * NEVENTS_SCALE.get(particle, 1)
                    jobs.append((SHIELDING_GAS, particle, energy, n, al_mm))

    print(f"Total jobs to submit: {len(jobs)}")
    print(f"  Gases    : {args.gases}")
    print(f"  Particles: {args.particles}")
    print(f"  Events   : {args.nevents} per job (x{NEVENTS_SCALE} scale by particle)")
    print(f"  Output   : {outdir}")
    print(f"  Exe      : {exe}")

    if args.dry_run:
        print("\n--- DRY RUN: first 10 jobs ---")
        for (gas, particle, energy, nevents, al_mm) in jobs[:10]:
            print(f"  {make_job_tag(gas, particle, energy, al_mm)}  ({nevents} events)")
        if len(jobs) > 10:
            print(f"  ... and {len(jobs)-10} more")
        return

    # Write files
    wrapper = write_wrapper_script(job_dir, exe, setup_script)
    sub_file = write_condor_submit(job_dir, wrapper, jobs,
                                   outdir, args.flavour)

    print(f"\nSubmit file: {sub_file}")
    print("Submitting to HTCondor...")

    ret = os.system(f"condor_submit {sub_file}")
    if ret != 0:
        print("ERROR: condor_submit failed")
        sys.exit(1)

    print("\nDone! Monitor with:")
    print("  condor_q")
    print(f"\nAfter completion, collect results:")
    print(f"  python3 scripts/collect_results.py --indir {outdir}")


if __name__ == "__main__":
    main()
