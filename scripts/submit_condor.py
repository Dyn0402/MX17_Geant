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

GASES = ["ArCF4", "HeEth", "ArCO2", "ArCF4Iso", "NeIso"]

# Energies in MeV for each particle type
PARTICLE_ENERGIES = {
    # Gammas: very wide range from near-IR equivalent (eV) to 10 MeV
    # Sparse sampling -- we want cross-section shape, not fine resolution
    "gamma": [
        1e-5,   # 10 eV     -- below K-edges, mostly Rayleigh
        1e-4,   # 100 eV
        1e-3,   # 1 keV     -- below Ar K-edge (3.2 keV), photoelectric dominates
        3e-3,   # 3 keV
        5e-3,   # 5 keV     -- Fe-55 region
        1e-2,   # 10 keV
        3e-2,   # 30 keV
        6e-2,   # 60 keV    -- Co-57 region
        1.22e-1,# 122 keV   -- Co-57 gamma
        5e-1,   # 500 keV   -- Compton starts dominating
        1.0,    # 1 MeV
        2.0,
        5.0,
        10.0,
    ],
    # Electrons: 1–10 MeV (as requested -- relativistic MIP regime)
    "electron": [
        1.0,
        2.0,
        3.0,
        5.0,
        7.0,
        10.0,
    ],
    # Neutrons: thermal to fast (1 eV to 20 MeV)
    "neutron": [
        1e-8,   # 10 neV  -- cold neutron
        2.5e-8, # 25 neV  -- thermal (room temperature)
        1e-6,   # 1 ueV
        1e-4,   # 0.1 meV
        1e-3,   # 1 meV
        1e-2,   # 10 meV  -- epithermal
        0.1,    # 100 keV -- fast
        1.0,    # 1 MeV
        2.0,
        5.0,
        14.0,   # 14 MeV  -- D-T fusion neutrons
        20.0,
    ],
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
    p.add_argument("--outdir",   default=None,
                   help="Output directory (default: $HOME/mm_results)")
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


def make_job_tag(gas, particle, energy_mev):
    """Create a short unique tag for a job."""
    e_str = f"{energy_mev:.6g}MeV".replace(".", "p")
    return f"{gas}_{particle}_{e_str}"


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

        echo "=== mm_sim job ==="
        echo "  Node     : $(hostname)"
        echo "  Gas      : $GAS"
        echo "  Particle : $PARTICLE"
        echo "  Energy   : $ENERGY MeV"
        echo "  Events   : $NEVENTS"
        echo "  Output   : $OUTFILE"
        echo "  Seed     : $SEED"
        echo "=================="

        "{exe}" \\
            -g "$GAS"      \\
            -p "$PARTICLE" \\
            -e "$ENERGY"   \\
            -n "$NEVENTS"  \\
            -o "$OUTFILE"  \\
            -s "$SEED"

        echo "Job done: $(date)"
    """)
    wrapper.write_text(content)
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def write_condor_submit(job_dir: Path, wrapper: Path, jobs: list,
                        outdir: Path, nevents: int, flavour: str) -> Path:
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
        "# Arguments: gas particle energy_MeV nevents outfile_base seed",
        "arguments             = $(gas) $(particle) $(energy) $(nevents) $(outfile) $(seed)",
        "",
        "queue gas,particle,energy,nevents,outfile,seed,tag from (",
    ]

    import random
    rng = random.Random(42)

    for (gas, particle, energy_mev) in jobs:
        tag     = make_job_tag(gas, particle, energy_mev)
        outfile = str(outdir / tag)
        seed    = rng.randint(1, 2**31 - 1)
        lines.append(f"  {gas}, {particle}, {energy_mev}, {nevents}, {outfile}, {seed}, {tag}")

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

    outdir = Path(args.outdir) if args.outdir else Path.home() / "mm_results"
    outdir.mkdir(parents=True, exist_ok=True)

    # Job submission directory (logs, sub file, wrapper)
    job_dir = outdir / "condor"
    job_dir.mkdir(exist_ok=True)

    # Setup script (must be accessible from worker nodes -- use EOS or AFS path)
    setup_script = str(Path(__file__).parent.resolve() / "setup_lxplus.sh")

    # Build job list
    jobs = []
    for gas in args.gases:
        for particle in args.particles:
            if particle not in PARTICLE_ENERGIES:
                print(f"WARNING: Unknown particle '{particle}', skipping")
                continue
            for energy in PARTICLE_ENERGIES[particle]:
                jobs.append((gas, particle, energy))

    print(f"Total jobs to submit: {len(jobs)}")
    print(f"  Gases    : {args.gases}")
    print(f"  Particles: {args.particles}")
    print(f"  Events   : {args.nevents} per job")
    print(f"  Output   : {outdir}")
    print(f"  Exe      : {exe}")

    if args.dry_run:
        print("\n--- DRY RUN: first 10 jobs ---")
        for j in jobs[:10]:
            print(f"  {make_job_tag(*j)}")
        if len(jobs) > 10:
            print(f"  ... and {len(jobs)-10} more")
        return

    # Write files
    wrapper = write_wrapper_script(job_dir, exe, setup_script)
    sub_file = write_condor_submit(job_dir, wrapper, jobs,
                                   outdir, args.nevents, args.flavour)

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
