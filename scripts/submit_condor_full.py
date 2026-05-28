#!/usr/bin/env python3
"""
submit_condor_full.py
Submits full-experiment stack simulation jobs to HTCondor.

Primary scan
  Electrons : 200 log-spaced points, 0.1–12 MeV
  Triton    : He-3(n,p)T thermal capture product (KE_T ≈ 0.191 MeV) + range
  Proton    : He-3(n,p)T thermal capture product (KE_p ≈ 0.573 MeV) + range

Particles originate at the centre of the He-3 gas volume (-m full mode).
The gun fires along +z through the full downstream material stack.

Usage
    python3 scripts/submit_condor_full.py [options]

After jobs finish
    python3 scripts/collect_results.py --indir <outdir>
"""

import argparse
import math
import os
import random
import stat
import sys
import textwrap
from pathlib import Path

# ============================================================
# SCAN GRID
# ============================================================

def _logspace(lo, hi, n):
    """n log-spaced values from lo to hi (no numpy required)."""
    la, lb = math.log10(lo), math.log10(hi)
    return [round(10 ** (la + i * (lb - la) / (n - 1)), 8) for i in range(n)]


# 200 log-spaced lepton energies, 0.1–18 MeV, plus exact round numbers.
# The log grid already covers the range densely; the round-number extras
# guarantee there are exact integer/half-integer points for clean plot labels.
_log_grid   = _logspace(0.1, 18.0, 200)
_round_grid = [0.5] + list(range(1, 19))          # 0.5, 1, 2, …, 18
# Merge: drop a round number if a log point is already within 0.1 % of it
# (only catches the 18.0 endpoint which is exact in both grids).
_tol = 0.001
_extras = [r for r in _round_grid
           if not any(abs(r - x) / r < _tol for x in _log_grid)]
LEPTON_ENERGIES = sorted(_log_grid + _extras)

PARTICLE_ENERGIES = {
    # Main scan — electron and positron transmission/calorimetry study
    "electron": LEPTON_ENERGIES,
    "positron": LEPTON_ENERGIES,  # same range; positrons annihilate at end of range

    # He-3(n,p)T capture products for thermal neutrons (E_n ≈ 25 meV):
    #   KE_triton ≈ Q × M_p/(M_p+M_T) = 0.764 × 1/4 = 0.191 MeV
    #   KE_proton ≈ Q × M_T/(M_p+M_T) = 0.764 × 3/4 = 0.573 MeV
    # Extended range to cover higher-energy neutron interactions
    "triton": [0.191, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0],
    "proton":  [0.573, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0],
}

# Default gas list
GASES_DEFAULT = ["ArIso"]

# Events per job (no multiplier for full-mode particles)
NEVENTS_DEFAULT = 50000
NEVENTS_SCALE = {
    "electron": 1,
    "positron": 1,
    "triton":   1,
    "proton":   1,
}

# Condor
JOB_FLAVOUR = "workday"   # ~8 h; use "longlunch" (~2 h) for testing


# ============================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Submit full-experiment stack scan to HTCondor",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--dry-run", action="store_true",
                   help="Print job list without submitting")
    p.add_argument("--outdir",
                   default="/eos/user/d/dneff/mx17_geant_sim_results/full",
                   help="Output directory for ROOT files (EOS)")
    p.add_argument("--jobdir",
                   default="/afs/cern.ch/user/d/dneff/condor/mx17_geant_full",
                   help="Directory for condor files and logs (AFS)")
    p.add_argument("--nevents", type=int, default=NEVENTS_DEFAULT,
                   help="Events per job")
    p.add_argument("--exe", default=None,
                   help="Path to mm_sim executable (auto-detected if omitted)")
    p.add_argument("--flavour", default=JOB_FLAVOUR,
                   help="HTCondor job flavour")
    p.add_argument("--gases", nargs="+", default=GASES_DEFAULT,
                   help="Gas mixtures to scan")
    p.add_argument("--particles", nargs="+",
                   default=list(PARTICLE_ENERGIES.keys()),
                   help="Particle types to include")
    p.add_argument("--cfrp", type=float, default=1.5,
                   help="CFRP wall thickness for LS cells [mm]")
    return p.parse_args()


def find_exe():
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


def make_tag(gas, particle, energy_mev):
    e_str = f"{energy_mev:.8g}MeV".replace(".", "p")
    return f"full_{gas}_{particle}_{e_str}"


def write_wrapper(job_dir: Path, exe: str, setup_script: str,
                  cfrp_mm: float) -> Path:
    """
    Bash wrapper re-sources Geant4 on the worker node and runs mm_sim in
    full-experiment mode. Arguments are positional: gas particle energy
    nevents outfile seed.
    """
    wrapper = job_dir / "run_full_job.sh"
    content = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -e
        source "{setup_script}"

        GAS="$1"
        PARTICLE="$2"
        ENERGY="$3"
        NEVENTS="$4"
        OUTFILE="$5"
        SEED="$6"

        echo "=== mm_sim full-experiment job ==="
        echo "  Node     : $(hostname)"
        echo "  Gas      : $GAS"
        echo "  Particle : $PARTICLE"
        echo "  Energy   : $ENERGY MeV"
        echo "  Events   : $NEVENTS"
        echo "  Output   : $OUTFILE"
        echo "  Seed     : $SEED"
        echo "=================================="

        "{exe}" \\
            -m full      \\
            -g "$GAS"    \\
            -p "$PARTICLE" \\
            -e "$ENERGY" \\
            -n "$NEVENTS" \\
            -o "$OUTFILE" \\
            -s "$SEED"   \\
            -c {cfrp_mm}

        echo "Job done: $(date)"
    """)
    wrapper.write_text(content)
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def write_submit(job_dir: Path, wrapper: Path, jobs: list,
                 outdir: Path, flavour: str) -> Path:
    """
    Single Condor submit file; one queue entry per job using
    'queue ... from' (Condor 8.8+ syntax).
    """
    submit_file = job_dir / "mm_full_scan.sub"
    log_dir = job_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    lines = [
        f"executable            = {wrapper}",
        f"output                = {log_dir}/$(tag).out",
        f"error                 = {log_dir}/$(tag).err",
        f"log                   = {log_dir}/condor.log",
        f'+JobFlavour           = "{flavour}"',
        "request_cpus          = 1",
        "request_memory        = 2048",
        "request_disk          = 2048",
        "should_transfer_files = YES",
        "when_to_transfer_output = ON_EXIT",
        "",
        'requirements          = (OpSysAndVer =?= "AlmaLinux9")',
        "",
        "# Arguments: gas particle energy_MeV nevents outfile_base seed",
        "arguments             = $(gas) $(particle) $(energy) $(nevents) $(outfile) $(seed)",
        "",
        "queue gas,particle,energy,nevents,outfile,seed,tag from (",
    ]

    rng = random.Random(42)
    for (gas, particle, energy, nevents) in jobs:
        tag     = make_tag(gas, particle, energy)
        outfile = str(outdir / tag)
        seed    = rng.randint(1, 2**31 - 1)
        lines.append(f"  {gas}, {particle}, {energy}, {nevents}, {outfile}, {seed}, {tag}")

    lines.append(")")
    submit_file.write_text("\n".join(lines) + "\n")
    return submit_file


def main():
    args = parse_args()

    exe = args.exe or find_exe()
    if not exe or not Path(exe).is_file():
        print("ERROR: mm_sim executable not found.")
        print("  Build first:  bash scripts/build.sh")
        print("  Or specify:   --exe /path/to/mm_sim")
        sys.exit(1)
    exe = str(Path(exe).resolve())

    outdir  = Path(args.outdir)
    job_dir = Path(args.jobdir)
    outdir.mkdir(parents=True, exist_ok=True)
    job_dir.mkdir(parents=True, exist_ok=True)

    setup_script = str(Path(__file__).parent.resolve() / "setup_lxplus.sh")

    # Build job list
    jobs = []
    for gas in args.gases:
        for particle in args.particles:
            if particle not in PARTICLE_ENERGIES:
                print(f"WARNING: unknown particle '{particle}', skipping")
                continue
            for energy in PARTICLE_ENERGIES[particle]:
                n = args.nevents * NEVENTS_SCALE.get(particle, 1)
                jobs.append((gas, particle, energy, n))

    n_electron = sum(1 for j in jobs if j[1] == "electron")
    n_other    = len(jobs) - n_electron

    print(f"Jobs to submit : {len(jobs)}")
    print(f"  Electrons    : {n_electron} energy points")
    print(f"  Other        : {n_other} points (triton + proton)")
    print(f"  Gases        : {args.gases}")
    print(f"  Events/job   : {args.nevents}")
    print(f"  CFRP (LS)    : {args.cfrp} mm")
    print(f"  Output       : {outdir}")
    print(f"  Exe          : {exe}")

    if args.dry_run:
        print("\n--- DRY RUN: first 10 jobs ---")
        for (gas, particle, energy, nevents) in jobs[:10]:
            print(f"  {make_tag(gas, particle, energy)}  ({nevents} events)")
        if len(jobs) > 10:
            print(f"  ... and {len(jobs)-10} more")
        print(f"\nLepton energies: {LEPTON_ENERGIES[0]:.4f} – {LEPTON_ENERGIES[-1]:.4f} MeV "
              f"({len(LEPTON_ENERGIES)} points: 200 log-spaced + {len(_extras)} round-number extras)")
        return

    wrapper  = write_wrapper(job_dir, exe, setup_script, args.cfrp)
    sub_file = write_submit(job_dir, wrapper, jobs, outdir, args.flavour)

    print(f"\nSubmit file : {sub_file}")
    print("Submitting to HTCondor...")

    ret = os.system(f"condor_submit {sub_file}")
    if ret != 0:
        print("ERROR: condor_submit failed")
        sys.exit(1)

    print("\nDone! Monitor with:")
    print("  condor_q")
    print(f"\nAfter completion, collect results:")
    print(f"  python3 scripts/collect_results.py --indir {args.outdir}")


if __name__ == "__main__":
    main()
