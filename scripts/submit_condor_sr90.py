#!/usr/bin/env python3
"""
submit_condor_sr90.py
Submits Sr-90/Y-90 calibration simulation jobs to HTCondor.

Geometry: full downstream stack (MM → PCB → scint wall → LS calorimeter)
          but WITHOUT the He-3 pressurised target or its Al/CFRP capsule.
          Source is placed in air at the position of the former He-3 centre,
          with 226.5 mm of air between source and the MM entrance.
Mode flag: -m sr90

Energy range: 0.1 – 3.0 MeV covering the Sr-90 (0.546 MeV) and Y-90
              (2.28 MeV) beta endpoints with fine resolution.
              Exact endpoints and round numbers are added on top of the
              log-spaced grid.

Usage
    python3 scripts/submit_condor_sr90.py [options]
    python3 scripts/submit_condor_sr90.py --dry-run   # inspect without submitting

After jobs finish
    python3 scripts/analyze_full_experiment.py \\
        --indir  <outdir> --prefix sr90 \\
        --gas ArIso --particles electron positron
"""

import argparse
import math
import os
import random
import stat
import sys
import textwrap
from pathlib import Path

# ── Energy grid ──────────────────────────────────────────────────────────────

def _logspace(lo, hi, n):
    la, lb = math.log10(lo), math.log10(hi)
    return [round(10 ** (la + i * (lb - la) / (n - 1)), 8) for i in range(n)]


# 150 log-spaced points from 0.1 to 3.0 MeV
_log_grid = _logspace(0.1, 3.0, 150)

# Exact round numbers and physics endpoints
_round_grid = (
    [0.2, 0.3, 0.4, 0.5]          # Sr-90 region
    + [0.546]                      # Sr-90 beta endpoint
    + [0.6, 0.7, 0.8, 0.9]
    + list(range(1, 4))            # 1, 2, 3 MeV
    + [1.5, 2.28, 2.5]            # Y-90 endpoint and neighbours
)

_tol = 0.001   # drop round number if a log point is already within 0.1 %
_extras = [r for r in _round_grid
           if not any(abs(r - x) / r < _tol for x in _log_grid)]
LEPTON_ENERGIES = sorted(_log_grid + _extras)

PARTICLE_ENERGIES = {
    "electron": LEPTON_ENERGIES,
    "positron": LEPTON_ENERGIES,
}

GASES_DEFAULT  = ["ArIso"]
NEVENTS_DEFAULT = 50000
NEVENTS_SCALE  = {"electron": 1, "positron": 1}
JOB_FLAVOUR    = "longlunch"   # ~2 h; jobs are fast (low energy, short range)


# ── Condor helpers ────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Submit Sr-90 calibration simulation to HTCondor",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--mode", default="sr90", choices=["sr90", "sr90nomm"],
                   help="sr90: full stack with MM  |  sr90nomm: no MM or PCB")
    p.add_argument("--outdir",
                   default="/eos/user/d/dneff/mx17_geant_sim_results/sr90",
                   help="Output directory for ROOT files (EOS); sr90nomm jobs "
                        "go into a 'nomm' subdirectory automatically")
    p.add_argument("--jobdir",
                   default="/afs/cern.ch/user/d/dneff/condor/mx17_geant_sr90",
                   help="Condor submit files and logs (AFS)")
    p.add_argument("--nevents", type=int, default=NEVENTS_DEFAULT)
    p.add_argument("--exe",     default=None,
                   help="Path to mm_sim (auto-detected if omitted)")
    p.add_argument("--flavour", default=JOB_FLAVOUR)
    p.add_argument("--gases",   nargs="+", default=GASES_DEFAULT)
    p.add_argument("--particles", nargs="+",
                   default=list(PARTICLE_ENERGIES.keys()))
    p.add_argument("--cfrp",   type=float, default=1.5,
                   help="CFRP wall thickness for LS cells [mm]")
    return p.parse_args()


def find_exe():
    script_dir = Path(__file__).parent
    for c in [script_dir.parent / "build" / "mm_sim",
              script_dir.parent / "install" / "bin" / "mm_sim",
              Path(os.environ.get("MM_SIM_EXE", ""))]:
        if c.is_file():
            return str(c.resolve())
    return None


def make_tag(gas, particle, energy_mev, mode="sr90"):
    e_str = f"{energy_mev:.8g}MeV".replace(".", "p")
    return f"{mode}_{gas}_{particle}_{e_str}"


def eos_xrd(path: str) -> str:
    """Convert /eos/... path to root://eosuser.cern.ch//eos/... for xrdcp."""
    if path.startswith("/eos/"):
        return "root://eosuser.cern.ch/" + path
    return path


def write_wrapper(job_dir: Path, exe: str, setup_script: str,
                  cfrp_mm: float, mode: str = "sr90") -> Path:
    """
    Wrapper that xrdcp's the executable from EOS, runs locally, then
    xrdcp's output back.  Batch nodes don't POSIX-mount /eos.
    Arguments: gas particle energy nevents outfile seed
      outfile is the desired EOS destination (used to derive local name + xrdcp target).
    """
    wrapper = job_dir / f"run_{mode}_job.sh"
    exe_xrd = eos_xrd(exe)

    content = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -e
        source "{setup_script}"

        GAS="$1"; PARTICLE="$2"; ENERGY="$3"
        NEVENTS="$4"; OUTFILE="$5"; SEED="$6"
        LOCAL=$(basename "$OUTFILE")

        echo "=== mm_sim {mode} job ==="
        echo "  Node     : $(hostname)"
        echo "  Gas      : $GAS  Particle : $PARTICLE"
        echo "  Energy   : $ENERGY MeV   Events : $NEVENTS"
        echo "=========================="

        # Copy executable from EOS to local scratch
        xrdcp -s "{exe_xrd}" ./mm_sim
        chmod +x ./mm_sim

        ./mm_sim \\
            -m {mode}      \\
            -g "$GAS"      \\
            -p "$PARTICLE" \\
            -e "$ENERGY"   \\
            -n "$NEVENTS"  \\
            -o "./$LOCAL"  \\
            -s "$SEED"     \\
            -c {cfrp_mm}

        # Copy output back to EOS (dirname gives /eos/..., prepend root://server)
        for f in "./$LOCAL"*; do
            xrdcp -s "$f" "root://eosuser.cern.ch/$(dirname $OUTFILE)/$(basename $f)"
        done

        echo "Job done: $(date)"
    """)
    wrapper.write_text(content)
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def write_submit(job_dir: Path, wrapper: Path, jobs: list,
                 outdir: Path, flavour: str, mode: str = "sr90") -> Path:
    submit_file = job_dir / f"mm_{mode}_scan.sub"
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
        "arguments             = $(gas) $(particle) $(energy) $(nevents) $(outfile) $(seed)",
        "",
        "queue gas,particle,energy,nevents,outfile,seed,tag from (",
    ]
    rng = random.Random(42)
    for (gas, particle, energy, nevents) in jobs:
        tag     = make_tag(gas, particle, energy, mode)
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
        print("ERROR: mm_sim not found. Build first: bash scripts/build.sh")
        sys.exit(1)

    mode    = args.mode
    # sr90nomm output goes in a separate subdirectory automatically
    base_outdir = Path(args.outdir)
    outdir  = base_outdir if mode == "sr90" else base_outdir.parent / (base_outdir.name + "_nomm")
    job_dir = Path(args.jobdir) if mode == "sr90" else Path(args.jobdir + "_nomm")
    outdir.mkdir(parents=True, exist_ok=True)
    job_dir.mkdir(parents=True, exist_ok=True)

    setup_script = str(Path(__file__).parent.resolve() / "setup_lxplus.sh")

    jobs = []
    for gas in args.gases:
        for particle in args.particles:
            if particle not in PARTICLE_ENERGIES:
                print(f"WARNING: unknown particle '{particle}', skipping")
                continue
            for energy in PARTICLE_ENERGIES[particle]:
                n = args.nevents * NEVENTS_SCALE.get(particle, 1)
                jobs.append((gas, particle, energy, n))

    n_e = sum(1 for j in jobs if j[1] == "electron")
    n_p = sum(1 for j in jobs if j[1] == "positron")
    print(f"Mode         : {mode}")
    print(f"Jobs to submit : {len(jobs)}")
    print(f"  Electrons    : {n_e}  Positrons : {n_p}")
    print(f"  Energy range : {LEPTON_ENERGIES[0]:.4f} – {LEPTON_ENERGIES[-1]:.4f} MeV"
          f"  ({len(LEPTON_ENERGIES)} points)")
    print(f"  Gases        : {args.gases}")
    print(f"  Events/job   : {args.nevents}")
    print(f"  Output       : {outdir}")

    if args.dry_run:
        print(f"\n--- DRY RUN ({mode}): first 8 jobs ---")
        for g, p, e, n in jobs[:8]:
            print(f"  {make_tag(g, p, e, mode)}  ({n} events)")
        if len(jobs) > 8:
            print(f"  ... and {len(jobs)-8} more")
        return

    wrapper  = write_wrapper(job_dir, exe, setup_script, args.cfrp, mode)
    sub_file = write_submit(job_dir, wrapper, jobs, outdir, args.flavour, mode)

    print(f"\nSubmit file : {sub_file}")
    ret = os.system(f"condor_submit {sub_file}")
    if ret != 0:
        print("ERROR: condor_submit failed")
        sys.exit(1)

    print("\nDone! Monitor with: condor_q")
    print(f"\nAfter completion, analyse with:")
    print(f"  python3 scripts/analyze_full_experiment.py \\")
    print(f"      --indir {outdir} --prefix {mode} \\")
    print(f"      --gas {args.gases[0]} --particles electron positron")
    print(f"\nThen run the Sr-90 calibration analysis:")
    print(f"  python3 sr90_calibration/analyze_sr90_calibration.py \\")
    print(f"      --spectrum sr90_calibration/Sr90_Y90_Beta_Spectrum.csv \\")
    print(f"      --summary  <outfile>_{mode}_electron.csv")


if __name__ == "__main__":
    main()
