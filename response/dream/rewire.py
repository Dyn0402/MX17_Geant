#!/usr/bin/env python3
"""
rewire.py — put simulated channels on det3's real connector wiring (2026-10-18).

Stage B writes pad column k into decoded channel k (run.py, "the map is the identity"). det3's run_config says
every connector on both views is "inverted": inside each 64-channel connector the channel order runs opposite
to the strips. wft's strip map (nTof_x17 `wft.io.strip_position_map`) applies that, so for FEU 7/8 it gives

    position(channel j) = 0.78 mm * (64 * (j // 64) + 63 - j % 64)

and a simulated track is read mirrored inside its connector and torn apart across a connector boundary: a track
over pad columns 315..325 shows up as two pieces ~99 mm apart. run.py's comment called the identity a pure
relabelling that only matters channel by channel; it also matters to anything that reads strip positions from
the det3 map (wft on sim, the extent/span comparisons) whenever a track crosses a boundary. The sim beam sits
on columns ~262..298 at 0-10 deg (3 % of tracks cross column 256) but moves to ~320 at 20 deg, where more than
a quarter of the tracks cross.

`connector_inverted_perm()` is its own inverse. Writing pad column k into channel perm[k] makes the det3 map
read position 0.78 * k, the true layout.

    python -m response.dream.rewire IN_DIR OUT_DIR      # every *.root decoded file in IN_DIR
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

N_CHAN, CONNECTOR = 512, 64


def connector_inverted_perm(n_chan: int = N_CHAN, width: int = CONNECTOR) -> np.ndarray:
    k = np.arange(n_chan)
    return width * (k // width) + (width - 1) - (k % width)


def rewire_file(src: str, dst: str) -> int:
    """Rewrite one decoded file with its channels permuted. Returns the number of events."""
    import uproot
    from .daq import write_decoded, verify
    t = uproot.open(src)["nt"]
    d = t.arrays(["eventId", "ftst", "amplitude"], library="np")
    n_ev = len(d["eventId"])
    first = np.asarray(d["amplitude"][0])
    n_sample = first.size // N_CHAN
    adc = np.stack([np.asarray(a).reshape(n_sample, N_CHAN) for a in d["amplitude"]])
    perm = connector_inverted_perm()
    out = np.empty_like(adc)
    out[:, :, perm] = adc            # pad column k -> channel perm[k]
    write_decoded(dst, out, event_ids=d["eventId"], ftst=d["ftst"])
    bad = verify(dst, out)
    if bad:
        raise RuntimeError(f"{dst}: {bad} events do not read back")
    return n_ev


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    src_dir, dst_dir = argv[0], argv[1]
    os.makedirs(dst_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(src_dir, "*.root")))
    if not files:
        raise SystemExit(f"no .root files in {src_dir}")
    for f in files:
        n = rewire_file(f, os.path.join(dst_dir, os.path.basename(f)))
        print(f"  {os.path.basename(f)}: {n} events rewired")


if __name__ == "__main__":
    main()
