#!/usr/bin/env python3
"""Survey the W2 family slabs: health, wall time, cross-family consistency.

Reads only the tiny `meta` member of each slab, so it is cheap enough to run
against the EOS fuse mount; the arrays are never touched.

WHAT `prompt_capture` IS, since it invites exactly one wrong reading. It is
`v0.mean()` from `solve_prompt_cg(vpad)` in `w2_production.solve_family` — a
FULL-GRID CG solve that does not depend on `ifam`. So it is a per-BOX, per-DRIVE
constant that every family of that box recomputes independently, not that
family's share of anything. Summing it over families is meaningless (and gives
a number that looks plausible: 64 x the Y value lands near the expected total by
coincidence of the row count).

That redundancy is useful, though: because every family of a box solves the same
full-grid problem, the values MUST agree across families. A family that
disagrees has a corrupted assembly or a bad drive mask, and this is the cheapest
place to see it — before the ~1 h combine, and before four products are built on
top of it. That is what this script checks.

    python3 w2_slab_survey.py <slabdir>
"""
import glob
import io
import json
import os
import sys
import zipfile

import numpy as np

TOL = 1e-9          # same full-grid CG solve; agreement should be ~machine


def meta_of(path):
    z = zipfile.ZipFile(path)
    return json.loads(str(np.load(io.BytesIO(z.read("meta.npy")),
                                  allow_pickle=True)))


def main():
    d = sys.argv[1]
    rows = []
    for p in sorted(glob.glob(os.path.join(d, "w2slab_*.npz"))):
        try:
            rows.append((os.path.basename(p), meta_of(p)))
        except Exception as e:                      # truncated / mid-upload
            print(f"UNREADABLE {os.path.basename(p)}: {e}")
    if not rows:
        print("no slabs yet")
        return 0

    walls, ny, seen = [], set(), {}
    fams = {}
    bad = []
    for name, m in rows:
        box = m["box"].replace("box=", "")
        walls.append(m["wall_s"])
        ny.add(m["ny"])
        fams.setdefault(box, set()).add(m["ifam"])
        for drive, v in m["prompt_capture"].items():
            key = (box, drive)
            if key not in seen:
                seen[key] = (v, name)
            elif abs(seen[key][0] - v) > TOL:
                bad.append((name, box, drive, v, seen[key][0], seen[key][1]))
        if m["boundary"][:2] != "W2":
            bad.append((name, box, "BOUNDARY", m["boundary"], "W2...", "-"))

    print(f"slabs present : {len(rows)}/66")
    print(f"ny            : {sorted(ny)}  (production must be [512])")
    print(f"wall_s        : min {min(walls):.0f}  med {np.median(walls):.0f}  "
          f"max {max(walls):.0f}")
    for box in sorted(fams):
        want = next(m["n_fam"] for _, m in rows
                    if m["box"].replace("box=", "") == box)
        have = fams[box]
        missing = sorted(set(range(want)) - have)
        print(f"box {box}: {len(have)}/{want} families"
              + (f", missing {missing[:12]}{'...' if len(missing) > 12 else ''}"
                 if missing else " — COMPLETE"))
    print("\nper-drive full-grid prompt capture (must be identical across "
          "families of a box):")
    for (box, drive), (v, src) in sorted(seen.items()):
        print(f"  {box} {drive:>4s} = {v!r}   (first seen in {src})")

    if bad:
        print("\n*** CROSS-FAMILY DISAGREEMENT — do NOT combine:")
        for r in bad:
            print("   ", r)
        return 1
    print("\nall families agree to <1e-9 — assembly consistent so far")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
