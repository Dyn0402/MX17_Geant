#!/usr/bin/env python3
"""
check_grid_active.py — refuse to run on a Garfield build that ignores the field map's inactive regions.

An unpatched Garfield 927e5c21 parses the region flag of a 3D xyz map but never stores it
(setup_garfield.sh, "LOCAL PATCH"), so electrons drift through the mesh wires and the transparency
comes out wrong without any error. The patched build sees ~1.6 % of meshfield_production.txt's
bounding box as inactive (wires + below the anode); the unpatched one sees 0 (2026-10-16, lxplus).

    python3 check_grid_active.py meshfield_production.txt     # exit 1 if nothing is inactive
"""
import ctypes
import itertools
import os
import sys

import ROOT
import Garfield  # noqa: F401  (loads the Garfield dictionary)

G = ROOT.Garfield
grid = G.ComponentGrid()
if not grid.LoadElectricField(sys.argv[1], "xyz", True, True):
    sys.exit("[check_grid_active] failed to load the map")
grid.SetMedium(G.MediumMagboltz())
box = [ctypes.c_double(0) for _ in range(6)]
grid.GetBoundingBox(*box)
box = [v.value for v in box]
N = 40
inactive = sum(
    not grid.GetMedium(*[box[a] + (box[a + 3] - box[a]) * (c + .5) / N for a, c in enumerate(ijk)])
    for ijk in itertools.product(range(N), repeat=3))
frac = inactive / N ** 3
print(f"[check_grid_active] inactive fraction {frac:.4f}", flush=True)
# os._exit skips the ROOT/Garfield teardown, which sometimes segfaults after the answer is known.
os._exit(0 if frac > 0.005 else 1)
