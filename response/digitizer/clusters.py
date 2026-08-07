#!/usr/bin/env python3
"""
clusters.py — read a Geant4 ClusterTree into the digitizer's frame (plan §6).

Stage A writes cluster positions in the WORLD frame and, alongside them, a
one-entry `Meta` tree carrying the world -> active-area transform. Using the
world coordinates directly is wrong and would fail quietly: for the 2026-08-07
production the drift gas sits at world z in [-12.5, +17.5] mm, so a digitizer
that assumed "z is the drift distance from the ESL" would hand negative depths
to the diffusion model for half the gap and put the other half at the wrong
depth. Always go through the Meta transform:

    x_act = sign_x * (x_world - origin_x)
    y_act = sign_y * (y_world - origin_y)
    z_act = sign_z * (z_world - origin_z)

which for that file maps the gas to z_act in [0, 30.2] mm, i.e. 0 at the ESL
surface and 30 mm at the cathode — what §7's drift model wants.

Only DriftGas and AmpGas clusters are ionisation the response chain should see.
AmpGas ones are already inside the amplification gap: they skip the drift and
are amplified where they sit, and NEEDED_INPUTS records that the ESL inter-strip
grooves are scored as AmpGas too. Everything else (PCB, mesh, windows) is
material budget, not signal.

`nPrimary` is the electron count to use, NOT edep/W. SteppingAction already
applied W and carried the sub-W remainder probabilistically, so re-deriving
from edep would double-count the conversion and lose the remainder.
"""

from __future__ import annotations

import numpy as np

GAS_VOLUMES = ("DriftGas", "AmpGas")


class ClusterFile:
    """A Geant4 ClusterTree, transformed into the active-area frame."""

    def __init__(self, path, volumes=GAS_VOLUMES):
        import uproot
        f = uproot.open(path)
        meta = f["Meta"].arrays(library="np")
        self.origin = np.array([meta["origin_x_mm"][0], meta["origin_y_mm"][0],
                                meta["origin_z_mm"][0]])
        self.sign = np.array([meta["sign_x"][0], meta["sign_y"][0],
                              meta["sign_z"][0]])
        self.active_width_mm = float(meta["active_width_mm"][0])
        if not bool(meta["valid"][0]):
            raise ValueError(f"{path}: Meta.valid is False — the geometry did "
                             "not expose an active area, so cluster positions "
                             "cannot be placed. Re-run Stage A.")

        t = f["ClusterTree"]
        d = t.arrays(["eventID", "x", "y", "z", "time", "edep", "nPrimary",
                      "volume"], library="np")
        vol = np.array([v.decode() if isinstance(v, bytes) else str(v)
                        for v in d["volume"]])
        keep = np.isin(vol, list(volumes))

        self.event = d["eventID"][keep]
        self.x = self.sign[0] * (d["x"][keep] - self.origin[0])
        self.y = self.sign[1] * (d["y"][keep] - self.origin[1])
        self.z = self.sign[2] * (d["z"][keep] - self.origin[2])
        self.t = d["time"][keep]
        self.edep = d["edep"][keep]
        self.n_e = d["nPrimary"][keep].astype(int)
        self.volume = vol[keep]
        self.path = path

    def events(self):
        """Yield (eventID, slice) per event, in file order."""
        order = np.argsort(self.event, kind="stable")
        ev = self.event[order]
        bounds = np.flatnonzero(np.diff(ev)) + 1
        for grp in np.split(order, bounds):
            yield int(self.event[grp[0]]), grp

    def describe(self):
        return {
            "file": self.path.split("/")[-1],
            "n_clusters": int(len(self.z)),
            "n_events": int(len(np.unique(self.event))),
            "n_primary_total": int(self.n_e.sum()),
            "n_primary_per_event": float(self.n_e.sum()
                                         / max(len(np.unique(self.event)), 1)),
            "z_mm_range": [float(self.z.min()), float(self.z.max())],
            "t_ns_range": [float(self.t.min()), float(self.t.max())],
            "active_width_mm": self.active_width_mm,
        }


if __name__ == "__main__":
    import argparse
    import pprint
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    a = ap.parse_args()
    pprint.pprint(ClusterFile(a.path).describe())
