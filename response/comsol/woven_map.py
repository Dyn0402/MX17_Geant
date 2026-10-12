#!/usr/bin/env python3
"""
woven_map.py -- the Garfield field map of the true woven mesh, from WovenMesh.java.

COMSOL's sandbox forbids file I/O from model code, so `MAP=1 WovenMesh` computes the
grid, the flags and the shell projections itself and prints them as CSV lines, which
run_lxplus.sh collects into $OUT. This script turns that CSV into the map:

  python3 woven_map.py assemble <comsol.csv> <meshfield_woven.txt>

It re-derives every flag with classify() below and exits nonzero on any disagreement.

The grid and the flag semantics are those of meshcell/solve_fieldmap.py, so the map
drops into the same Garfield code: 1 um step, x, y in [-33.5, 33.5], z in [-168, 218],
the written flag eroded by one step and the in-wire shell nodes carrying the field of
their projection 0.4 um outside the nearest wire surface (see _sample_plane there).
The difference is the wire geometry: four sinusoidal half-wires on the cell faces
(WovenMesh.java header). Load with mirror periodicity:

    grid.LoadElectricField("meshfield_woven.txt", "xyz", true, true);
    grid.EnableMirrorPeriodicityX(); grid.EnableMirrorPeriodicityY();
"""
import sys
import numpy as np

P, RW, OVL = 67.0, 9.5, 1.0
A = RW - OVL / 2
STEP, STANDOFF = 1.0, 0.4
XS = np.arange(-P / 2, P / 2 + 1e-9, STEP)            # 68
ZS = np.arange(-168.0, 218.0 + 1e-9, STEP)            # 387

# (axis along, transverse coordinate, sign): z = sign * A sin(pi s / p)
WIRES = [("x", -P / 2, -1.0), ("x", P / 2, 1.0), ("y", -P / 2, 1.0), ("y", P / 2, -1.0)]


def _nearest_on_axis(s0, t, z, c, sg):
    """
    Nearest point of the axis (s, c, sg*A sin(pi s/p)) to the point (s0, t, z), by
    Newton on d/ds |.|^2. The axis curvature radius >= p^2/(pi^2 A) = 50 um exceeds the
    wire radius, so near the surface the minimum is unique and starts from s = s0.
    """
    k = np.pi / P
    s = s0.copy()
    for _ in range(8):
        zs = sg * A * np.sin(k * s)
        dz = sg * A * k * np.cos(k * s)
        d2z = -sg * A * k * k * np.sin(k * s)
        f = (s - s0) + (zs - z) * dz
        fp = 1.0 + dz * dz + (zs - z) * d2z
        s = s - f / fp
    zs = sg * A * np.sin(k * s)
    d = np.sqrt((s - s0) ** 2 + (t - c) ** 2 + (z - zs) ** 2)
    return s, zs, d


def wire_distances(x, y, z):
    """Distance from each point to each wire axis, and the nearest axis points."""
    out = []
    for ax, c, sg in WIRES:
        if ax == "x":
            s, zs, d = _nearest_on_axis(x, y, z, c, sg)
            out.append((d, np.stack([s, np.full_like(s, c), zs])))
        else:
            s, zs, d = _nearest_on_axis(y, x, z, c, sg)
            out.append((d, np.stack([np.full_like(s, c), s, zs])))
    return out


def classify(x, y, z):
    """flag: 0 deep in a wire, 1 gas, 2 in-wire shell (within one step of the surface)."""
    wd = wire_distances(x, y, z)
    D = np.stack([d for d, _ in wd])
    dmin = D.min(0)
    flag = np.where(dmin >= RW, 1, np.where(dmin >= RW - STEP, 2, 0))
    # projection of shell points radially out of the nearest wire, STANDOFF outside
    near = D.argmin(0)
    C = np.stack([c for _, c in wd])                   # [wire, xyz, n]
    cpt = C[near, :, np.arange(len(x))].T              # [xyz, n]
    pts = np.stack([x, y, z])
    v = pts - cpt
    sc = (RW + STANDOFF) / np.maximum(dmin, 1e-6)
    proj = cpt + v * sc
    return flag, proj


def grid():
    Z, X, Y = np.meshgrid(ZS, XS, XS, indexing="ij")   # z outer, x, y inner (as solve_fieldmap)
    return X.ravel(), Y.ravel(), Z.ravel()


def cmd_assemble(csv, out):
    """CSV from `MAP=1 WovenMesh` -> ComponentGrid xyz map at `out`."""
    d = np.loadtxt(csv, comments="#")
    ix, iy, iz = d[:, 0].astype(int), d[:, 1].astype(int), d[:, 2].astype(int)
    flag = d[:, 3].astype(int)
    n = len(XS)
    assert len(d) == n * n * len(ZS), (len(d), n * n * len(ZS))
    x, y, z = XS[ix], XS[iy], ZS[iz]
    # cross-check COMSOL-side flags against the numpy classifier
    f_np, _ = classify(x, y, z)
    nmis = int(np.sum(f_np != flag))
    print(f"flags: gas {np.sum(flag == 1)}, shell {np.sum(flag == 2)}, deep {np.sum(flag == 0)}; "
          f"disagree with numpy classify: {nmis}")
    vals = d[:, 4:8]
    bad = (flag > 0) & ~np.all(np.isfinite(vals), 1)
    if bad.any():
        print(f"WARNING {bad.sum()} evaluation points outside the COMSOL domain -> flagged dead")
    alive = (flag > 0) & ~bad
    vals[~alive] = 0.0
    cm = 1e-4
    with open(out, "w") as fh:
        fh.write(f"# XMIN = {XS[0]*cm:.5f}, XMAX = {XS[-1]*cm:.5f}, NX = {n}\n")
        fh.write(f"# YMIN = {XS[0]*cm:.5f}, YMAX = {XS[-1]*cm:.5f}, NY = {n}\n")
        fh.write(f"# ZMIN = {ZS[0]*cm:.4f}, ZMAX = {ZS[-1]*cm:.4f}, NZ = {len(ZS)}\n")
        rows = np.column_stack([x * cm, y * cm, z * cm, vals[:, 1], vals[:, 2], vals[:, 3],
                                vals[:, 0], alive.astype(int)])
        np.savetxt(fh, rows, fmt="%.6e %.6e %.6e %.6e %.6e %.6e %.6e %d")
    print(f"wrote {out}: {alive.sum()} active of {len(d)}")
    return nmis


if __name__ == "__main__":
    # python3 woven_map.py assemble <comsol.csv> <meshfield_woven.txt>
    assert sys.argv[1] == "assemble"
    sys.exit(1 if cmd_assemble(sys.argv[2], sys.argv[3]) else 0)
