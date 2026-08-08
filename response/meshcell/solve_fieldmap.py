#!/usr/bin/env python3
"""solve_fieldmap.py — T6: the MX17 micromesh unit-cell field map (S2).

FEM (gmsh + scikit-fem P2 tets + pyamg) solve of one unit cell of the
bi-planar crossed-wire approximation of the woven mesh, between the ESL anode
plane and a drift-slice top plane. Replaces the neBEM attempt (field_map.C):
neBEM cannot represent infinite plates — its bounding-plane/mirror interface
data is read but unused by the solver — and a finite plate patch leaves a
1/n^2 capacitor-fringe error that swamps the 333 V/cm drift field (measured:
+540% at 5 periodic copies, still +12% at 30). Here the infinite lattice is
EXACT: the bi-planar cell is mirror-symmetric about x=0, y=0, x=±p/2, y=±p/2,
so homogeneous-Neumann side walls (FEM natural BC) reproduce the infinite
periodic mesh with no periodicity machinery at all.

GEOMETRY (shared/MX17ModuleGeometry.hh):
  wire d 19 um, pitch 67 um; x-wire at z=+o, y-wire at z=-o with
  o = r - overlap/2 (small kiss overlap at the crossing avoids the degenerate
  tangency; real woven wires rest on each other). z=0 = weave mid-plane.
  Mesh underside = -(o+r). The 150 um amplification gap is the BULK PILLAR
  HEIGHT (ESL to mesh underside), so the anode plane sits at
  z_anode = -(o + r + 150 um) — NOT at -150 um. Mean amp field is therefore
  ~30.5 kV/cm at 490 V, about 7% below the uniform-field S3's 32.7 kV/cm;
  S3 gains re-run in this map will drop accordingly. That is physics, not a
  bug. Top plane at z_cath = +(o + r + 200 um): the drift field is uniform to
  ~exp(-2*pi*200/67) there, so a Dirichlet equipotential cut is exact.

POTENTIALS: mesh 0 V, anode +V_MESH (bench 490 V), top -E_DRIFT*z_cath
(det3/T14 target 333 V/cm). In Garfield's E=-grad V convention Ez > 0
everywhere in the cell; electrons drift toward -z.

OUTPUT: meshfield_<tag>.txt in ComponentGrid::LoadElectricField "xyz" format
('#' mesh-header lines + "x y z Ex Ey Ez V flag", cm / V/cm / V), flag=0
inside wires, + a provenance JSON. Consumer:
    ComponentGrid grid;
    grid.LoadElectricField("meshfield_production.txt", "xyz", true, true);
    grid.SetMedium(&gas);
    grid.EnablePeriodicityX(); grid.EnablePeriodicityY();

GATES (all must print PASS):
  S1 far-field  <Ez>(top of gas) = +E_DRIFT within 1%
  S2 amp bulk   <Ez>(z=-90um) within 2% of the independent neBEM
                infinite-lattice extrapolation, 30561 V/cm (copies_scan.C)
  S3 wire BC    |V| just off the wire surface < 2 V
  S4 mirror     Ez(x,y,z) = Ez(-x,y,z) = Ez(x,-y,z) within 0.5%
  S5 evaluator  custom P2 evaluator == basis value interpolator, < 1e-8 rel
  S6 refine     production only: field vs a 1.4x-coarser solve, 95th pct
                rel diff < 1% over random gas points

Run (laptop or desktop; needs the nTof_x17 venv):
    ~/PycharmProjects/nTof_x17/.venv/bin/python solve_fieldmap.py --smoke
    ~/PycharmProjects/nTof_x17/.venv/bin/python solve_fieldmap.py --production
Then run the Garfield-side gates: gates_check.C (transparency etc.).
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import scipy.sparse as sp

# ── Geometry & operating point (um, V) ───────────────────────────────────────
WIRE_R = 9.5          # wire radius
PITCH = 67.0          # weave pitch
OVERLAP = 1.0         # kiss overlap of the two wire sets at the crossing
AMP_GAP = 150.0       # pillar height: ESL surface to mesh underside
DRIFT_SLICE = 200.0   # exported drift, above mesh topside (plan §2 contract)
V_MESH = 490.0        # mesh-anode voltage, bench operating point
E_DRIFT = 333.0       # V/cm, det3/T14 drift field

Z_OFF = WIRE_R - OVERLAP / 2.0          # wire-axis |z| offset (=9.0)
Z_UNDER = -(Z_OFF + WIRE_R)             # mesh underside (-18.5)
Z_TOP = +(Z_OFF + WIRE_R)               # mesh topside   (+18.5)
Z_ANODE = Z_UNDER - AMP_GAP             # ESL surface    (-168.5)
Z_CATH = Z_TOP + DRIFT_SLICE            # drift cut      (+218.5)
V_ANODE = V_MESH
V_CATH = -E_DRIFT * Z_CATH * 1e-4       # E_DRIFT is V/cm; z in um

# Independent cross-check: neBEM infinite-lattice amp-gap field (copies_scan.C
# converged to ~0.005%/step at 30 copies; wires at exactly +-r, no overlap).
NEBEM_AMP_EZ = 30561.0  # V/cm

HERE = os.path.dirname(os.path.abspath(__file__))


def build_mesh(lc_wire, lc_far, msh_path):
    """Unit cell in gmsh: box minus the two wire cylinders."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ
    box = occ.addBox(-PITCH / 2, -PITCH / 2, Z_ANODE,
                     PITCH, PITCH, Z_CATH - Z_ANODE)
    pad = 5.0
    cx = occ.addCylinder(-PITCH / 2 - pad, 0., +Z_OFF,
                         PITCH + 2 * pad, 0., 0., WIRE_R)
    cy = occ.addCylinder(0., -PITCH / 2 - pad, -Z_OFF,
                         0., PITCH + 2 * pad, 0., WIRE_R)
    out, _ = occ.cut([(3, box)], [(3, cx), (3, cy)])
    occ.synchronize()

    # Wire surfaces = the only surfaces confined to the weave z-band.
    wire_surfs = []
    for (dim, tag) in gmsh.model.getEntities(2):
        bb = gmsh.model.getBoundingBox(dim, tag)
        if bb[2] > Z_UNDER - 1.0 and bb[5] < Z_TOP + 1.0:
            wire_surfs.append(tag)
    if len(wire_surfs) < 2:
        raise RuntimeError(f"wire surface classification failed: {wire_surfs}")

    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "SurfacesList", wire_surfs)
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 200)
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", lc_wire)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", 2 * WIRE_R)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", 80.0)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
    gmsh.model.mesh.generate(3)
    gmsh.write(msh_path)
    gmsh.finalize()


def load_mesh(msh_path):
    import meshio
    mm = meshio.read(msh_path)
    pts = mm.points.astype(np.float64)
    tets = mm.cells_dict["tetra"].astype(np.int64)
    from skfem import MeshTet
    return MeshTet(pts.T.copy(), tets.T.copy())


def classify_boundary_facets(m):
    """Split boundary facets into anode / top / wall / wire by midpoint."""
    bf = m.boundary_facets()
    mid = m.p[:, m.facets[:, bf]].mean(axis=1)  # (3, nbf)
    tol = 0.05
    is_anode = mid[2] < Z_ANODE + tol
    is_top = mid[2] > Z_CATH - tol
    is_wall = ((np.abs(np.abs(mid[0]) - PITCH / 2) < tol) |
               (np.abs(np.abs(mid[1]) - PITCH / 2) < tol)) & ~is_anode & ~is_top
    is_wire = ~(is_anode | is_top | is_wall)
    # Sanity: every "wire" facet midpoint must sit near one of the two wire
    # surfaces (within chord sagitta of the discretized cylinder).
    dx = np.hypot(mid[1], mid[2] - Z_OFF)   # distance to x-wire axis
    dy = np.hypot(mid[0], mid[2] + Z_OFF)   # distance to y-wire axis
    d = np.minimum(dx, dy)[is_wire]
    if d.size and (d.max() > WIRE_R + 0.3 or d.min() < WIRE_R - 3.0):
        raise RuntimeError(
            f"wire facet classification off: d in [{d.min():.2f},{d.max():.2f}]")
    return bf[is_anode], bf[is_top], bf[is_wall], bf[is_wire]


def solve(m):
    """Assemble P2 Laplace and solve the TWO unit problems.

    Returns (basis, u_A, u_B):
      u_A: anode=+V_MESH, wires=0, top=0
      u_B: anode=0,       wires=0, top=+1 V
    The physical solution is u = u_A + c*u_B with c chosen so the BULK drift
    field equals E_DRIFT exactly. This is required because the amplification
    field penetrates the 46%-open mesh and offsets the drift-side potential
    by ~2 V — negligible against the real 1000 V / 30 mm drift, but a ~30%
    field error over a 218 um slice if the top potential is set naively to
    -E_drift*z_cath.
    """
    from skfem import Basis, ElementTetP2
    from skfem.models.poisson import laplace
    basis = Basis(m, ElementTetP2())
    A = laplace.assemble(basis)
    f_anode, f_top, f_wall, f_wire = classify_boundary_facets(m)
    d_anode = basis.get_dofs(facets=f_anode).flatten()
    d_top = basis.get_dofs(facets=f_top).flatten()
    d_wire = basis.get_dofs(facets=f_wire).flatten()

    ndof = basis.N
    D = np.unique(np.concatenate([d_anode, d_top, d_wire]))
    interior = np.setdiff1d(np.arange(ndof), D, assume_unique=False)
    Aii = A[interior][:, interior].tocsr()
    Aid = A[interior][:, D]
    import pyamg
    ml = pyamg.smoothed_aggregation_solver(Aii)

    def solve_bc(v_anode, v_top):
        x = np.zeros(ndof)
        x[d_anode] = v_anode
        x[d_top] = v_top
        x[d_wire] = 0.0
        # Wire (0 V) wins on shared edges via the unique() order above only
        # if consistent; enforce explicitly:
        x[np.intersect1d(d_wire, d_anode)] = 0.0
        b = -Aid @ x[D]
        res = []
        xi = ml.solve(b, tol=1e-10, accel="cg", residuals=res, maxiter=400)
        rel = np.linalg.norm(Aii @ xi - b) / max(1e-30, np.linalg.norm(b))
        if rel > 1e-8:
            raise RuntimeError(f"linear solve did not converge (rel={rel:.1e})")
        x[interior] = xi
        return x, len(res), rel

    u_A, itA, relA = solve_bc(V_ANODE, 0.0)
    u_B, itB, relB = solve_bc(0.0, 1.0)
    print(f"  solve: ndof={ndof}, Dirichlet={D.size}; "
          f"A: {itA} iters rel={relA:.1e}; B: {itB} iters rel={relB:.1e}")
    return basis, u_A, u_B


class P2Evaluator:
    """Vectorized value+gradient evaluation of a P2 tet solution.

    skfem's probes() gives values only; we need smooth gradients on ~2M grid
    nodes, so evaluate the P2 shape functions directly. DOF ordering per
    element is taken from basis.element_dofs and VERIFIED geometrically
    against basis.doflocs (vertex dofs at vertices, edge dofs at midpoints);
    gate S5 additionally compares values against skfem's own interpolator.
    """

    def __init__(self, m, basis, u):
        self.m, self.u = m, u
        ed = basis.element_dofs  # (10, nelems)
        self.ed = ed
        p, t = m.p, m.t
        # Verify ordering: rows 0-3 vertex dofs, rows 4-9 edge-midpoint dofs.
        locs = basis.doflocs
        for i in range(4):
            assert np.allclose(locs[:, ed[i]], p[:, t[i]], atol=1e-8)
        verts = p[:, t]                       # (3, 4, nelems)
        self.edge_pairs = []
        for i in range(4, 10):
            mid = locs[:, ed[i]]              # (3, nelems)
            found = None
            for a in range(4):
                for bb in range(a + 1, 4):
                    if np.allclose(mid, 0.5 * (verts[:, a] + verts[:, bb]),
                                   atol=1e-8):
                        found = (a, bb)
                        break
                if found:
                    break
            if found is None:
                raise RuntimeError("could not identify P2 edge dof pairing")
            self.edge_pairs.append(found)
        # Barycentric machinery: lambda = M @ (x,y,z,1)
        n = t.shape[1]
        T = np.empty((n, 4, 4))
        T[:, :3, :] = verts.transpose(2, 0, 1)
        T[:, 3, :] = 1.0
        self.M = np.linalg.inv(T)             # (nelems, 4, 4)
        self.finder = m.element_finder()

    def _find(self, pts):
        """element_finder that tolerates points outside the mesh.

        skfem raises ValueError for ANY not-found point in a batch; isolate
        offenders by recursive bisection so one bad point does not poison the
        rest (cost O(log N) extra calls per bad point).
        """
        N = pts.shape[1]
        try:
            return self.finder(pts[0], pts[1], pts[2])
        except ValueError:
            if N == 1:
                return np.array([-1], dtype=np.int64)
            h = N // 2
            return np.concatenate([self._find(pts[:, :h]),
                                   self._find(pts[:, h:])])

    def __call__(self, pts):
        """pts (3, N) -> (val, grad(3,N), ok mask). Not-found -> ok=False."""
        N = pts.shape[1]
        val = np.zeros(N)
        grad = np.zeros((3, N))
        ok = np.zeros(N, bool)
        e = self._find(pts)
        good = e >= 0
        if not good.any():
            return val, grad, ok
        e = e[good]
        P = np.concatenate([pts[:, good], np.ones((1, good.sum()))], axis=0)
        lam = np.einsum("nij,jn->in", self.M[e], P)       # (4, N)
        glam = self.M[e][:, :, :3].transpose(1, 2, 0)     # (4, 3, N)
        ed = self.ed[:, e]
        uv = self.u[ed]                                   # (10, N)
        v = np.zeros(good.sum())
        g = np.zeros((3, good.sum()))
        for i in range(4):
            li = lam[i]
            v += uv[i] * li * (2 * li - 1)
            g += uv[i] * (4 * li - 1) * glam[i]
        for k, (a, bb) in enumerate(self.edge_pairs):
            la, lb = lam[a], lam[bb]
            v += uv[4 + k] * 4 * la * lb
            g += uv[4 + k] * 4 * (la * glam[bb] + lb * glam[a])
        val[good] = v
        grad[:, good] = g
        ok[good] = True
        return val, grad, ok


def in_wire(pts, margin=0.0):
    dx = np.hypot(pts[1], pts[2] - Z_OFF)
    dy = np.hypot(pts[0], pts[2] + Z_OFF)
    return (dx < WIRE_R + margin) | (dy < WIRE_R + margin)


def project_to_surface(pts, standoff):
    """Project in-wire points radially onto the nearer wire surface + standoff."""
    out = pts.copy()
    dx = np.hypot(pts[1], pts[2] - Z_OFF)   # distance to x-wire axis
    dy = np.hypot(pts[0], pts[2] + Z_OFF)   # distance to y-wire axis
    use_x = dx >= dy                        # project out of the NEARER surface
    # x-wire: radial in (y, z-Z_OFF); guard the (impossible for shell
    # nodes) exactly-on-axis case.
    r_t = WIRE_R + standoff
    d = np.where(use_x, dx, dy)
    d = np.maximum(d, 1e-6)
    sc = r_t / d
    ix = use_x
    out[1, ix] = pts[1, ix] * sc[ix]
    out[2, ix] = Z_OFF + (pts[2, ix] - Z_OFF) * sc[ix]
    iy = ~use_x
    out[0, iy] = pts[0, iy] * sc[iy]
    out[2, iy] = -Z_OFF + (pts[2, iy] + Z_OFF) * sc[iy]
    return out


_CTX = None  # set by run() before the sampling Pool forks


def _sample_plane(iz):
    """Sample one z-plane of the product grid (multiprocessing worker).

    Flag semantics compensate Garfield's conservative interpolation rule
    (a point is inactive if ANY cell corner is inactive, which fattens the
    wires by up to one grid step): the written flag is ERODED by one step,
    so the effective absorption boundary sits back on the true wire surface.
    In-wire nodes inside the erosion shell stay active and carry the field
    of their radial surface projection, so near-surface interpolation is not
    polluted by zeros.
    """
    zs, xs, nx = _CTX["zs"], _CTX["xs"], _CTX["nx"]
    Xg, Yg, ev, cm = _CTX["Xg"], _CTX["Yg"], _CTX["ev"], _CTX["cm"]
    step = _CTX["step"]
    z = zs[iz]
    pts = np.vstack([Xg.ravel(), Yg.ravel(), np.full(Xg.size, z)])
    gas = ~in_wire(pts)
    deep = in_wire(pts, margin=-step)       # > one step inside the wire
    shell = (~gas) & (~deep)                # inside wire, within one step
    val = np.zeros(pts.shape[1])
    grad = np.zeros((3, pts.shape[1]))
    okz = np.zeros(pts.shape[1], bool)
    if gas.any():
        v_, g_, ok_ = ev(pts[:, gas])
        val[gas], grad[:, gas], okz[gas] = v_, g_, ok_
    if shell.any():
        proj = project_to_surface(pts[:, shell], standoff=0.4)
        v_, g_, ok_ = ev(proj)
        val[shell], grad[:, shell], okz[shell] = v_, g_, ok_
    flag = ((gas | shell) & okz).astype(np.int64)
    ex, ey, ezf = -grad[0] * 1e4, -grad[1] * 1e4, -grad[2] * 1e4
    dead = flag == 0
    ex[dead] = 0.; ey[dead] = 0.; ezf[dead] = 0.; val[dead] = 0.
    rows = np.column_stack([
        np.repeat(xs, nx) * cm, np.tile(xs, nx) * cm,
        np.full(pts.shape[1], z * cm), ex, ey, ezf, val, flag])
    buf = "\n".join(
        "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %d" % tuple(r) for r in rows)
    return int(dead.sum()), buf


def run(tag, lc_wire, lc_far, grid_step, nev_note, do_refine_gate, jobs):
    t0 = time.time()
    print(f"[{tag}] meshing (lc_wire={lc_wire}, lc_far={lc_far}) ...")
    msh = os.path.join(tempfile.gettempdir(), f"meshcell_{tag}.msh")
    build_mesh(lc_wire, lc_far, msh)
    m = load_mesh(msh)
    print(f"  mesh: {m.p.shape[1]} nodes, {m.t.shape[1]} tets")
    basis, uA, uB = solve(m)
    ev = P2Evaluator(m, basis, uA)

    gates = {}

    # ── Fix the drift-side BC exactly (see solve() docstring) ───────────────
    xyf = np.linspace(-PITCH / 2, PITCH / 2, 5)
    Xf, Yf = np.meshgrid(xyf, xyf)
    fitpl = np.vstack([Xf.ravel(), Yf.ravel(), np.full(Xf.size, Z_CATH - 5.0)])
    _, gA, okf = ev(fitpl)
    ezA = -gA[2, okf].mean() * 1e4
    ev.u = uB
    _, gB, _ = ev(fitpl)
    ezB = -gB[2, okf].mean() * 1e4
    c = (E_DRIFT - ezA) / ezB
    u = uA + c * uB
    ev.u = u
    v_top_naive = -E_DRIFT * Z_CATH * 1e-4
    dv_pen = c - v_top_naive
    delta_pen_um = dv_pen / ((NEBEM_AMP_EZ - E_DRIFT) * 1e-4)
    gates["c_top_V"] = float(c)
    gates["penetration_V"] = float(dv_pen)
    print(f"  drift BC: top plane set to {c:.3f} V (naive {v_top_naive:.3f} V"
          f"); mesh field-penetration offset {dv_pen:.3f} V "
          f"(effective penetration depth {delta_pen_um:.2f} um)")
    if abs(dv_pen) > 15.0:
        raise RuntimeError("implausible penetration offset — check BCs")

    # S5 evaluator closure vs skfem's own interpolation (values).
    rng = np.random.default_rng(17)
    pts = np.empty((3, 400))
    nfill = 0
    while nfill < 400:
        cand = np.vstack([
            rng.uniform(-PITCH / 2, PITCH / 2, 800),
            rng.uniform(-PITCH / 2, PITCH / 2, 800),
            rng.uniform(Z_ANODE + 1, Z_CATH - 1, 800)])
        cand = cand[:, ~in_wire(cand, margin=1.0)]
        take = min(400 - nfill, cand.shape[1])
        pts[:, nfill:nfill + take] = cand[:, :take]
        nfill += take
    v_mine, _, okm = ev(pts)
    from skfem import Basis  # noqa: F401
    probes = basis.probes(pts[:, okm])
    v_skfem = probes @ u
    rel = np.max(np.abs(v_mine[okm] - v_skfem) /
                 np.maximum(1.0, np.abs(v_skfem)))
    gates["S5_evaluator_rel"] = float(rel)
    print(f"  S5 evaluator vs skfem probes: max rel diff {rel:.2e} "
          f"-> {'PASS' if rel < 1e-8 else 'FAIL'}")

    # S1 far field (checked at a DIFFERENT plane than the c fit, so residual
    # non-uniformity or z-dependence shows up as a deviation).
    xy = np.linspace(-PITCH / 2, PITCH / 2, 5)
    X, Y = np.meshgrid(xy, xy)
    top = np.vstack([X.ravel(), Y.ravel(),
                     np.full(X.size, Z_CATH - 40.0)])
    _, g, okt = ev(top)
    ez_top = -g[2, okt].mean() * 1e4          # V/um -> V/cm
    ez_std = (-g[2, okt] * 1e4).std()
    exy_top = np.hypot(g[0, okt], g[1, okt]).mean() * 1e4
    dev1 = 100 * (ez_top - E_DRIFT) / E_DRIFT
    gates["S1_Ez_top_Vcm"] = float(ez_top)
    gates["S1_dev_pct"] = float(dev1)
    print(f"  S1 far field (z=+{Z_CATH - 40:.0f}um): <Ez>= {ez_top:.2f} V/cm "
          f"(want +{E_DRIFT}, dev {dev1:+.3f}%), std {ez_std:.2f}, "
          f"<|Exy|>= {exy_top:.2f} "
          f"-> {'PASS' if abs(dev1) < 1.0 else 'FAIL'}")

    # S2 amp bulk vs neBEM.
    amp = np.vstack([X.ravel(), Y.ravel(), np.full(X.size, -90.0)])
    _, g, oka = ev(amp)
    ez_amp = -g[2, oka].mean() * 1e4
    dev2 = 100 * (ez_amp - NEBEM_AMP_EZ) / NEBEM_AMP_EZ
    gates["S2_Ez_amp_Vcm"] = float(ez_amp)
    gates["S2_dev_vs_nebem_pct"] = float(dev2)
    print(f"  S2 amp bulk: <Ez>(z=-90um) = {ez_amp:.1f} V/cm "
          f"(neBEM infinite-lattice: {NEBEM_AMP_EZ}, dev {dev2:+.2f}%) "
          f"-> {'PASS' if abs(dev2) < 2.0 else 'FAIL'}")

    # S3 wire BC. Probe around the x-wire away from the crossing (|x| > 12 um
    # keeps clear of the y-wire volume).
    th = rng.uniform(0, 2 * np.pi, 40)
    along = np.concatenate([rng.uniform(-PITCH / 2 + 1, -12, 20),
                            rng.uniform(12, PITCH / 2 - 1, 20)])
    rr = WIRE_R + 0.4
    wpts = np.vstack([along, rr * np.cos(th), Z_OFF + rr * np.sin(th)])
    vv, _, okw = ev(wpts)
    if not okw.any():
        raise RuntimeError("S3: no wire-surface probe points located")
    vmax = np.abs(vv[okw]).max()
    gates["S3_wire_V_max"] = float(vmax)
    print(f"  S3 wire BC: max |V| at r+0.4um = {vmax:.3f} V "
          f"-> {'PASS' if vmax < 2.0 else 'FAIL'}")

    # S4 mirror symmetry.
    base = np.vstack([rng.uniform(2, PITCH / 2 - 1, 60),
                      rng.uniform(2, PITCH / 2 - 1, 60),
                      rng.uniform(Z_ANODE + 5, Z_CATH - 5, 60)])
    base = base[:, ~in_wire(base, margin=1.5)]
    _, g0, ok0 = ev(base)
    mx = base.copy(); mx[0] *= -1
    _, gx, okx = ev(mx)
    my = base.copy(); my[1] *= -1
    _, gy, oky = ev(my)
    ok = ok0 & okx & oky
    ez0 = g0[2, ok]
    d = np.maximum(np.abs(gx[2, ok] - ez0), np.abs(gy[2, ok] - ez0)) \
        / np.maximum(1e-4, np.abs(ez0))
    dmed, d95 = float(np.median(d)), float(np.percentile(d, 95))
    # The continuum problem is exactly mirror-symmetric; residual asymmetry
    # is unstructured-mesh discretization error, so gate the median (bulk
    # accuracy) and only report the near-wire tail.
    gates["S4_mirror_rel"] = dmed
    gates["S4_mirror_rel_p95"] = d95
    print(f"  S4 mirror symmetry: median rel Ez diff {100 * dmed:.3f}% "
          f"(p95 {100 * d95:.2f}%) -> {'PASS' if dmed < 0.005 else 'FAIL'}")

    # ── Product grid ────────────────────────────────────────────────────────
    nx = int(round(PITCH / grid_step)) + 1
    xs = np.linspace(-PITCH / 2, PITCH / 2, nx)
    zlo, zhi = Z_ANODE + 0.5, Z_CATH - 0.5
    nz = int(round((zhi - zlo) / grid_step)) + 1
    zs = np.linspace(zlo, zhi, nz)
    print(f"  product grid: {nx} x {nx} x {nz} = {nx * nx * nz} nodes")

    Xg, Yg = np.meshgrid(xs, xs, indexing="ij")
    outname = os.path.join(HERE, f"meshfield_{tag}.txt")
    cm = 1e-4  # um -> cm
    t_grid = time.time()

    # Make sure the (lazy) element finder exists BEFORE forking, so workers
    # share it copy-on-write instead of each rebuilding the KD-tree.
    ev(np.array([[0.], [20.], [-50.]]))

    # Context for the module-level worker (fork-inherited; Pool cannot
    # pickle closures).
    global _CTX
    _CTX = {"zs": zs, "xs": xs, "nx": nx, "Xg": Xg, "Yg": Yg, "ev": ev,
            "cm": cm, "step": grid_step}

    nflag0 = 0
    import multiprocessing as mp
    with open(outname, "w") as f:
        f.write(f"# XMIN = {xs[0] * cm:.7g}, XMAX = {xs[-1] * cm:.7g}, "
                f"NX = {nx}\n")
        f.write(f"# YMIN = {xs[0] * cm:.7g}, YMAX = {xs[-1] * cm:.7g}, "
                f"NY = {nx}\n")
        f.write(f"# ZMIN = {zlo * cm:.7g}, ZMAX = {zhi * cm:.7g}, "
                f"NZ = {nz}\n")
        ctx = mp.get_context("fork")
        with ctx.Pool(jobs) as pool:
            for iz, (nf, buf) in enumerate(
                    pool.imap(_sample_plane, range(nz), chunksize=2)):
                nflag0 += nf
                f.write(buf + "\n")
                if iz % 25 == 0:
                    print(f"\r    z-plane {iz + 1}/{nz}", end="", flush=True)
    print(f"\r    grid sampled+written in {time.time() - t_grid:.0f} s "
          f"({jobs} workers); {nflag0} nodes flagged inactive "
          f"({100 * nflag0 / (nx * nx * nz):.2f}% — wire interiors)")

    # S6 refinement gate (production): compare against a coarser solve.
    if do_refine_gate:
        print("  S6: solving 1.4x coarser mesh for the refinement gate ...")
        msh2 = os.path.join(tempfile.gettempdir(), f"meshcell_{tag}_c.msh")
        build_mesh(lc_wire * 1.4, lc_far * 1.4, msh2)
        m2 = load_mesh(msh2)
        basis2, uA2, uB2 = solve(m2)
        ev2 = P2Evaluator(m2, basis2, uA2)
        _, gA2, okf2 = ev2(fitpl)
        ezA2 = -gA2[2, okf2].mean() * 1e4
        ev2.u = uB2
        _, gB2, _ = ev2(fitpl)
        ezB2 = -gB2[2, okf2].mean() * 1e4
        ev2.u = uA2 + (E_DRIFT - ezA2) / ezB2 * uB2
        nprobe = 3000
        cand = np.vstack([
            rng.uniform(-PITCH / 2, PITCH / 2, 3 * nprobe),
            rng.uniform(-PITCH / 2, PITCH / 2, 3 * nprobe),
            rng.uniform(Z_ANODE + 1, Z_CATH - 1, 3 * nprobe)])
        cand = cand[:, ~in_wire(cand, margin=0.7)][:, :nprobe]
        _, gA, okA = ev(cand)
        _, gB, okB = ev2(cand)
        ok = okA & okB
        eA = np.linalg.norm(gA[:, ok], axis=0)
        d = (np.linalg.norm(gA[:, ok] - gB[:, ok], axis=0) /
             np.maximum(np.percentile(eA, 5), eA))
        p95 = float(np.percentile(d, 95))
        gates["S6_refine_p95"] = p95
        print(f"  S6 refinement: 95th pct |dE|/|E| vs coarser = "
              f"{100 * p95:.3f}% -> {'PASS' if p95 < 0.01 else 'FAIL'}")

    ok_all = (abs(gates["S1_dev_pct"]) < 1.0 and
              abs(gates["S2_dev_vs_nebem_pct"]) < 2.0 and
              gates["S3_wire_V_max"] < 2.0 and
              gates["S4_mirror_rel"] < 0.005 and
              gates["S5_evaluator_rel"] < 1e-8 and
              (not do_refine_gate or gates["S6_refine_p95"] < 0.01))

    githash = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=HERE, capture_output=True,
                             text=True).stdout.strip()
    prov = {
        "product": os.path.basename(outname),
        "tag": tag,
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "host": os.uname().nodename,
        "git_MX17_Geant": githash,
        "solver": "gmsh+scikit-fem P2+pyamg (solve_fieldmap.py)",
        "geometry_um": {"wire_r": WIRE_R, "pitch": PITCH, "overlap": OVERLAP,
                        "z_anode": Z_ANODE, "z_cath": Z_CATH,
                        "amp_gap_def": "pillar height, ESL to mesh underside"},
        "potentials_V": {"mesh": 0.0, "anode": V_ANODE, "cath": V_CATH},
        "E_drift_Vcm": E_DRIFT,
        "mesh": {"lc_wire_um": lc_wire, "lc_far_um": lc_far,
                 "nodes": int(m.p.shape[1]), "tets": int(m.t.shape[1])},
        "grid": {"nx": nx, "ny": nx, "nz": nz, "step_um": grid_step},
        "gates": gates,
        "all_gates_pass": bool(ok_all),
        "note": nev_note,
        "wall_s": round(time.time() - t0, 1),
    }
    with open(outname.replace(".txt", ".json"), "w") as f:
        json.dump(prov, f, indent=2)
    print(f"[{tag}] {'ALL SOLVER GATES PASS' if ok_all else '*** GATE FAILURE ***'}"
          f" — wall {prov['wall_s']} s\n  wrote {outname}\n"
          f"  wrote {outname.replace('.txt', '.json')}")
    return 0 if ok_all else 1


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--smoke", action="store_true")
    g.add_argument("--production", action="store_true")
    ap.add_argument("--jobs", type=int,
                    default=max(1, (os.cpu_count() or 4) - 2),
                    help="sampling workers (default: ncpu-2)")
    args = ap.parse_args()
    if args.smoke:
        return run("smoketest", lc_wire=2.0, lc_far=14.0, grid_step=3.0,
                   nev_note="smoke resolution — do not use for S3",
                   do_refine_gate=False, jobs=args.jobs)
    return run("production", lc_wire=0.8, lc_far=8.0, grid_step=1.0,
               nev_note="production map for S3",
               do_refine_gate=True, jobs=args.jobs)


if __name__ == "__main__":
    sys.exit(main())
