"""Prompt weighting potential of ONE readout channel including its buried one-sided bus (2026-10-19).

Why: det3's X view shows a one-sided in-phase copy (NB_ASYM.md). The gerbers put every bus line midway between two pad
rows, joined to its own pads by 390 um stubs that all point to the same side, so each channel's bus runs under the
inter-pad gap toward ONE neighbour. The response kernels (pads only, no buses) are symmetric by construction.
This solve asks how big a one-sided copy that geometry makes, against the bus depth and the gap width.

Geometry (z up, um; pad top at z = 0): mesh (V=0) at 69 + 150 = 219; gas 69..219; kapton (3.5) 19..69 with the resistive
layer on top at 69 (a dielectric surface at prompt time); glue (3.2) 0..19 and in the inter-pad gaps; pads -26..0;
FR4 (4.4) below; bus layer L5 at depth h1 below the pads, L6 a further h2 below; ground at the bottom (BOTTOM_UM).
In-plane: the channel is row 0 along x (the bus direction), its pads every other column (x = 0.78 i, i odd);
the even-column pads of row 0 and all pads of other rows belong to other channels (V = 0). Every pad has a 0.25 mm via
from the pad down to L7 with a 0.5 mm ring on L5 and L6; the own channel's vias, rings, stubs (0.39 mm toward -y) and
bus (0.1 mm wide at y = -0.39) are at V = 1, everything else at 0. The other view's buses (along y, midway between
columns, stubs +x) sit on the other bus layer at V = 0. Periodic in x (1.56 mm), V = 0 at |y| = NROW rows.

Output: phi(x, y) on the resistive-layer plane; folded with a Gaussian cloud (sigma S) whose centre is uniform over a
row: I_k = mean induced on channel 0 for clouds centred on row k. The copy (as measured) = (I_{-1} - I_{+1}) / I_0.
Usage: bus_wpot.py [--h1 264] [--h2 264] [--layer 5|6] [--gap 100] [--dxy 20] [--sigma 300,400] [--no-bus]"""
import argparse, time, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

ap = argparse.ArgumentParser()
ap.add_argument('--h1', type=float, default=264.); ap.add_argument('--h2', type=float, default=264.)
ap.add_argument('--layer', type=int, default=5); ap.add_argument('--gap', type=float, default=100.)
ap.add_argument('--dxy', type=float, default=20.); ap.add_argument('--sigma', default='300,450')
ap.add_argument('--no-bus', action='store_true'); ap.add_argument('--nrow', type=int, default=5)
ap.add_argument('--bottom', type=float, default=1400.); ap.add_argument('--out')
a = ap.parse_args()
P = 780.; PAD = P - a.gap; CU = 26.; GLUE = 19.; KAP = 50.; AMP = 150.; Z_ESL = GLUE + KAP; Z_MESH = Z_ESL + AMP
z5 = (-CU - a.h1 - CU, -CU - a.h1); z6 = (z5[0] - a.h2 - CU, z5[0] - a.h2)
own_z, oth_z = (z5, z6) if a.layer == 5 else (z6, z5)

# grid: uniform in x, y; z nodes chosen to land on every interface
d = a.dxy; nx = int(round(2 * P / d)); ny = int(round(2 * a.nrow * P / d))
xc = (np.arange(nx) + .5) * d - P / 2                      # pads at x = 0.78 i -> own pad at x = +390 (i = 1 -> shift)
yc = (np.arange(ny) + .5) * d - a.nrow * P
def zgrid():
    knots = sorted({Z_MESH, Z_ESL, GLUE, 0., -CU, *z5, *z6, -a.bottom})[::-1]
    zs = [knots[0]]
    for hi, lo in zip(knots[:-1], knots[1:]):
        n = max(1, int(np.ceil((hi - lo) / (13. if lo > -700 else 60.)))); zs += list(np.linspace(hi, lo, n + 1)[1:])
    return np.array(zs)
zn = zgrid(); zc = .5 * (zn[:-1] + zn[1:]); dz = zn[:-1] - zn[1:]; nz = len(zc)

eps = np.empty((nx, ny, nz)); V = np.full((nx, ny, nz), np.nan)          # NaN = free
for k, z in enumerate(zc):
    eps[:, :, k] = 1. if z > Z_ESL else 3.5 if z > GLUE else 3.2 if z > -CU else 4.4
X, Y = np.meshgrid(xc, yc, indexing='ij')
# pads: centres at (x = -390 + 780 i, y = 780 j); own pads: row 0, i odd -> x = +390 (mod 1560)
def pad_mask(cx, cy): return (np.abs(X - cx) <= PAD / 2) & (np.abs(Y - cy) <= PAD / 2)
def circ(cx, cy, r): return (X - cx) ** 2 + (Y - cy) ** 2 <= r ** 2
def wrapx(c): return [c - 2 * P, c, c + 2 * P]
own_pad = np.zeros((nx, ny), bool); oth_pad = np.zeros((nx, ny), bool); own_via = np.zeros_like(own_pad); oth_via = np.zeros_like(own_pad)
own_ring = np.zeros_like(own_pad); oth_ring = np.zeros_like(own_pad)
for j in range(-a.nrow, a.nrow + 1):
    for cx0 in (-P / 2, P / 2):
        for cx in wrapx(cx0):
            own = (j == 0 and cx0 == P / 2)
            (own_pad if own else oth_pad)[...] |= pad_mask(cx, j * P)
            (own_via if own else oth_via)[...] |= circ(cx, j * P, 125.)
            (own_ring if own else oth_ring)[...] |= circ(cx, j * P, 250.)
# buses on the own layer: every row j has a bus at y = j P - P/2 (stubs -y); the own one is row 0
bus_own = np.zeros_like(own_pad); bus_oth = np.zeros_like(own_pad)
for j in range(-a.nrow, a.nrow + 1):
    line = np.abs(Y - (j * P - P / 2)) <= 50.
    (bus_own if j == 0 else bus_oth)[...] |= line
# stubs: from each pad of that bus's channel (row j, every other column: own row -> x = +390) toward -y by 390
for j in range(-a.nrow, a.nrow + 1):
    cx0 = P / 2 if j % 2 == 0 else -P / 2                      # checkerboard: next row's channel owns the other columns
    for cx in wrapx(cx0):
        stub = (np.abs(X - cx) <= 50.) & (Y <= j * P) & (Y >= j * P - P / 2)
        (bus_own if j == 0 else bus_oth)[...] |= stub
if a.no_bus: bus_own[:] = False
# other view's buses on the other layer: along y at x = cx + 390 (stubs +x), all at 0
oth_lines = np.zeros_like(own_pad)
for cx0 in (-P / 2, P / 2):
    for cx in wrapx(cx0):
        oth_lines |= np.abs(X - (cx + P / 2)) <= 50.
for k, z in enumerate(zc):
    if -CU < z < 0:
        V[:, :, k][oth_pad] = 0.; V[:, :, k][own_pad] = 1.
    if z < -CU:                                                  # vias run from the pads down to L7 (below L6)
        if z > z6[0] - 100:
            V[:, :, k][oth_via] = 0.; V[:, :, k][own_via] = 1.
    if own_z[0] < z < own_z[1]:
        V[:, :, k][oth_ring | bus_oth] = 0.; V[:, :, k][own_ring | bus_own] = 1.
        V[:, :, k][own_via] = 1.
    if oth_z[0] < z < oth_z[1]:
        V[:, :, k][oth_ring | oth_lines] = 0.; V[:, :, k][own_ring] = 1.

# finite-volume Laplacian, harmonic eps on faces; Dirichlet: mesh above (0) at z = Z_MESH, ground below, 0 at |y| edges
free = np.isnan(V); idx = -np.ones(V.shape, np.int64); idx[free] = np.arange(free.sum()); N = int(free.sum())
rows, cols, vals = [], [], []; b = np.zeros(N); diag = np.zeros(N)
hm = lambda e1, e2: 2 * e1 * e2 / (e1 + e2)
DZ = dz[None, None, :]
g = hm(eps, np.roll(eps, -1, 0)) * d * DZ / d                    # x faces (periodic)
def couple_arr(ia, ib, va, vb, g):
    fa = ia >= 0; fb = ib >= 0
    m = fa & fb; rows.extend([ia[m], ib[m]]); cols.extend([ib[m], ia[m]]); vals.extend([-g[m], -g[m]])
    np.add.at(diag, ia[m], g[m]); np.add.at(diag, ib[m], g[m])
    m = fa & ~fb; np.add.at(diag, ia[m], g[m]); np.add.at(b, ia[m], g[m] * vb[m])
    m = ~fa & fb; np.add.at(diag, ib[m], g[m]); np.add.at(b, ib[m], g[m] * va[m])
couple_arr(idx, np.roll(idx, -1, 0), V, np.roll(V, -1, 0), g)
g = hm(eps[:, :-1], eps[:, 1:]) * d * DZ / d
couple_arr(idx[:, :-1], idx[:, 1:], V[:, :-1], V[:, 1:], g)
for side in (0, -1):                                              # y edges: V = 0 half a cell outside
    gi = eps[:, side] * d * dz[None, :] / (d / 2); ii = idx[:, side]; m = ii >= 0; np.add.at(diag, ii[m], gi[m])
dzf = .5 * (dz[:-1] + dz[1:])
g = hm(eps[:, :, :-1], eps[:, :, 1:]) * d * d / dzf[None, None, :]
couple_arr(idx[:, :, :-1], idx[:, :, 1:], V[:, :, :-1], V[:, :, 1:], g)
for k, vv in ((0, 0.), (-1, 0.)):                                 # mesh on top, ground at bottom (half cell away)
    gi = eps[:, :, k] * d * d / (dz[k] / 2); ii = idx[:, :, k]; m = ii >= 0; np.add.at(diag, ii[m], gi[m]); np.add.at(b, ii[m], gi[m] * vv)
A = sp.csr_matrix((np.concatenate(vals) if vals else [], (np.concatenate(rows), np.concatenate(cols))), shape=(N, N)) + sp.diags(diag)
t0 = time.time(); M = sp.diags(1 / diag)
x, info = spl.cg(A, b, rtol=1e-8, maxiter=40000, M=M)
phi = V.copy(); phi[free] = x
kE = int(np.flatnonzero(zc > Z_ESL)[-1])                         # cell just above the kapton top = the resistive layer
S = phi[:, :, kE].mean(0)                                         # averaged along x (cloud centres uniform in x)
# fold: cloud centre uniform over row k, Gaussian sigma in y; x already averaged
res = dict(args=vars(a), cg_info=int(info), seconds=round(time.time() - t0, 1), N=N,
           pad_centre=float(phi[np.argmin(np.abs(xc - P / 2)), np.argmin(np.abs(yc)), kE]))
for s in [float(v) for v in a.sigma.split(',')]:
    I = {}
    for k in (-2, -1, 0, 1, 2):
        cs = np.linspace(k * P - P / 2, k * P + P / 2, 40)
        I[k] = float(np.mean([np.sum(S * np.exp(-.5 * ((yc - c) / s) ** 2)) / np.sum(np.exp(-.5 * ((yc - c) / s) ** 2)) for c in cs]))
    res[f'sigma{int(s)}'] = dict(I=I, copy=(I[-1] - I[1]) / I[0], copy2=(I[-2] - I[2]) / I[0])
# profile across rows at the ESL (x-averaged) for the record
res['profile_y_um'] = yc[::5].round(1).tolist(); res['profile_phi'] = S[::5].round(6).tolist()
print(json.dumps({k: v for k, v in res.items() if not k.startswith('profile')}, indent=1))
if a.out: json.dump(res, open(a.out, 'w'))
