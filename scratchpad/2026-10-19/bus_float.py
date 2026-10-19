"""Weighting potential of one X channel with the resistive stripes as FLOATING conductors (2026-10-19 evening).

Why: the det3 June-15 X copy lives in events whose charge sits over a Y pad, is the same at every Y pad, and the static
solve (bus_wpot.py, ESL as a plain dielectric surface) gives ~0 on the data's own estimator (fwd_parity.py). At 2 MOhm/sq
a 550 um stripe equalises across its width in ~26 ns, so at signal times each stripe acts as a conductor with zero net
charge. If a stripe over row j also covers the gap toward row j+1, where row j+1's bus runs, charge on that stripe
couples into line j+1 (copy into the next line, in phase, X only, large only when the X central is small = Y-pad events).

Geometry as bus_wpot.py (channel = row 0 along x, own pads at x = +390 mod 1560, own bus at y = -390 on L5), plus stripes
along x (parallel to the X rows, as the tau_Y > tau_X anisotropy says) of width W at pitch 780 (no beat: the data show no
31.2 mm period), centred at y = OFF + 780 m on the cell layer just above the kapton. Floating: phi = phi0 + sum_t V_t psi_t
with V from zero net charge on every stripe (capacitance matrix from the per-stripe solves). Stripes are infinite along x
(periodic, 1.56 mm), i.e. the fast-spreading limit along the stripe too.
--stripes none|float|ground. Output: the ESL map (charge on a stripe sees the stripe's potential; in a gap, phi) as
<out>_esl.npy for fwd_parity.py, and the old-style row copy in <out>.json.
Usage: bus_float.py --out bus_par/f_off200.json [--off 200] [--w 550] [--stripes float] [--dxy 20] [--no-bus] [--h1 264]"""
import argparse, time, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

ap = argparse.ArgumentParser()
ap.add_argument('--h1', type=float, default=264.); ap.add_argument('--h2', type=float, default=264.)
ap.add_argument('--gap', type=float, default=100.); ap.add_argument('--dxy', type=float, default=20.)
ap.add_argument('--no-bus', action='store_true'); ap.add_argument('--nrow', type=int, default=4)
ap.add_argument('--bottom', type=float, default=1400.); ap.add_argument('--out', required=True)
ap.add_argument('--stripes', default='float', choices=['none', 'float', 'ground'])
ap.add_argument('--w', type=float, default=550.); ap.add_argument('--off', type=float, default=0.)
a = ap.parse_args()
P = 780.; PAD = P - a.gap; CU = 26.; GLUE = 19.; KAP = 50.; AMP = 150.; Z_ESL = GLUE + KAP; Z_MESH = Z_ESL + AMP
z5 = (-CU - a.h1 - CU, -CU - a.h1); z6 = (z5[0] - a.h2 - CU, z5[0] - a.h2)
d = a.dxy; nx = int(round(2 * P / d)); ny = int(round(2 * a.nrow * P / d))
xc = (np.arange(nx) + .5) * d - P / 2; yc = (np.arange(ny) + .5) * d - a.nrow * P
def zgrid():
    knots = sorted({Z_MESH, Z_ESL, Z_ESL + 13., GLUE, 0., -CU, *z5, *z6, -a.bottom})[::-1]
    zs = [knots[0]]
    for hi, lo in zip(knots[:-1], knots[1:]):
        n = max(1, int(np.ceil((hi - lo) / (13. if lo > -700 else 60.)))); zs += list(np.linspace(hi, lo, n + 1)[1:])
    return np.array(zs)
zn = zgrid(); zc = .5 * (zn[:-1] + zn[1:]); dz = zn[:-1] - zn[1:]; nz = len(zc)
kE = int(np.flatnonzero(zc > Z_ESL)[-1])                       # the stripe / ESL cell layer
eps = np.empty((nx, ny, nz))
for k, z in enumerate(zc):
    eps[:, :, k] = 1. if z > Z_ESL else 3.5 if z > GLUE else 3.2 if z > -CU else 4.4
X, Y = np.meshgrid(xc, yc, indexing='ij')
def pad_mask(cx, cy): return (np.abs(X - cx) <= PAD / 2) & (np.abs(Y - cy) <= PAD / 2)
def circ(cx, cy, r): return (X - cx) ** 2 + (Y - cy) ** 2 <= r ** 2
def wrapx(c): return [c - 2 * P, c, c + 2 * P]
own_pad = np.zeros((nx, ny), bool); oth_pad = np.zeros_like(own_pad); own_via = np.zeros_like(own_pad); oth_via = np.zeros_like(own_pad)
own_ring = np.zeros_like(own_pad); oth_ring = np.zeros_like(own_pad)
for j in range(-a.nrow, a.nrow + 1):
    for cx0 in (-P / 2, P / 2):
        for cx in wrapx(cx0):
            own = (j == 0 and cx0 == P / 2)
            (own_pad if own else oth_pad)[...] |= pad_mask(cx, j * P)
            (own_via if own else oth_via)[...] |= circ(cx, j * P, 125.)
            (own_ring if own else oth_ring)[...] |= circ(cx, j * P, 250.)
bus_own = np.zeros_like(own_pad); bus_oth = np.zeros_like(own_pad)
for j in range(-a.nrow, a.nrow + 1):
    (bus_own if j == 0 else bus_oth)[...] |= np.abs(Y - (j * P - P / 2)) <= 50.
for j in range(-a.nrow, a.nrow + 1):
    cx0 = P / 2 if j % 2 == 0 else -P / 2
    for cx in wrapx(cx0):
        (bus_own if j == 0 else bus_oth)[...] |= (np.abs(X - cx) <= 50.) & (Y <= j * P) & (Y >= j * P - P / 2)
if a.no_bus: bus_own[:] = False
oth_lines = np.zeros_like(own_pad)
for cx0 in (-P / 2, P / 2):
    for cx in wrapx(cx0):
        oth_lines |= np.abs(X - (cx + P / 2)) <= 50.
# stripes: one per row j, centred at y = OFF + j P, width W, kept inside the domain
stripes = []
if a.stripes != 'none':
    for j in range(-a.nrow, a.nrow + 1):
        m = np.abs(Y - (a.off + j * P)) <= a.w / 2
        if m[0].any(): stripes.append(m)
# Dirichlet classes: 0 = grounded electrodes, 1 = own channel, 2 + t = stripe t
D = np.full((nx, ny, nz), -1, np.int16)
for k, z in enumerate(zc):
    if -CU < z < 0:
        D[:, :, k][oth_pad] = 0; D[:, :, k][own_pad] = 1
    if z < -CU and z > z6[0] - 100:
        D[:, :, k][oth_via] = 0; D[:, :, k][own_via] = 1
    if z5[0] < z < z5[1]:
        D[:, :, k][oth_ring | bus_oth] = 0; D[:, :, k][own_ring | bus_own] = 1; D[:, :, k][own_via] = 1
    if z6[0] < z < z6[1]:
        D[:, :, k][oth_ring | oth_lines] = 0; D[:, :, k][own_ring] = 1
for t, m in enumerate(stripes):
    D[:, :, kE][m] = 0 if a.stripes == 'ground' else 2 + t
free = D < 0; idx = -np.ones(D.shape, np.int64); idx[free] = np.arange(free.sum()); N = int(free.sum())
hm = lambda e1, e2: 2 * e1 * e2 / (e1 + e2)
DZ = dz[None, None, :]; dzf = .5 * (dz[:-1] + dz[1:])
faces = [  # (index arrays a, b, class arrays, conductance) for x (periodic), y, z neighbours
    (idx, np.roll(idx, -1, 0), D, np.roll(D, -1, 0), hm(eps, np.roll(eps, -1, 0)) * DZ),
    (idx[:, :-1], idx[:, 1:], D[:, :-1], D[:, 1:], hm(eps[:, :-1], eps[:, 1:]) * DZ),
    (idx[:, :, :-1], idx[:, :, 1:], D[:, :, :-1], D[:, :, 1:], hm(eps[:, :, :-1], eps[:, :, 1:]) * d * d / dzf[None, None, :])]
rows, cols, vals = [], [], []; diag = np.zeros(N)
bpart = []          # (free index, class, conductance) for every free-Dirichlet face
for ia, ib, ca, cb, g in faces:
    fa = ia >= 0; fb = ib >= 0; m = fa & fb
    rows.extend([ia[m], ib[m]]); cols.extend([ib[m], ia[m]]); vals.extend([-g[m], -g[m]])
    np.add.at(diag, ia[m], g[m]); np.add.at(diag, ib[m], g[m])
    m = fa & ~fb; np.add.at(diag, ia[m], g[m]); bpart.append((ia[m], cb[m], g[m]))
    m = ~fa & fb; np.add.at(diag, ib[m], g[m]); bpart.append((ib[m], ca[m], g[m]))
for side in (0, -1):                                            # y edges at 0, half a cell outside
    gi = eps[:, side] * d * dz[None, :] / (d / 2); ii = idx[:, side]; m = ii >= 0; np.add.at(diag, ii[m], gi[m])
for k in (0, -1):                                               # mesh on top, ground at the bottom
    gi = eps[:, :, k] * d * d / (dz[k] / 2); ii = idx[:, :, k]; m = ii >= 0; np.add.at(diag, ii[m], gi[m])
A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N)) + sp.diags(diag)
fi = np.concatenate([p[0] for p in bpart]); fc = np.concatenate([p[1] for p in bpart]); fg = np.concatenate([p[2] for p in bpart])
ncls = 2 + len(stripes); M = sp.diags(1 / diag)
def solve(cls):
    """Field with Dirichlet class `cls` at 1, every other electrode at 0."""
    b = np.zeros(N); m = fc == cls; np.add.at(b, fi[m], fg[m])
    x, info = spl.cg(A, b, rtol=1e-9, maxiter=60000, M=M)
    phi = np.where(D == cls, 1., 0.).astype(float); phi[free] = x
    return phi, info
def charge(phi, cls):
    """Net (scaled) charge on the cells of class cls: sum over faces to non-member cells of g (phi_c - phi_n)."""
    q = 0.; mem = D == cls
    for ax, g in ((0, faces[0][4]), (1, faces[1][4]), (2, faces[2][4])):
        if ax == 0:
            ma, mb, pa, pb = mem, np.roll(mem, -1, 0), phi, np.roll(phi, -1, 0)
        elif ax == 1:
            ma, mb, pa, pb = mem[:, :-1], mem[:, 1:], phi[:, :-1], phi[:, 1:]
        else:
            ma, mb, pa, pb = mem[:, :, :-1], mem[:, :, 1:], phi[:, :, :-1], phi[:, :, 1:]
        q += np.sum(g[ma & ~mb] * (pa - pb)[ma & ~mb]) + np.sum(g[mb & ~ma] * (pb - pa)[mb & ~ma])
    return q
t0 = time.time()
phi0, info0 = solve(1)
res = dict(args=vars(a), n_stripes=len(stripes), N=N)
if a.stripes == 'float' and stripes:
    psi = [solve(2 + t)[0] for t in range(len(stripes))]
    Cm = np.array([[charge(psi[t], 2 + s) for t in range(len(stripes))] for s in range(len(stripes))])
    q0 = np.array([charge(phi0, 2 + s) for s in range(len(stripes))])
    Vs = np.linalg.solve(Cm, -q0)
    phi = phi0 + sum(v * p for v, p in zip(Vs, psi))
    res['stripe_V'] = Vs.round(5).tolist(); res['stripe_y'] = [float(a.off + j * P) for j in range(-a.nrow, a.nrow + 1)]
else:
    phi = phi0
res['seconds'] = round(time.time() - t0, 1); res['cg_info'] = int(info0)
E = phi[:, :, kE]
for s in (300., 450.):
    S = E.mean(0); I = {}
    for k in (-1, 0, 1):
        cs = np.linspace(k * P - P / 2, k * P + P / 2, 40)
        I[k] = float(np.mean([np.sum(S * np.exp(-.5 * ((yc - c) / s) ** 2)) / np.sum(np.exp(-.5 * ((yc - c) / s) ** 2)) for c in cs]))
    res[f'sigma{int(s)}'] = dict(I=I, copy=(I[-1] - I[1]) / I[0])
np.save(a.out.replace('.json', '_esl.npy'), E.astype(np.float32))
json.dump(res, open(a.out, 'w'), indent=1); print(json.dumps({k: v for k, v in res.items() if k != 'args'}))
