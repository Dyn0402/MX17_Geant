"""2-D weighting potentials with the ESL strips as FLOATING conductors (2026-10-21 evening), for ion_path_toy.py --esl float.
At 5 MOhm/sq a 550 um strip equalises across its width in ~65 ns, so during the ~200 ns ion drift the strips act as
conductors with zero net charge (route 9's limit, here in the 2-D cross-section). Same geometry and frame as
ion_path_toy.py (31.2 mm periodic box, 40 pads 680/780 at z = 0 with linear gap ramps, dielectric 0..69 um eps 3.39, gas to
the mesh at 219 um); strips of width W at k*800 um on the z = 69 um interface.
psi = psi_0 + sum_t V_t psi_t, psi_0: pad drive with all strips at 0; psi_t: strip t at 1; V from zero net charge on every
strip (C V = -q_0, charges = residual of the full operator on the strip nodes). One sparse LU for all 79 solves.
Output: float_psi_w<W>.npz with PSI (40, NZ, NX) on the toy's ion heights (zs) and grid (x).
Usage: float_psi.py [--w 550] [--zb 10] [--dz 3]"""
import argparse, os, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl

ap = argparse.ArgumentParser()
ap.add_argument('--w', type=float, default=550.); ap.add_argument('--zb', type=float, default=10.)
ap.add_argument('--dz', type=float, default=3.); ap.add_argument('--dx', type=float, default=5.)
a = ap.parse_args()
P, SP, BOX, PADW = 780., 800., 31200., 680.
A_D, B_M = 69., 219.; EPS1 = (3.2 * 19. + 3.5 * 50.) / 69.
NX = int(round(BOX / a.dx)); x = (np.arange(NX) + .5) * a.dx
NZ = int(round(B_M / a.dz)); zn = np.arange(NZ + 1) * a.dz; jE = int(round(A_D / a.dz))
assert abs(zn[jE] - A_D) < 1e-9, 'ESL must sit on a z node'
ORG = -(512 - 1) / 2 * P; padc = np.mod(ORG + np.arange(40) * P, BOX)

# node (i, j), j = 0..NZ; index j-major inside each column i -> narrow band in z, periodic in x
idx = lambda i, j: (i % NX) * (NZ + 1) + j
N = NX * (NZ + 1)
ev = np.where(zn[1:] <= A_D + 1e-9, EPS1, 1.)               # vertical link j -> j+1 permittivity
eh = np.where(zn < A_D - 1e-9, EPS1, np.where(np.abs(zn - A_D) < 1e-9, (EPS1 + 1) / 2, 1.))   # horizontal links per row
r, c, v = [], [], []
I = np.arange(NX)
for j in range(NZ + 1):
    w_e = eh[j] / a.dx ** 2
    for di in (-1, 1):
        r.append(idx(I, j)); c.append(idx(I + di, j)); v.append(np.full(NX, -w_e))
    diag = 2 * w_e * np.ones(NX)
    if j > 0:
        w = ev[j - 1] / a.dz ** 2; r.append(idx(I, j)); c.append(idx(I, j - 1)); v.append(np.full(NX, -w)); diag += w
    if j < NZ:
        w = ev[j] / a.dz ** 2; r.append(idx(I, j)); c.append(idx(I, j + 1)); v.append(np.full(NX, -w)); diag += w
    r.append(idx(I, j)); c.append(idx(I, j)); v.append(diag)
L = sp.csr_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))), shape=(N, N))

ph = np.mod(x, SP); onstrip = (ph < a.w / 2) | (ph >= SP - a.w / 2)
sid = np.rint(x / SP).astype(int) % int(round(BOX / SP))    # strip id 0..38
dirich = np.zeros(N, bool); dirich[idx(I, 0)] = True; dirich[idx(I, NZ)] = True
snodes = idx(I[onstrip], jE); dirich[snodes] = True
F = np.flatnonzero(~dirich); D = np.flatnonzero(dirich)
LFF = L[F][:, F].tocsc(); LFD = L[F][:, D]
t0 = time.time(); lu = spl.splu(LFF, permc_spec='COLAMD'); print(f'LU {time.time() - t0:.0f} s, N free {len(F)}', flush=True)

def solve(phiD):
    psi = np.zeros(N); psi[D] = phiD[D]; psi[F] = lu.solve(-(LFD @ phiD[D])); return psi
def strip_charges(psi):
    q = (L @ psi)[snodes]; return np.bincount(sid[onstrip], weights=q, minlength=39)

def pad_drive(xc):
    d = ((x - xc + BOX / 2) % BOX) - BOX / 2; ad = np.abs(d)
    return np.where(ad <= PADW / 2, 1., np.clip((P / 2 + (P - PADW) / 2 - ad) / (P - PADW), 0, 1))

NS = 39; Pst = []
for t in range(NS):
    phi = np.zeros(N); phi[snodes[sid[onstrip] == t]] = 1.; Pst.append(solve(phi))
Cm = np.stack([strip_charges(p) for p in Pst], 1)            # C[t, t'] = charge on t with t' at 1
print(f'strip solves done {time.time() - t0:.0f} s; C diag {Cm[0, 0]:.3g}, off {Cm[1, 0]:.3g}', flush=True)
Pst = np.stack(Pst)                                            # (39, N)

Z0 = A_D + a.zb; zs = np.linspace(Z0, B_M, 120)               # the toy's ion heights
PSI = np.empty((40, len(zs), NX))
for k in range(40):
    phi = np.zeros(N); phi[idx(I, 0)] = pad_drive(padc[k]); p0 = solve(phi)
    V = np.linalg.solve(Cm, -strip_charges(p0)); psi = (p0 + V @ Pst).reshape(NX, NZ + 1)
    for iz, z in enumerate(zs):                                # linear interpolation in z
        jf = z / a.dz; j0 = min(int(jf), NZ - 1); f = jf - j0
        PSI[k, iz] = (1 - f) * psi[:, j0] + f * psi[:, j0 + 1]
    if k == 0:
        print(f'pad 0: strip potentials {np.round(np.sort(V)[-3:], 3)}; psi over own pad at ESL '
              f'{psi[int(padc[0] / a.dx), jE]:.3f}', flush=True)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'float_psi_w{a.w:g}.npz')
np.savez_compressed(out, PSI=PSI.astype(np.float32), zs=zs, x=x, w=a.w)
print('wrote', out, f'{time.time() - t0:.0f} s')
