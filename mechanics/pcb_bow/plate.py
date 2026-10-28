#!/usr/bin/env python3
"""plate.py -- large-deflection two-layer plate model of the MX17 readout board over the support-plate aperture.

The board (layer 1: PCB, 1.2-2.4 mm FR4) and its backing (layer 2: 5 mm Rohacell, optionally with the 25 um mylar
skin bonded to it) span the 402 x 402 mm through-aperture of the 8 mm Al support plate.  Both layers share one
transverse deflection w(x, y) (they are in contact: the gas overpressure pushes the board onto the foam), and each has
its own in-plane displacement (u_i, v_i) at its own neutral surface.  They are coupled by an interface shear spring of
stiffness k(x, y) [N/m^3] acting on the slip between the PCB's lower face and the foam's upper face:

    slip_x = (u1 - u2) + d * w_x,     d = distance between the two layers' neutral surfaces.

k = 0 is a loose backing (two plates bending independently, D = D1 + D2), k -> inf a perfect bond (full composite
section), and a finite or patchy k a taped or partly debonded backing.

Energy (z up = toward the cathode, w > 0 toward the cathode; overpressure and, on the bench, gravity push w < 0):
  bending   1/2 sum_i D_i [w_xx^2 + w_yy^2 + 2 nu w_xx w_yy + 2 (1 - nu) w_xy^2]       (Kirchhoff, own axes)
  membrane  1/2 sum_i eps_i^T A_i eps_i,  eps = sym grad u_i + 1/2 grad w grad w^T       (von Karman)
  interface 1/2 k (slip_x^2 + slip_y^2)
  loads     - q w  (q < 0 pushes the board away from the cathode)  - m (w_xx + w_yy)     (mesh-tension moment m)
Edges: w = 0; rotation clamped, simply supported, or an elastic rotational spring K_rot [N m / m / rad]; in-plane
immovable (bolted, friction) or free.

Discretisation: square grid of N x N nodes.  Bending by finite differences (w_xx, w_yy at nodes, ghost nodes for the
edge condition; w_xy at cell centres, which carries the edge twist); membrane and interface terms on bilinear Q4 cells
with 2 x 2 Gauss points.  Newton iteration with load stepping.  Validated in validate.py.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

NU = 0.3


def layer_props(layers, nu=NU):
    """layers: list of (E [Pa], t [m]) from top to bottom, bonded perfectly to each other.
    Returns dict: A (membrane stiffness E t / (1 - nu^2) summed), D about the composite neutral surface, the neutral
    surface depth below the layer group's top face, and total thickness."""
    z, top = [], 0.0
    for E, t in layers:
        z.append(top + t / 2); top += t
    Et = np.array([E * t for E, t in layers]); z = np.array(z)
    zn = float((Et * z).sum() / Et.sum())
    D = sum(E * (t ** 3 / 12 + t * (zi - zn) ** 2) for (E, t), zi in zip(layers, z)) / (1 - nu ** 2)
    return dict(A=float(Et.sum() / (1 - nu ** 2)), D=float(D), zn=zn, t=top)


class Plate:
    def __init__(self, a=0.402, N=61, layer1=None, layer2=None, k=0.0, edge='clamped', K_rot=None,
                 inplane='immovable', nu=NU, theta0=0.0):
        """layer1, layer2: dicts from layer_props (layer2 may be None).  k: scalar or callable k(x, y) [N/m^3]
        (x, y in m from the centre).  edge: 'clamped' | 'ss' | 'spring' (with K_rot)."""
        self.a, self.N, self.nu = a, N, nu
        self.h = h = a / (N - 1)
        self.L1, self.L2 = layer1, layer2
        self.edge, self.K_rot, self.inplane = edge, K_rot, inplane
        self.nl = 1 if layer2 is None else 2
        n = N * N
        self.n = n
        xs = np.linspace(-a / 2, a / 2, N)
        self.xs = xs
        X, Y = np.meshgrid(xs, xs, indexing='ij')       # node (i, j): x = xs[i], y = xs[j]; index i * N + j
        self.X, self.Y = X, Y
        idx = np.arange(n).reshape(N, N)
        self.idx = idx
        bnd = np.zeros((N, N), bool); bnd[0, :] = bnd[-1, :] = bnd[:, 0] = bnd[:, -1] = True
        self.bnd = bnd.ravel()

        # --- bending: curvature operators ---------------------------------------------------------------------
        Dsum = layer1['D'] + (layer2['D'] if layer2 else 0.0)
        self.Dsum = Dsum
        rows_xx, cols_xx, vals_xx = [], [], []
        rows_yy, cols_yy, vals_yy = [], [], []
        wt_node = []
        r = 0
        clamped = edge == 'clamped'
        for i in range(N):
            for j in range(N):
                onb_x = i in (0, N - 1); onb_y = j in (0, N - 1)
                if onb_x and onb_y:
                    continue                                    # corner: w = 0, no bending energy
                if (onb_x or onb_y) and not clamped:
                    continue                                    # SS / spring: w_xx = w_yy = 0 on the edge
                # w_xx
                def add(rr, cc, vv, R, C, V):
                    R.append(rr); C.append(cc); V.append(vv)
                if onb_x:                                       # clamped, edge normal to x: ghost w_-1 = w_1
                    ii = 1 if i == 0 else N - 2
                    add(r, idx[ii, j], 2 / h ** 2, rows_xx, cols_xx, vals_xx)
                    wn = 0.5 * h ** 2
                elif onb_y:                                     # along the edge: w = 0 -> w_xx = 0
                    wn = 0.5 * h ** 2
                else:
                    for di, c in ((-1, 1), (0, -2), (1, 1)):
                        add(r, idx[i + di, j], c / h ** 2, rows_xx, cols_xx, vals_xx)
                    wn = h ** 2
                if onb_y:
                    jj = 1 if j == 0 else N - 2
                    add(r, idx[i, jj], 2 / h ** 2, rows_yy, cols_yy, vals_yy)
                elif not onb_x:
                    for dj, c in ((-1, 1), (0, -2), (1, 1)):
                        add(r, idx[i, j + dj], c / h ** 2, rows_yy, cols_yy, vals_yy)
                wt_node.append(wn); r += 1
        nr = r
        Kxx = sp.csr_matrix((vals_xx, (rows_xx, cols_xx)), shape=(nr, n))
        Kyy = sp.csr_matrix((vals_yy, (rows_yy, cols_yy)), shape=(nr, n))
        Wn = sp.diags(np.array(wt_node))
        # twist at cell centres
        rr, cc, vv = [], [], []
        c = 0
        for i in range(N - 1):
            for j in range(N - 1):
                for (di, dj, s) in ((0, 0, 1), (1, 0, -1), (0, 1, -1), (1, 1, 1)):
                    rr.append(c); cc.append(idx[i + di, j + dj]); vv.append(s / h ** 2)
                c += 1
        Kxy = sp.csr_matrix((vv, (rr, cc)), shape=((N - 1) ** 2, n))
        Wc = sp.identity((N - 1) ** 2) * h ** 2
        self.Kb = Dsum * (Kxx.T @ Wn @ Kxx + Kyy.T @ Wn @ Kyy + nu * (Kxx.T @ Wn @ Kyy + Kyy.T @ Wn @ Kxx)
                          + 2 * (1 - nu) * Kxy.T @ Wc @ Kxy)
        # mesh-tension moment: Pi_m = -m * integral (w_xx + w_yy)
        self.fm_unit = np.asarray((Kxx.T @ Wn @ np.ones(nr)) + (Kyy.T @ Wn @ np.ones(nr))).ravel()
        # elastic rotational spring on the edges (edge = 'spring'), via a ghost node g per edge node eliminated
        # exactly: the edge node carries w_nn = (w1 + g)/h^2 (weight h^2/2) and the spring 1/2 K ((w1 - g)/2h)^2 h;
        # minimising over g leaves a diagonal stiffness 8AB/(A+B) on w1, A = D/(4h^2), B = K/(8h).  K = 0 -> SS,
        # K -> inf -> the clamped ghost (g = w1).
        Ks = sp.csr_matrix((n, n))
        if edge == 'spring':
            A_ = Dsum / (4 * h ** 2); B_ = K_rot / (8 * h)
            kd = 8 * A_ * B_ / (A_ + B_) if K_rot > 0 else 0.0
            diag = np.zeros(n)
            for t in range(1, N - 1):
                for p1 in (idx[1, t], idx[N - 2, t], idx[t, 1], idx[t, N - 2]):
                    diag[p1] += kd
            Ks = sp.diags(diag)
        self.Kb = (self.Kb + Ks).tocsr()
        # clamped edge with an imposed inward slope theta0 (rad; < 0 = the board dips into the aperture from the
        # clamp, i.e. a mounting rotation): ghost g = w1 - 2 h theta0 -> edge-node energy (2D/h^2)(w1 - h theta0)^2/...
        # i.e. a load 2 D theta0 / h on each first-interior node next to an edge.
        self.f_theta = np.zeros(n)
        if theta0 != 0.0 and edge == 'clamped':
            for t in range(1, N - 1):
                for p1 in (idx[1, t], idx[N - 2, t], idx[t, 1], idx[t, N - 2]):
                    self.f_theta[p1] += 2 * Dsum * theta0 / h

        # --- Q4 Gauss-point operators ----------------------------------------------------------------------------
        g = 1 / np.sqrt(3)
        gps = [(-g, -g), (g, -g), (-g, g), (g, g)]
        rN, cN, vN, rX, cX, vX, rY, cY, vY = [], [], [], [], [], [], [], [], []
        gx, gy = [], []
        q = 0
        for i in range(N - 1):
            for j in range(N - 1):
                corners = [(i, j, -1, -1), (i + 1, j, 1, -1), (i, j + 1, -1, 1), (i + 1, j + 1, 1, 1)]
                for (xi, eta) in gps:
                    for (ci, cj, sx, sy) in corners:
                        node = idx[ci, cj]
                        Nv = 0.25 * (1 + sx * xi) * (1 + sy * eta)
                        dNx = 0.25 * sx * (1 + sy * eta) * 2 / h
                        dNy = 0.25 * sy * (1 + sx * xi) * 2 / h
                        rN.append(q); cN.append(node); vN.append(Nv)
                        rX.append(q); cX.append(node); vX.append(dNx)
                        rY.append(q); cY.append(node); vY.append(dNy)
                    gx.append(xs[i] + (1 + xi) * h / 2); gy.append(xs[j] + (1 + eta) * h / 2)
                    q += 1
        ng = q
        self.ng = ng
        self.Ng = sp.csr_matrix((vN, (rN, cN)), shape=(ng, n))
        self.Gx = sp.csr_matrix((vX, (rX, cX)), shape=(ng, n))
        self.Gy = sp.csr_matrix((vY, (rY, cY)), shape=(ng, n))
        self.wg = np.full(ng, h * h / 4)
        self.gx, self.gy = np.array(gx), np.array(gy)
        # interface slip: one-point quadrature at cell centres (the 2 x 2 rule locks as k -> inf)
        rN, cN, vN, rX, cX, vX, rY, cY, vY = [], [], [], [], [], [], [], [], []
        cxs, cys = [], []
        q = 0
        for i in range(N - 1):
            for j in range(N - 1):
                for (ci, cj, sx, sy) in ((i, j, -1, -1), (i + 1, j, 1, -1), (i, j + 1, -1, 1), (i + 1, j + 1, 1, 1)):
                    node = idx[ci, cj]
                    rN.append(q); cN.append(node); vN.append(0.25)
                    rX.append(q); cX.append(node); vX.append(0.5 * sx / h)
                    rY.append(q); cY.append(node); vY.append(0.5 * sy / h)
                cxs.append(xs[i] + h / 2); cys.append(xs[j] + h / 2); q += 1
        nc = q
        self.Nc = sp.csr_matrix((vN, (rN, cN)), shape=(nc, n))
        self.Cx = sp.csr_matrix((vX, (rX, cX)), shape=(nc, n))
        self.Cy = sp.csr_matrix((vY, (rY, cY)), shape=(nc, n))
        self.cx, self.cy = np.array(cxs), np.array(cys)
        self.wc = np.full(nc, h * h)
        kk = k(self.cx, self.cy) if callable(k) else np.full(nc, float(k))
        self.kc = np.asarray(kk, float)
        # neutral surface of layer 1 to its lower face, plus layer 2's top face to its neutral surface
        self.d = ((layer1['t'] - layer1['zn']) + layer2['zn']) if layer2 else 0.0

        # dof layout: w, (u1, v1), (u2, v2)
        self.ndof = n * (1 + 2 * self.nl)
        fixed = np.zeros(self.ndof, bool)
        fixed[:n] = self.bnd
        if inplane == 'immovable':
            for b in range(2 * self.nl):
                fixed[n * (1 + b):n * (2 + b)] = self.bnd
        else:   # free in-plane: remove rigid-body modes by pinning the centre node's (u, v) of each layer (symmetry)
            ctr = idx[N // 2, N // 2]
            for b in range(2 * self.nl):
                fixed[n * (1 + b) + ctr] = True
        self.free = ~fixed
        self.Dlist = [layer1] + ([layer2] if layer2 else [])

    # ---------------------------------------------------------------------------------------------------------------
    def _split(self, x):
        n = self.n
        w = x[:n]
        uv = [(x[n * (1 + 2 * l):n * (2 + 2 * l)], x[n * (2 + 2 * l):n * (3 + 2 * l)]) for l in range(self.nl)]
        return w, uv

    def residual_tangent(self, x, f):
        n, nu = self.n, self.nu
        Gx, Gy, Ng, W = self.Gx, self.Gy, self.Ng, self.wg
        w, uv = self._split(x)
        wx, wy = Gx @ w, Gy @ w
        R = np.zeros(self.ndof)
        R[:n] = self.Kb @ w - f
        blocks = {}
        def addb(key, M):
            blocks[key] = blocks[key] + M if key in blocks else M
        addb((0, 0), self.Kb)
        Ntot = [np.zeros(self.ng), np.zeros(self.ng), np.zeros(self.ng)]
        dWx, dWy = sp.diags(wx) @ Gx, sp.diags(wy) @ Gy
        dWxy = sp.diags(wy) @ Gx + sp.diags(wx) @ Gy
        for l, (u, v) in enumerate(uv):
            A11 = self.Dlist[l]['A']; A12 = nu * A11; A66 = (1 - nu) / 2 * A11
            ex = Gx @ u + 0.5 * wx ** 2; ey = Gy @ v + 0.5 * wy ** 2; gxy = Gy @ u + Gx @ v + wx * wy
            Nx = A11 * ex + A12 * ey; Ny = A12 * ex + A11 * ey; Nxy = A66 * gxy
            for t, NN in zip(Ntot, (Nx, Ny, Nxy)):
                t += NN
            iu, iv = 1 + 2 * l, 2 + 2 * l
            R[n * iu:n * (iu + 1)] += Gx.T @ (W * Nx) + Gy.T @ (W * Nxy)
            R[n * iv:n * (iv + 1)] += Gy.T @ (W * Ny) + Gx.T @ (W * Nxy)
            R[:n] += Gx.T @ (W * (Nx * wx + Nxy * wy)) + Gy.T @ (W * (Ny * wy + Nxy * wx))
            # material tangent J^T C J with J rows: ex -> [dWx | Gx | 0], ey -> [dWy | 0 | Gy], gxy -> [dWxy | Gy | Gx]
            Jw = (dWx, dWy, dWxy); Ju = (Gx, None, Gy); Jv = (None, Gy, Gx)
            C = np.array([[A11, A12, 0], [A12, A11, 0], [0, 0, A66]])
            for a_, (Jwa, Jua, Jva) in enumerate(zip(Jw, Ju, Jv)):
                for b_, (Jwb, Jub, Jvb) in enumerate(zip(Jw, Ju, Jv)):
                    cab = C[a_, b_]
                    if cab == 0: continue
                    S = sp.diags(W * cab)
                    for (ra, Ja) in ((0, Jwa), (iu, Jua), (iv, Jva)):
                        if Ja is None: continue
                        for (rb, Jb) in ((0, Jwb), (iu, Jub), (iv, Jvb)):
                            if Jb is None: continue
                            addb((ra, rb), Ja.T @ S @ Jb)
        # geometric stiffness (ww block)
        Nx, Ny, Nxy = Ntot
        addb((0, 0), Gx.T @ sp.diags(W * Nx) @ Gx + Gy.T @ sp.diags(W * Ny) @ Gy
             + Gx.T @ sp.diags(W * Nxy) @ Gy + Gy.T @ sp.diags(W * Nxy) @ Gx)
        # interface slip
        if self.nl == 2 and np.any(self.kc > 0):
            (u1, v1), (u2, v2) = uv
            d, Nc, Cx, Cy = self.d, self.Nc, self.Cx, self.Cy
            sx = Nc @ (u1 - u2) + d * (Cx @ w)
            sy = Nc @ (v1 - v2) + d * (Cy @ w)
            Kw = self.wc * self.kc
            R[n:2 * n] += Nc.T @ (Kw * sx); R[3 * n:4 * n] -= Nc.T @ (Kw * sx)
            R[2 * n:3 * n] += Nc.T @ (Kw * sy); R[4 * n:5 * n] -= Nc.T @ (Kw * sy)
            R[:n] += d * (Cx.T @ (Kw * sx) + Cy.T @ (Kw * sy))
            S = sp.diags(Kw)
            Jx = {0: d * Cx, 1: Nc, 3: -Nc}
            Jy = {0: d * Cy, 2: Nc, 4: -Nc}
            for Jd in (Jx, Jy):
                for ra, Ja in Jd.items():
                    for rb, Jb in Jd.items():
                        addb((ra, rb), Ja.T @ S @ Jb)
        nb = 1 + 2 * self.nl
        Bm = [[blocks.get((i, j)) for j in range(nb)] for i in range(nb)]
        for i in range(nb):
            if Bm[i][i] is None:
                Bm[i][i] = sp.csr_matrix((n, n))
        K = sp.bmat(Bm, format='csr')
        return R, K

    def solve(self, q, m=0.0, steps=6, tol=1e-9, maxit=40, x0=None, linear=False, verbose=False):
        """q: uniform pressure [Pa], signed (negative = away from the cathode).  m: mesh-tension moment [N m/m].
        Returns the nodal deflection w (N x N, m)."""
        h = self.h
        wt = np.full((self.N, self.N), h * h)
        wt[0, :] *= 0.5; wt[-1, :] *= 0.5; wt[:, 0] *= 0.5; wt[:, -1] *= 0.5
        fq = wt.ravel()
        x = np.zeros(self.ndof) if x0 is None else x0.copy()
        fr = self.free
        if linear:
            f = q * fq + m * self.fm_unit + self.f_theta
            _, K = self.residual_tangent(np.zeros(self.ndof), f)
            x[fr] = spla.spsolve(K[fr][:, fr].tocsc(), np.r_[f, np.zeros(self.ndof - self.n)][fr])
            self.x = x
            return x[:self.n].reshape(self.N, self.N)
        for s in range(1, steps + 1):
            lam = s / steps
            f = lam * (q * fq + m * self.fm_unit + self.f_theta)
            for it in range(maxit):
                R, K = self.residual_tangent(x, f)
                Rf = R[fr]
                nrm = np.linalg.norm(Rf)
                if it == 0: nrm0 = max(nrm, 1e-30)
                dx = spla.spsolve(K[fr][:, fr].tocsc(), -Rf)
                x[fr] += dx
                if verbose: print(f'  step {s} it {it} |R| {nrm:.3e} |dx| {np.abs(dx).max():.3e}')
                if np.abs(dx[:self.n - self.bnd.sum()]).max() < tol * max(1e-6, np.abs(x[:self.n]).max()) \
                        and it > 0:
                    break
        self.x = x
        return x[:self.n].reshape(self.N, self.N)

    def sweep(self, qs, m=0.0, substeps=2):
        """continuation over a list of pressures (Pa, signed, monotonic): returns a list of nodal w arrays."""
        out, x, qprev = [], None, 0.0
        for qv in qs:
            ws = None
            for s in range(1, substeps + 1):
                qq = qprev + (qv - qprev) * s / substeps
                ws = self.solve(qq, m=m, steps=1, x0=x)
                x = self.x
            qprev = qv
            out.append(ws)
        return out

    # ---------------------------------------------------------------------------------------------------------------
    def interp(self, w, xq, yq):
        """bilinear interpolation of a nodal field at points (m, centre frame)."""
        from scipy.interpolate import RegularGridInterpolator
        return RegularGridInterpolator((self.xs, self.xs), w)(np.c_[np.ravel(xq), np.ravel(yq)]).reshape(np.shape(xq))


def fit_depth(plate, w, umax=2.0, wmax=1.8, ngrid=41):
    """Emulate the bench measurement: fit h = c + gu U + gw W + kuu U^2/2 + kww W^2/2 + kuw U W (U, W in units of
    100 mm from the centre) uniformly over |U| <= umax, |W| <= wmax (the efficient area: 399 x 360 mm) and return
    the fitted depth at 150 mm, 0.5 kbar 1.5^2 (FB summary.py), plus the true dip w(150 mm) - w(0) on axis.
    Heights are positive toward the cathode: a board pushed away from the cathode has w(0) < w(150) -> depth > 0."""
    U, W = np.meshgrid(np.linspace(-umax, umax, ngrid), np.linspace(-wmax, wmax, ngrid), indexing='ij')
    hh = plate.interp(w, U * 0.1, W * 0.1) * 1e3   # mm
    B = np.c_[np.ones(U.size), U.ravel(), W.ravel(), 0.5 * U.ravel() ** 2, 0.5 * W.ravel() ** 2, (U * W).ravel()]
    beta, *_ = np.linalg.lstsq(B, hh.ravel(), rcond=None)
    kuu, kww = beta[3], beta[4]
    kbar = 0.5 * (kuu + kww)
    c = plate.interp(w, np.array([0.0]), np.array([0.0]))[0] * 1e3
    e = plate.interp(w, np.array([0.150, 0.0]), np.array([0.0, 0.150])) * 1e3
    return dict(fit_depth_mm=float(0.5 * kbar * 1.5 ** 2), kuu=float(kuu), kww=float(kww),
                dip150_mm=float(np.mean(e) - c), centre_mm=float(-c))
