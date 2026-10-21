"""Strips with a resistivity gradient across their width: can uneven paste make the X asymmetry? (2026-10-22, Dylan's idea)
Idea: more paste on one side of each strip -> lower resistivity there; charge landing on the strip spreads faster toward
that side. A strip settles to one potential whatever its profile (the end state depends only on its shape and charge),
so only the transient can matter: the centroid of a deposit drifts toward higher conductivity (d<x>/dt = <dD/dx>) and
then relaxes to the strip centre.
1-D model across one strip (x, 5 um cells, insulating gaps = no-flux edges): dq/dt = d/dx (D(x) dq/dx), D = sigma(x)/c,
sigma(x) = sigma_mean (1 + g x / (w/2)), c = eps0 eps_r / 69 um (glue + kapton). At 5 MOhm/sq a 550 um strip equalises in
~67 ns; --fast multiplies sigma by 100. Gap landings stay put (or --funnel moves them to the nearer strip edge).
Signal on column k: S_k(t) = int p(x, t) psi_k(x, ESL) dx - psi_k(x0, z_ion(t)) (electron charge spreading on the ESL,
avalanche ion drifting straight up to the mesh in T ns), transparent-ESL layered weighting potentials over the 31.2 mm
superperiod (as ion_path_toy.py), the sim's DreamShaper, a1 = alpha+1 - alpha-1 (asym_lib). Every strip has the same
gradient direction (+x for g > 0).
Scores: A = board-mean a1 (a one-direction skew), O = <sgn(delta) a1> over |delta| >= 50 um (registration-correlated,
paired_odd.py convention), s1, and r(f) = (W+1 + W-1)/W_c at 25 % rise / peak / 50 % fall (time structure).
Usage: grad_strip_toy.py [--g 0,0.5,0.9] [--sig 150,300] [--ncen 12000] -> grad_strip_toy.json"""
import argparse, json, os, sys, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..')); sys.path.insert(0, os.path.join(HERE, '..', 'prompt'))
from response.dream.shaper import DreamShaper

ap = argparse.ArgumentParser()
ap.add_argument('--g', default='0,0.5,0.9'); ap.add_argument('--sig', default='150,300')
ap.add_argument('--ncen', type=int, default=12000); ap.add_argument('--T', type=float, default=200.)
ap.add_argument('--zb', type=float, default=10.); ap.add_argument('--w', type=float, default=550.)
ap.add_argument('--out', default='grad_strip_toy.json')
a = ap.parse_args()

P, SP, BOX, PADW = 780., 800., 31200., 680.
A_D, B_M = 69., 219.; EPS1 = (3.2 * 19. + 3.5 * 50.) / 69.
DX = 5.; NX = int(round(BOX / DX)); x = (np.arange(NX) + .5) * DX
ORG = -(512 - 1) / 2 * P; padc = np.mod(ORG + np.arange(40) * P, BOX); NPAD = 40
delta = ((padc + 400.) % 800.) - 400.

def pad_drive(xc):
    d = ((x - xc + BOX / 2) % BOX) - BOX / 2; ad = np.abs(d)
    return np.where(ad <= PADW / 2, 1., np.clip((P / 2 + (P - PADW) / 2 - ad) / (P - PADW), 0, 1))

def layered(V, zs):
    k = 2 * np.pi * np.fft.rfftfreq(NX, DX); Vh = np.fft.rfft(V); out = np.empty((len(zs), NX))
    ka, kc = k * A_D, k * (B_M - A_D)
    with np.errstate(all='ignore'):
        ch, sh, s2, c2 = np.cosh(ka), np.sinh(ka), np.sinh(kc), np.cosh(kc)
        B = -Vh * (EPS1 * sh * s2 + ch * c2) / (EPS1 * ch * s2 + sh * c2); C = (Vh * ch + B * sh) / s2
    e = A_D / EPS1 + (B_M - A_D)
    for i, z in enumerate(zs):
        f = Vh * np.cosh(k * z) + B * np.sinh(k * z) if z <= A_D else C * np.sinh(k * (B_M - z))
        f[0] = Vh[0] * (1 - (z / EPS1) / e) if z <= A_D else Vh[0] * ((B_M - z) / e)
        out[i] = np.fft.irfft(np.where(np.isfinite(f), f, 0), NX)
    return out

Z0 = A_D + a.zb; zs = np.concatenate([[A_D], np.linspace(Z0, B_M, 120)])
PSI = np.stack([layered(pad_drive(padc[c]), zs) for c in range(NPAD)])      # (40, 121, NX); zs[0] = the ESL plane
PSI_E, PSI_I = PSI[:, 0], PSI[:, 1:]
NT = 900; t = np.arange(NT) * 1.0
iz = np.clip(np.searchsorted(zs[1:], Z0 + (B_M - Z0) * np.clip(t / a.T, 0, 1)), 0, 119)

# strips: centres at m*SP in the box frame (39 per superperiod), cells of the strip-local grid
half = a.w / 2; NS = int(round(BOX / SP)); NC = int(round(a.w / DX)); xc_loc = (np.arange(NC) + .5) * DX - half
C_F = 8.854e-12 * EPS1 / (A_D * 1e-6)                                      # F/m^2
def spread_ops(g, rho):
    """p(x, t) for a unit deposit in each strip cell: (NC from, NC to, NT), dq/dt = d/dx (D dq/dx), D = sigma/c."""
    sig_face = (1 + g * (np.arange(1, NC) * DX - half) / half) / rho            # interior faces, S
    Dface = sig_face / C_F * 1e12 / 1e9                                         # m^2/s -> um^2/ns
    L = np.zeros((NC, NC)); i = np.arange(NC - 1); w_ = Dface / DX ** 2
    L[i, i] -= w_; L[i + 1, i + 1] -= w_; L[i, i + 1] += w_; L[i + 1, i] += w_
    lam, V = np.linalg.eigh(L)
    E = np.exp(np.clip(lam[:, None] * t[None, :], -700, 0))                     # (NC, NT)
    return np.einsum('am,mt,bm->abt', V, E, V)                                   # K[a -> b, t]

sh = DreamShaper(t_peak_ns=166.)
RNG = np.random.default_rng(5); res = {}
for fast in (False, True):
    rho = 5e6 / (100. if fast else 1.)
    for g in map(float, a.g.split(',')):
        t0 = time.time(); K = spread_ops(g, rho)
        # per strip m: electron-part signal on all 40 columns for a deposit in each cell: (NS, NC, 40, NT)
        Gs = np.empty((NS, NC, NPAD, NT), np.float32)
        for m in range(NS):
            jx = np.clip(((m * SP + xc_loc) % BOX / DX).astype(int), 0, NX - 1)  # strip cells on the box grid
            Gs[m] = np.einsum('abt,kb->akt', K, PSI_E[:, jx])
        for funnel in (0, 1):
            for sig in map(float, a.sig.split(',')):
                rows = []
                for c in range(NPAD):
                    for x0 in padc[c] + RNG.uniform(-P / 2, P / 2, a.ncen // NPAD):
                        xo = np.mod(x0 + sig * RNG.standard_normal(64), BOX)
                        ph = ((xo + half) % SP) - half                            # position relative to the nearest strip's
                        m = np.rint((xo - ph) / SP).astype(int) % NS              # centre and that strip's index
                        on = np.abs(ph) < half
                        if funnel:
                            # gap landing -> nearer strip edge (inside by half a cell)
                            # (ph is measured from the NEAREST strip centre, so the nearer edge is always strip m's)
                            ph = np.where(on, ph, np.where(ph >= 0, half - DX / 2, -half + DX / 2))
                            on = np.ones_like(on)
                        S = np.zeros((NPAD, NT))
                        jx0 = np.clip((xo / DX).astype(int), 0, NX - 1)
                        for i in range(64):
                            if on[i]:
                                ci = int(np.clip((ph[i] + half) // DX, 0, NC - 1)); S += Gs[m[i], ci]
                                jx0[i] = int(((m[i] * SP + xc_loc[ci]) % BOX) / DX)
                            else:
                                S += PSI_E[:, jx0[i]][:, None]                     # static charge in the gap
                            S -= PSI_I[:, iz, jx0[i]]                             # the ion, straight up from the landing
                        I = np.diff(S, axis=1, prepend=0.) / 1e-9
                        W = np.stack([sh.apply(I[k]) for k in range(NPAD)])
                        pk = int(W.max(1).argmax()); Wc = W[pk]; n2 = Wc @ Wc
                        a1 = (W[(pk + 1) % NPAD] - W[pk - 1]) @ Wc / n2; s1 = (W[(pk + 1) % NPAD] + W[pk - 1]) @ Wc / n2
                        tp = int(Wc.argmax()); nb = W[(pk + 1) % NPAD] + W[pk - 1]
                        r25 = np.flatnonzero(Wc[:tp] < .25 * Wc[tp]); f50 = tp + np.flatnonzero(Wc[tp:] < .5 * Wc[tp])
                        rows.append((pk, a1, s1, nb[r25[-1]] if len(r25) else np.nan, Wc[r25[-1]] if len(r25) else np.nan,
                                     nb[tp], Wc[tp], nb[f50[0]] if len(f50) else np.nan, Wc[f50[0]] if len(f50) else np.nan))
                R = np.array(rows); dl = delta[R[:, 0].astype(int)]; mk = np.abs(dl) >= 50
                rr = lambda i: float(np.nansum(R[:, i]) / np.nansum(np.where(np.isfinite(R[:, i]), R[:, i + 1], np.nan)))
                key = f"{'fast' if fast else '5M'} g{g:g} f{funnel} s{sig:g}"
                res[key] = dict(A=float(R[:, 1].mean()), A_err=float(R[:, 1].std() / np.sqrt(len(R))),
                                O=float(np.mean(np.sign(dl[mk]) * R[mk, 1])), s1=float(R[:, 2].mean()),
                                r25=rr(3), peak=rr(5), f50=rr(7), n=len(R))
                d = res[key]
                print(f"{key:22s} A {d['A']:+.4f}±{d['A_err']:.4f}  O {d['O']:+.4f}  s1 {d['s1']:.3f}  "
                      f"r(25%/peak/fall50) {d['r25']:.3f}/{d['peak']:.3f}/{d['f50']:.3f}   ({time.time() - t0:.0f} s)", flush=True)
res['args'] = vars(a)
json.dump(res, open(os.path.join(HERE, a.out), 'w'), indent=1)
