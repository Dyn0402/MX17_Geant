"""Ion return paths under a saturated charged-gap funnel: can they skew the X response? (2026-10-21 evening)
The skewed-induction family (FOLLOWUP 'What this leaves'): its resistive-strip member is already in every full-sim arm
(the W2 X box is the 31.2 mm superperiod, 800 um strips over 780 um pads at every offset). What the sim cannot make:
ions drift straight up from the landing point, and every channel gets the same ion time profile. With a saturated
funnel an electron from above a gap lands on the strip edge, and its avalanche ions drift back along the same field
line, ending at the mesh above where the electron came from. The integrated charge is path-independent (mesh at
psi = 0); the neighbours' pulse SHAPES are not, and a1 is an in-phase projection. Is that a sum-conserving
n-1 -> n+1 transfer of the data's size?

2-D (x across the strips, z up), exact layered Fourier weighting potentials over the 31.2 mm superperiod: pads 680/780
at z = 0 (linear ramp across the pad gaps), glue 19 um (eps 3.2) + kapton 50 um (3.5), ESL at z = 69 treated as
transparent (prompt limit; the resistive late-time part is the W2 kernel's and is the same in all cases), 150 um gas
gap, mesh at psi = 0. X pads as full columns (no Y-pad comb): overstates X coupling by <~2.
Ions: born ZB above the ESL at the landing x, constant speed, drift time T to the mesh; lateral path
x(z) = x_land + L (1 - exp(-(z - z0)/LAM)) / (1 - exp(-(b - z0)/LAM)), L = x_orig - x_land (field-line return).
Pulses: the sim's DreamShaper (166 ns); a1 = alpha_+1 - alpha_-1, alpha_k = <W_k, W_c> / <W_c, W_c> as asym_lib.
Clouds: Gaussian sigma (no resistive spreading), centres uniform over the board, peak column = max shaped amplitude.
Cases: none (no funnel), sim (funnel, ions straight up from the edge: what the digitizer does), path (funnel, ions
back along their field line). Score O = <sgn(delta) a1> over |delta| >= 50 um (paired_odd.py's definition and
delta convention), all clouds, plus s1 = alpha_+1 + alpha_-1.
--esl float: the strips as floating conductors during the ion drift (float_psi.py must have been run for that --w/--zb).
Usage: ion_path_toy.py [--T 200] [--lam 60] [--zb 10] [--w 550] [--sig 150,250,350] [--esl transparent|float] [--out f.json]"""
import argparse, json, os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from response.dream.shaper import DreamShaper

ap = argparse.ArgumentParser()
ap.add_argument('--T', type=float, default=200.); ap.add_argument('--lam', type=float, default=60.)
ap.add_argument('--zb', type=float, default=10.); ap.add_argument('--w', type=float, default=550.)
ap.add_argument('--sig', default='150,250,350'); ap.add_argument('--out', default='ion_path_toy.json')
ap.add_argument('--ncen', type=int, default=4000)
ap.add_argument('--esl', default='transparent', choices=['transparent', 'float'])
a = ap.parse_args()

P, SP, BOX = 780., 800., 31200.                     # um
NPAD = int(round(BOX / P)); PADW = 680.
A_D, B_M = 19. + 50., 69. + 150.                    # ESL height, mesh height
DX = 5.; NX = int(round(BOX / DX)); x = (np.arange(NX) + .5) * DX
ORG = -(512 - 1) / 2 * P                             # board frame: pad c centre ORG + c P, strip centres at k SP
padc = np.mod(ORG + np.arange(NPAD) * P, BOX)        # 40 pad centres in the box (columns c = 0..39 mod 40)

def pad_drive(xc):
    d = ((x - xc + BOX / 2) % BOX) - BOX / 2; ad = np.abs(d)
    return np.clip((P / 2 + (P - PADW) / 2 - ad) / (P - PADW), 0, 1) * (ad > PADW / 2) + (ad <= PADW / 2)

def layered(V, zs):
    """psi(x, z) for pad-plane potential V(x), mesh at B_M = 0; dielectric slab eps_d (0..A_D), gas above."""
    eps1 = (3.2 * 19. + 3.5 * 50.) / 69.             # thickness-weighted (series would be ~3.43; same to 1 %)
    k = 2 * np.pi * np.fft.rfftfreq(NX, DX); Vh = np.fft.rfft(V); out = np.empty((len(zs), NX))
    ka, kc = k * A_D, k * (B_M - A_D)
    with np.errstate(all='ignore'):
        # layer1: A cosh(kz) + B sinh(kz), A = Vh; layer2: C sinh(k(b-z))
        # continuity at a: A ch + B sh = C s2 ; eps1 (A sh + B ch) = -C c2
        ch, sh, s2, c2 = np.cosh(ka), np.sinh(ka), np.sinh(kc), np.cosh(kc)
        B = -Vh * (eps1 * sh * s2 + ch * c2) / (eps1 * ch * s2 + sh * c2); C = (Vh * ch + B * sh) / s2
    for i, z in enumerate(zs):
        if z <= A_D: f = Vh * np.cosh(k * z) + B * np.sinh(k * z)
        else:        f = C * np.sinh(k * (B_M - z))
        # k = 0: series capacitor, linear in each layer
        e = A_D / eps1 + (B_M - A_D); v0 = Vh[0]
        f[0] = v0 * (1 - (z / eps1) / e) if z <= A_D else v0 * ((B_M - z) / e)
        out[i] = np.fft.irfft(np.where(np.isfinite(f), f, 0), NX)
    return out

Z0 = A_D + a.zb; NZ = 120; zs = np.linspace(Z0, B_M, NZ)            # ion heights along the path
NB = 7; ds = np.arange(-3, 4)                                         # channels c-3..c+3
PSI = np.stack([layered(pad_drive(padc[c]), zs) for c in range(NPAD)])   # (40, NZ, NX)
if a.esl == 'float':                                                  # strips as floating conductors (float_psi.py)
    zf = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'float_psi_w{a.w:g}.npz'))
    assert np.allclose(zf['zs'], zs) and len(zf['x']) == NX, 'float_psi grid mismatch: rerun float_psi.py with the same --zb/--w'
    PSI = zf['PSI'].astype(float)
print('weighting potentials done; psi_c at the ESL over its own pad centre', float(PSI[0, 0, int(padc[0] / DX)]))

# landing (funnel) map, strips 'w' wide centred at k*SP (board frame), board x = box x here (box frame = board mod 31.2)
half = a.w / 2; ph = np.mod(x, SP)
gap = (ph >= half) & (ph < SP - half)
xland_f = np.where(gap, np.where(ph < SP / 2, x - ph + half, x - ph + SP - half), x)

sh = DreamShaper(t_peak_ns=166.); NT = 900; t = np.arange(NT) * 1.0
s = np.clip(t / a.T, 0, 1)                                               # path fraction vs time
zt = Z0 + (B_M - Z0) * s
iz = np.clip(np.searchsorted(zs, zt), 0, NZ - 1)
prof = (1 - np.exp(-(zt - Z0) / a.lam)) / (1 - np.exp(-(B_M - Z0) / a.lam))

def waves(case, xi):
    """Induced charge Q_c(t) on all 40 columns for a unit avalanche from electron origin x[xi] (box grid index)."""
    xl = x[xi] if case == 'none' else xland_f[xi]
    L = 0. if case != 'path' else (x[xi] - xl)
    xt = np.mod(xl + L * prof, BOX); jx = np.clip(np.rint(xt / DX - .5).astype(int), 0, NX - 1)
    psi_t = PSI[:, iz, jx]                                               # (40, NT) psi along the path
    return psi_t[:, :1] - psi_t                                          # Q(t) = psi(start) - psi(now) (ion leaving)

delta = ((padc + 400.) % 800.) - 400.                                    # sim_delta.py convention, per column
rng = np.random.default_rng(1)
res = {}
for case in ('none', 'sim', 'path'):
    # per-origin Q waveforms on a 5 um grid is 6240 x 40 x 900 = too big: build clouds by sampling origins instead
    for sig in map(float, a.sig.split(',')):
        out = []
        for c in range(NPAD):
            cen = padc[c] + rng.uniform(-P / 2, P / 2, a.ncen // NPAD)
            for x0 in cen:
                xo = np.mod(x0 + sig * rng.standard_normal(64), BOX); xi = np.clip((xo / DX).astype(int), 0, NX - 1)
                Q = sum(waves(case, i) for i in xi)                       # (40, NT)
                I = np.diff(Q, axis=1, prepend=0.) / 1e-9
                W = np.stack([sh.apply(I[k]) for k in range(NPAD)])
                amp = W.max(1); pk = int(amp.argmax())
                Wc = W[pk]; n2 = Wc @ Wc
                al = {d: W[(pk + d) % NPAD] @ Wc / n2 for d in (-2, -1, 1, 2)}
                # +x order: columns increase with x in the box frame (padc increasing in c within a period)
                out.append((pk, al[1] - al[-1], al[1] + al[-1], al[2] - al[-2]))
        out = np.array(out); dl = delta[out[:, 0].astype(int)]; m = np.abs(dl) >= 50
        O = float(np.mean(np.sign(dl[m]) * out[m, 1])); s1 = float(out[:, 2].mean())
        O2 = float(np.mean(np.sign(dl[m]) * out[m, 3]))
        # a1 vs delta in 100 um bins for the record
        bins = np.arange(-400, 401, 100); prof_d = [float(out[(dl >= lo) & (dl < lo + 100), 1].mean()) for lo in bins[:-1]]
        res.setdefault(case, {})[f'{sig:g}'] = dict(O=O, O2=O2, s1=s1, n=int(len(out)), a1_vs_delta=prof_d)
        print(f'{case:5s} sigma {sig:4.0f}  O {O:+.4f}  O(+-2) {O2:+.4f}  s1 {s1:.3f}  a1 vs delta {np.round(prof_d, 3)}', flush=True)
res['args'] = vars(a)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), a.out), 'w'), indent=1)
