#!/usr/bin/env python3
"""
compare_kernels3d.py -- the chain's S1/W2 kernels G against a direct COMSOL solve.

Until 2026-10-12 the S1 kernels had only ever been checked against the chain's own
solvers (sum rules, a real-space FD lift, V6). Induce3D.java solves the same board
directly in 3D -- a point avalanche on the centre of an ESL strip and of a Y pad, the
real 550/800 um strips, the real checkerboard comb channels -- with no weighting-field
theorem. Its `deposit` case (electrons only, landed by ~1 ns) is the channel kernel G
for a charge on the ESL, i.e. exactly what the product stores.

Matching point in the chain's coordinates (solver/kernels.py conventions):
  strip centres at x = 0 mod 800 um; pad column c at kernels.pad_x(c);
  row 0 at y = 0 in both the Y and the X kernel frames;
  Y row j     -> G_Y[j mod 2](t, y = 780 j um, x0)
  X column i  -> G_X[c0 + i](t, y = 0, x0)
with c0 an even column (its row-0 pad is a Y pad) whose centre is within 10 um of a
strip centre (780/800 um never align exactly).

Prints G(t) per channel at a few times, then the DREAM-shaped responses (peak ratio,
peak time, and the neighbour peak delay, which the chain gets ~2x too long vs data).

Usage: python3 -m response.comsol.compare_kernels3d <induce3d_deposit.csv> <greens_comb_w2_*.npz>
"""
import sys
import numpy as np

from ..common import constants as C
from ..solver import kernels as K
from ..dream.shaper import DreamShaper


def load_comsol(path):
    with open(path) as fh:
        head = fh.readline().strip()
        cols = fh.readline().strip().split(",")[1:]
    d = np.loadtxt(path, delimiter=",", skiprows=2)
    t, I = d[:, 0], d[:, 1:]
    q = np.vstack([np.zeros(I.shape[1]), np.cumsum(0.5 * (I[1:] + I[:-1]) * np.diff(t)[:, None], 0)])
    return head, cols, t, q / 1e-15        # induced charge per unit avalanche charge


def chain_G(npz, cols):
    d = np.load(npz)
    t = d["t"]; x = d["x"]; yY = d["y_Y"]; yX = d["y_X"]
    # an even column within 10 um of a strip centre
    cand = np.arange(0, C.N_PAD_PER_SUPER, 2)
    px = K.pad_x(cand)
    ph = (K.esl_phase(px) + C.ESL_PITCH_M / 2) % C.ESL_PITCH_M - C.ESL_PITCH_M / 2
    c0 = int(cand[np.argmin(np.abs(ph))])
    x0 = float(K.pad_x(c0))
    ix = int(np.argmin(np.abs(x - x0)))
    print(f"# chain source: column {c0} at x0 = {x0*1e3:.3f} mm, strip phase "
          f"{float(ph[np.argmin(np.abs(ph))])*1e6:+.0f} um, x sample {x[ix]*1e3:.3f} mm")
    out = {}
    xcol = d["x_cols"]
    for c in cols:
        n = int(c[1:])
        if c[0] == "Y":
            G = d["G_Y_even"] if n % 2 == 0 else d["G_Y_odd"]
            iy = int(np.argmin(np.abs(yY - n * C.PAD_PITCH_M)))
            out[c] = np.asarray(G[:, iy, ix], float)
        else:
            # the column n pads to the +x side of c0; G_X is indexed by column phase
            want = K.pad_x((c0 + n) % C.N_PAD_PER_SUPER)
            kc = int(np.argmin(np.abs(((xcol - want) + C.SUPERPERIOD_M / 2) % C.SUPERPERIOD_M - C.SUPERPERIOD_M / 2)))
            iy = int(np.argmin(np.abs(yX - 0.0)))
            out[c] = np.asarray(d["G_X"][kc, :, iy, ix], float)
    return t, out


def shaped_from_charge(tq, q, tg, sh, dt):
    """charge on a non-uniform grid -> bin-averaged current on tg -> DREAM shaper"""
    edges = np.append(tg, tg[-1] + dt) * 1e-9
    qi = np.interp(edges, tq, q)
    qi[0] = 0.0          # the chain's G(0) is the prompt step: it arrives in the first bin
    cur = np.diff(qi) / (dt * 1e-9)
    return sh.apply(cur, dt_ns=dt)


def main():
    head, cols, tc, qc = load_comsol(sys.argv[1])
    print("#", head)
    tk, Gk = chain_G(sys.argv[2], cols)
    sgn = np.sign(qc[np.argmin(np.abs(tc - 20e-9)), cols.index("Y0")])
    qc *= sgn
    times = [2e-9, 10e-9, 30e-9, 100e-9, 300e-9, 1e-6]
    times = [x for x in times if x <= tc[-1]]
    print("\nG(t): induced charge per unit charge on the ESL   COMSOL / chain")
    print("chan " + "".join(f"{x*1e9:>16.0f}ns" for x in times))
    for k, c in enumerate(cols):
        row = []
        for x in times:
            a = np.interp(x, tc, qc[:, k]); b = np.interp(x, tk, Gk[c])
            row.append(f"{a:+.4f}/{b:+.4f}")
        print(f"{c:4s} " + "".join(f"{r:>18s}" for r in row))

    dt = 1.0
    tg = np.arange(0.0, tc[-1] * 1e9, dt)
    sh = DreamShaper(dt_ns=dt)
    print("\nDREAM-shaped (180 ns): peak COMSOL/chain, peak time COMSOL / chain")
    res = {}
    for k, c in enumerate(cols):
        sc = shaped_from_charge(tc, qc[:, k], tg, sh, dt)
        sk = shaped_from_charge(tk, Gk[c], tg, sh, dt)
        ic, ik = np.argmax(np.abs(sc)), np.argmax(np.abs(sk))
        res[c] = (sc[ic], tg[ic], sk[ik], tg[ik])
        print(f"{c:4s} peak {sc[ic]:+.4f} / {sk[ik]:+.4f}  ratio {sc[ic]/sk[ik]:6.3f}   "
              f"t {tg[ic]:5.0f} / {tg[ik]:5.0f} ns")
    for nb in ("X0", "X1", "Y1"):
        if nb in res:
            print(f"neighbour delay {nb} - Y0:  COMSOL {res[nb][1]-res['Y0'][1]:+.0f} ns   "
                  f"chain {res[nb][3]-res['Y0'][3]:+.0f} ns")


if __name__ == "__main__":
    main()
