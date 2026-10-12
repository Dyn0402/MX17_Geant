#!/usr/bin/env python3
"""
compare_induce2d.py -- does the chain's ion operator survive a direct COMSOL solve?

Induce2D.java gives two COMSOL answers in the same 2D stack:

  deposit  electrons born at the mesh, no ions. Their current on a pad is the
           chain's surface kernel G(t): charge landing on the ESL, seen through
           the resistive layer.
  pair     electron/ion pairs born on exp(-z/Lz) at t = 0, ions drifting to the
           mesh over T_ion. This is the truth, with no weighting-field theorem.

The chain (kernel_lut.apply_ion_transit) builds the pair signal from G alone:

    I_chain = (1 - f) G  +  f (rect_T (x) G),      f = <1 - z/g> = f_ion

i.e. the ion is coupled through the same surface kernel as a landed charge.
If the resistive layer screens the in-gas ion differently, I_pair and I_chain
part company, and the fitted f_eff on the central pad says by how much.

Usage:  python3 -m response.comsol.compare_induce2d deposit.csv pair.csv [T_ion_ns]

Works for both Induce2D and InduceAxi output (the central terminal is "p0").
"""
import sys
import numpy as np

from ..dream.shaper import DreamShaper

GAP_UM, LZ_UM = 150.0, 14.8


def load(path):
    with open(path) as fh:
        head = fh.readline().strip()
        cols = fh.readline().strip().split(",")
    d = np.loadtxt(path, delimiter=",", skiprows=2)
    return head, cols, d[:, 0], d[:, 1:]


def resample(t_s, cur, dt_ns, t_max_ns):
    """
    Bin-average the current onto a uniform grid by differencing the cumulative
    charge. Interpolating the current itself loses the ~1 ns prompt spike of the
    pair case (COMSOL stores 0.1 ns steps only to 2 ns) and with it a few % of
    the charge.
    """
    tg = np.arange(0.0, t_max_ns, dt_ns)
    edges = np.append(tg, tg[-1] + dt_ns) * 1e-9
    out = []
    for c in cur.T:
        qc = np.concatenate([[0.0], np.cumsum(0.5 * (c[1:] + c[:-1]) * np.diff(t_s))])
        out.append(np.diff(np.interp(edges, t_s, qc)) / (dt_ns * 1e-9))
    return tg, np.stack(out)


def f_ion_profile():
    """<1 - z/g> over the truncated exp(-z/Lz) birth profile (the chain's f_ion in 2D)."""
    z = np.linspace(0, GAP_UM, 20001)
    w = np.exp(-z / LZ_UM)
    return float(np.sum(w * (1 - z / GAP_UM)) / np.sum(w))


def main():
    dep_p, pair_p = sys.argv[1], sys.argv[2]
    t_ion = float(sys.argv[3]) if len(sys.argv) > 3 else 340.0
    dt = 1.0
    hd, cols, td, Id = load(dep_p)
    hp, colsp, tp, Ip = load(pair_p)
    assert cols == colsp, (cols, colsp)
    t_max = min(td[-1], tp[-1]) * 1e9
    tg, G = resample(td, Id, dt, t_max)
    _, P = resample(tp, Ip, dt, t_max)
    q = 1e-12                                   # C/m, Induce2D's Q
    # Normalise both to the collected charge on the central pad at the end,
    # sign-fixed so the central pad's net charge is positive.
    ic = cols.index("p0") - 1
    # bins are averages (resample), so the charge is a plain sum; trapz would
    # halve the first bin, which holds the whole prompt spike
    QG = G[ic].sum() * dt * 1e-9
    QP = P[ic].sum() * dt * 1e-9
    print(f"# {hd}\n# {hp}")
    print(f"net charge on p0 / Q:  deposit {QG / q:+.4f}   pair {QP / q:+.4f}")
    for k, name in enumerate(cols[1:]):
        print(f"  {name:4s} deposit {G[k].sum()*dt*1e-9/q:+.4f}   pair {P[k].sum()*dt*1e-9/q:+.4f}")

    # The ion current of the whole cloud in the parallel-plate weighting field:
    # an ion born at z drifts (g - z) at v = g/T_ion and induces v/g while it
    # moves, so  i(t) = (v/g) P(z < g - v t),  whose integral is <1 - z/g> = f_ion.
    # Normalised here to unit area; the chain couples it through G.
    v = GAP_UM / t_ion
    zc = np.clip(GAP_UM - v * (tg + 0.5 * dt), 0, GAP_UM)
    cdf = (1 - np.exp(-zc / LZ_UM)) / (1 - np.exp(-GAP_UM / LZ_UM))
    ion = cdf / cdf.sum()
    conv = np.stack([np.convolve(g, ion)[: len(tg)] for g in G])
    f0 = f_ion_profile()
    sh = DreamShaper(dt_ns=dt)

    def shaped(x):
        return np.stack([sh.apply(r, dt_ns=dt) for r in x])

    SP, SG, SC = shaped(P), shaped(G), shaped(conv)
    # least-squares f_eff on p0: SP ~ a*((1-f) SG + f SC)  =  A SG + B SC
    M = np.stack([SG[ic], SC[ic]], 1)
    (A, B), *_ = np.linalg.lstsq(M, SP[ic], rcond=None)
    f_eff = B / (A + B)
    chain = (1 - f0) * SG + f0 * SC
    sc = SP[ic].max() / chain[ic].max()
    # deposit and pair carry the same Q, so G is already per unit avalanche charge and
    # the chain prediction is absolute: no rescaling (a final-charge rescale is
    # meaningless once the sheet drains the pad, as in the axisymmetric model)
    print(f"\nchain p0 shaped peak, absolute (no free scale): pair/chain = {SP[ic].max() / chain[ic].max():.3f}")
    print(f"\nf_ion (profile <1-z/g>) = {f0:.3f}")
    print(f"fitted f_eff on p0      = {f_eff:.3f}   (data want 0.17-0.31; chain uses ~0.90)")
    for k, name in enumerate(cols[1:]):
        def rise(y):
            m = y.max(); i10 = np.argmax(y >= 0.1 * m); i90 = np.argmax(y >= 0.9 * m)
            return (i90 - i10) * dt, np.argmax(y) * dt
        rp, pp = rise(SP[k]) if SP[k].max() > 0 else (np.nan, np.nan)
        rc, pc = rise(chain[k]) if chain[k].max() > 0 else (np.nan, np.nan)
        print(f"  {name:4s} shaped peak pair/chain = {SP[k].max() / (sc * chain[k].max()):.3f}   "
              f"rise pair {rp:.0f} ns chain {rc:.0f} ns   peak t pair {pp:.0f} chain {pc:.0f} ns")
    np.savez(pair_p.replace(".csv", "_cmp.npz"), t_ns=tg, cols=np.array(cols[1:]),
             pair=P, deposit=G, conv=conv, shaped_pair=SP, shaped_chain=chain, f_eff=f_eff, f0=f0)


if __name__ == "__main__":
    main()
