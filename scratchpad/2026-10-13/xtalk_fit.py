# Can prompt channel-to-channel crosstalk explain why data neighbours (+-1) peak earlier than
# the chain's? Event-aligned averages (central strip's leading 50 % crossing, per-event
# normalised to the central peak, railed events dropped) of the central strip and the +-1
# neighbours, data and chain, vertical tracks (waveforms_pm3.npz). Then fit
#     data_nb = a * chain_nb + k * data_c + k' * d(data_c)/dt
# and report a, k, k' and the neighbour peak delay before/after.
# Run with the nTof_x17 .venv in the dir holding the npz (FRESH_LOOK §10.3).
import sys, numpy as np
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_shape import cross
DT = 60.0
grid = np.arange(-4, 20, 0.25)

def aligned(W):
    keep = W[:, 3].max(1) < 3500
    W = W[keep]
    C, N = [], []
    for w in W:
        l = cross(w[3], .5, -1)
        if np.isnan(l): continue
        pk = w[3].max()
        x = np.arange(w.shape[1])
        C.append(np.interp(grid + l, x, w[3] / pk))
        N.append(0.5 * (np.interp(grid + l, x, w[2] / pk) + np.interp(grid + l, x, w[4] / pk)))
    return np.mean(C, 0), np.mean(N, 0)

z = np.load('waveforms_pm3.npz')
for view in ('x', 'y'):
    dc, dn = aligned(z[f'data_{view}_wf'])
    sc, sn = aligned(z[f'sim_{view}_wf'])
    ddc = np.gradient(dc, grid)
    tpk = lambda y: grid[np.argmax(y)] * DT
    M = np.stack([sn, dc, ddc], 1)
    (a, k, kp), *_ = np.linalg.lstsq(M, dn, rcond=None)
    fit = M @ [a, k, kp]
    M0 = sn[:, None]
    (a0,), *_ = np.linalg.lstsq(M0, dn, rcond=None)
    rms = lambda r: float(np.sqrt(np.mean(r ** 2)))
    print(f"view {view}: central peak t data {tpk(dc):.0f} / chain {tpk(sc):.0f} ns (from the leading 50%)")
    print(f"   +-1 neighbour peak t: data {tpk(dn):.0f}, chain {tpk(sn):.0f}  -> delay vs central "
          f"data {tpk(dn)-tpk(dc):+.0f}, chain {tpk(sn)-tpk(sc):+.0f} ns")
    print(f"   neighbour/central peak: data {dn.max():.3f}, chain {sn.max():.3f}")
    print(f"   fit data_nb = a chain_nb: a = {a0:.3f}, rms {rms(dn - a0*sn):.4f}")
    print(f"   fit + k data_c + k' d/dt data_c: a = {a:.3f}, k = {k:+.4f}, k' = {kp:+.4f} [/sample], "
          f"rms {rms(dn - fit):.4f}; model nb peak t {tpk(fit):.0f} ns")
    print("   t[ns]     " + " ".join(f"{g*DT:5.0f}" for g in grid[::8]))
    for name, y in (("data c", dc), ("data nb", dn), ("chain nb", sn), ("fit nb", fit)):
        print(f"   {name:9s} " + " ".join(f"{v:5.3f}" for v in y[::8]))
