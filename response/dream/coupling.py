"""
coupling.py — inter-strip coupling between the readout strips and the FEU inputs.

FRESH_LOOK §21.5–21.6: the chain gets the central strip's timing and the ±1 neighbours'
amplitude right, but the neighbours peak +120..+150 ns after the central strip, while the
data show +45/+60 ns. A fit of the data neighbour against the chain neighbour needs
a prompt copy of the central pulse (k) plus a derivative term (k' ≈ 0.15–0.2 per 60 ns
sample). Resistive-layer spreading cannot produce either, because it is slow by
construction. This is the signature of coupling capacitance between adjacent channels:

  * DERIVATIVE term. Each input node has a finite series impedance R (protection
    network, trace). The central current moves its node voltage, v_c = R i_c, and the
    coupling capacitance C_c injects tau di_c/dt into each neighbour, with tau = R C_c.
    That current leaves the central node, so the central keeps i_c - 2 tau di_c/dt,
    roughly i_c delayed by 2 tau. The term is bipolar and carries no net charge.
  * IN-PHASE term. With a finite-gain CSA the node follows the integrated charge, so
    C_c forms a divider: each neighbour receives k i_c and the central keeps
    (1 - 2k) i_c.

    di_n = k (i_{n-1} + i_{n+1} - 2 i_n) + tau d/dt (i_{n-1} + i_{n+1} - 2 i_n)

(Kirchhoff at each node: whatever a neighbour gains, the central loses.)

DIAGNOSTIC ONLY: an optional third value KFREE adds kf i_c to each neighbour WITHOUT taking it from the
central. That does not conserve charge, so it is not a circuit model. It asks whether the data's prompt
neighbour share (FRESH_LOOK §21.8 and §21.11: the residual in-phase k ≈ 0.05 with no central loss) closes
the rest of the ±1 delay.

This acts on the induced currents, per view and in channel order (channel = strip
index, run.py), before the shaper. Both steps are linear, so the order is only a
matter of convention. The truth sidecar keeps the uncoupled induced charge.
"""

from __future__ import annotations

import numpy as np


def parse_coupling(spec):
    """'K,TAU_NS[,KFREE]' or 'X:K,TAU[,KFREE]/Y:...' -> {'X': (k, tau_ns, kfree), 'Y': (...)}."""
    if spec is None:
        return None
    out = {}
    parts = spec.split("/")
    for p in parts:
        if ":" in p:
            view, vals = p.split(":", 1)
            views = (view.strip().upper(),)
        else:
            if len(parts) > 1:
                raise ValueError(f"mix of per-view and global coupling: {spec!r}")
            views, vals = ("X", "Y"), p
        f = [float(v) for v in vals.split(",")]
        if len(f) not in (2, 3):
            raise ValueError(f"coupling wants K,TAU or K,TAU,KFREE, got {vals!r}")
        k, tau, kf = f[0], f[1], (f[2] if len(f) == 3 else 0.0)
        if not 0.0 <= k < 0.5 or tau < 0 or not 0.0 <= kf < 0.5:
            raise ValueError(f"coupling needs 0 <= k, kfree < 0.5 and tau >= 0, got {k}, {tau}, {kf}")
        for v in views:
            out[v] = (k, tau, kf)
    return out


def couple(cur, coupling, dt_ns):
    """Apply inter-strip coupling to {(view, idx): current waveform}. Returns a new dict.
    Neighbours that carried no current gain entries."""
    if not coupling or not cur:
        return cur
    out = {key: np.asarray(w, dtype=float).copy() for key, w in cur.items()}
    n = max(len(w) for w in cur.values())
    for (view, idx), w in cur.items():
        if view not in coupling:
            continue
        k, tau, kf = coupling[view]
        w = np.asarray(w, dtype=float)
        # Central-difference derivative. np.gradient keeps the length and the sum is
        # ~0, so the derivative term carries no net charge.
        give = k * w + (tau / dt_ns) * np.gradient(w) if tau > 0 else k * w
        out[(view, idx)] -= 2.0 * give
        give = give + kf * w          # the diagnostic term is not taken from the central
        for j in (idx - 1, idx + 1):
            tgt = out.get((view, j))
            if tgt is None:
                tgt = out[(view, j)] = np.zeros(n)
            tgt[:len(give)] += give
    return out
