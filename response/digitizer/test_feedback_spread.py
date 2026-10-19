#!/usr/bin/env python3
"""
test_feedback_spread.py — Digitizer.add_feedback_spread (2026-10-19).

Checks, on a stand-in object (no kernel LUT needed):
  1. mean total charge per primary = gbar / (1 - bG), as add_feedback;
  2. one generation's displacement projected on one axis is a Cauchy of
     half-width H (median |dx| = H, P(|dx| > 3H) = 1 - 2/pi atan 3, less the
     ~1 % the 10 mm cap removes);
  3. no displacement exceeds rmax;
  4. secondaries carry their parent's time;
  5. --feedback-spread-um 0 keeps the historical path (transport returns the
     primaries only; add_feedback is what runs).

    python3 -m response.digitizer.test_feedback_spread
"""
import types
import numpy as np
from .digitize import Digitizer


def stand_in(bg=0.45, H_um=150.0, rmax_mm=10.0, seed=1):
    d = types.SimpleNamespace(rng=np.random.default_rng(seed), feedback_bg=bg, gbar=2.0e4, theta=1.15,
                              feedback_spread_um=H_um, feedback_rmax_mm=rmax_mm)
    return d


def main():
    n = 200000
    d = stand_in()
    x = np.zeros(n); y = np.zeros(n); t = np.arange(n, dtype=float)
    g = d.rng.gamma(1 + d.theta, d.gbar / (1 + d.theta), n)
    X, Y, T, Q = Digitizer.add_feedback_spread(d, x, y, t, g)
    want = g.sum() / (1 - d.feedback_bg)
    ratio = Q.sum() / want
    print(f"1. total charge / (sum g / (1 - bG)) = {ratio:.4f}")
    assert abs(ratio - 1) < 0.01

    # one generation only: secondaries get zero gain, so they cannot branch (stub the module's polya_sample)
    from . import digitize as D
    real = D.polya_sample
    D.polya_sample = lambda rng, gbar, theta, m: np.zeros(m)
    try:
        d2 = stand_in(bg=0.45)
        Xg, Yg, Tg, _ = Digitizer.add_feedback_spread(d2, x, y, t, g * 4)
    finally:
        D.polya_sample = real
    dx = np.abs(Xg[n:]) * 1e6
    H = d2.feedback_spread_um
    med = np.median(dx); tail = (dx > 3 * H).mean(); want_tail = 1 - 2 / np.pi * np.arctan(3)
    print(f"2. median |dx| = {med:.1f} um (H = {H}), P(|dx| > 3H) = {tail:.4f} (Cauchy {want_tail:.4f}), "
          f"n secondaries {len(dx)}")
    assert abs(med / H - 1) < 0.03 and abs(tail - want_tail) < 0.015   # the 10 mm cap trims the tail ~1 %
    rr = np.hypot(Xg[n:], Yg[n:])
    print(f"3. max r = {rr.max() * 1e3:.2f} mm (cap {d2.feedback_rmax_mm})")
    assert rr.max() <= d2.feedback_rmax_mm * 1e-3 + 1e-12
    # 4. time: every secondary's time equals some primary time; with t = index and x = 0 the parent is not
    # recoverable from position, so check the multiset: secondary times are integers within the primaries' range
    assert np.all(np.isin(Tg[n:], t))
    print("4. secondary times are parent times")
    # 5. spread 0: transport's branch is `if spread > 0` -> add_feedback; check the method exists and the branch text
    import inspect
    src = inspect.getsource(Digitizer.transport)
    assert "if self.feedback_spread_um > 0:" in src and "gain = self.add_feedback(gain)" in src
    print("5. spread 0 keeps add_feedback (historical path)")
    print("PASS")


if __name__ == "__main__":
    main()
