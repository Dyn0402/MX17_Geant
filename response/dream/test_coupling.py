#!/usr/bin/env python3
"""
test_coupling.py — checks on the inter-strip coupling (coupling.py).

  * charge: every channel's integral is conserved in sum; the derivative term moves none
    between strips, the in-phase term moves 2k of the central to the neighbours.
  * timing: the derivative term makes each neighbour lead the central and delays the
    central by ~2 tau.
  * the shaper commutes with it (both linear).

    python3 -m response.dream.test_coupling
"""

import numpy as np

from .coupling import couple, parse_coupling
from .shaper import DreamShaper


def main():
    t = np.arange(1500.0)
    w = np.where(t > 100, (t - 100) / 60.0 * np.exp(-(t - 100) / 60.0), 0.0)
    cur = {("X", 10): w, ("Y", 20): 0.5 * w}

    o = couple(cur, parse_coupling("0,10"), 1.0)
    s = {k: v.sum() / w.sum() for k, v in o.items()}
    assert abs(s[("X", 10)] - 1) < 1e-6 and abs(s[("X", 11)]) < 1e-6, s
    assert o[("X", 11)].argmax() < w.argmax() < o[("X", 10)].argmax()
    print(f"tau 10 ns: central peak {w.argmax()} -> {o[('X', 10)].argmax()} ns, "
          f"neighbour leads at {o[('X', 11)].argmax()} ns, charge moved 0")

    o = couple(cur, parse_coupling("X:0.05,0"), 1.0)
    assert abs(o[("X", 10)].sum() / w.sum() - 0.90) < 1e-9
    assert abs(o[("X", 9)].sum() / w.sum() - 0.05) < 1e-9
    assert ("Y", 21) not in o and np.array_equal(o[("Y", 20)], 0.5 * w)
    print("k 0.05 on X only: central keeps 0.90, each neighbour 0.05; Y untouched")

    sh = DreamShaper(dt_ns=1.0)
    cpl = parse_coupling("0.03,15")
    a = {k: sh.apply(v) for k, v in couple(cur, cpl, 1.0).items()}
    b = couple({k: sh.apply(v) for k, v in cur.items()}, cpl, 1.0)
    err = max(np.abs(a[k] - b[k]).max() for k in a) / np.abs(a[("X", 10)]).max()
    assert err < 1e-3, err
    print(f"shaper commutes with coupling (max rel. diff {err:.1e})")
    print("OK")


if __name__ == "__main__":
    main()
