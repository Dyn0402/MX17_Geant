#!/usr/bin/env python3
"""
check_ion_template_repro.py — does the FIXED S3 producer reproduce v2's ion template?

    python3 scripts/check_ion_template_repro.py <v2.json> <repro_slice.json> [more slices...]

WHY THIS EXISTS. The first S3 campaign wrote `i_elec`/`i_ion` arrays that were
identically zero in all 56 slices: `ComponentConstant` has no weighting
potential until `SetWeightingPotential()` is called, and Garfield computes the
induced current from psi rather than the weighting field by default, so every
signal came out zero while the avalanche and ion drift ran perfectly normally
(gain and sigma0 were unaffected, which is why it went unnoticed). Commit
42390d1 fixed it, and `aval_calib_v2.json` came from a re-run afterwards —
whose raw slices were never archived. So v2's templates are trustworthy but
were not, until this, reproducible from anything on EOS.

WHAT IT CHECKS. f_ion = Qi/(Qe+Qi), the fraction of the induced charge carried
by the ion tail. It is the right observable because it is a RATIO: insensitive
to gain, to slice count, and to the ~1 % run-to-run spread that makes v2 a
different realization from the archived raw (v2's 470 V gain is 22854.3 against
the archived 23107.4). It is also the number with an independent internal
cross-check inside v2 — 0.9079 from the current integrals against 0.9077 from
the alpha_z histogram, agreeing to four digits, which is itself the proof that
no ion charge fell off the end of the 400 ns window.

A PASS here means the fixed producer reproduces the physics of v2's template
from scratch, which closes the provenance gap as far as it can be closed
without re-running all 56 slices.
"""

from __future__ import annotations

import json
import sys

import numpy as np

TOL = 0.005          # absolute on f_ion; run-to-run spread is well under this


def f_ion(i_elec, i_ion, dt):
    qe = float(np.asarray(i_elec, dtype=float).sum()) * dt
    qi = float(np.asarray(i_ion, dtype=float).sum()) * dt
    tot = qe + qi
    return (qi / tot if tot else float("nan")), qe, qi


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    v2 = json.load(open(argv[0]))
    key = [k for k in v2["points"] if "490" in k]
    if not key:
        print("FAIL: no 490 V point in the v2 calib")
        return 3
    p = v2["points"][key[0]]
    ref, qe_r, qi_r = f_ion(p["i_elec"], p["i_ion"], p["signal_dt_ns"])
    print(f"v2 reference  {key[0]}")
    print(f"  nev {p['nev_total']}  dt {p['signal_dt_ns']} ns")
    print(f"  Qe {qe_r:.5g}  Qi {qi_r:.5g}  f_ion {ref:.6f}\n")

    # Slices are charge-weighted by their event count, exactly as collect.py
    # merges them — a plain mean over slices would misweight unequal nev.
    tot_e = tot_i = 0.0
    nev = 0
    for path in argv[1:]:
        s = json.load(open(path))
        n = s.get("nev", s.get("nev_total", 0))
        dt = s.get("signal_dt_ns", p["signal_dt_ns"])
        fi, qe, qi = f_ion(s["i_elec"], s["i_ion"], dt)
        nz = int(np.count_nonzero(np.asarray(s["i_elec"], dtype=float)))
        print(f"  {path.split('/')[-1]}: nev {n}  nonzero i_elec bins {nz}  "
              f"f_ion {fi:.6f}")
        if nz == 0:
            print("    ^ ALL ZERO — this producer still has the psi bug")
        tot_e += qe * n
        tot_i += qi * n
        nev += n

    if not nev:
        print("\nFAIL: no events")
        return 3
    got = tot_i / (tot_e + tot_i)
    d = got - ref
    print(f"\nmerged over {nev} events:  f_ion {got:.6f}")
    print(f"v2:                        f_ion {ref:.6f}")
    print(f"difference                 {d:+.6f}  (tol {TOL})")
    ok = abs(d) <= TOL
    print("\n" + ("PASS — the fixed producer reproduces v2's ion template; the "
                  "provenance gap is closed."
                  if ok else
                  "FAIL — the fixed producer does NOT reproduce v2. v2's "
                  "templates need re-deriving, not just re-archiving."))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
