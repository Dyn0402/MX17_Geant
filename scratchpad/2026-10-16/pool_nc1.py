#!/usr/bin/env python3
"""
Pool the 56-slice null-collision (nc1) meshfield campaign into one 490 V point, the way
aval_calib_meshfield_pooled.json was built for nc0: every slice ran on the SAME map
(meshfield_production.txt, solved at 490 V), so the --voltage label is ignored and
collect.py's per-voltage grouping does not apply.

Run from the MX17_Geant src root (for `response.avalanche.collect`):
    python3 pool_nc1.py <raw_dir> <out.json> [<nc0_pooled.json>]
"""
import json
import sys

from response.avalanche.collect import load, merge

raw, out = sys.argv[1], sys.argv[2]
slices = [s for sl in load(raw).values() for s in sl]
fm = {s["field_model"] for s in slices}
assert len(fm) == 1 and fm.pop().startswith("meshfield"), f"mixed/unexpected field models: {fm}"
m = merge(slices)
m.update(voltage_V=490.0, gas_file=slices[0]["gas"], null_collision=True)
json.dump({"schema": "aval_calib_meshfield_pooled/2",
           "note": ("Null-collision steps (--null-collision), patched Garfield 927e5c21p with our headers first "
                    "on ROOT_INCLUDE_PATH (condor 4474926, 2026-10-16). All slices used the same 490 V map "
                    "regardless of their --voltage label; pooled as one 490 V point, as the nc0 file."),
           "point": m}, open(out, "w"))
p = m["polya"]
print(f"nc1: {m['n_slices']} slices, nev {m['nev_total']}, survival {p['survival']:.4f} +- {p['survival_err']:.4f}, "
      f"gain|surv {p['gain_mean']:.0f}, theta {p['theta']:.2f}, sigma0 {m['sigma0_um']:.1f} um, "
      f"t_arr {m['t_arrival_mean_ns']:.2f} ns")
if len(sys.argv) > 3:
    ref = json.load(open(sys.argv[3]))["point"]
    q = ref["polya"]
    print(f"nc0: {ref['n_slices']} slices, nev {ref['nev_total']}, survival {q['survival']:.4f} +- {q['survival_err']:.4f}, "
          f"gain|surv {q['gain_mean']:.0f}, theta {q['theta']:.2f}, sigma0 {ref['sigma0_um']:.1f} um, "
          f"t_arr {ref['t_arrival_mean_ns']:.2f} ns")
    print(f"nc1/nc0: survival x{p['survival'] / q['survival']:.4f}, gain|surv x{p['gain_mean'] / q['gain_mean']:.4f}, "
          f"mean charge x{p['gain_mean_uncond'] / q['gain_mean_uncond']:.4f}")
