"""aval_calib_v3 = v2's measured ion templates + the survival block.

NOT a from-raw regeneration. `collect.py` run over the archived raw slices
reproduces everything EXCEPT i_elec/i_ion, because the raw slices' own
i_elec/i_ion are all zeros (checked 2026-08-08, every voltage) — so v2's
templates came from a different S3 pass and are not reproducible from what is
archived. Regenerating from raw therefore produces a file that looks newer and
is strictly worse: it would be refused by digitize.py's template guard.
"""
import json, sys
v2, v3raw, out = sys.argv[1], sys.argv[2], sys.argv[3]
a = json.load(open(v2)); b = json.load(open(v3raw))
n = 0
for key, pt in a["points"].items():
    src = b["points"].get(key)
    if not src:
        raise SystemExit(f"{key} missing from the from-raw pass — cannot merge")
    for f in ("survival", "survival_err", "n_seeds", "gain_mean_uncond"):
        pt["polya"][f] = src["polya"][f]
    n += 1
a["schema"] = "aval_calib/3"
a["provenance_note"] = (
    "v2 measured i_elec/i_ion + the survival block recovered from the raw "
    "slices (audit A7/Fix 7). The raw slices' own i_elec/i_ion are all zero, "
    "so the templates are NOT reproducible from the archived raw; they are "
    "carried forward from v2 unchanged.")
json.dump(a, open(out, "w"))
print(f"merged {n} points -> {out}")
