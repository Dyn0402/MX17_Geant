#!/bin/bash
# transparency_curve.sh — S2 deliverable (a): epsilon(E_amp/E_drift), per the
# runbook's "Cheap follow-ups". Loops E_DRIFT through solve_fieldmap.py at
# smoke resolution (curve SHAPE only — the bench point at 333 V/cm is quoted
# from the already-accepted production map, git 5763342, not re-solved here),
# then runs gates_check.C's G7 machinery on each map.
#
# Usage: ./transparency_curve.sh [jobs] [nev]
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

JOBS="${1:-4}"     # solve_fieldmap.py sampling workers; kept low so this can
                   # run alongside the T7 avalanche campaign on the same box
NEV="${2:-400}"    # G7 events per point — smoke-resolution convention

POINTS=(50 100 200 500 700 1000)   # 333 comes from the production map instead
OUT_CSV=transparency_curve.csv

source ~/PycharmProjects/nTof_x17/garfield_sim/setup_garfield.sh
ROOT_BIN="$(command -v root || echo /home/dylan/Software/root_6_36_06/bin/root)"

echo "E_drift_Vcm,E_amp_Vcm,ratio_amp_drift,eps,eps_err,funnel_rms_um,source" \
    > "$OUT_CSV"

# ── Bench point: already-accepted production map, no re-solve ──────────────
python3 - "$OUT_CSV" <<'PYEOF'
import json, sys
d = json.load(open("meshfield_production.json"))
g = d["gates"]
# eps/err/funnel_rms are the runbook's accepted, hand-verified numbers
# (production G7 run, git 5763342) — not re-derivable from the JSON, which
# only stores the solver-side S-gates, not the Garfield-side G7 measurement.
with open(sys.argv[1], "a") as f:
    f.write(f"333.0,{g['S2_Ez_amp_Vcm']:.1f},"
            f"{g['S2_Ez_amp_Vcm']/333.0:.3f},0.955,0.005,,production_accepted\n")
PYEOF

for E in "${POINTS[@]}"; do
  tag=$(printf "edrift%04d" "$E")
  map="meshfield_${tag}.txt"
  echo "[curve] E_drift=$E V/cm: solving ($tag) ..."
  python3 solve_fieldmap.py --smoke --e-drift "$E" --tag "$tag" --jobs "$JOBS" \
      > "solve_${tag}.log" 2>&1

  echo "[curve] E_drift=$E V/cm: G7 via gates_check.C ..."
  "$ROOT_BIN" -b -q -e "gSystem->Load(\"libGarfield\");
    gSystem->CompileMacro(\"gates_check.C\",\"k\");
    gates_check(\"$map\", $NEV, $E);" \
      > "gates_${tag}.log" 2>&1 || true

  python3 - "$OUT_CSV" "gates_${tag}.log" "$tag" "$E" <<'PYEOF'
import json, re, sys
out_csv, gatelog, tag, E = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
txt = open(gatelog).read()
m_eps = re.search(r"G7 transparency:\s*([\d.]+)\s*\+-\s*([\d.]+).*funnel rms\s*([\d.]+)\s*um", txt)
d = json.load(open(f"meshfield_{tag}.json"))
e_amp = d["gates"]["S2_Ez_amp_Vcm"]
if m_eps:
    eps, err, rms = m_eps.group(1), m_eps.group(2), m_eps.group(3)
else:
    eps, err, rms = "NaN", "NaN", "NaN"
    print(f"[curve] WARNING: could not parse G7 from {gatelog}", file=sys.stderr)
with open(out_csv, "a") as f:
    f.write(f"{E:.1f},{e_amp:.1f},{e_amp/E:.3f},{eps},{err},{rms},smoke\n")
PYEOF
done

echo "[curve] wrote $OUT_CSV:"
cat "$OUT_CSV"
