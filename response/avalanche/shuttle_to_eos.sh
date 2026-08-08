#!/bin/bash
# shuttle_to_eos.sh — push the overnight field-map ladder + HV-scan campaign
# from /media/ucla (desktop, parked there because root has no headroom) to
# the canonical EOS bulk store, via lxplus (desktop has no local EOS mount;
# rsync-over-ssh to lxplus, which does, is the established "push directly"
# path -- design/RESPONSE_SIM_PLAN.md section 2, verified working
# 2026-08-08 with a round-trip test file).
#
# Requires a live Kerberos ticket (klist) for the lxplus GSSAPI hop.
#
# Usage: ./shuttle_to_eos.sh [what]
#   what: ladder | archive490 | hvscan | all   (default: all)
set -e

EOS=/eos/experiment/ntof/data/x17/response_sim
UCLA=/media/ucla/mx17_response_sim
WHAT="${1:-all}"

check_ticket() {
  if ! klist -s 2>/dev/null; then
    echo "[shuttle] no live Kerberos ticket -- run kinit first. Aborting." >&2
    exit 1
  fi
}
check_ticket

push() {
  local src="$1" dst="$2"
  if [ ! -d "$src" ]; then
    echo "[shuttle] skip (missing) $src"
    return 0
  fi
  echo "[shuttle] $src -> lxplus:$dst"
  ssh lxplus "mkdir -p $dst"
  rsync -av --partial "$src/" "lxplus:$dst/"
}

case "$WHAT" in
  ladder|all)
    push "$UCLA/meshfield_ladder" "$EOS/s2/meshfield_ladder" ;;
esac
case "$WHAT" in
  archive490|all)
    push "$UCLA/avalanche/results_meshfield_490V_20260808" \
         "$EOS/avalanche/raw_meshfield_490V_20260808" ;;
esac
case "$WHAT" in
  hvscan|all)
    push "$UCLA/avalanche/results_meshfield_hvscan" \
         "$EOS/avalanche/raw_meshfield_hvscan_20260808" ;;
esac

# Small merged/pooled JSONs travel with git normally, but ship a copy to EOS
# too since that is where every other calib JSON in this tree already lives.
for f in aval_calib_meshfield_pooled.json aval_calib_meshfield_hvscan.json; do
  if [ -f "$f" ]; then
    echo "[shuttle] $f -> lxplus:$EOS/avalanche/$f"
    rsync -av "$f" "lxplus:$EOS/avalanche/$f"
  fi
done
echo "[shuttle] done ($WHAT)"
