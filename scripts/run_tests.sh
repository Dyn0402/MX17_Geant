#!/bin/bash
# The response-chain test suite, in one place so it is not re-typed per machine.
#   scripts/run_tests.sh <s1_product.npz> <aval_calib.json> [dense_ny.npz]
set -uo pipefail
K="${1:?usage: run_tests.sh <kernel.npz> <calib.json> [ny1024.npz]}"
C="${2:?}"; NY="${3:-}"
cd "$(dirname "$0")/.."
fail=0
run() { echo "### $1"; shift; "$@" 2>&1 | tail -2; [ "${PIPESTATUS[0]}" -eq 0 ] || fail=1; }
run test_longitudinal       python3 -m response.digitizer.test_longitudinal
run test_lut_vs_solver      python3 -m response.digitizer.test_lut_vs_solver --kernel "$K"
run test_charge_audit       python3 -m response.digitizer.test_charge_audit --kernel "$K" --calib "$C" --n-rep 16
run test_induce_equivalence python3 -m response.digitizer.test_induce_equivalence --kernel "$K" --calib "$C"
run test_daq_wft            python3 -m response.dream.test_daq_wft
# test_ny_grid is a CERTIFICATE, not a gate: its verdict is a statement about
# whichever pair you hand it, and the production answer is already known — at
# ny=512 it FAILS by design (0.452 % vs a 0.3 % bar, audit C6), which is why
# production moved to ny=1024. Run it for the record, do not fail the suite on
# it.
if [ -n "$NY" ]; then
  echo "### test_ny_grid (informational — see audit C6)"
  python3 -m response.digitizer.test_ny_grid "$K" "$NY" 2>&1 | tail -2
fi
echo; [ "$fail" -eq 0 ] && echo "SUITE PASS" || echo "SUITE FAIL"
exit "$fail"
