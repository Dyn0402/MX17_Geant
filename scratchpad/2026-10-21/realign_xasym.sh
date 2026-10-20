#!/bin/bash
# Re-alignment check of nTof_x17 wft.xasym (2026-10-21): per bench run, the production alignment (mx_june_wft/01_alignment.py)
# on the production events.parquet and on a copy with x_p0 corrected, both into scratch outputs (production untouched).
# Prints the X/Y resolution lines of each. Usage: realign_xasym.sh <scratch dir>
set -u
S="$1"; mkdir -p "$S"; cd ~/PycharmProjects/nTof_x17
PY=.venv/bin/python
while read -r key build scale; do
  [ -z "$key" ] && continue
  T=$($PY -c "import sys; sys.path[:0]=['mx_june_cosmic_qa','cosmic_bench_analysis','.']; from qa_config import get_config; import os; print(os.path.join(get_config('$key').OUT_BASE,'wft','events.parquet'))" 2>/dev/null | tail -1)
  mkdir -p "$S/$key"
  $PY - "$T" "$S/$key" "$build" "$scale" <<'PYEOF'
import sys, shutil, pandas as pd
sys.path.insert(0, '.')
from wft import xasym
t, d, b, sc = sys.argv[1], sys.argv[2], sys.argv[3].replace('_', ' '), float(sys.argv[4])
df = pd.read_parquet(t); shutil.copy(t.replace('.parquet', '.meta.json'), f'{d}/corr.meta.json')
df['x_p0'] = xasym.correct_x(df.x_p0.values, df.y_p0.values, b, sc); df.to_parquet(f'{d}/corr.parquet')
print('corrected', len(df), b, sc)
PYEOF
  for arm in orig corr; do
    tab=$([ $arm = orig ] && echo "$T" || echo "$S/$key/corr.parquet")
    echo "== $key $arm"
    $PY mx_june_wft/01_alignment.py "$key" --table "$tab" --out "$S/$key/align_$arm" 2>&1 | grep -E "resolution|clean events" 
  done
done <<'LIST'
g_det3 det3_June-15 0.96
g_det3_wknd det3_June-15 0.44
sat_det3 det3_June-15 0.44
g_det2 det2 1.0
o22_long_det2 det2 1.0
LIST
