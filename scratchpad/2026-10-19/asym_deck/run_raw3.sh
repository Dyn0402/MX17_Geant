#!/bin/bash
# yoga driver, third batch: the raw extractor over every cosmic-bench (run, sub-run, detector) of ../survey_jobs.tsv,
# both views, 2 files each -> $O/survey/s<NNN>_<x|y>.npz, index $O/survey/index.tsv (file, run_dir, sub, det, view, feu)
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python; D=$(cd $(dirname $0) && pwd); S=$D/raw_asym.py
O=/media/dylan/data/claude_tmp/tmp/nb/raw/survey; mkdir -p $O; : > $O/index.tsv; i=0
while IFS=$'\t' read -r rd sub det; do
  i=$((i+1)); [[ $rd == *beam_july* ]] && continue
  for v in x y; do
    feu=$(python3 -c "import json;rc=json.load(open('$rd/run_config.json'));print([d for d in rc['detectors'] if d['name']=='$det'][0]['dream_feus']['${v}_1'][0])")
    f=$O/s$(printf %03d $i)_$v.npz
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$rd" "$sub" "$det" "$v" "$feu" >> $O/index.tsv
    [ -s $f ] || $P $S "$rd/$sub/decoded_root" $feu $f --files 2 > ${f%.npz}.log 2>&1 &
    while [ $(jobs -r | wc -l) -ge 5 ]; do sleep 3; done
  done
done < $D/../survey_jobs.tsv
wait; echo RAW3 DONE
