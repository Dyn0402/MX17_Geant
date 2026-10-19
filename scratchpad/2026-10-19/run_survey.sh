#!/bin/bash
# Yoga driver (2026-10-19): dump_det.py over survey_jobs.tsv (run_dir, subrun, detector) -> survey/<n>.npz + index.tsv.
D=$(cd $(dirname $0) && pwd); P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python; O=/media/dylan/data/claude_tmp/tmp/nb/survey
cd ~; i=0
while IFS=$'\t' read -r rd sub det; do
  i=$((i+1)); f=$O/s$(printf %03d $i).npz
  [ -s $f ] || $P $D/dump_det.py "$rd" "$sub" "$det" $f --files 3 --jobs 6 > $O/s$(printf %03d $i).log 2>&1
  printf '%s\t%s\t%s\t%s\n' "$f" "$rd" "$sub" "$det" >> $O/index.tsv
  echo "$(date +%T) $i $det $sub $(tail -1 $O/s$(printf %03d $i).log)"
done < $D/survey_jobs.tsv
echo "$(date +%T) SURVEY DONE"
