#!/bin/bash
# yoga driver, second batch: det3 X over time / drift field / beam, and the matching Y -> $O/<tag>.npz (6 files each)
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python; S=$(cd $(dirname $0) && pwd)/raw_asym.py
O=/media/dylan/data/claude_tmp/tmp/nb/raw; B=/home/dylan/x17/cosmic_bench; BJ=/media/dylan/data/x17/beam_july/runs
while read tag dir feu; do
  [ -s $O/$tag.npz ] || $P $S $dir/decoded_root $feu $O/$tag.npz --files 6 > $O/$tag.log 2>&1 &
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 5; done
done <<L
t623s_x $B/det3_det4/mx17_det3_det4_overnight_6-23-26/short_run 3
t623l_x $B/det3_det4/mx17_det3_det4_overnight_6-23-26/longer_run 3
t623L_x $B/det3_det4/mx17_det3_det4_overnight_6-23-26/long_run 3
t622L_x $B/det2_det3/mx17_det2_det3_overnight_6-22-26/long_run 3
t622s_x $B/det2_det3/mx17_det2_det3_overnight_6-22-26/short_run 3
t627d100_x $B/det3/mx17_det3_saturday_scan_6-27-26/drift_scan_resist_490V_drift_100V 7
t627d1100_x $B/det3/mx17_det3_saturday_scan_6-27-26/drift_scan_resist_490V_drift_1100V 7
t627L_x $B/det3/mx17_det3_saturday_scan_6-27-26/long_run_resist_490V_drift_1000V 7
t628wk_x $B/det3/mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check 7
b79_x $BJ/run_79/stat090_0000 3
b79_y $BJ/run_79/stat090_0000 4
L
wait; echo RAW2 DONE
