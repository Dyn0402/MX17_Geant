#!/bin/bash
# yoga driver: raw_asym.py over the deck's run list -> $O/<tag>_<view>.npz
P=/home/dylan/PycharmProjects/nTof_x17/.venv/bin/python; S=$(cd $(dirname $0) && pwd)/raw_asym.py
O=/media/dylan/data/claude_tmp/tmp/nb/raw; mkdir -p $O; B=/home/dylan/x17/cosmic_bench
while read tag dir feu; do
  [ -s $O/$tag.npz ] || $P $S $B/$dir/decoded_root $feu $O/$tag.npz --files ${NF:-3} > $O/$tag.log 2>&1 &
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 5; done
done <<'L'
d3jun_x det3_det4/mx17_det3_det4_overnight_6-23-26/long_run 3
d3jun_y det3_det4/mx17_det3_det4_overnight_6-23-26/long_run 4
d4_x det3_det4/mx17_det3_det4_overnight_6-23-26/long_run 6
d4_y det3_det4/mx17_det3_det4_overnight_6-23-26/long_run 8
d3wk_x det3/mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check 7
d3wk_y det3/mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check 8
d3may_x det_3/mx17_det3_HV_Scan_5-5-26/resist_470V_drift_900V 3
d3may_y det_3/mx17_det3_HV_Scan_5-5-26/resist_470V_drift_900V 4
d2_x det2_det3/mx17_det2_det3_overnight_6-22-26/long_run 6
d2_y det2_det3/mx17_det2_det3_overnight_6-22-26/long_run 8
d6_x det6_det7/mx17_det6_det7_overnight_6-26-26/longer_run 3
d6_y det6_det7/mx17_det6_det7_overnight_6-26-26/longer_run 4
d7_x det6_det7/mx17_det6_det7_overnight_6-26-26/longer_run 6
d7_y det6_det7/mx17_det6_det7_overnight_6-26-26/longer_run 8
d1_x det_1/mx17_det1_daytime_run_1-28-26/overnight_run 4
d1_y det_1/mx17_det1_daytime_run_1-28-26/overnight_run 6
L
wait; echo RAW DONE
