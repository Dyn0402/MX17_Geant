source $HOME/Software/root_6_30/bin/thisroot.sh
cd /home/dylan/PycharmProjects/nTof_x17/garfield_sim && source setup_garfield.sh >/dev/null 2>&1; cd /home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-13c/magboltz
python3 drift_table.py he_c2h6_96p5_3p5.json he 96.5 c2h6 3.5 > he.log 2>&1 &
python3 drift_table.py ar_cf4_90_10.json ar 90 cf4 10 > cf4.log 2>&1 &
python3 drift_table.py ar_ic4h10_95_5.json ar 95 ic4h10 5 > ariso.log 2>&1 &
wait
