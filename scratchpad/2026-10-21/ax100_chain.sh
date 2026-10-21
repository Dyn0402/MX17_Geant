#!/bin/bash
# Runs on the desktop (2026-10-21 night, FOLLOWUP §9): wait for the sigma_xx = 100 sigma_yy X kernel (condor 4526340), build
# the 5 MOhm hybrid on lxplus, submit the paired arms glax100 / glax100f1 (th00 + thy*, seeds 7-12, paired with glg0), wait,
# pull their npz to ~/x17/nb_asym_raw/simother/. Waits on cluster ids and counts only when condor_q succeeds (an empty or
# failed condor_q never ends a wait); held jobs are listed and do not block. Markers: KERNEL_DONE, ARMS_SUBMITTED <id>,
# AX100_DONE; FAIL lines on errors.
set -u
D=/afs/cern.ch/work/d/dneff/mx17_s1; KCL=${1:-4526340}; T=gap0_ax100
live() { ssh lxplus "q=\$(condor_q $1 -af JobStatus 2>/dev/null) && echo -n \"\$q\" | grep -c [^5]"; }
until [ "$(live $KCL)" = "0" ]; do sleep 300; done
ssh lxplus "condor_q $KCL -af ProcId HoldReason 2>/dev/null" | sed 's/^/HELD /'
if ! ssh lxplus "test -s /eos/user/d/dneff/tmp/gx/${T}_nx1560/gx_w2_rho5M_${T}_nx1560.npz"; then
    echo "FAIL: no kernel product for $T (see $D/logs/gxaniso_-_100.err)"; exit 1; fi
echo KERNEL_DONE
ssh lxplus "set -e; source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh >/dev/null 2>&1; cd /tmp && xrdcp -s -f root://eosuser.cern.ch//eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/products/greens_comb_w2_rho5M_dk50um_g19um.npz /tmp/dn_prod.npz && xrdcp -s -f root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/${T}_nx1560/gx_w2_rho5M_${T}_nx1560.npz /tmp/dn_gx.npz && python3 $D/build_hybrid.py /tmp/dn_prod.npz /tmp/dn_gx.npz /tmp/dn_hyb.npz && xrdcp -s -f /tmp/dn_hyb.npz root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/hybrid_rho5M_${T}_nx1560.npz; rm -f /tmp/dn_prod.npz /tmp/dn_gx.npz /tmp/dn_hyb.npz" || { echo "FAIL: hybrid build"; exit 1; }
SUB=$(ssh lxplus "cd $D && : > ax100_list.txt && for arm in 'ax100:hybrid_rho5M_${T}_nx1560.npz:0' 'ax100f1:hybrid_rho5M_${T}_nx1560.npz:1'; do IFS=: read t k fu <<<\"\$arm\"; for a in th00 thy05 thy10 thy15 thy20 thy30; do for s in 7 8 9 10 11 12; do echo \"\$t, \$k, \$fu, \$a, \$s\" >> ax100_list.txt; done; done; done; sed -e 's|aniso_list.txt|ax100_list.txt|' -e 's|^# Anisotropic-ESL arms.*|# Fast-equalising strips (2026-10-21 night): hybrid kernel sigma_xx = 100 sigma_yy, with and without the funnel;|' stageBC_aniso.sub > stageBC_ax100.sub && condor_submit stageBC_ax100.sub")
echo "$SUB" | tail -1
ACL=$(echo "$SUB" | grep -o "submitted to cluster [0-9]*" | grep -o "[0-9]*$")
[ -n "$ACL" ] || { echo "FAIL: arms not submitted"; exit 1; }
echo ARMS_SUBMITTED $ACL
sleep 600
until [ "$(live $ACL)" = "0" ]; do sleep 600; done
ssh lxplus "condor_q $ACL -af ProcId HoldReason 2>/dev/null" | sed 's/^/HELD /'
mkdir -p $HOME/x17/nb_asym_raw/simother
rsync -a --include='glax100_*' --include='glax100f1_*' --exclude='*' lxplus:$D/simgapleak_npz/ $HOME/x17/nb_asym_raw/simother/
ls $HOME/x17/nb_asym_raw/simother | grep -c '^glax100'
echo AX100_DONE
