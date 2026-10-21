#!/bin/bash
# Runs on the desktop (2026-10-21 evening): wait for the narrow-strip X kernels (condor 4526243), build the 5 MOhm hybrids on
# lxplus, submit the paired full-sim arms w190, w190f1, w290, w290f1 (funnel width matched; th00 + thy*, seeds 7-12, paired
# with glg0). Then wait for the arms and pull their npz to the desktop job tmp.
set -u
D=/afs/cern.ch/work/d/dneff/mx17_s1
until [ "$(ssh lxplus 'condor_q 4526243 -af JobStatus 2>/dev/null | wc -l')" = "0" ]; do sleep 300; done
ssh lxplus "ls /eos/user/d/dneff/tmp/gx/ | grep -E 'w190|w290'"
ssh lxplus "set -e; source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh >/dev/null 2>&1; cd /tmp && xrdcp -s -f root://eosuser.cern.ch//eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/products/greens_comb_w2_rho5M_dk50um_g19um.npz /tmp/dn_prod.npz && for t in gap0_w190 gap0_w290; do xrdcp -s -f root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/\${t}_nx1560/gx_w2_rho5M_\${t}_nx1560.npz /tmp/dn_gx.npz; python3 $D/build_hybrid.py /tmp/dn_prod.npz /tmp/dn_gx.npz /tmp/dn_hyb.npz; xrdcp -s -f /tmp/dn_hyb.npz root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/hybrid_rho5M_\${t}_nx1560.npz; done; rm -f /tmp/dn_prod.npz /tmp/dn_gx.npz /tmp/dn_hyb.npz
cd $D && condor_submit stageBC_narrow.sub"
echo CHAIN_SUBMITTED
sleep 600
until [ "$(ssh lxplus 'condor_q -af Cmd Args 2>/dev/null | grep -c "run_funnelw_job.sh glw[12]90"')" = "0" ]; do sleep 600; done
mkdir -p $HOME/x17/nb_asym_raw/simnarrow
rsync -a --include='glw190*' --include='glw290*' --exclude='*' lxplus:$D/simgapleak_npz/ $HOME/x17/nb_asym_raw/simnarrow/
ls $HOME/x17/nb_asym_raw/simnarrow | wc -l
echo ARMS_DONE
