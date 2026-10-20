#!/bin/bash
# Runs on the desktop/laptop: wait for the anisotropic X products (condor 4520948), build the 5 MOhm hybrids on lxplus,
# submit the paired full-sim arms (glax3, glax10, glax10w690, glax10w690f1; th00 + thy*, seeds 7-12). (2026-10-21)
set -u
D=/afs/cern.ch/work/d/dneff/mx17_s1
until [ "$(ssh lxplus 'condor_q 4520948 -af JobStatus 2>/dev/null | wc -l')" = "0" ]; do sleep 300; done
ssh lxplus "ls /eos/user/d/dneff/tmp/gx/ | grep ax"
ssh lxplus "set -e; source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh >/dev/null 2>&1; cd /tmp && xrdcp -s -f root://eosuser.cern.ch//eos/experiment/ntof/data/x17/response_sim/s1_w2_ny512/products/greens_comb_w2_rho5M_dk50um_g19um.npz /tmp/dn_prod.npz && for t in gap0_ax3 gap0_ax10 gap0_w690_ax10; do xrdcp -s -f root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/\${t}_nx1560/gx_w2_rho5M_\${t}_nx1560.npz /tmp/dn_gx.npz; python3 $D/build_hybrid.py /tmp/dn_prod.npz /tmp/dn_gx.npz /tmp/dn_hyb.npz; xrdcp -s -f /tmp/dn_hyb.npz root://eosuser.cern.ch//eos/user/d/dneff/tmp/gx/hybrid_rho5M_\${t}_nx1560.npz; done; rm -f /tmp/dn_prod.npz /tmp/dn_gx.npz /tmp/dn_hyb.npz
cd $D && : > aniso_list.txt; for arm in 'ax3:hybrid_rho5M_gap0_ax3_nx1560.npz:0' 'ax10:hybrid_rho5M_gap0_ax10_nx1560.npz:0' 'ax10w690:hybrid_rho5M_gap0_w690_ax10_nx1560.npz:0' 'ax10w690f1:hybrid_rho5M_gap0_w690_ax10_nx1560.npz:1'; do IFS=: read t k fu <<<\"\$arm\"; for a in th00 thy05 thy10 thy15 thy20 thy30; do for s in 7 8 9 10 11 12; do echo \"\$t, \$k, \$fu, \$a, \$s\" >> aniso_list.txt; done; done; done; sed -e 's|gapwidth_list.txt|aniso_list.txt|' -e 's|^# Strip-width arms.*|# Anisotropic-ESL arms (2026-10-21): hybrid kernels (production 5M Y + sigma_xx = A sigma_yy X), A 3, 10, 10 + 690 um, + funnel;|' stageBC_gapwidth.sub > stageBC_aniso.sub && condor_submit stageBC_aniso.sub"
echo CHAIN_SUBMITTED
