#!/bin/bash
# Finishing half of narrow_chain.sh for cluster 4526326 (2026-10-21): the chain's original wait loop grepped condor_q
# Args, which is undefined for this sub, and exited at once. Wait on the cluster, then pull the glw190/glw290 npz.
set -u
D=/afs/cern.ch/work/d/dneff/mx17_s1; CL=${1:-4526326}
until [ "$(ssh lxplus "q=\$(condor_q $CL -af JobStatus 2>/dev/null) && echo -n \"\$q\" | grep -c [^5]")" = "0" ]; do sleep 600; done
ssh lxplus "condor_q $CL -af ProcId JobStatus HoldReason 2>/dev/null" | sed 's/^/HELD /'
mkdir -p $HOME/x17/nb_asym_raw/simnarrow
rsync -a --include='glw190*' --include='glw290*' --exclude='*' lxplus:$D/simgapleak_npz/ $HOME/x17/nb_asym_raw/simnarrow/
ls $HOME/x17/nb_asym_raw/simnarrow | wc -l
echo FINISH_DONE
