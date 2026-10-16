#!/bin/bash
# Copies each finished nc1 slice off the failing /media/ucla disk as soon as it appears:
# to yoga (LAN) and to EOS (canonical). Only files untouched for >60 s are copied, so a slice still being written is skipped.
# It keeps going until the campaign has finished and a final pass has run.
SRC=/media/ucla/mx17_response_sim/avalanche/results_meshfield_nc1_20261016
YOGA=yoga:/media/dylan/data/mx17_response_sim/avalanche/results_meshfield_nc1_20261016
EOS=lxplus:/eos/experiment/ntof/data/x17/response_sim/avalanche/raw_meshfield_nc1_20261016
pass() {
  cd $SRC || return
  files=$(find . -maxdepth 2 -type f -mmin +1 | sed 's|^\./||')
  [ -z "$files" ] && return
  echo "$files" | rsync -a --files-from=- . $YOGA/ || echo "$(date) YOGA RSYNC FAILED"
  echo "$files" | rsync -a --files-from=- . $EOS/  || echo "$(date) EOS RSYNC FAILED"
  echo "$(date) mirrored $(ls $SRC/*.json 2>/dev/null | wc -l) json"
}
while ! { grep -q "all slices submitted" $SRC/campaign.log && ! pgrep -f "mx17_aval_calib.py.*nc1_20261016" >/dev/null; }; do
  pass; sleep 120
done
sleep 70; pass
# verify: compare sizes on yoga to the source
ssh -o BatchMode=yes yoga "cd ${YOGA#yoga:} && find . -type f -printf '%P %s\n' | sort" > /tmp/mirror_yoga.$$
(cd $SRC && find . -type f -printf '%P %s\n' | sort) > /tmp/mirror_src.$$
if cmp -s /tmp/mirror_src.$$ /tmp/mirror_yoga.$$; then echo "$(date) MIRROR DONE, yoga matches ($(wc -l < /tmp/mirror_src.$$) files)"; else echo "$(date) MIRROR MISMATCH"; diff /tmp/mirror_src.$$ /tmp/mirror_yoga.$$ | head; fi
rm -f /tmp/mirror_*.$$
