#!/bin/bash
# Waits for the first relay pass (pid 684857) to end, then reruns it once: it skips finished points,
# redoes r5verdict1_thx20 (push to yoga was killed at 14:11 on the first pass) and reruns the yoga analysis.
while kill -0 684857 2>/dev/null; do sleep 60; done
echo "$(date) first pass ended; rerun"
bash ~/x17/relay_1016.sh
