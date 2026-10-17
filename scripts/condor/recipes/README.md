# Stage B/C baseline recipes

Each `.args` file holds the flags a .sub passes as `run_stageBC_w2.sh`'s 5th argument. Pair it with the rho tag, noise file
and clusters named below. `run.py`'s own defaults are left bare (no ions scaling, no feedback, 30 mm) so that earlier runs
stay reproducible; the physics baseline lives here.

| recipe | rho tag | noise | clusters | from |
|---|---|---|---|---|
| `det3_v3.args` | rho5M | `calib/noise_det3_satscan_feu07_cmq.json` | `clusters_ang_e17/mx17_muons_3k_ang_<ang>_e17_t0.root` | FRESH_LOOK §21.7 (2026-10-18) |

det3 v3: ion reconstructed K0 2.05 (systematic 1.80 / 2.25), drift column 26.5 mm (systematic 27.9 = GAP_STUDY), bG 0.45
(joint fit at that column; s1.75 instead of K0 2.05 → bG ~0.443), nc1 avalanche calib, DREAM tp 166 ns.
