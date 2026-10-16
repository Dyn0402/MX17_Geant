# Avalanche slices with null-collision steps — 490 V meshfield (straight map), 2026-10-16

Run: `aval_nullcoll_490.sh` (desktop, local Garfield on ROOT 6.40), seeds 49001–49004, 140 ev each,
the same seeds with (`nc1`) and without (`nc0`) `--null-collision`. Raw JSONs: desktop `~/x17/aval_nullcoll_2026-10-16/`.

| seed  | surv nc0 | surv nc1 | <g> nc0 | <g> nc1 |
|-------|---------:|---------:|--------:|--------:|
| 49001 | 0.950 | 0.993 | 23 850 | 22 108 |
| 49002 | 0.957 | 0.986 | 22 830 | 24 631 |
| 49003 | 0.921 | 1.000 | 20 694 | 21 545 |
| 49004 | 0.971 | 1.000 | 21 485 | 24 580 |
| **pooled (560 seeds)** | **0.950** | **0.995** | | |

- Gain given survival (the Polya the digitizer applies): nc0 23 384 ± 679, nc1 23 341 ± 654, **ratio 0.998 ± 0.040**.
- Unconditional mean gain: ratio 1.045 ± 0.043, which is just the survival change.
- Matches the mesh-cell G7 result (0.951 → 0.990). The production pooled calib's survival of 0.9559 is a
  stepping artefact, so the chain's mesh transparency is ~4 % low. The gain per surviving electron does not move.

## Follow-up launched 16:18

Full 56-slice campaign with `NULLCOLL=1`, 12-way, out `/media/ucla/mx17_response_sim/avalanche/results_meshfield_nc1_20261016/`
(root and /media/dylan/data are both <4 GB free).

**Prediction for the pooled nc1 point** (vs `aval_calib_meshfield_pooled.json`, nc0: survival 0.9559, gain|surv 24 094, theta 1.14):
survival 0.990–0.997; gain|surv within ±3 % of 24 094; theta within ±0.1; sigma0 and t_arrival within 3 %.
Falsifier: gain|surv shifts >5 %, or survival <0.985.
