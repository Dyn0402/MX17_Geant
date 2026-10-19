# Neighbour delay vs strip distance, X vs Y (2026-10-19, desktop)

Vertical tracks (data: M3 |ax|,|ay| < 3°, 287 X / 305 Y events after dropping railed centrals; sim: th00, ~2000 events).
`dump_pmN_m3.py` (±5 strips, run on yoga; npz in yoga `/media/dylan/data/claude_tmp/tmp/nb/{v3rw,v4c}/waveforms_pm5.npz`),
`nb_far.py` → `nb_far.log` (bootstrap errors), `nb_far_plot.py` → `nb_far_shapes.png`. t = 0 at the central's leading 50 %.
±1 reproduces 2026-10-18/nb_compare.py exactly.

Peak delay vs central [ns] (folded ±k), and the 50 % leading-edge delay (more robust on the far strips):

| | data X | v3rw X | v4c X | data Y | v3rw Y | v4c Y |
|---|---|---|---|---|---|---|
| ±1 peak | +45 ± 14 | +135 | +75 | +60 ± 10 | +135 | +60 |
| ±2 peak | +315 ± 53 | +225 | +90 | +420 ± 21 | +330 | +300 |
| ±3 peak | +420 ± 110 | +330 | +300 | +600 ± 29 | +600 | +810 (flat top) |
| ±4 peak | noisy (±140) | no signal | no signal | after end of record | after end | after end |
| ±1 edge | −25 ± 7 | +22 | +2 | −11 ± 5 | +34 | +1 |
| ±2 edge | −7 ± 30 | +91 | +16 | +88 ± 17 | +178 | +109 |
| ±3 edge | +48 ± 31 | +240 | +164 | +279 ± 31 | +375 | +395 |
| ±4 edge | +74 ± 35 | – | – | +579 ± 42 | +863 | +886 |

Peak / central: X ±1 .368/.372/.385, ±2 .069/.103/.095, ±3 .012/.010/.006, ±4 .003/.001/.000;
Y ±1 .408/.435/.456, ±2 .156/.157/.140, ±3 .046/.031/.032, ±4 .019/.012/.012 (still rising at record end).

Reading:
- X and Y differ beyond ±1: X far strips are a short pulse (gone by ~+900 ns); Y far strips are a slow rise lasting past the
  record (±4 still rising at +1.4 µs), in data and sim alike: the resistive spread runs along Y.
- The data's far strips (X ±3/±4, Y ±3/±4) carry a prompt component that rises with the central's leading edge
  (X ±3 plateau ~0.010, X ±4 ~0.0025, Y ±4 step ~0.006) and then the slow/late part. The sim has only the late part.
  That is why the data's far edges are early (Y ±4 by ~300 ns, X ±3 by ~120 ns) and the far amplitudes 1.5–2× the sim's.
- v4c's prompt anode spread (σ0 420 µm) fixes ±1 but makes X ±2 a prompt copy (peak +90 vs data +315); the data ±2 X
  starts prompt and peaks late. The σ0 spread is not the data's mechanism at ±2.
- Data X is side-asymmetric: −1/+1 peak 0.416/0.321, −2/+2 0.090/0.051 (5σ). Y is symmetric (0.416/0.401), sim symmetric.
  Not a track lean (selected mean ax 0.01°, ay 0.13°). Open.
