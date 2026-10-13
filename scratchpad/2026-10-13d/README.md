# 2026-10-13d: protected-FEU input network (FRESH_LOOK §17)

- `protection_rc_scan.py` puts a protection-network filter (a high-pass with series C into R to
  ground, or a low-pass with series R into the strip capacitance) in front of the DREAM shaper,
  in the §7.5/§10.2 toy. It fits the inclined det3 X shapes and reports the rise and the
  undershoot. Run it with the nTof_x17 `.venv`, from
  `yoga:~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/` (it needs both npz files).
  Output: `protection_rc_scan.log`.
