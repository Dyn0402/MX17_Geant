# Is the det3 X asymmetry only over the Y pads? (2026-10-21)

The old pad-class labels, (p_x + p_y) % 2, use the X peak, which the asymmetry itself moves (peak migration), so a
class split can't locate the effect. Two X-independent labels were tried.

1. **M3 column** (`m3_align.py`, tagged dumps from yoga `/media/dylan/data/claude_tmp/tmp/nb/tagged`): affine fit of the
   3-channel centroid on M3 (x, y, slope). Residual core sigma 0.9-2.0 mm on det2/3/4, both views (det6/det7 don't align).
   Pads are 0.78 mm, so M3 can't place an event on a pad. Dead end.
2. **Y-view share alone** (`ysh_class.py`; raw survey, NB_RAW=~/x17/nb_asym_raw). s = alpha_Y(-1) + alpha_Y(+1): over a
   Y pad the charge sits on one Y row (small s), over an X pad it splits between two rows (large s). On det3 May
   (quiet, labels can't migrate) the Y-pad fraction falls 0.80 (s < 0.3) -> 0.35 (1.25-1.5) -> 0.20 (> 1.5).
   - **det3 June-15 a1X over all events is flat in s:** +0.17/+0.19 (s < 0.3), +0.11..+0.14 (0.5-1.0), +0.15/+0.18 (s > 1.25).
     An effect only over Y pads would fall 3-4x. **The asymmetry is the same over X and Y pads.**
   - Within one s bin the old labels split it +0.26 / -0.48 (narrow): that is migration, not location.
   - det2's positive band (rows 128-255): +0.17..+0.19 at small s, +0.10/+0.03 (±0.04) at s > 1.1: a mild fall, weak stats.
   - Fits an across-strip mechanism (the strips run along y, so the pad type under the charge, which alternates along y,
     shouldn't matter) and excludes anything tied to the Y pads (vias, Y-pad copper, bus beside the Y pads).
3. **Sim arms** (`ysh_class_sim.py`, O = <sgn(delta) a1>; glg0, glax10, glw690f1, glax10w690f1 in `ysh_class_sim.json`,
   xo0/xom200 in `ysh_class_sim_xo.json`): ~1.6k matched events per arm, ±0.03 per bin. Roughly flat too, so not a
   discriminator at this statistics. The sim reproduces the label artefact (X-labelled narrow bins -0.2..-0.4).
