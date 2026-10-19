# Raw-ADC check of the X asymmetry (2026-10-19)

Deck: https://dylan-neff.web.cern.ch/notes/det3-x-asymmetry.html, built by `make_asym_deck.py` (observables in
`asym_lib.py`). Rebuild: `~/PycharmProjects/nTof_x17/.venv/bin/python make_asym_deck.py` then
`python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py det3-x-asymmetry.html --slug det3-x-asymmetry --force --deploy`.
Raw survey done (151 detector-views x sub-runs): det3 June-15 X median 0.137 (0.04-0.21, 57 sub-runs), its Y 0.00 +- 0.04;
det3 May X -0.008; det4 X 0.05, det7 X 0.06; det1/2/6 ~0. 2D map (6-24, 4495 events): X blocks 1-2 ~0.25, X block 0 ~0.1.

Extractor `raw_asym.py` (yoga): decoded_root ADC minus per-channel pedestal only. No strip map, no common-mode
subtraction, no clustering, no M3; DAQ channel numbers. Outputs: yoga `/media/dylan/data/claude_tmp/tmp/nb/raw/`,
copied to desktop `~/x17/nb_asym_raw/`. Batch 3 (`run_raw3.sh`, every cosmic sub-run in ../survey_jobs.tsv, both views)
was started on yoga at wrap time and writes `raw/survey/` + `index.tsv`. Beam files give 0 events with this selection (flash).

Good-pulse cut used below: 200 < A_c < 3500, central peak sample 4..22, |mean(samples 0-2)| < 0.1 A_c.
P = fraction of events with A(p+1) > A(p-1); a1 = in-phase (p+1 minus p-1) projection on the central pulse.

| sample | P | a1 | per chip (0..7) |
|---|---|---|---|
| det3 June-15 X, 6-24 long (FEU 3) | 0.85 | +0.177 | .08 .10 .22 .29 .28 .22 .17 .10 |
| det3 June-15 X, 6-23 long (FEU 3) | 0.62 | +0.158 | .07 .15 .19 .21 .18 .16 .12 .09 |
| det3 June-15 X, 6-28 weekend (FEU 7) | 0.56 | +0.092 | .04 .11 .10 .09 .13 .08 .10 .02 |
| det3 June-15 X, 6-27 drift 100 V | 0.74 | +0.159 | -.03 .06 .09 .33 .23 .18 .14 .08 |
| det3 June-15 X, 6-27 drift 1100 V | 0.54 | +0.064 | |
| det3 June-15 Y | 0.49 | -0.008 | flat ~0 |
| det3 May build X (same FEU 3) | 0.49 | -0.012 | flat ~0 |
| det1 X / det2 X / det6 X (FEU 3) | 0.50/0.48/0.50 | +.004/-.009/+.012 | flat |
| det4 X / det7 X (Stephan June-22 batch) | 0.56/0.53 | +0.054/+0.034 | |

Findings:
1. The asymmetry is in the raw ADC: on det3 June-15 X the p+1 neighbour is bigger in 85 % of events; every
   other board/view through the same code sits at 50 %. (R-L)/c is 0.20-0.23 at every sample: an in-phase copy.
2. It is NOT uniform across the board: high in chips 2-5 (board centre), ~0.1 at the edges, same shape on 6-23, 6-24, 6-27.
   A PCB trace geometry would be uniform. Next: 2D map (match X/Y by eventId, d3jun_x + d3jun_y).
3. TRAP: the survey's low 6-23 short/longer values (0.035 -> 0.057, looked like a turn-on over the evening) are a hot
   channel (ch 492, 3259/4493 events of the short run). With a hot-channel veto the longer run has the long run's chip profile.
   The survey (copy_survey.csv) needs a hot-channel veto before any time-dependence claim.
4. Noise (quiet events, common mode removed): corr(d=1)-corr(d=2) is -0.03 on det3 June X, same as det3 Y and May X;
   injecting a 0.15 copy gives +0.13. Not copied in noise -> not the FEU chain after the input.
5. Open: whether the one-sided bus side predicts the observed direction (gerber -> map orientation unpinned), and
   which data view sits on L5 vs L6 (bus solve gives 0.010 on L5, 0.001 on L6). May and June-15 det3 may be different PCBs.

Original deck plan (slide-note, `make_asym_deck.py` to write): setup cartoon; suspects map; raw event; 3-step processing with
struck steps; coin-flip bars; L/c vs R/c hists; mean waveforms; copy shape; copy model + uncopy; hot-channel trap;
same-FEU grid; noise; board map; amplitude; survey; time/drift; bus geometry; verdicts; what it doesn't rule out.
