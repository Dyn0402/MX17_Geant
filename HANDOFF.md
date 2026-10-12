# Handoff

## T14 fresh look (sim vs data waveforms) — updated 2026-10-12 (dylan-MS-7C84)

**Resume:** shaper confirmed 180 ns, so the ion-suppression puzzle is back. Next: FEU pulser shot, or a resistive-MM ion mechanism.

**Goal:** find out why the end-to-end MX17 response simulation (Geant4 → Garfield → DREAM)
gives slower and wider pulses than the det3 cosmic-bench data. The rise-time / f_ion
contradiction of T14 is the problem being worked here.

**Done (2026-10-12, details in `design/report/FRESH_LOOK_2026-10-12.md` §6–7.6):**
- §6: the quality cut hid a third of the amplitude deficit. The rise time is unaffected.
- §7.1: an independent toy model agrees with the chain, so the chain code has no bug.
- §7.2: vertical tracks can't measure the ion fraction. The rise there mostly reflects where
  the peak lands on a flat top. At matched peak delay, data agrees with full ions.
- §7.3: at 10° and 20°, the data pins only the convolution shaper ⊗ ion. The pulse is about as
  narrow as a 180 ns shaper alone.
- §7.5: an ion-mobility rescale is an exact time rescale, so no Garfield rerun is needed. At
  180 ns, full ions need T_ion ≤ 110 ns (×3–4 Ar⁺), which is unphysical. August's Blanc's-law
  "species 4 % slower" check is invalid: it mixed two different ions
  (`response/meshcell/ion_template_check.py:193`).
- §7.6 + `design/DREAM_SETTINGS.md`: **the shaper is 180 ns**, read back from the target run's
  own archived cfg. The word order is proven (high word first). All 3574 archived cfgs were
  swept, and every data FEU in the cosmic, n_TOF and SPS runs is 180 ns. The 123 ns branch is
  closed.

- Published the state of play to the site: 12-slide note
  https://dylan-neff.web.cern.ch/notes/mx17-response-sim-status.html and an animated chain
  walk-through https://dylan-neff.web.cern.ch/notes/mx17-response-chain.html (toy pulse vs the
  real det3 average; the ion toggle shows the rise problem). It supersedes the August
  `t14-sim-vs-data-waveforms` note.

**In progress / where it stopped:** nothing half-done. The conclusion is that with the true
180 ns shaper the inclined data needs the ion signal largely absent (f_eff ≈ 0.2). That is
August's reading, and it still has **no mechanism**.

**Next steps:**
1. Measure the shaper's *actual* response: fire the FEU pulser once through the test capacitor
   (bench). Chip spread and PZC trim aren't in any register.
2. Look for an ion-signal suppression mechanism specific to the **resistive** bulk Micromegas.
   Sheet screening (≤ 8 %) is the only one sized so far.
3. Still unexplained: the early "foot" before data pulses, and the −0.2 negative lobe on the
   sim's Y channels (§7.2).
4. Amplitude (§6) and gain slope remain separate open threads (FRESH_LOOK §1, §4 step 3).

**Gotchas / decisions:**
- Don't use vertical tracks for the ion fraction. Use 10°/20°, matched by angle, with no
  quality cut.
- Archived cfgs live at `<area>/Run|runs/<run>/<subrun>/raw_daq_data/*.cfg` on EOS
  `/eos/experiment/ntof/data/x17/`. A recursive `find` over x17 is too slow; use globs.
- The CEA FEU source (`~/Feu/Firmware5/...`) is on the DAQ PC only. `daq_lxplus` was not
  answering on 2026-10-12.

**Key files & commands:**
- `design/report/FRESH_LOOK_2026-10-12.md` — the running write-up (§6–7.6).
- `design/DREAM_SETTINGS.md`, `design/dream_cfg/` — DREAM settings record, per-run CSV,
  target cfg.
- `scratchpad/2026-10-12/` (README inside) — toy model, waveform analyses,
  `wf_ion_T_scan.py`, `dream_settings_by_run.py`.
- Waveform npz files: laptop `~/x17/response_sim/stageB_w2/fresh_look_2026-10-12/`. Run the
  `wf_*.py` scripts with both npz files in the cwd, using the nTof_x17 `.venv` python (needs
  pandas).
- `python ~/PycharmProjects/nTof_x17_DAQ/dream_scripts/decode_dream_cfg.py [--sweep] <cfg|dir>`
- `design/report/deck/make_sim_status_deck.py`, `make_chain_animation.py`: site generators. Run with
  `--wf-dir <dir with both npz> --out <html>`, then publish with
  `python3 ~/PycharmProjects/dylan-cern-site/scripts/add-note.py <html> --slug <slug> --force --deploy`.
  Update the LEDGER in the deck script when FRESH_LOOK numbers change.
