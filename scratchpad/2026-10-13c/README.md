# scratchpad 2026-10-13c — every chamber x gas in the archive (FRESH_LOOK §13)

Run with `~/PycharmProjects/nTof_x17/.venv/bin/python` from this directory (nTof_x17 imported read-only).
Per-event products (parquet + `_wf.npy`, ~1 GB) were written to the job tmp dir, not kept;
regenerate with `extract_bench.py` (minutes per subrun, run 10 in parallel).

Data on the desktop: det0 April runs pulled from EOS `cosmic_bench/det_0_old/` to
`/home/dylan/det0_old/` (NVMe; `~/x17` -> /media/dylan/data is 100 % full). Old det3 in
`~/x17/cosmic_bench/det_3_old/`. Two det0 He final_505 files are 0-byte / unreadable on EOS itself.

- `extract_bench.py BASE RUN SUBRUN DET_Z OUT [ALIGN]` — threshold-free peak-strip max on M3-good
  rays, +-2-strip waveforms, M3 slopes; decoded_root or filtered_root ([512][32] channel-major).
  det0 needs ALIGN `x:215:1:ref_y,y:191:-1:ref_x` (track window + noisy-channel mask); old det3
  alignment is `x:229:1:ref_y,y:153:-1:ref_x` (not needed there, gain is high).
- `ladders_all.py` -> `ladders_all.csv`, `fits_all.csv`, `fits_all.log`: quantile ladders, exp /
  Garfield-slope feedback / free-slope feedback fits, all sets incl. June det3 (§9's csv).
- `width_vs_v.py` -> `width_vs_v.log`: relative width of the amplitude distribution vs V.
- `rise_all.py` -> `rise_all.log`: peak-strip 10-90 % rise / FWHM / delay by track inclination.
- `point_rise.py` -> `point_rise.log`: toy point/short-column rise vs T_ion, f_ion.
- `he_template_fit.py` -> `he_template_fit.log`: det0 Ar/iso = He template (x) ion kernel (shaper-free).
- `peak_vs_tion.py` -> `peak_vs_tion.log`: peak per charge vs T_ion (toy, T14 geometry).
- `plot_s13.py` -> `design/report/figs/s13_ladders_and_ion_2026-10-13.png`.
- `magboltz/drift_table.py` — drift v/D for He/C2H6, Ar/CF4, Ar/iso; running on lxplus
  `~/work/mx17_gas_oct13/` (the desktop's PyROOT is broken: ROOT 6.30 built for Python 3.10).
