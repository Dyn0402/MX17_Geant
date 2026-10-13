#!/usr/bin/env python3
"""Gain ladders of the May 2026 "det3" (bulked 5-5-26; NOT the June det3 of T14)
in two gases, same chamber, one day apart:

  mx17_det3_HV_Scan_5-5-26       Ar/iC4H10 95/5   resist 460-530 V, drift 900 V
  mx17_det3_ArCF4_HV_Scan_5-6-26 Ar/CF4 90/10     resist 530-670 V, drift 900 V

Same observable and selection as nTof_x17 mx17_sim_wft/hv_slope/extract.py
(threshold-free peak-strip max of the raw waveform, M3-good rays only, so the
event population is voltage independent); its view_rows() is reused unchanged.
nTof_x17 is imported read-only.

    ~/PycharmProjects/nTof_x17/.venv/bin/python extract_may.py RUN SUBRUN OUT.parquet
"""
import os, re, sys
import numpy as np, pandas as pd
NT = os.path.expanduser('~/PycharmProjects/nTof_x17')
for p in (NT, f'{NT}/mx_june_cosmic_qa', f'{NT}/cosmic_bench_analysis', f'{NT}/mx17_sim_wft/hv_slope'):
    sys.path.insert(0, p)
from qa_config import _Config, setup_paths
setup_paths()
import extract as X
from wft.io import strip_position_map

BASE = '/home/dylan/x17/cosmic_bench/det_3_old/'
X.DET_Z = 232.0                      # mx17_1 det_center z in these runs' config

run, sub, out = sys.argv[1:4]
volt = int(re.search(r'resist_(\d+)V', sub).group(1))
cfg = _Config('may', run, sub, feus=[3, 4], det_z=232.0, base_path=BASE, det_name='mx17_1')
pos = strip_position_map(cfg)
refs = X.ray_positions(cfg)
print(f'[{sub}] {volt} V  {len(refs):,} M3-good rays', flush=True)
rows = []
for view, feu in (('x', 3), ('y', 4)):
    r = X.view_rows(os.path.join(BASE, run, sub, 'decoded_root'), feu, pos[feu], refs, sub, volt, run, view)
    print(f'   {view}: {len(r):,} events', flush=True)
    rows += r
pd.DataFrame(rows).to_parquet(out)
print('wrote', out)
