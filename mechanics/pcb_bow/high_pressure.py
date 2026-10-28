#!/usr/bin/env python3
"""high_pressure.py -- how much overpressure would the measured dish need?  Four bracketing scenarios of
scenarios.py swept to 100 mbar.  Writes results/high_pressure.json."""
import json
import os
import sys
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scenarios as S

P = [0, 5, 10, 15, 20, 30, 40, 60, 80, 100]
SC = [dict(name='H_t1.3_none_ss', t=1.3e-3, backing='none', edge='ss', p_mbar=P),
      dict(name='H_t1.6_loose_ss', t=1.6e-3, backing='loose', edge='ss', p_mbar=P),
      dict(name='H_t1.6_loose', t=1.6e-3, backing='loose', p_mbar=P),
      dict(name='H_t1.6_bonded_sandwich', t=1.6e-3, backing='bonded_sandwich', p_mbar=P)]

if __name__ == '__main__':
    with Pool(4) as pool:
        res = pool.map(S.run, SC, chunksize=1)
    json.dump(res, open(os.path.join(S.RES, 'high_pressure.json'), 'w'), indent=1)
    for r in res:
        print(r['name'], ' '.join(f"{x['p_mbar']}:{x['fit_depth_mm']:.2f}/{x['dip150_mm']:.2f}" for x in r['rows']))
