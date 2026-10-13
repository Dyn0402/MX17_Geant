#!/usr/bin/env python3
"""Magboltz drift v, D_L, D_T over 50-1000 V/cm for one mixture at Saclay pressure
(same settings as nTof_x17 garfield_sim/mm_drift_velocity_table.py).
  python3 drift_table.py OUT.json gas1 frac1 gas2 frac2 [...]"""
import sys, json, ctypes, numpy as np
import ROOT; ROOT.gROOT.SetBatch(True); import Garfield  # noqa
out, comp = sys.argv[1], sys.argv[2:]
g = ROOT.Garfield.MediumMagboltz()
args = []
for i in range(0, len(comp), 2): args += [comp[i], float(comp[i + 1])]
g.SetComposition(*args); g.SetTemperature(293.15); g.SetPressure(745.83)
g.SetFieldGrid(50., 1000., 16, True); g.GenerateGasTable(5)
rows = []
for e in np.logspace(np.log10(50), np.log10(1000), 40):
    vx, vy, vz = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    g.ElectronVelocity(0., 0., -e, 0., 0., 0., vx, vy, vz)
    dl, dt = ctypes.c_double(), ctypes.c_double()
    g.ElectronDiffusion(0., 0., -e, 0., 0., 0., dl, dt)
    rows.append(dict(E_Vcm=float(e), v_um_per_ns=abs(vz.value) * 1e4, dL_sqrtcm=dl.value, dT_sqrtcm=dt.value))
json.dump(dict(composition=comp, pressure_torr=745.83, rows=rows), open(out, 'w'), indent=1)
print('wrote', out)
