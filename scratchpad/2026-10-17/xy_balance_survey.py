#!/usr/bin/env python3
"""X/Y charge balance per (chamber, FEU pair, run), 2026-10-17. Separates an FEU effect from a
chamber effect on the sim's Y/X offset (FRESH_LOOK §20.2): det3 ran on FEU 7/8 (6-27) and on 3/4 (6-22).
Per event and plane: q = sum of combined_hits `amplitude` (saturation-corrected) over strips >= 100 ADC;
both planes fired, < 50 hits in total (spark veto). No M3 match: every cosmic, same mix in every run.
Also the median peak-strip amplitude per plane (the ledger observable's analogue).
Usage: xy_balance_survey.py   (reads qa_config from nTof_x17/mx_june_cosmic_qa)"""
import os, sys, glob
import numpy as np, uproot
QA = os.path.expanduser('~/PycharmProjects/nTof_x17/mx_june_cosmic_qa')
sys.path.insert(0, QA)
import qa_config
THR, VETO = 100.0, 50
seen = set()
print(f'{"key":18s} {"det":7s} {"run/subrun":58s} {"X/Y feu":7s} {"n":>6s} {"f med":>6s} {"ln qX/qY":>8s} {"pkX/pkY":>7s} {"sat frac":>8s}')
for key, c in qa_config.RUNS.items():
    files = sorted(glob.glob(os.path.join(c.combined_hits_dir, '*.root')))
    sig = (c.combined_hits_dir, tuple(c.MX17_FEUS))
    if not files or sig in seen: continue
    seen.add(sig)
    fx, fy = c.MX17_FEUS
    qx, qy, px, py, sat = [], [], [], [], []
    for fn in files:
        try:
            a = uproot.open(fn)['hits'].arrays(['eventId', 'feu', 'amplitude', 'saturated'], library='np')
        except Exception as e:
            continue
        ev, feu, amp, s = a['eventId'], a['feu'], a['amplitude'], a['saturated']
        if not (np.isin(fx, feu) and np.isin(fy, feu)): break
        nall = np.bincount(np.unique(ev, return_inverse=True)[1])
        ue, inv = np.unique(ev, return_inverse=True)
        ok_ev = nall < VETO
        m = amp >= THR
        for f_, qs, ps in ((fx, qx, px), (fy, qy, py)):
            sel = m & (feu == f_)
            q = np.bincount(inv[sel], weights=amp[sel], minlength=len(ue))
            p = np.zeros(len(ue)); np.maximum.at(p, inv[sel], amp[sel])
            qs.append(q[ok_ev]); ps.append(p[ok_ev])
        ss = np.zeros(len(ue)); np.maximum.at(ss, inv[(feu == fx) | (feu == fy)], s[(feu == fx) | (feu == fy)].astype(float))
        sat.append(ss[ok_ev])
    if not qx: continue
    qx, qy, px, py, sat = map(np.concatenate, (qx, qy, px, py, sat))
    both = (qx > 0) & (qy > 0)
    if both.sum() < 200: continue
    f = qx[both] / (qx[both] + qy[both])
    print(f'{key:18s} {c.DET_NAME:7s} {c.RUN + "/" + c.SUB_RUN:58.58s} {fx}/{fy:<5d} {both.sum():6d} {np.median(f):6.3f} '
          f'{np.median(np.log(qx[both] / qy[both])):+8.3f} {np.median(px[both]) / np.median(py[both]):7.3f} {sat[both].mean():8.3f}')
