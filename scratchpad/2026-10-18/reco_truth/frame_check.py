"""Truth-frame check for sim_calib.py: ref-pinned chi2 vs a shift of the pinned mesh position and vs the tan sign.
The minimum must sit at zero shift with the stated sign, on every view and angle."""
import sys, numpy as np
sys.argv = [sys.argv[0]] + sys.argv[1:]
from sim_calib import *
from wft import model as wm
pts = [(f'/media/dylan/data/claude_tmp/v3rw/v3k265b45_{a}', f'/media/dylan/data/claude_tmp/stageA_e17/mx17_muons_3k_ang_{a}_e17_t0.root', *ANG[a])
       for a in ('th00', 'thx10', 'thx20', 'thy10', 'thy20')]
ev = build_cache(pts, 12)
cal = CalibrationBundle.load(R06); cal.v_drift = V_SIM
wm.use_calibration(cal); hyper = cal.hyper
def chi(e, plane, dp, sgn):
    P = e[plane]; W, noise, pos, sat = wm.prep_plane(P, plane)
    w = sgn * e[f'tan_{plane}'] * V_SIM * 1e-3; p0 = e[f'ref_mesh_{plane}'] + dp
    return min(wm.chi2_plane(plane, W, noise, pos, sat, p0, w, float(t), hyper)[0] for t in np.arange(150, 900, 15.0))
for ip in range(len(pts)):
    ks = [k for k in ev if k // 100000 == ip]
    for plane in 'xy':
        shifts = np.array([-1.0, -0.5, -0.25, 0, 0.25, 0.5, 1.0])
        c = np.array([[chi(ev[k], plane, d, 1) for d in shifts] for k in ks]).sum(0)
        cf = sum(chi(ev[k], plane, 0, -1) for k in ks)
        print(os.path.basename(pts[ip][0]), plane, 'best shift', shifts[np.argmin(c)], ' chi2(rel to 0):', np.round(c / c[3], 3), ' flipped-sign/right:', round(cf / c[3], 3), flush=True)
