"""Peak-strip +-3 waveforms for the angled-ladder legs and the matching data,
no quality cut, data window = sim median +- HW. Fresh look 2026-10-12."""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/PycharmProjects/nTof_x17/mx17_sim_wft"))
import t14_compare as T
from wft.io import FeuReader, strip_position_map
B = os.path.expanduser("~/x17/response_sim/stageB_w2"); O = f"{B}/fresh_look_2026-10-12"
HW = 3.0
def grab(d, ids, feu, pos):
    valid = ~np.isnan(pos); order = np.flatnonzero(valid)[np.argsort(pos[valid])]
    rank = {c: i for i, c in enumerate(order)}; W, E = [], []
    for path in sorted(glob.glob(os.path.join(d, "decoded_root", f"*_{feu:02d}.root"))):
        rdr = FeuReader(path); noise = np.where(rdr.noise > 0, rdr.noise, np.inf)
        want = ids & set(int(e) for e in rdr.event_ids)
        if not want: continue
        for eid, _f, wfm in rdr.iter_events(want):
            amp = wfm.max(axis=1); over = valid & (amp / noise >= T.SIGMA)
            if not over.any(): continue
            c = np.flatnonzero(over); pk = int(c[np.argmax(amp[c])])
            if pk not in rank: continue
            r = rank[pk]; blk = np.full((7, wfm.shape[1]), np.nan)
            for k in range(-3, 4):
                if 0 <= r + k < len(order): blk[k + 3] = wfm[order[r + k]]
            W.append(blk); E.append(int(eid))
    return np.array(W, np.float32), np.array(E)
ev_d = pd.read_parquet(f"{O}/events_data_noq.parquet")
out = {}
for tag in ("th00", "thx10", "thx20"):
    sd = f"{B}/w2_rho2M_ang/{tag}"
    ev_s = pd.read_parquet(glob.glob(f"{sd}/events_ang_{tag}.parquet")[0])
    for v in "xy": ev_s[v + "_quality_ok"] = True
    cfg = T._leg_cfg(sd, "sat_det3"); pm = strip_position_map(cfg)
    v, feu = "x", cfg.MX17_FEU_X
    m = float(ev_s[f"{v}_theta_deg"].dropna().median()); win = (m - HW, m + HW)
    ids_s = T._cap(T._sel_ids(ev_s, v)); ids_d = T._cap(T._sel_ids(ev_d, v, win))
    for leg, d, ids in (("sim", sd, ids_s), ("data", f"{B}/data_leg_satdet3", ids_d)):
        W, E = grab(d, ids, feu, pm[feu])
        out[f"{tag}_{leg}_wf"] = W; out[f"{tag}_{leg}_eid"] = E
        print(tag, leg, "window", [round(x, 1) for x in win], len(W), flush=True)
np.savez_compressed(f"{O}/waveforms_angled_pm3.npz", **out); print("done")
