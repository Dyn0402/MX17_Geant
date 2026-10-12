"""Dump peak-strip +-3 neighbour waveforms (position-ordered) for the T14 legs.
Same selection, seeding threshold and peak-strip choice as t14_compare.extract_view.
Fresh look 2026-10-12."""
import os, sys, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/PycharmProjects/nTof_x17/mx17_sim_wft"))
import t14_compare as T
from wft.io import FeuReader, strip_position_map

B = os.path.expanduser("~/x17/response_sim/stageB_w2")
O = os.path.join(B, "fresh_look_2026-10-12")
legs = {"sim": (f"{B}/w2_rho2M_TRANSPFIX/default", f"{O}/events_transpfix_noq.parquet"),
        "data": (f"{B}/data_leg_satdet3", f"{O}/events_data_noq.parquet")}
summ = json.load(open(f"{O}/t14_compare_transpfix_noq/t14_summary.json"))
wins = summ["theta_windows_deg"]
cfg = T._leg_cfg(legs["sim"][0], "sat_det3")
pos_maps = strip_position_map(cfg)
out = {}
for leg, (d, evp) in legs.items():
    ev = pd.read_parquet(evp)
    for view, feu in (("x", cfg.MX17_FEU_X), ("y", cfg.MX17_FEU_Y)):
        ids = T._cap(T._sel_ids(ev, view, None if leg == "sim" else wins[view]))
        pos = pos_maps[feu]; valid = ~np.isnan(pos)
        order = np.flatnonzero(valid)[np.argsort(pos[valid])]
        rank = {c: i for i, c in enumerate(order)}
        W, E, N = [], [], []
        for path in sorted(glob.glob(os.path.join(d, "decoded_root", f"*_{feu:02d}.root"))):
            rdr = FeuReader(path)
            noise = np.where(rdr.noise > 0, rdr.noise, np.inf)
            want = ids & set(int(e) for e in rdr.event_ids)
            if not want:
                continue
            for eid, _f, wfm in rdr.iter_events(want):
                amp = wfm.max(axis=1); over = valid & (amp / noise >= T.SIGMA)
                if not over.any():
                    continue
                c = np.flatnonzero(over); pk = int(c[np.argmax(amp[c])])
                if pk not in rank:
                    continue
                r = rank[pk]
                blk = np.full((7, wfm.shape[1]), np.nan)
                for k in range(-3, 4):
                    if 0 <= r + k < len(order):
                        blk[k + 3] = wfm[order[r + k]]
                W.append(blk); E.append(int(eid)); N.append(float(np.median(rdr.noise)))
        out[f"{leg}_{view}_wf"] = np.array(W, np.float32)
        out[f"{leg}_{view}_eid"] = np.array(E)
        print(leg, view, len(W), flush=True)
np.savez_compressed(f"{O}/waveforms_pm3.npz", **out)
print("done")
