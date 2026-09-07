#!/usr/bin/env python3
"""Close the amplitude ledger's last open multiplicative row: primary ionisation.

The ledger (OVERNIGHT_2026-08-10.md section 7) eliminated every candidate for
the x0.63 charge deficit except the avalanche gain, leaving two rows owing
numbers. Diffusion dilution is closed by the q_sum argument (an integral cannot
be diluted by spreading). This closes the other: does Geant4 produce the right
number of drift-gap electrons per unit track length for Ar/iC4H10 95/5?

A x1.6 error here would be glaring and would move the whole ledger, because
primary ionisation multiplies every downstream stage.

Literature for an Ar-dominant MIP mixture at NTP (Sauli, PDG gas tables):
  * total ionisation ~90-100 e-/cm  (Ar 94, iC4H10 ~195 at 5 % -> ~99 for 95/5)
  * primary clusters ~25-30 /cm     (Ar 25, iC4H10 ~84 at 5 % -> ~28 for 95/5)
  * W-value ~26 eV/pair             (Ar 26, iC4H10 23)

Read-only on the frozen Stage A cluster file. Path length is taken from each
event's own drift-gap track extent rather than assumed to be the 30 mm gap,
since the gun is inclined for some points and tracks can clip the edge.
"""
import json

import numpy as np
import uproot

SRC = "/media/dylan/data/x17/response_sim/clusters/mx17_muons_500_t0.root"
GAP_MM = 30.0

f = uproot.open(SRC)
ev = f["EventTree"].arrays(library="np")
cl = f["ClusterTree"].arrays(library="np")

n_ev = len(ev["eventID"])
print(f"source     : {SRC}")
print(f"events     : {n_ev}")

# ---- per-event totals from the EventTree summary branches -------------------
nprim = ev["nPrimDrift"].astype(float)
nclus = ev["nClusDrift"].astype(float)
edep = ev["edepDrift"].astype(float)          # eV (verified: /npc gives W ~26)

# ---- path length in the drift gap, per event, from the cluster hits ---------
eid = cl["eventID"]
z = cl["z"]
vol = cl["volume"]
# drift-gap hits only: use the volume label the Stage A writer uses
vols, counts = np.unique(vol, return_counts=True)
print(f"volumes    : {dict(zip([str(v) for v in vols], counts.tolist()))}")

drift_lab = max(zip(counts, vols))[1]         # the most populated = drift gap
m = vol == drift_lab
order = np.argsort(eid[m])
e_s, z_s = eid[m][order], z[m][order]
x_s, y_s = cl["x"][m][order], cl["y"][m][order]

path = {}
for e in np.unique(e_s):
    k = e_s == e
    if k.sum() < 2:
        continue
    dz = z_s[k].max() - z_s[k].min()
    dx = x_s[k].max() - x_s[k].min()
    dy = y_s[k].max() - y_s[k].min()
    path[int(e)] = float(np.sqrt(dx * dx + dy * dy + dz * dz))

ids = np.array(sorted(path))
L = np.array([path[i] for i in ids])          # mm
sel = L > 0.8 * GAP_MM                        # full-gap crossers only
idx = np.searchsorted(ev["eventID"], ids[sel])

Lc = L[sel] / 10.0                            # cm
npc = nprim[idx] / Lc
ncc = nclus[idx] / Lc
edc = edep[idx] / Lc

print(f"full-gap crossers: {sel.sum()} of {n_ev} "
      f"(path {np.median(L[sel]):.2f} mm median)")

print("\n=== Geant4 Stage A, per cm of track in the drift gap ===")
for lab, a, lit in (("total ionisation e-/cm", npc, "90-100"),
                    ("primary clusters /cm", ncc, "25-30"),
                    ("energy deposit keV/cm", edc / 1000.0, "~2.4-2.6")):
    print(f"  {lab:<26} median {np.median(a):8.2f}   mean {a.mean():8.2f}"
          f"   literature {lit}")

w_ev = edc / npc                              # eV per electron
print(f"  {'implied W [eV/pair]':<26} median {np.median(w_ev):8.2f}"
      f"   mean {w_ev.mean():8.2f}   literature ~26")

e_per_clus = npc / ncc
print(f"  {'electrons per cluster':<26} median {np.median(e_per_clus):8.2f}"
      f"   mean {e_per_clus.mean():8.2f}   literature ~3.5")

print("\n=== ledger verdict ===")
tot = float(np.median(npc))
lo, hi = 90.0, 100.0
if lo <= tot <= hi:
    print(f"  total ionisation {tot:.1f} e-/cm is INSIDE the literature "
          f"{lo:.0f}-{hi:.0f} -> this row CANNOT carry the x1.6 deficit")
else:
    f_need = np.mean([lo, hi]) / tot
    print(f"  total ionisation {tot:.1f} e-/cm is OUTSIDE {lo:.0f}-{hi:.0f}; "
          f"a factor {f_need:.3f} would be needed to reach the middle")

out = dict(source=SRC, n_events=int(n_ev), n_full_gap=int(sel.sum()),
           median_path_mm=float(np.median(L[sel])),
           total_ionisation_per_cm=float(np.median(npc)),
           clusters_per_cm=float(np.median(ncc)),
           edep_keV_per_cm=float(np.median(edc)),
           implied_W_eV=float(np.median(w_ev)),
           electrons_per_cluster=float(np.median(e_per_clus)))
with open("primary_ionisation_check.json", "w") as fh:
    json.dump(out, fh, indent=1)
print("\nwrote primary_ionisation_check.json")
