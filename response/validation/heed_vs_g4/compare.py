"""Heed vs Geant4 bare-slab primaries (2026-10-15).
    python3 compare.py <heed run dir> <slab run dir>
Run files: heed_mip -> heed_<gas>[_30mm].txt; slab -> g4_<gas>[_mx17]_<gap>_<phys>.txt (see README.md)."""
import sys
import numpy as np
HD, GD = sys.argv[1], sys.argv[2]
def st(n): return f"mean {n.mean():7.2f} median {np.median(n):6.1f}"
def mpv(x, bw):
    h,b=np.histogram(x,bins=np.arange(0,np.quantile(x,.9),bw)); return b[h.argmax()]+bw/2
rows=[]
for g,gap,hf in [("ArIso","4mm","heed_ArIso"),("NeIso","4mm","heed_NeIso"),("ArCF4Iso","4mm","heed_ArCF4Iso"),("ArIso","30mm","heed_ArIso_30mm")]:
    h=np.loadtxt(f"{HD}/{hf}.txt"); hn,he=h[:,1],h[:,2]; Wh=he.sum()/hn.sum()
    bw=25 if gap=="4mm" else 100
    print(f"== {g} {gap}  Heed W_eff {Wh:.2f} eV")
    print(f"  {'Heed':16s} n {st(hn)}  | edep {st(he)} mpv {mpv(he,bw):5.0f}")
    for v in ["opt4","pai","paiphot","mx17_opt4" ,"mx17_pai"]:
        tag = f"g4_{g}_{gap}_{v}" if not v.startswith("mx17") else f"g4_{g}_mx17_{gap}_{v[5:]}"
        try: d=np.loadtxt(f"{GD}/{tag}.txt")
        except OSError: continue
        e,n=d[:,0],d[:,1]; W=open(f"{GD}/{tag}.txt").readline().split("W=")[1].strip()
        print(f"  G4 {v:10s}W={W:6s} n {st(n)} ({n.mean()/hn.mean():.3f}x mean, {np.median(n)/np.median(hn):.3f}x med) | edep {st(e)} mpv {mpv(e,bw):5.0f} ({np.median(e)/np.median(he):.3f}x med)")
