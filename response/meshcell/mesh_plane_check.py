#!/usr/bin/env python3
"""
mesh_plane_check.py — why the T6 map's amp-gap field is 31.0 kV/cm and not
V/g = 490 V / 150 um = 32.7 kV/cm.

The 150 um gap is measured from the ESL to the mesh UNDERSIDE (the pillar
height), but a woven mesh with 48 um openings between 19 um wires is not a
solid plate: the amplification field leaks up into the openings, so the
gas-averaged potential at the underside is ~25 V, not 0. Extrapolating the
(perfectly linear) bulk potential, it reaches the mesh's 0 V about 8 um
INSIDE the weave, ~1.5 um above the lower wire layer's axis. That plane is
the mesh's equivalent electrode, the effective gap is ~158 um, and
490 / 158 um = 31.0 kV/cm. A single wire grid's textbook equivalent-plane
shift, (a / 2 pi) ln(a / 2 pi r) = 1.2 um past the wire axes for a = 67 um,
r = 9.5 um, says the same thing. neBEM gives the bulk field to 1.5 %
(meshfield_production.json, gate S2).

(The "38.3 kV/cm" quoted in the 2026-10-15 ion review was 490 V / 128 um,
the generic bulk-Micromegas gap; MX17's is 150 um, NEEDED_INPUTS §4.)

    python3 mesh_plane_check.py meshfield_production.txt
"""
import numpy as np, sys
a=np.loadtxt(sys.argv[1],comments='#')
z=a[:,2]*1e4; V=a[:,6]; Ez=a[:,5]; f=a[:,7]
zs=np.unique(z)
rows=[]
for zz in zs:
    s=(z==zz)&(f==1)
    if s.sum()==0: rows.append((zz,np.nan,np.nan,0)); continue
    rows.append((zz,V[s].mean(),Ez[s].mean(),s.sum()/(z==zz).sum()))
r=np.array(rows)
for zz in (-168,-150,-100,-60,-40,-30,-25,-20,-18,-15,-10,-5,0,5,10,15,18,20,30,50,100,200):
    i=np.argmin(abs(r[:,0]-zz)); print(f"z={r[i,0]:7.1f} um  <V>={r[i,1]:8.3f} V  <Ez>={r[i,2]/1e3:7.3f} kV/cm  gasfrac={r[i,3]:.2f}")
b=(r[:,0]>-160)&(r[:,0]<-60)
p=np.polyfit(r[b,0],r[b,1],1)
print("amp-bulk fit: E=%.1f V/cm, V(anode=-168.5)=%.2f, V=0 at z=%.2f um (mesh underside -18.5, wire axes -9/+9)"%(-p[0]*1e4,np.polyval(p,-168.5),-p[1]/p[0]))
print("effective gap = %.1f um ; 490/150 = %.1f kV/cm"%(-p[1]/p[0]+168.5, 490/150*10))
