"""Ion birth-height distribution: how much does it move the ion-current timing?
Uses the meshfield S3 alpha_z_hist, T6 map E_z(z), psi_readout and funnel endpoints
(same inputs as response/meshcell/ion_template_check.py). Run with the nTof_x17 venv."""
import json,numpy as np,pandas as pd
R='/home/dylan/CLionProjects/MX17_Geant/response/'
df=pd.read_csv(R+'meshcell/meshfield_production.txt',sep=r"\s+",comment="#",header=None,names="x y z ex ey ez v flag".split())
g=df[df.flag==1].groupby("z")["ez"].mean(); ez=np.column_stack([g.index.values*1e4,g.values])
psi=json.load(open(R+'meshcell/psi_readout.json')); zr=np.array(psi['psi_profile']['z_above_anode_um']); pv=np.array(psi['psi_profile']['psi']); za=psi['geometry']['z_anode_um']
print("E_z(h) [kV/cm]:", [(h, round(np.interp(za+h,ez[:,0],ez[:,1])/1e3,1)) for h in (0,2,5,10,20,40,75,110,130,140,150)])
cal=json.load(open(R+'avalanche/aval_calib_meshfield_pooled.json'))['point']; zh=cal['alpha_z_hist']
c=np.array(zh['counts'],float); e=np.array(zh['edges'],float); zc=.5*(e[1:]+e[:-1]); w=c/c.sum()
fun=json.load(open(R+'meshcell/funnel_ion_endpoints.json'))['ion_endpoints']
zend=np.median([r['ze_above_anode_um'] for r in fun['rows'] if r['outcome']=='absorbed_on_mesh'])
n=101325/(1.380649e-23*293)/1e6; sc=2.6868e19/n
tab=np.loadtxt('/home/dylan/Software/garfield/Data/IonMobility_Ar+_Ar.txt',comments='#')[:,:2]
grid=np.linspace(0,zend,3001); dz=grid[1]-grid[0]; psig=np.interp(grid,zr,pv)
E=np.maximum(np.interp(za+grid,ez[:,0],ez[:,1]),1)
def run(K0fun,births,wts):
    v=K0fun(E/n/1e-17)*sc*E*1e4/1e9
    t=np.concatenate([[0],np.cumsum(dz/(.5*(v[1:]+v[:-1])))])
    tb=np.linspace(0,700,3501); q=np.zeros(3500)
    for h0,wi in zip(births,wts):
        s=grid>=h0
        tt=t[s]-t[s][0]; dq=-np.diff(psig[s])
        np.add.at(q,np.clip(np.searchsorted(tb,tt[1:])-1,0,3499),dq*wi)
    cc=np.cumsum(q)/q.sum(); return [round(float(np.interp(f,cc,tb[1:])),0) for f in (.1,.5,.9,.99)]
ar=lambda td: np.interp(td,tab[:,0],tab[:,1]); flat=lambda K: (lambda td: K+0*td)
zm=(w*zc).sum()
for name,K in (("Ar+ table",ar),("flat K0=2.05",flat(2.05)),("flat K0=2.25",flat(2.25))):
    print(f"{name:14s} q10/50/90/99 ns:  all@0 {run(K,[0],[1])}  all@mean {run(K,[zm],[1])}  full hist {run(K,zc,w)}")
# what if gain differed: shift decay length
for L in (10,15,20,30):
    ww=np.exp(-zc/L); ww/=ww.sum()
    print(f"exp decay {L} um, flat 2.05: {run(flat(2.05),zc,ww)}, f_ion(pp)={1-(ww*zc).sum()/150:.3f}")
