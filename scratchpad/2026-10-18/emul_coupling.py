import sys, numpy as np
sys.path.insert(0, '/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_shape import cross
DT=60.; grid=np.arange(-4,20,.25)
def aligned(W):
    W=W[W[:,3].max(1)<3500]; C,N,F=[],[],[]
    for w in W:
        l=cross(w[3],.5,-1)
        if np.isnan(l): continue
        pk=w[3].max(); x=np.arange(w.shape[1])
        C.append(np.interp(grid+l,x,w[3]/pk)); N.append(.5*(np.interp(grid+l,x,w[2]/pk)+np.interp(grid+l,x,w[4]/pk)))
        a=np.clip(w.max(1),0,None); F.append(a[3]/a.sum())
    return np.mean(C,0),np.mean(N,0),np.mean(F)
def cpl(W,k,tau):
    L=np.zeros_like(W,dtype=float); W=W.astype(float)
    L[:,1:-1]=W[:,:-2]+W[:,2:]-2*W[:,1:-1]
    return W+k*L+(tau/DT)*np.gradient(L,axis=2)
z=np.load(sys.argv[1]); tpk=lambda y:grid[np.argmax(y)]*DT
for v in 'xy':
    dc,dn,df=aligned(z[f'data_{v}_wf'])
    print(f"{v}: data delay {tpk(dn)-tpk(dc):+.0f} nb/c {dn.max():.3f} c_t {tpk(dc):.0f} pk/sum7 {df:.3f}  nb(t=0) {dn[16]:.3f}")
    for k in (0,.03,.06):
        for tau in (0,10,20,30,45):
            sc,sn,sf=aligned(cpl(z[f'sim_{v}_wf'],k,tau))
            print(f"   k {k:.2f} tau {tau:2d}: delay {tpk(sn)-tpk(sc):+4.0f} nb/c {sn.max():.3f} c_t {tpk(sc):4.0f} pk/sum7 {sf:.3f} nb(t=0) {sn[16]:.3f} rms {np.sqrt(np.mean((sn-dn)**2)):.4f}")
print("residual curves (data - sim) at 60 ns steps from -240 ns")
for v in 'xy':
    dc,dn,_=aligned(z[f'data_{v}_wf'])
    for k,tau in ((0,0),(.05,15)):
        sc,sn,_=aligned(cpl(z[f'sim_{v}_wf'],k,tau))
        print(f"{v} k{k} t{tau} nb  "+" ".join(f"{x:+.3f}" for x in (dn-sn)[::4][:18]))
        print(f"{v} k{k} t{tau} c   "+" ".join(f"{x:+.3f}" for x in (dc-sc)[::4][:18]))
