import sys, numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,'.'); import toy_pulse as T
DT=60; Q=(.1,.25,.5,.75,.9)
def cross(w,frac,side):
    ipk=int(np.argmax(w)); lvl=frac*w[ipk]
    b=np.flatnonzero(w[:ipk+1]<lvl)
    if not len(b): return np.nan
    i=b[-1]; j=i+1
    return i+(lvl-w[i])/(w[j]-w[i]) if w[j]!=w[i] else float(i)
def job(over):
    p=dict(T.BASE); p.update(over)
    B=T.sampled_peak_strips(p,1500,seed=3,target_peak=1650)
    W=B[:,3]; W=W[W.max(1)<3500]
    l10=np.array([cross(w,.1,-1) for w in W]); l50=np.array([cross(w,.5,-1) for w in W]); l90=np.array([cross(w,.9,-1) for w in W]); pk=np.argmax(W,1)
    r=lambda x:[int(np.nanquantile(x,q)) for q in Q]
    return over, r((l50-l10)*DT), r((l90-l50)*DT), r((pk-l90)*DT), r((l90-l10)*DT)
grid=[{'f_ion':0.908},{'f_ion':0.5},{'f_ion':0.2},{'f_ion':0.0},{'t_ion':170},
      {'c_spread':0.6,'tau_s':150},{'f_ion':0.0,'c_spread':0.6,'tau_s':150},{'sigma_p':0.6},{'sigma_p':0.6,'f_ion':0.0}]
with ProcessPoolExecutor(10) as ex:
    for over,a,b,c,d in ex.map(job,grid):
        print(over); print('   10->50',a,' 50->90',b,' 90->pk',c,' rise',d)
