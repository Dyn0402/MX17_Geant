import sys, numpy as np
sys.path.insert(0,'/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_shape import cross
import toy_pulse as T
DT=60; grid=np.arange(-4,16,1.0)
def avg(W):
    W=W[(W.max(1)<3500)]; acc=np.zeros_like(grid); n=0
    for w in W:
        l=cross(w,.5,-1)
        if np.isnan(l): continue
        acc+=np.interp(grid+l,np.arange(len(w)),w/w.max()); n+=1
    return acc/n
def _m():
  pass
if __name__=="__main__":
    za=np.load('waveforms_angled_pm3.npz')
    for tag,th in (('thx20',20.0),('thx10',10.7)):
        print(f'== {tag}  t[ns] from leading 50%:',' '.join(f'{int(g*DT):5d}' for g in grid))
        rows=[('data',za[f'{tag}_data_wf'][:,3]),('chain',za[f'{tag}_sim_wf'][:,3])]
        for name,over in (('toy p1',{}),('toy p5',{'psi_p':5.0}),('toy p8',{'psi_p':8.0}),('toy f0',{'f_ion':0.0})):
            rows.append((name,T.sampled_peak_strips({**T.BASE,**over,'theta':th},1500,seed=2,target_peak=2000)[:,3]))
        for name,W in rows: print(f'  {name:6s}',' '.join(f'{x:5.2f}' for x in avg(W)))
