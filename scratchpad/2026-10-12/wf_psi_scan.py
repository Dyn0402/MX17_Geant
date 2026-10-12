import sys, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,'/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_matched import table
import toy_pulse as T
Q=(.1,.25,.5,.75,.9); bins=[0,250,350,450,550,700,900]
def job(a):
    over,th=a
    B=T.sampled_peak_strips({**T.BASE,**over,'theta':th},2000,seed=11,target_peak=2000)
    return a, table(B[:,3])
z0=np.load('waveforms_pm3.npz'); za=np.load('waveforms_angled_pm3.npz')
data={0:table(za['th00_data_wf'][:,3]),10.7:table(za['thx10_data_wf'][:,3]),20.0:table(za['thx20_data_wf'][:,3])}
cfgs=[{'psi_p':1.0},{'psi_p':2.0},{'psi_p':3.0},{'psi_p':5.0},{'psi_p':8.0},{'f_ion':0.0}]
jobs=[(c,th) for c in cfgs for th in (0,10.7,20.0)]
with ProcessPoolExecutor(12) as ex: res=list(ex.map(job,jobs))
def summ(d,th):
    s=f'rise q {[int(d.rise.quantile(q)) for q in Q]}'
    if th==0:
        c=pd.cut(d.delay,bins); s+=f'  early(<350) {np.mean(d.delay<350):.2f}  matched edge {[int(x) if x==x else -1 for x in d.groupby(c,observed=False).edge.median()]}'
    return s
for th in (0,10.7,20.0):
    print(f'== theta {th}:  data  {summ(data[th],th)}')
    for (c,t),d in res:
        if t==th: print(f'   {str(c):18s}  {summ(d,th)}')
