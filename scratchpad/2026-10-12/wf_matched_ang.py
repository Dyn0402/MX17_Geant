import sys, numpy as np, pandas as pd
sys.path.insert(0,'/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_matched import table
import toy_pulse as T
z=np.load('waveforms_angled_pm3.npz')
bins=[0,150,250,350,450,550,700,900]
Q=(.1,.25,.5,.75,.9)
for tag,th in (('th00',0),('thx10',10.7),('thx20',20.0)):
    print(f'===== {tag}   delay bins', [f'{a}-{b}' for a,b in zip(bins[:-1],bins[1:])])
    sets=[('sim',table(z[f'{tag}_sim_wf'][:,3])),('data',table(z[f'{tag}_data_wf'][:,3]))]
    for name,f in (('toy f0.91',0.908),('toy f0',0.0)):
        B=T.sampled_peak_strips({**T.BASE,'f_ion':f,'theta':th},2000,seed=11,target_peak=2000)
        sets.append((name,table(B[:,3])))
    for name,d in sets:
        c=pd.cut(d.delay,bins)
        print(f'  {name:9s} n={len(d):4d} rise q {[int(d.rise.quantile(q)) for q in Q]}  edge q {[int(d.edge.quantile(q)) for q in Q]}')
        print(f'  {"":9s} frac  ',' '.join(f'{x:5.2f}' for x in d.groupby(c,observed=False).size()/len(d)))
        print(f'  {"":9s} edge  ',' '.join(f'{x:5.0f}' for x in d.groupby(c,observed=False).edge.median()))
