import sys, numpy as np, pandas as pd
sys.path.insert(0,'/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_shape import cross
import toy_pulse as T
DT=60
def table(W):
    W=W[W.max(1)<3500]
    l10=np.array([cross(w,.1,-1) for w in W]); l50=np.array([cross(w,.5,-1) for w in W]); l90=np.array([cross(w,.9,-1) for w in W])
    pk=np.argmax(W,1).astype(float)
    return pd.DataFrame(dict(delay=(pk-l10)*DT, edge=(l50-l10)*DT, rise=(l90-l10)*DT))
def _main():
    z=np.load("waveforms_pm3.npz")
    bins=[0,150,250,350,450,550,700,900,1500]
    rows={}
    for v in 'xy':
        for leg in ('sim','data'): rows[(v,leg)]=table(z[f'{leg}_{v}_wf'][:,3])
    for name,over in (('toy f0.91',{}),('toy f0',{'f_ion':0.0})):
        B=T.sampled_peak_strips({**T.BASE,**over},2500,seed=5,target_peak=2000)
        rows[('x',name)]=table(B[:,3])
    for v in 'xy':
        print(f'== view {v}: delay bin (ns, peak after 10% onset) ->', [f'{a}-{b}' for a,b in zip(bins[:-1],bins[1:])])
        for leg in (['sim','data','toy f0.91','toy f0'] if v=='x' else ['sim','data']):
            d=rows[(v,leg)]; c=pd.cut(d.delay,bins)
            frac=d.groupby(c,observed=False).size()/len(d); e=d.groupby(c,observed=False).edge.median()
            print(f'  {leg:9s} frac  ',' '.join(f'{x:5.2f}' for x in frac))
            print(f'  {"":9s} edge  ',' '.join(f'{x:5.0f}' for x in e))
if __name__=="__main__":
    _main()
