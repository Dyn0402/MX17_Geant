import numpy as np
from wf_shape import cross
z=np.load('waveforms_pm3.npz'); DT=60
grid=np.arange(-6,22,1.0)
def avg(W,ref):
    L=np.array([cross(r,.5,-1) for r in ref]); acc=np.zeros_like(grid); n=0
    for w,l,r in zip(W,L,ref):
        if np.isnan(l): continue
        acc+=np.interp(grid+l,np.arange(len(w)),w)/r.max(); n+=1
    return acc/n
for v in 'xy':
    print('== view',v,'   t[ns]:',' '.join(f'{int(g*DT):5d}' for g in grid[::2]))
    for leg in ('sim','data'):
        B=z[f'{leg}_{v}_wf']; pk=B[:,3]; full=np.isfinite(B).all((1,2))
        ok=(np.nanmax(B,(1,2))<3500)&full
        S=np.nansum(B[ok],1)
        a_pk=avg(pk[ok],pk[ok]); a_sum=avg(S,S); a_nb=avg(B[ok][:,[2,4]].sum(1),pk[ok])
        for name,a in (('peak strip',a_pk),('sum +-3',a_sum),('+-1 nbrs/pk',a_nb)):
            print(f'  {leg:4s} {name:12s} n={ok.sum():4d}',' '.join(f'{x:5.2f}' for x in a[::2]))
