import numpy as np
from wf_shape import cross
z=np.load('waveforms_pm3.npz'); DT=60
Q=(.1,.25,.5,.75,.9)
for v in 'xy':
    print('== view',v)
    for leg in ('sim','data'):
        W=z[f'{leg}_{v}_wf'][:,3]; W=W[W.max(1)<3500]
        c=lambda f:np.array([cross(w,f,-1) for w in W])
        l10,l50,l90=c(.1),c(.5),c(.9); pk=np.argmax(W,1)
        for name,x in (('10->50%',(l50-l10)*DT),('50->90%',(l90-l50)*DT),('90%->peak',(pk-l90)*DT),('10->90 (rise)',(l90-l10)*DT)):
            print(f'  {leg:4s} {name:14s}',[int(np.nanquantile(x,q)) for q in Q])
