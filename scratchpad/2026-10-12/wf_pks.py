import numpy as np, pandas as pd
from wf_shape import cross
z=np.load('waveforms_pm3.npz'); DT=60
for leg,pq in (('sim','events_transpfix_noq.parquet'),('data','events_data_noq.parquet')):
    ev=pd.read_parquet(pq).set_index('event_id')
    for v in 'xy':
        W=z[f'{leg}_{v}_wf'][:,3]; eid=z[f'{leg}_{v}_eid']; ok=W.max(1)<3500
        W=W[ok]; eid=eid[ok]
        l10=np.array([cross(w,.1,-1) for w in W]); l50=np.array([cross(w,.5,-1) for w in W])
        d=pd.DataFrame(dict(edge=(l50-l10)*DT,pks=np.argmax(W,1),l10=l10,t0=ev.loc[eid,f'{v}_t0'].values,w0=W[:,0],w1=W[:,1]))
        g=d.groupby(d.pks.clip(5,20))
        print(f'{leg} {v}: pks  : ' + ' '.join(f'{k:4d}' for k in g.groups))
        print(f'       n    : ' + ' '.join(f'{x:4d}' for x in g.size()))
        print(f'       edge : ' + ' '.join(f'{x:4.0f}' for x in g.edge.median()))
        print(f'       l10  : ' + ' '.join(f'{x:4.1f}' for x in g.l10.median()))
        print(f'       t0   : ' + ' '.join(f'{x:4.0f}' for x in g.t0.median()))
