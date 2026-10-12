import numpy as np, pandas as pd
from wf_shape import cross
z=np.load('waveforms_pm3.npz'); DT=60
ev=pd.read_parquet('events_data_noq.parquet').set_index('event_id')
for v in 'xy':
    W=z[f'data_{v}_wf'][:,3]; eid=z[f'data_{v}_eid']; a=W.max(1); ok=a<3500
    l10=np.array([cross(w,.1,-1) for w in W]); l50=np.array([cross(w,.5,-1) for w in W])
    e=(l50-l10)*DT
    df=pd.DataFrame(dict(edge=e,amp=a,ok=ok,pks=np.argmax(W,1),
        th=ev.loc[eid,f'{v}_theta_deg'].values,p0=ev.loc[eid,f'{v}_p0'].values,t0=ev.loc[eid,f'{v}_t0'].values,
        nst=ev.loc[eid,f'{v}_n_strips'].values,q=ev.loc[eid,f'{v}_q_sum'].values))
    d=df[df.ok]
    print('== data',v,'n',len(d))
    print('  edge hist (ns):',np.histogram(d.edge.dropna(),bins=[0,60,90,120,150,180,210,240,300,400,600,2000])[0])
    for col in ('amp','th','p0','pks','nst'):
        b=pd.qcut(d[col].rank(method='first'),4,labels=False)
        print(f'  median edge by {col:4s} quartile:',[int(x) for x in d.groupby(b).edge.median()],
              ' (bin medians', [round(x,1) for x in d.groupby(b)[col].median()],')')
