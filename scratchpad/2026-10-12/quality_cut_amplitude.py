import pandas as pd, numpy as np
d=pd.read_parquet('events_data.parquet'); s=pd.read_parquet('events_transpfix.parquet')
for v in 'xy':
    base={}
    th0=s.loc[s[v+'_ok']&s[v+'_quality_ok'],v+'_theta_deg'].median()
    def sel(e):
        return e[e[v+'_ok']&~e.spark&((e[v+'_theta_deg']-th0).abs()<3)].copy()
    S,D=sel(s),sel(d)
    for g in (S,D): g['nchi']=g[v+'_chi2']/g[v+'_dof']/(g[v+'_q_sum']/1000)**2
    print(f'== {v}  n sim {len(S)} data {len(D)}')
    print('   norm-chi2 quantiles p50/p90/p99 sim', np.round(S.nchi.quantile([.5,.9,.99]).values,2),' data',np.round(D.nchi.quantile([.5,.9,.99]).values,2))
    print('   n_candidates>1 frac sim %.3f data %.3f'%((S[v+'_n_candidates']>1).mean(),(D[v+'_n_candidates']>1).mean()))
    cuts={'quality_ok (as frozen)':lambda g:g[v+'_quality_ok'],
          'no quality cut':lambda g:np.ones(len(g),bool),
          'norm-chi2 < data p95':lambda g:g.nchi<D.nchi.quantile(.95),
          'no cut, n_candidates==1':lambda g:g[v+'_n_candidates']==1}
    for k,f in cuts.items():
        ms,md=S[f(S)][v+'_q_sum'],D[f(D)][v+'_q_sum']
        r=ms.median()/md.median()
        # bootstrap error
        rng=np.random.default_rng(1); bs=[np.median(rng.choice(ms.values,len(ms)))/np.median(rng.choice(md.values,len(md))) for _ in range(300)]
        print(f'   {k:28s} pass sim {f(S).mean():.3f} data {np.mean(f(D)):.3f}   q_sum ratio {r:.4f} ± {np.std(bs):.4f}  deficit x{1/r:.3f}')
