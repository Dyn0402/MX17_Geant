import pandas as pd, numpy as np, sys
Q=[.05,.1,.25,.5,.75,.9]
def load(dirn,leg,v): return pd.read_parquet(f'{dirn}/wf_{leg}_{v}.parquet')
for v in 'xy':
    print(f'== view {v}')
    for dirn in sys.argv[1:]:
        s,d=load(dirn,'sim',v),load(dirn,'data',v)
        line=f'  {dirn:8s} n sim {len(s)} data {len(d)}  sat sim {s.saturated.mean():.3f} data {d.saturated.mean():.3f}'
        rs,rd=s.rise_ns.dropna(),d.rise_ns.dropna()
        line+=f'\n     rise median sim {rs.median():.1f} data {rd.median():.1f} ratio {rs.median()/rd.median():.3f}; offsets sim-data at q='+' '.join(f'{(rs.quantile(q)-rd.quantile(q)):.0f}' for q in Q)
        fs,fd=s.fwhm_ns.dropna(),d.fwhm_ns.dropna()
        line+=f'\n     fwhm ratio {fs.median()/fd.median():.3f}   peak_amp median ratio {s.peak_amp.median()/d.peak_amp.median():.4f}  q_event {s.q_event.median()/d.q_event.median():.4f}  undefined-rise data {d.rise_ns.isna().mean()-d.saturated.mean():.3f}'
        print(line)
