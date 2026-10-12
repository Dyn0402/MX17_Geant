import numpy as np
z=np.load('waveforms_pm3.npz'); DT=60
def cross(w,frac,side):
    ipk=int(np.argmax(w)); a=w[ipk]; lvl=frac*a
    if side<0:
        b=np.flatnonzero(w[:ipk+1]<lvl)
        if not len(b): return np.nan
        i=b[-1]; j=i+1
    else:
        b=np.flatnonzero(w[ipk:]<lvl)
        if not len(b): return np.nan
        j=ipk+b[0]; i=j-1
    return i+(lvl-w[i])/(w[j]-w[i]) if w[j]!=w[i] else float(i)
if __name__=="__main__":
  pass
for v in ([] if __name__!="__main__" else "xy"):
    print('== view',v)
    for leg in ('sim','data'):
        W=z[f'{leg}_{v}_wf'][:,3]; a=W.max(1); ok=a<3500
        Wu=W[ok]
        L10=np.array([cross(w,.1,-1) for w in Wu]); R10=np.array([cross(w,.1,1) for w in Wu])
        L50=np.array([cross(w,.5,-1) for w in Wu]); R50=np.array([cross(w,.5,1) for w in Wu])
        pk=np.argmax(Wu,1)
        dur10=(R10-L10)*DT; fwhm=(R50-L50)*DT; tpk=(pk-L10)*DT
        q=lambda x:[int(np.nanquantile(x,p)) for p in (.1,.25,.5,.75,.9)]
        print(f' {leg:4s} unsat {ok.mean():.2f}  dur(10%->10%) {q(dur10)}  fwhm {q(fwhm)}  peak-after-10% {q(tpk)}')
        # average aligned on leading-edge 50% crossing (sub-sample via interpolation)
        grid=np.arange(-10,25,0.5); acc=np.zeros_like(grid); n=0
        for w,l in zip(Wu,L50):
            if np.isnan(l): continue
            acc+=np.interp(grid+l,np.arange(len(w)),w/w.max(),left=np.nan,right=np.nan) if False else np.interp(grid+l,np.arange(len(w)),w/w.max()); n+=1
        acc/=n
        print('   avg (aligned at leading 50%), t[ns]:', ' '.join(f'{int(g*DT)}:{x:.2f}' for g,x in zip(grid[::4],acc[::4])))
