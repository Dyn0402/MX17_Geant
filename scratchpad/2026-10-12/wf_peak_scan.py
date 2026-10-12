import sys, numpy as np, itertools
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,'/home/dylan/CLionProjects/MX17_Geant/scratchpad/2026-10-12')
from wf_tail20 import avg, grid
import toy_pulse as T
DT=60
za=np.load('waveforms_angled_pm3.npz')
D={th:avg(za[f'{tag}_data_wf'][:,3]) for tag,th in (('thx20',20.0),('thx10',10.7))}
from wf_matched import table
R={th:table(za[f'{tag}_data_wf'][:,3]).rise.median() for tag,th in (('thx20',20.0),('thx10',10.7))}
_ORIG=T.shaper
def job(c):
    fam,tp,f,Tion=c
    T.shaper=lambda family,t_peak=180.0: _ORIG(family,tp)
    out=[]
    for th in (20.0,10.7):
        B=T.sampled_peak_strips({**T.BASE,'shaper':fam,'f_ion':f,'t_ion':Tion,'theta':th},800,seed=4,target_peak=2000)[:,3]
        a=avg(B); m=grid>=-2
        out.append((float(np.sqrt(np.mean((a[m]-D[th][m])**2))), float(table(B).rise.median())))
    return c,out
cfgs=list(itertools.product(['crrc2','manual_1'],[76,123,180,228],[0.0,0.5,0.908],[340.0]))+[(fam,tp,0.908,170.0) for fam in ('crrc2','manual_1') for tp in (76,123)]
with ProcessPoolExecutor(14) as ex: res=list(ex.map(job,cfgs))
print('data rise median 20deg %.0f  10deg %.0f'%(R[20.0],R[10.7]))
print('shaper   tpk  f_ion T_ion | shape-rms20 rise20 | shape-rms10 rise10')
for c,o in sorted(res,key=lambda r:r[1][0][0]+r[1][1][0]):
    print(f'{c[0]:8s} {c[1]:4d} {c[2]:5.2f} {c[3]:5.0f} |   {o[0][0]:.3f}    {o[0][1]:4.0f}  |   {o[1][0]:.3f}    {o[1][1]:4.0f}')
