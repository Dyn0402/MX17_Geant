import sys, json, itertools
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, '.')
import toy_pulse as T

def job(args):
    over, tp = args
    p = dict(T.BASE); p.update(over)
    r = T.run(p, n_events=2000, seed=7, target_peak=tp)
    r['target_peak'] = tp
    return r

grid = []
for tp in (1650, 3000):
    for f in (0.0, 0.2, 0.4, 0.6, 0.908):
        grid.append(({'f_ion': f}, tp))
    for T_ion in (120, 200, 260):
        grid.append(({'t_ion': T_ion}, tp))
    for c in (0.3, 0.6):
        for tau in (100, 230, 410):
            grid.append(({'c_spread': c, 'tau_s': tau}, tp))
    grid.append(({'v': 0.0366}, tp))
    grid.append(({'v': 0.0281}, tp))
    grid.append(({'f_ion': 0.0, 'c_spread': 0.3, 'tau_s': 230}, tp))
with ProcessPoolExecutor(12) as ex:
    res = list(ex.map(job, grid))
json.dump(res, open(sys.argv[1], 'w'), indent=0)
for r in res:
    p = r['params']; d = {k: p[k] for k in ('f_ion','t_ion','c_spread','tau_s','v')}
    print(r['target_peak'], d, 'rise q', [round(x) for x in r['rise_q']], 'fwhm', round(r['fwhm_med']), 'sat %.2f' % r['sat'])
