import sys, numpy as np
import partial_m3 as P
for name in sys.argv[1:]:
    for view in ('x', 'y'):
        r = P.session(name, view); e = r['_ev']; M = e['M']; g = np.isfinite(M[:, 2]); p = P.L.load(f"{P.RAW}/{P.SESS[name][0]}_{view}.npz")['p']
        cx = [np.corrcoef(p[g & (p // 64 == c)], M[g & (p // 64 == c), 4])[0, 1] for c in range(8) if (g & (p // 64 == c)).sum() > 20]
        cy = [np.corrcoef(p[g & (p // 64 == c)], M[g & (p // 64 == c), 5])[0, 1] for c in range(8) if (g & (p // 64 == c)).sum() > 20]
        mx = [np.median(M[g & (p // 64 == c), 4]) for c in range(8) if (g & (p // 64 == c)).sum() > 20]
        my = [np.median(M[g & (p // 64 == c), 5]) for c in range(8) if (g & (p // 64 == c)).sum() > 20]
        print(f'{name} {view}: chip median M3 x', np.round(mx), ' y', np.round(my))
        a = np.abs(e['dtn']); 
        for nm, th in (('thx', M[:, 2]), ('thy', M[:, 3])):
            s = [f'{lo}-{hi}°: {np.mean(a[g & (np.abs(th) >= lo) & (np.abs(th) < hi)] < 80):.2f} fwhm {np.nanmedian(e["fw"][g & (np.abs(th) >= lo) & (np.abs(th) < hi)]):.0f}' for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 20), (20, 40))]
            print('   ', nm, ' P(|dtn|<80), FWHM by |angle|:', ' | '.join(s))
