"""Score the mixed-response grid (run_mix.sh) and the pure offset arms against det3 June-15 sessions (2026-10-20).
Also builds a "plain copy" family (x(n+1) += a x(n) on the xo0 arm) as the comparison model.
Data fingerprints cached to fp_data.npz-like pickle (fp_data.pkl). Writes eval_mix.json.
Usage: eval_mix.py"""
import os, re, glob, json, pickle, numpy as np
import fingerprint as FP, r5lib as R
HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.expanduser('~/x17/nb_asym_raw/simxoff'); MIX = '/media/dylan/data/claude_tmp/tmp/nb/xmix/npz'
SESS = {'6-27/28 1000V': lambda r: r['build'] == 'det3 June-15' and r['feu'] in (7, 8) and r['drift'] == '1000',
        '6-22/23 1000V': lambda r: r['build'] == 'det3 June-15' and 'det2_det3_overnight_6-22' in r['run'],
        '6-23/24 600V': lambda r: r['build'] == 'det3 June-15' and '6-23-26' in r['run'],
        'all June-15': lambda r: r['build'] == 'det3 June-15', 'det2': lambda r: r['build'] == 'det2', 'det3 May': lambda r: r['build'] == 'det3 May'}
pk = os.path.join(HERE, 'fp_data.pkl')
if os.path.exists(pk): DATA = pickle.load(open(pk, 'rb'))
else:
    DATA = {k: FP.fp_data(s) for k, s in SESS.items()}; pickle.dump(DATA, open(pk, 'wb'))
if __name__ == "__main__":
    for k, o in DATA.items(): FP.show('DATA ' + k, o)
ANGS = ('th00', 'thy05', 'thy10', 'thy15', 'thy20', 'thy30')
def load_pair(fx, fy):
    zx, zy = np.load(fx), np.load(fy)
    if zx['wf'].ndim != 3 or zy['wf'].ndim != 3: return None
    dx, dy = R.L.load(fx), R.L.load(fy)
    kx = dx['fi'].astype(np.int64) * 10**9 + dx['ev']; ky = dy['fi'].astype(np.int64) * 10**9 + dy['ev']
    _, ix, iy = np.intersect1d(kx, ky, return_indices=True)
    return dx['W'][ix], dy['W'][iy], dx['p'][ix], dy['p'][iy]
def arm(xpat):
    parts = []
    for ang in ANGS:
        for s in (5, 6):
            fx = xpat.format(ang=ang, s=s); fy = f'{SIM}/xo0_{ang}_s{s}_y.npz'
            if os.path.exists(fx) and (p := load_pair(fx, fy)) is not None: parts.append(p)
    return [np.concatenate([q[i] for q in parts]) for i in range(4)] if parts else None
def main():
    MODELS = {}
    for off in (0, 100, 200, 300):
        a = arm(f'{SIM}/xo{"" if off == 0 else "m"}{off}_{{ang}}_s{{s}}_x.npz')
        MODELS[f'D-{off} w1'] = FP.fp_arrays(*a)
        if off == 0:
            for c in (.1, .2, .3):
                Wx = a[0].astype(float).copy(); Wx[:, 1:] += c * a[0][:, :-1]
                MODELS[f'copy {c}'] = FP.fp_arrays(Wx, *a[1:])
    for f in sorted(glob.glob(f'{MIX}/m*_th00_s5_x.npz')):
        m = re.search(r'm(\d+)_w([\d.]+)_th00', f); off, w = m.group(1), m.group(2)
        a = arm(f'{MIX}/m{off}_w{w}_{{ang}}_s{{s}}_x.npz')
        if a is not None: MODELS[f'D-{off} w{w}'] = FP.fp_arrays(*a)
    OUT = {}
    M0 = MODELS['D-0 w1']
    refs = ('6-22/23 1000V', '6-23/24 600V', '6-27/28 1000V', 'all June-15')
    print(f'\n{"model":14s} a1     coin  | ' + ' | '.join(f'{r:>14s}' for r in refs) + '   (chi2/n: A_Y, A_X absolute; fY, m1, p1 as effects vs det2 / sim D0)')
    for nm, o in MODELS.items():
        sc = {r: FP.score_delta(o, M0, DATA[r], DATA['det2']) for r in refs}; OUT[nm] = dict(marg=o['marg'], score=sc)
        print(f'{nm:14s} {o["marg"]["a1"]:+.3f} {o["marg"]["coin"]:.3f} | ' + ' | '.join(f'{c:7.0f}/{n:3d}     ' for c, n in sc.values()))
    for k in ('6-22/23 1000V', '6-27/28 1000V'):
        print(f'\nbest vs {k}:')
        for nm in sorted(MODELS, key=lambda n: FP.score_delta(MODELS[n], M0, DATA[k], DATA['det2'])[0])[:5]:
            FP.show(nm, MODELS[nm])
            print('   per-quantity chi2:', {q: round(FP.score_delta(MODELS[nm], M0, DATA[k], DATA['det2'], tuple(x for x in (q,) if x in ('A_Y', 'A_X')), tuple(x for x in (q,) if x in ('fY', 'm1', 'p1')))[0]) for q in ('A_Y', 'A_X', 'fY', 'm1', 'p1')})
    json.dump(OUT, open(os.path.join(HERE, 'eval_mix.json'), 'w'), indent=1)


if __name__ == "__main__":
    main()
