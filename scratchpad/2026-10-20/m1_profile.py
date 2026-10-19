"""Is the det3 X copy a readout-side (M1 card / connector / cable) effect? Data handles (2026-10-20):
 1. profile across the board in 16-column bins (true columns): a per-card effect steps at the 64-column connector edges,
    a detector effect varies smoothly through them;
 2. DAQ neighbours across a connector boundary (DAQ 64k-1 -> 64k: different connector + DREAM chip, columns 127 apart):
    electronics after the connectors would copy there too; the detector cannot;
 3. distance from the readout edge along the X line (Y position of the hit): coupling along the line or in the fan-out.
det3 June-15, all matched sub-runs, Y-pad and all events. Writes m1_profile.json."""
import os, sys, json, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '2026-10-19', 'asym_deck'); sys.path.insert(0, D)
import make_asym_deck as M, asym_lib as L, asym_more as A
rows = M.survey_rows(); OUT = {}
col = lambda ch: 64 * (ch % 512 // 64) + 63 - ch % 64
for bu in ('det3 June-15', 'det2', 'det3 May'):
    E = {k: [] for k in ('a1', 'px', 'py', 'yp', 'W')}
    for rx, ry in A.pairs(rows, lambda r: r['build'] == bu):
        dx, dy, ix, iy = A.matched(rx['f'], ry['f'])
        if len(ix) < 20: continue
        par = (dx['p'][ix] + dy['p'][iy]) % 2; ysh = (dy['al'][1] + dy['al'][-1])[iy]
        E['yp'].append(par == int(ysh[par == 1].mean() < ysh[par == 0].mean()))
        E['a1'].append(dx['a1'][ix]); E['px'].append(dx['p'][ix] % 512); E['py'].append(dy['p'][iy] % 512); E['W'].append(dx['W'][ix])
    E = {k: np.concatenate(v) for k, v in E.items()}
    o = {}
    # 1. profile in 16-column bins, interior of each connector vs its edges
    c = col(E['px']); inner = (E['px'] % 64 >= 2) & (E['px'] % 64 <= 61)
    prof = []
    for b in range(32):
        m = inner & (c // 16 == b) & E['yp']
        prof.append([b * 16 + 8, float(E['a1'][m].mean()) if m.sum() > 40 else None, float(E['a1'][m].std() / np.sqrt(max(m.sum(), 1))), int(m.sum())])
    o['profile_ypad'] = prof
    # step test: within each connector, quarter means; jump between last quarter of connector k and first of k+1
    q = [p_[1] for p_ in prof]; jumps = []
    for k in range(7):
        a_, b_ = prof[4 * k + 3], prof[4 * k + 4]
        if a_[1] is not None and b_[1] is not None:
            jumps.append([k, b_[1] - a_[1], float(np.hypot(a_[2], b_[2]))])
    inq = []
    for k in range(8):
        for qq in range(3):
            a_, b_ = prof[4 * k + qq], prof[4 * k + qq + 1]
            if a_[1] is not None and b_[1] is not None: inq.append([k, b_[1] - a_[1], float(np.hypot(a_[2], b_[2]))])
    o['jumps_at_connector'] = jumps; o['steps_inside'] = inq
    # 2. DAQ boundary: peak at DAQ 64k-1 (local 63); its DAQ n+1 = 64k is window index 5 -> another connector
    W = E['W']; cc = W[:, 4]; n2 = np.einsum('ij,ij->i', cc, cc)
    al = lambda k: np.einsum('ij,ij->i', W[:, 4 + k], cc) / n2
    for nm, m in (('peak local 63 (n+1 across boundary)', E['px'] % 64 == 63), ('peak local 0 (n-1 across boundary)', E['px'] % 64 == 0),
                  ('interior', inner)):
        m = m & (E['px'] < 511) & (E['px'] > 0)
        o[nm] = dict(n=int(m.sum()), a_p1=float(al(1)[m].mean()) if m.sum() else None, a_m1=float(al(-1)[m].mean()) if m.sum() else None,
                     se=float(al(1)[m].std() / np.sqrt(max(m.sum(), 1))))
    # 3. along the X line: Y chip of the hit (distance from the X connectors at one board edge)
    o['vs_ychip_ypad'] = [[k, float(E['a1'][E['yp'] & inner & (E['py'] // 64 == k)].mean()) if (E['yp'] & inner & (E['py'] // 64 == k)).sum() > 40 else None,
                           float(E['a1'][E['yp'] & inner & (E['py'] // 64 == k)].std() / np.sqrt(max((E['yp'] & inner & (E['py'] // 64 == k)).sum(), 1)))] for k in range(8)]
    OUT[bu] = o
    print(f'== {bu}')
    print('  profile (Y-pad, 16-col bins):', ' '.join(f'{p_[1]:+.2f}' if p_[1] is not None else '  -  ' for p_ in prof))
    print('  jumps at connector edges   :', ' '.join(f'{j[1]:+.3f}±{j[2]:.3f}' for j in jumps))
    if inq: print(f'  rms step inside connectors : {np.sqrt(np.mean([s[1]**2 for s in inq])):.3f} (mean err {np.mean([s[2] for s in inq]):.3f});'
                  f' at edges {np.sqrt(np.mean([j[1]**2 for j in jumps])):.3f} (mean err {np.mean([j[2] for j in jumps]):.3f})')
    for nm in ('peak local 63 (n+1 across boundary)', 'peak local 0 (n-1 across boundary)', 'interior'):
        r = o[nm]; print(f'  {nm:36s} n={r["n"]:5d}  alpha+1 {r["a_p1"]:+.3f}  alpha-1 {r["a_m1"]:+.3f}  (se {r["se"]:.3f})' if r['n'] else f'  {nm}: none')
    print('  vs Y chip (along the X line):', ' '.join(f'{v[1]:+.2f}' if v[1] is not None else '  -  ' for v in o['vs_ychip_ypad']))
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'm1_profile.json'), 'w'))
