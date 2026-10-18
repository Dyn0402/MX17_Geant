"""Print the ladder compression table: data, ideal (no readout), and each sim point. C = 1 is unbiased."""
import json, sys
r = json.load(open(sys.argv[1])); idl = json.load(open(sys.argv[2]))
angs = ['thx05', 'thx10', 'thx15', 'thx20', 'thx30', 'thy05', 'thy10', 'thy15', 'thy20', 'thy30']
for est, iest in (('cfd', 'med'), ('cfd_xt', 'med_xt')):
    print(f'\n## estimator {est} (t on x)' if est == 'cfd' else f'\n## estimator {est} (x on t)')
    print('| | ' + ' | '.join(angs) + ' |'); print('|---' * (len(angs) + 1) + '|')
    f = lambda o, k: (f"{o[k]['C']:.3f}" if o and o.get(k) else '–')
    print('| data (M3) | ' + ' | '.join(f(r['data'].get(a), est) if r['data'].get(a, {}).get('n', 0) >= 30 else
                                      (f"({r['data'][a][est]['C']:.2f}, n={r['data'][a]['n']})" if r['data'].get(a, {}).get(est) else '–') for a in angs) + ' |')
    print('| ideal (no readout) | ' + ' | '.join(f(idl.get(a), iest) for a in angs) + ' |')
    for pt, d in r['sim'].items():
        print(f'| {pt} | ' + ' | '.join(f(d.get(a), est) for a in angs) + ' |')
