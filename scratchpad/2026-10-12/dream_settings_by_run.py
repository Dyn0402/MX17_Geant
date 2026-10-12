"""Per-run DREAM settings table from archived cfgs (pulled from EOS into CFG_ROOT).
Usage: python dream_settings_by_run.py <CFG_ROOT> <out.csv>
CFG_ROOT mirrors /eos/experiment/ntof/data/x17 (cosmic_bench/..., july_beam/runs/..., p2_sps_july/runs/...).
Uses nTof_x17_DAQ/dream_scripts/decode_dream_cfg.py. One row per (area, run, FEU role), with every
distinct value seen across the run's subruns joined by '/'."""
import csv, os, sys
sys.path.insert(0, os.path.expanduser('~/PycharmProjects/nTof_x17_DAQ/dream_scripts'))
from decode_dream_cfg import summarize
root, out = sys.argv[1], sys.argv[2]
F = ('Peaking_time_ns', 'Integration_time_constant', 'Polarity', 'External_input', 'input_range_fC',
     'Threshold_DAC', 'latency_samples')
rows = {}
for dp, _, names in os.walk(root):
    for n in names:
        if not n.endswith('.cfg') or '/raw_daq_data' not in dp:
            continue
        rel = os.path.relpath(os.path.join(dp, n), root).split('/')
        i = rel.index('raw_daq_data')
        area = '/'.join(rel[:2]) if rel[0] != 'cosmic_bench' else '/'.join(rel[:3])
        run = rel[i - 2] if rel[i - 2] not in ('runs', 'Run') else rel[i - 1]
        rep = summarize(os.path.join(dp, n))
        for feu, e in rep['feus'].items():
            role = 'Trg' if 'Trg' in e['roles'].values() else 'Dat'
            r = rows.setdefault((area, run, role), {'subruns': set(), 'feus': set(), 'samples': set(),
                                                    'T_sample_ns': set(), **{f: set() for f in F}})
            r['subruns'].add(rel[i - 1]); r['feus'].add(feu)
            r['samples'].add(e['samples']); r['T_sample_ns'].add(e['sample_period_ns'])
            for d in ([e['dream']] if e['dreams_identical'] else e['dream'].values()):
                for f in F:
                    r[f].add(d.get(f))
key = lambda s: sorted(s, key=str)
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['area', 'run', 'role', 'n_subruns', 'feus', 'samples', 'T_sample_ns'] + list(F))
    for (area, run, role), r in sorted(rows.items(), key=lambda kv: (kv[0][0], [(0, int(t), '') if t.isdigit() else (1, 0, t) for t in kv[0][1].replace('_', ' ').split()], kv[0][2])):
        w.writerow([area, run, role, len(r['subruns']), ' '.join(map(str, sorted(r['feus']))),
                    '/'.join(map(str, key(r['samples']))), '/'.join(map(str, key(r['T_sample_ns'])))] +
                   ['/'.join(map(str, key(r[f]))) for f in F])
print(len(rows), 'rows ->', out)
