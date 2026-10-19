import json, os, sys
base='/eos/experiment/ntof/data/x17/july_beam/runs'
with open('jobs.txt','w') as fh:
    for run in sys.argv[1:]:
        c=json.load(open(f'{base}/{run}/run_config.json'))
        for s in c['sub_runs']:
            if os.path.isdir(f"{base}/{run}/{s['sub_run_name']}/decoded_root"): fh.write(f"{run} {s['sub_run_name']}\n")
print(sum(1 for _ in open('jobs.txt')))
