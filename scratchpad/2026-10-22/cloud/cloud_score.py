"""Sim early-cloud width (2026-10-22; FOLLOWUP 2026-10-21 §10 next step 1). The data's neighbour/central ratio r(f) is
flat from the onset of the central pulse to its peak in both views on every board; the sim's (glg0 = v5a) rises 2-3x.
Arms (condor 4526808, run_cloud_job.sh, each one change on top of glg0, paired seeds/angles): cdt0s4 / cdt0s2 = no drift
diffusion + sigma0 400 / 250 (depth-independent clouds), cs4 = sigma0 400 with the wet diffusion, cdry = dry-gas D_T
(the bench's 'widens like dry gas'), cH4 = feedback spread 400 um. Same near-vertical cut and r(f) points as
2026-10-21/prompt/coupling_checks.py; X from every angle, Y from th00 only (thy* are inclined in Y).
Shape numbers per s1 bin: early/peak = r05/rpk and late/peak = rf50/rpk (data det3 May X 0.65 / 2.3).
Usage: cloud_score.py [out.json]   (NB_RAW defaults to ~/x17/nb_asym_raw; sim arms under simcloud/ and simgapleak/)"""
import os, sys, glob, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NB_RAW', os.path.expanduser('~/x17/nb_asym_raw'))
sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-10-21', 'prompt'))
import coupling_checks as CC
RAW = os.environ['NB_RAW']
DATA = {'det3 May': 'd3may', 'det3 June-15': 'd3jun', 'det2': 'd2'}
ARMS = [a for a in os.environ.get('ARMS', 'glg0,cdt0s4,cdt0s2,cs4,cdry,cH4').split(',') if a]

def sim_files(tag, view):
    pat = f'{tag}_th*_{view}.npz' if view == 'x' else f'{tag}_th00_*_{view}.npz'
    return sorted(glob.glob(f'{RAW}/simcloud/{pat}') + glob.glob(f'{RAW}/simgapleak/{pat}'))

def shape(pr):
    out = {}
    for b, row in pr.items():
        if 'peak' in row and 'r05' in row and 'f50' in row:
            out[b] = dict(early=row['r05'][0] / row['peak'][0], late=row['f50'][0] / row['peak'][0])
    return out

if __name__ == '__main__':
    R = {}
    for view in ('x', 'y'):
        sets = {f'{n} {view.upper()}': [f'{RAW}/{t}_{view}.npz'] for n, t in DATA.items()}
        sets.update({f'sim {a} {view.upper()}': sim_files(a, view) for a in ARMS})
        for name, files in sets.items():
            if not files: print(f'\n{name}: no files'); continue
            W, s1, n = CC.events(files); fl = CC.floor(s1); pr = CC.profile(W, s1)
            R[name] = dict(n=n, n_sel=len(s1), n_files=len(files), floor=fl, profile=pr, shape=shape(pr))
            CC.show(f'{name} [{len(files)} files]', n, len(s1), fl, pr)
            print('   shape early/peak, late/peak: ' + '  '.join(f'{b} {v["early"]:.2f} / {v["late"]:.2f}' for b, v in R[name]['shape'].items()))
    json.dump(R, open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'cloud_score.json'), 'w'), indent=1, default=float)
