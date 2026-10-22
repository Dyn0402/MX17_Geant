"""6-17 resist scan: does the bottom detector's gain follow its resist supply? (2026-10-22, §13.) Per sub-run (1 file,
raw_asym.py extraction in ~/x17/nb_asym_raw/scan617/): FEU 3 = mx17_1 (det2 per copy_survey; bottom, z 232, resist 3:3)
and FEU 7 = mx17_2 (top, z 702, resist 3:4). All isolated events: median 9-pad peak A9, fraction of 9-pad columns ending
before sample 14, central FWHM < 350 ns fraction, median s1. Set voltages from run_config. Usage: scan617.py"""
import glob, json, os, numpy as np
import session_profile as SP, event_ends as EE
RUN = os.path.expanduser('~/x17/cosmic_bench/det1_det2/mx17_det1_det2_overnight_6-17-26'); cfg = json.load(open(f'{RUN}/run_config.json'))
HV = {s['sub_run_name']: s['hvs'] for s in cfg['sub_runs']}; out = {}
print('sub-run                    set 3:3  3:4   | FEU 3 (bottom): n    A9   end<14 short  s1   | FEU 7 (top): n    A9   end<14 short  s1')
for f3 in sorted(glob.glob(os.path.expanduser('~/x17/nb_asym_raw/scan617/*_f3.npz'))):
    sub = os.path.basename(f3)[:-7]; row = [f"{sub:26s} {HV[sub]['3'].get('3')!s:>5s} {HV[sub]['3'].get('4')!s:>5s}  |"]; out[sub] = dict(hv=HV[sub])
    for feu in (3, 7):
        d = SP.L.load(f3.replace('_f3', f'_f{feu}')); W = d['W'].astype(float); S9 = W.sum(1); A9 = S9.max(1); e9 = EE.last_cross(S9, .3)
        c = W[:, 4]; fw = (SP.xcross(c, .5, False) - SP.xcross(c, .5, True)) * 60.; s1 = np.einsum('ij,ij->i', W[:, 3] + W[:, 5], c) / np.einsum('ij,ij->i', c, c)
        r = dict(n=len(A9), A9=float(np.median(A9)), early_end=float(np.mean(e9 < 14)), short=float(np.nanmean(fw < 350)), s1=float(np.median(s1)))
        out[sub][f'feu{feu}'] = r; row.append(f"{r['n']:5d} {r['A9']:6.0f}  {r['early_end']:.2f}  {r['short']:.2f}  {r['s1']:.2f}  |")
    print(' '.join(row))
json.dump(out, open('scan617.json', 'w'), indent=1)
