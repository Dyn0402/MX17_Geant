"""Copy survey table (2026-10-19): one row per (dump, view) from run_survey.sh's index.tsv, with the physical detector,
run date, FEU, mesh/drift HV and gas from each run's run_config.json, and the copy estimators of nb_copy_survey.py
(a_ch = channel p+1 minus channel p-1 content, wiring independent; a = mean(alpha_-1 - alpha_+1) by position over interior peaks; a2 = mean(alpha_-2 - alpha_+2) / mean(alpha_+1)).
Physical detector: mx17_3 and mx17_A are det3; mx17_1 is det3 in det_3/ runs, det4 in det_4/, det1 elsewhere;
mx17_B/C/D = det2/det6/det7; mx17_2/4/6/7 = det2/4/6/7; except the 6-17 run, where mx17_1 is the
June-12 board (later det2) and mx17_2 the June-15 board (later det3), per their build notes.
Usage: survey_table.py INDEX.tsv OUT.csv"""
import sys, os, json, numpy as np
HW = 5
def phys(rd, det):
    if det in ('mx17_3', 'mx17_A'): return 'det3'
    if 'det1_det2_overnight_6-17' in rd:             # the boards were renamed after this run (build notes in run_config)
        return {'mx17_1': 'det2', 'mx17_2': 'det3'}[det]
    if det == 'mx17_1': return 'det3' if '/det_3/' in rd else 'det4' if '/det_4/' in rd else 'det1'
    return {'mx17_B': 'det2', 'mx17_C': 'det6', 'mx17_D': 'det7'}.get(det, det.replace('mx17_', 'det'))
def copy_est(W, ch):
    amp = np.nanmax(W, 2); c = amp[:, HW]
    m = (c > 200) & (c < 3500) & np.isfinite(W[:, HW - 2:HW + 3, :16]).all((1, 2))
    if m.sum() < 50: return None
    W, ch = np.nan_to_num(W[m]), ch[m]; wc = W[:, HW]; n2 = np.einsum('ij,ij->i', wc, wc)
    al = {k: np.einsum('ij,ij->i', W[:, HW + k], wc) / n2 for k in (-2, -1, 1, 2)}
    inner = (ch[:, HW - 1] // 64 == ch[:, HW] // 64) & (ch[:, HW + 1] // 64 == ch[:, HW] // 64)
    d1 = (al[-1] - al[1])[inner]; d2 = al[-2] - al[2]
    # channel space (wiring independent): content of channel p+1 minus channel p-1, whichever positions they sit at
    p = ch[:, HW]; up = np.zeros(len(p)); dn = np.zeros(len(p)); ok = np.zeros(len(p), bool)
    for k in (-1, 1):
        up = np.where(ch[:, HW + k] == p + 1, al[k], up); dn = np.where(ch[:, HW + k] == p - 1, al[k], dn)
    ok = ((ch[:, HW - 1] == p + 1) | (ch[:, HW + 1] == p + 1)) & ((ch[:, HW - 1] == p - 1) | (ch[:, HW + 1] == p - 1))
    dch = (up - dn)[ok]
    orient = float(np.mean(ch[:, HW - 1] == p + 1))     # 1: position -1 is channel p+1 (connector-inverted)
    return (int(m.sum()), float(dch.mean()) if len(dch) else np.nan, float(dch.std() / np.sqrt(max(len(dch), 1))), orient, d1.mean(), d1.std() / np.sqrt(len(d1)), d2.mean() / al[1].mean(), d2.std() / np.sqrt(len(d2)) / al[1].mean(),
            al[1].mean(), float(np.median(c[m])))
rows = []
for line in open(sys.argv[1]):
    f, rd, sub, det = line.rstrip('\n').split('\t')
    if not os.path.exists(f): continue
    rc = json.load(open(f'{rd}/run_config.json')); d = [x for x in rc['detectors'] if x['name'] == det][0]
    s = [x for x in rc['sub_runs'] if x['sub_run_name'] == sub]; hvs = s[0].get('hvs', {}) if s else {}
    hv = lambda key: (hvs.get(str(d['hv_channels'][key][0]), {}) or {}).get(str(d['hv_channels'][key][1])) if key in d.get('hv_channels', {}) else None
    z = np.load(f)
    import glob, re
    fs = sorted(glob.glob(f'{rd}/{sub}/decoded_root/*.root')); mt = re.search(r'datrun_(\d{6})_(\d\d)H(\d\d)', fs[0]) if fs else None
    tstart = f'20{mt.group(1)[:2]}-{mt.group(1)[2:4]}-{mt.group(1)[4:]} {mt.group(2)}:{mt.group(3)}' if mt else ''
    for v in 'xy':
        if f'{v}_wf' not in z or len(z[f'{v}_wf']) < 50: continue
        r = copy_est(z[f'{v}_wf'], z[f'{v}_ch'].astype(int))
        if r is None: continue
        feu = d['dream_feus'][f'{v}_1'][0]
        rows.append([phys(rd, det), det, rd.split('/')[-1], sub, v.upper(), feu, hv('resist'), hv('drift'), rc.get('gas'),
                     tstart, *r])
hdr = 'phys,det,run,sub,view,feu,mesh_V,drift_V,gas,tstart,n,a_ch,ea_ch,frac_inverted,a,ea,a2,ea2,alpha_p1,median_c'
with open(sys.argv[2], 'w') as fh:
    fh.write(hdr + '\n')
    for r in rows: fh.write(','.join(str(x).replace(',', ';') for x in r) + '\n')
print(f'{len(rows)} rows -> {sys.argv[2]}')
