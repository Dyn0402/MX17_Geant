"""M3 reference-track join for the raw-ADC npz (2026-10-22, FOLLOWUP §12). rays files are matched to the npz's decoded
files by their '_NNN' file number; evn == decoded eventId within a file. Good track: Chi2X, Chi2Y < 5 (qa_config recipe)."""
import os, re, glob, numpy as np, uproot, awkward as ak
CHI2 = 5.0

def rays_for(decoded_file):
    run_dir = os.path.dirname(os.path.dirname(decoded_file)); m = re.search(r'(datrun_\d+_\d+H\d+_\d{3})_\d\d\.root$', decoded_file)
    for sub in ('m3_tracking_root_v2', 'm3_tracking_root'):
        f = glob.glob(f'{run_dir}/{sub}/*{m.group(1)}_rays.root')
        if f: return f[0]
    return None

def m3_table(decoded_file, zdet):
    """dict evn -> (n_all, n_good, thx, thy, x, y, chi2x, chi2y) using the single good track (nan if n_good != 1)."""
    rf = rays_for(decoded_file)
    if rf is None: return None
    r = uproot.open(rf)['T'].arrays(['evn', 'Z_Up', 'Z_Down', 'X_Up', 'Y_Up', 'X_Down', 'Y_Down', 'Chi2X', 'Chi2Y'], library='ak')
    good = (r['Chi2X'] < CHI2) & (r['Chi2Y'] < CHI2); n_all = ak.to_numpy(ak.num(r['Chi2X'])); n_good = ak.to_numpy(ak.sum(good, axis=1))
    one = n_good == 1; first = lambda v: ak.to_numpy(ak.firsts(r[v][good][one]))
    xu, xd, yu, yd, cx, cy = (first(v) for v in ('X_Up', 'X_Down', 'Y_Up', 'Y_Down', 'Chi2X', 'Chi2Y'))
    zu, zd = ak.to_numpy(r['Z_Up'])[one], ak.to_numpy(r['Z_Down'])[one]; dz = zd - zu; t = (zdet - zu) / dz
    tab = np.full((len(r), 8), np.nan); tab[:, 0], tab[:, 1] = n_all, n_good
    tab[one, 2:] = np.c_[np.degrees(np.arctan((xd - xu) / dz)), np.degrees(np.arctan((yd - yu) / dz)), xu + t * (xd - xu), yu + t * (yd - yu), cx, cy]
    return dict(zip(ak.to_numpy(r['evn']).astype(int).tolist(), tab.tolist()))
