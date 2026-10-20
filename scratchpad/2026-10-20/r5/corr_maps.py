"""Per-board X-centroid bias maps from the neighbour asymmetry (2026-10-21; next step 2 of FINDINGS.md).
What the field does to reconstruction: the 3-strip X centroid (in-phase content al[-1], 1, al[+1], DAQ order) is
pulled by dx = -(al[+1] - al[-1]) / (1 + al[+1] + al[-1]) columns along the TRUE column axis (DAQ n+1 = column - 1,
col = 64 (ch // 64) + 63 - ch % 64). Measured per 8 x 8 cell of (true column, true row) with matched X/Y events,
minus the board-wide symmetric reference (det3 May + det1: no field) so the estimator's own bias drops out.
det3 June-15 also per session, as a scale factor on its all-session map (amplitude follows cloud width, geometry
fixed). Output corr_maps.json: per build, cells [ix][iy] -> dx_um, err_um, n; plus the board mean and rms.
Usage: corr_maps.py"""
import os, json, numpy as np
import r5lib as R
from eval_mix import SESS
NC = 8; P = 780.
def dx_cells(E):
    px = E['px'] % 512; py = E['py'] % 512
    cx = 64 * (px // 64) + 63 - px % 64; cy = 64 * (py // 64) + 63 - py % 64
    ap, am = E['alx'][:, 5], E['alx'][:, 3]; dx = -(ap - am) / (1 + ap + am) * P
    i = cx * NC // 512; j = cy * NC // 512
    M = np.full((NC, NC), np.nan); S = np.full((NC, NC), np.nan); N = np.zeros((NC, NC), int)
    for a in range(NC):
        for b in range(NC):
            m = (i == a) & (j == b); N[a, b] = m.sum()
            if m.sum() >= 25: M[a, b] = dx[m].mean(); S[a, b] = dx[m].std() / np.sqrt(m.sum())
    return M, S, N, dx
REF = [R.ev(R.BUILDS[b]) for b in ('det3 May', 'det1')]
ref = float(np.concatenate([dx_cells(E)[3] for E in REF]).mean())
print(f'symmetric reference (det3 May + det1) mean dx {ref:+.1f} um')
OUT = dict(ref_um=ref, ncell=NC, axis='true column / true row, cell = 64 columns x 64 rows', builds={})
for b in ('det3 June-15', 'det2', 'det4', 'det6', 'det7', 'det3 May', 'det1'):
    E = R.ev(R.BUILDS[b]); M, S, N, dx = dx_cells(E); M = M - ref
    ok = np.isfinite(M); w = 1 / np.maximum(S[ok], 1.) ** 2
    chi = float((((M[ok]) / np.maximum(S[ok], 1.)) ** 2).sum())
    OUT['builds'][b] = dict(dx_um=np.where(ok, M, None).tolist(), err_um=np.where(ok, S, None).tolist(), n=N.tolist(),
                            mean_um=float((M[ok] * w).sum() / w.sum()), rms_um=float(np.sqrt(np.mean(M[ok] ** 2))),
                            chi2_zero=chi, ncells=int(ok.sum()), nev=int(len(dx)))
    o = OUT['builds'][b]
    print(f'{b:14s} n {o["nev"]:6d} cells {o["ncells"]:2d}  mean {o["mean_um"]:+6.1f} um  rms {o["rms_um"]:5.1f} um  '
          f'chi2(0) {chi:6.0f}/{o["ncells"]}')
# det3 June-15 per session: scale factor on the all-session map
E3 = R.ev(R.BUILDS['det3 June-15']); M3, S3, _, _ = dx_cells(E3); M3 = M3 - ref
OUT['det3_sessions'] = {}
for name in ('6-22/23 1000V', '6-23/24 600V', '6-27/28 1000V'):
    Es = R.ev(SESS[name]); m = np.ones(len(Es['px']), bool)
    Ms, Ss, _, _ = dx_cells(Es); Ms = Ms - ref; ok = np.isfinite(Ms) & np.isfinite(M3)
    w = 1 / (np.maximum(Ss[ok], 1.) ** 2); k = float((w * Ms[ok] * M3[ok]).sum() / (w * M3[ok] ** 2).sum())
    chi = float((w * (Ms[ok] - k * M3[ok]) ** 2).sum())
    OUT['det3_sessions'][name] = dict(scale=k, chi2=chi, ncells=int(ok.sum()), nev=int(m.sum()))
    print(f'  det3 {name:22s} n {int(m.sum()):6d}  scale {k:5.2f}  chi2 {chi:5.0f}/{int(ok.sum()) - 1}')
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'corr_maps.json'), 'w'), indent=1)
