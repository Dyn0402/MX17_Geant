"""Data: wft sign flips and core tan ratio against M3, by |tan_ref| (det3 production, r06). Same reference code as
data_line_res.py. Flip = sign(tan_fit) != sign(tan_ref) with |tan_fit/tan_ref - 1| > 0.25; core = |ratio - 1| <= 0.25."""
import sys, numpy as np
sys.argv = sys.argv[:1] + ['/dev/null']
exec(open('data_line_res.py').read().split("def main():")[0])
cfg = get_config('sat_det3')
params = cm.load_alignment(os.path.join(cfg.OUT_BASE, 'wft', 'alignment', 'alignment.json'))
df = compat.load_table(os.path.join(cfg.OUT_BASE, 'wft', 'events.parquet'))
results = compat.as_event_results(df)
rays = M3RefTracking(cfg.m3_tracking_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS)
xa, _ya, an = get_xy_angles(rays.ray_data); xa = params.ref_x_sign * np.array(xa)
cm.attach_reference_positions(results, rays, params, xa, an)
ref = {int(r.event_id): cm._rotate_ref_tangents(r, params) for r in results
       if not (np.isnan(r.ref_tan_theta_x) or np.isnan(r.ref_mesh_x_mm))}
df = df[~df['event_id'].duplicated(keep=False)].set_index('event_id')
for i, plane in enumerate('xy'):
    e = [k for k in df.index if k in ref and df.loc[k, f'{plane}_ok']]
    tr = np.array([ref[k][i] for k in e]); tf = df.loc[e, f'{plane}_tan_theta'].to_numpy()
    for lo, hi in ((0.06, 0.12), (0.12, 0.22), (0.22, 0.32), (0.32, 0.45), (0.45, 0.7)):
        m = (np.abs(tr) >= lo) & (np.abs(tr) < hi) & np.isfinite(tf)
        r = tf[m] / tr[m]; core = np.abs(r - 1) <= 0.25
        flip = (np.sign(tf[m]) != np.sign(tr[m])) & ~core
        print(f'{plane} |tan| {lo:.2f}-{hi:.2f}: n {m.sum():5d}  flip {flip.mean():.3f}  outside core {(~core).mean():.3f}  core ratio {np.median(r[core]):.4f}')
