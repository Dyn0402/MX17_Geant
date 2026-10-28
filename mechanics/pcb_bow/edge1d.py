#!/usr/bin/env python3
"""edge1d.py -- what rotational restraint does the real MX17 mounting give the board at the aperture edge?

The 2-D model (plate.py, 10 mm cells) cannot resolve the 10-20 mm wide clamping zone, so this is a fine (0.25 mm) 1-D
cross-section of the board in cylindrical bending, from the centre (symmetry) out past the frame, with the geometry
read off the STEP file (MX17 2026 v1.3, global placements via gmsh; distances from the active-area axis):

    0     - 201.0  free span over the 402 mm support-plate aperture (nothing behind the board)
    201.0 - 205.0  board on the Rohacell over the Al plate; gas above (inside the 410 mm frame aperture)
    205.0 - 206.55 frame inner lip above the board (if the frame bottoms on the board)
    206.55- 210.25 O-ring groove in the frame face; the O-ring ('joint', 3.7 mm wide) seals on the board, seal ~208.4
    210.25- 220.0  frame face above, Rohacell + plate below; M4-size screws at 214 mm on a 55 mm pitch
    220   - 235    board overhang (Rohacell + plate below, nothing above)

Supports:
  foam      Rohacell between board and plate outside the aperture: Winkler foundation k_f = E/t (1.4e10 N/m^3 for
            Rohacell 51), compression only (the tape cannot pull the board onto the plate), with a precompression
            delta0 from the screws (0 = just touching).
  frame     rigid, one-sided from above where its face meets the board (mode 'bottomed'), or absent (mode 'oring':
            the frame rides on the O-ring, which is then a preloaded line spring k_or from above at the seal).
  O-ring    sealing reaction F_or [N/m] pushing the board down at the seal (mode 'bottomed'), or the line spring.
  gas       uniform overpressure inside the seal line.
Outputs the effective rotational stiffness K_eff = M / theta at the aperture edge (x = 201 mm), the quantity the 2-D
model's 'spring' edge takes, plus the 1-D centre deflection against ideal clamped / pinned edges.
"""
import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AP, LIP0, GR0, SEAL, GR1, FR1, END = 201.0, 205.0, 206.55, 208.4, 210.25, 220.0, 235.0


def solve(D, q, kf=1.4e10, delta0=0.0, mode='bottomed', F_or=0.0, k_or=1e7, h_mm=0.25, kpen=1e14):
    """D [N m], q [Pa] (positive = toward the plate).
    Returns x [mm], w [m] (negative = toward the plate), K_eff [N m / m / rad]."""
    h = h_mm * 1e-3
    x = np.arange(0.0, END + 1e-9, h_mm)
    n = len(x)
    # curvature operator: node 0 symmetric ghost, free end (no curvature energy at the last node)
    rows, cols, vals, wt = [], [], [], []
    for i in range(n - 1):
        if i == 0:
            rows += [i, i]; cols += [0, 1]; vals += [-2 / h ** 2, 2 / h ** 2]; wt.append(h / 2)
        else:
            rows += [i] * 3; cols += [i - 1, i, i + 1]; vals += [1 / h ** 2, -2 / h ** 2, 1 / h ** 2]; wt.append(h)
    B = np.zeros((n - 1, n)); B[rows, cols] = vals
    Kb = D * B.T @ (np.array(wt)[:, None] * B)
    wn = np.full(n, h); wn[0] = wn[-1] = h / 2
    f = np.where(x < SEAL, -q, 0.0) * wn
    foam = x > AP
    frame = ((x >= LIP0) & (x < GR0)) | ((x >= GR1) & (x <= FR1)) if mode == 'bottomed' else np.zeros(n, bool)
    iseal = int(round(SEAL / h_mm))
    if mode == 'bottomed':
        f[iseal] -= F_or * 1.0                  # line load N/m on one node (node force per unit width)
    w = np.zeros(n)
    for it in range(200):
        Kt = Kb.copy(); ft = f.copy()
        a_foam = foam & (w < delta0 * 1e-3)    # foam compressed where w < delta0 (foam top surface at +delta0)
        Kt[np.diag_indices(n)] += np.where(a_foam, kf * wn, 0.0)
        ft += np.where(a_foam, kf * wn * delta0 * 1e-3, 0.0)
        a_fr = frame & (w > 0)
        Kt[np.diag_indices(n)] += np.where(a_fr, kpen * wn, 0.0)
        if mode == 'oring':
            Kt[iseal, iseal] += k_or            # preloaded line spring (bilateral while preloaded), N/m per m
            ft[iseal] -= F_or
        wnew = np.linalg.solve(Kt, ft)
        if it > 0 and np.array_equal(a_foam, foam & (wnew < delta0 * 1e-3)) and np.array_equal(a_fr, frame & (wnew > 0)):
            w = wnew; break
        w = wnew
    i = int(round(AP / h_mm))
    kap = (w[i - 1] - 2 * w[i] + w[i + 1]) / h ** 2
    th = (w[i + 1] - w[i - 1]) / (2 * h)
    # restraint moment the support region exerts on the span: M = D w'' at the edge; span rotation theta = -w'
    K = (D * kap) / th if th != 0 else np.inf
    return x, w, float(abs(K))


def main():
    out = dict(geometry=dict(aperture=AP, lip=LIP0, groove=[GR0, GR1], seal=SEAL, frame_outer=FR1, end=END), rows=[])
    E, nu = 20e9, 0.3
    for t_mm in (1.3, 1.6, 2.0):
        D = E * (t_mm * 1e-3) ** 3 / (12 * (1 - nu ** 2))
        for mode, kw in (('bottomed', dict(F_or=0.0)), ('bottomed', dict(F_or=2000.0)),
                         ('bottomed', dict(delta0=0.02)), ('bottomed', dict(kf=0.73e10)),
                         ('oring', dict(F_or=2000.0, k_or=1e7)), ('oring', dict(F_or=2000.0, k_or=1e6))):
            for q in (100.0, 500.0, 1000.0):
                x, w, K = solve(D, q, mode=mode, **kw)
                # reference strips of span 402 mm, uniform load: centre q L^4 / 384 D (clamped), 5 q L^4 / 384 D (pinned)
                L4 = (2 * AP * 1e-3) ** 4
                c_ref, p_ref = -q * L4 / (384 * D), -5 * q * L4 / (384 * D)
                frac = (w[0] - c_ref) / (p_ref - c_ref)      # 0 = clamped, 1 = pinned (centre deflection)
                row = dict(t_mm=t_mm, mode=mode, q_Pa=q, K_eff=K, K_a_over_D=K * 0.402 / D, frac_to_pinned=float(frac),
                           centre_mm=float(-w[0] * 1e3), clamped_mm=float(-c_ref * 1e3), pinned_mm=float(-p_ref * 1e3),
                           lift_max_um=float(w[x > AP].max() * 1e6), **{k: v for k, v in kw.items()})
                out['rows'].append(row)
                print(f"t={t_mm} {mode:8s} {str(kw):34s} q={q:5.0f}  K_eff={K:9.0f} N m/m  K a/D={row['K_a_over_D']:7.1f}"
                      f"  clamped->pinned {100 * frac:5.1f} %  max lift beyond aperture {row['lift_max_um']:6.2f} um")
    # profiles near the edge for the deck (1.6 mm PCB)
    D = E * 1.6e-3 ** 3 / (12 * (1 - nu ** 2))
    out['profiles'] = {}
    for lab, mode, kw, q in (('bottomed_5mbar', 'bottomed', dict(F_or=2000.0), 500.0),
                             ('bottomed_10mbar', 'bottomed', dict(F_or=2000.0), 1000.0),
                             ('oring_10mbar', 'oring', dict(F_or=2000.0, k_or=1e6), 1000.0)):
        x, w, K = solve(D, q, mode=mode, **kw)
        sel = x >= 150
        out['profiles'][lab] = dict(x_mm=x[sel][::8].tolist(), w_um=(w[sel][::8] * 1e6).tolist(), K_eff=K)
    json.dump(out, open(f'{HERE}/results/edge1d.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
