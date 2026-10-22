"""Level each macro photo on the ruler's top edge (2026-10-22, §14): per 200-px column block, the row of the strongest
dark->bright vertical step in the lower 2/3 (restricted to blocks where the ruler is present), a robust line fit, and a
rotation by its angle. Writes $LV/<tag>.png (default ./leveled; keep it out of the repo) (grey, full resolution) and a 1/4-scale preview with a 100-px grid.
Usage: level.py RAW_DIR"""
import os, sys, glob, json, numpy as np
from PIL import Image, ImageOps, ImageDraw
from scipy.ndimage import uniform_filter
LV = os.environ.get('LV', 'leveled'); os.makedirs(LV, exist_ok=True); info = {}
for i, f in enumerate(sorted(glob.glob(os.path.join(sys.argv[1], '*MACRO*.jpg')))):
    tag = f'p{i + 1}'; im = ImageOps.exif_transpose(Image.open(f)).convert('L'); g = np.asarray(im, float); H, W = g.shape
    gy = uniform_filter(np.diff(g, axis=0), (9, 31)); xs, ys, st = [], [], []
    for x in range(100, W - 100, 200):
        col = gy[H // 3:, x]; y = int(np.argmax(col)) + H // 3; xs.append(x); ys.append(y); st.append(col.max())
    xs, ys, st = map(np.array, (xs, ys, st)); ok = st > 0.5 * np.median(st)
    for _ in range(3):
        p = np.polyfit(xs[ok], ys[ok], 1); r = ys - np.polyval(p, xs); ok = ok & (np.abs(r) < max(3 * np.std(r[ok]), 6))
    ang = np.degrees(np.arctan(p[0])); rim = im.rotate(ang, resample=Image.BICUBIC, expand=False, fillcolor=0)
    rim.save(f'{LV}/{tag}.png'); info[tag] = dict(file=os.path.basename(f), angle_deg=float(ang), edge_row=float(np.polyval(p, W / 2)), n_blocks=int(ok.sum()))
    pv = rim.resize((W // 4, H // 4)).convert('RGB'); d = ImageDraw.Draw(pv)
    for y in range(0, H // 4, 25): d.line([(0, y), (W // 4, y)], fill=(255, 0, 0) if y % 100 == 0 else (90, 0, 0))
    pv.save(os.path.join(os.environ.get('PV', '.'), f'{tag}_lv.jpg'), quality=80)
    print(tag, info[tag])
json.dump(info, open('level.json', 'w'), indent=1)
