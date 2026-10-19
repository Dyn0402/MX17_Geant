"""Every draw (D01) and flash (D03) of a gerber layer inside a window, grouped by aperture, with shapes relative to the
0.78 mm pad grid (2026-10-19). Usage: gerber_layer_window.py GERBER [half-width mm, default 3]"""
import re, sys, numpy as np
from collections import Counter, defaultdict
pitch = 0.78; origin = -(512 - 1) / 2 * pitch; W = float(sys.argv[2]) if len(sys.argv) > 2 else 3.
ap_re = re.compile(r"%ADD(\d+)(\w+),([0-9.X]+)\*%"); fs_re = re.compile(r"%FSLAX(\d)(\d)Y(\d)(\d)\*%")
dn_re = re.compile(r"^(?:G54)?D(\d+)\*$"); op_re = re.compile(r"(?:X(-?\d+))?(?:Y(-?\d+))?D0?([123])\*")
aps, cur, x, y, dec = {}, None, 0, 0, 6; draws, flashes = defaultdict(list), defaultdict(list); regions = 0
for line in open(sys.argv[1], errors='replace'):
    line = line.strip()
    if (m := fs_re.match(line)): dec = int(m.group(2)); continue
    if (m := ap_re.match(line)): aps[int(m.group(1))] = f'{m.group(2)},{m.group(3)}'; continue
    if line.startswith('G36'): regions += 1
    if (m := dn_re.match(line)): cur = int(m.group(1)); continue
    for m in op_re.finditer(line):
        px, py = x, y
        if m.group(1) is not None: x = int(m.group(1))
        if m.group(2) is not None: y = int(m.group(2))
        s = 10. ** -dec; X0, Y0, X1, Y1 = px * s, py * s, x * s, y * s
        if m.group(3) == '1' and max(abs(X0), abs(X1), abs(Y0), abs(Y1)) < W: draws[cur].append((X0, Y0, X1, Y1))
        if m.group(3) == '3' and max(abs(X1), abs(Y1)) < W: flashes[cur].append((X1, Y1))
fr = lambda v: round(((v - origin) / pitch) % 1, 3)
print(sys.argv[1], f'regions (G36) in file: {regions}')
for k, v in draws.items():
    c = Counter((round(b[2] - b[0], 3), round(b[3] - b[1], 3), fr(b[0]), fr(b[1])) for b in v)
    print(f'  draw D{k} {aps.get(k)}: {len(v)}  top shapes (dx, dy, x0 frac, y0 frac):', c.most_common(5))
for k, v in flashes.items():
    c = Counter((fr(a), fr(b)) for a, b in v)
    print(f'  flash D{k} {aps.get(k)}: {len(v)}  grid fractions:', c.most_common(4))
