# step_extents.py STEP_FILE NAME... [--ring]: local extents of named products (and the O-ring radial extent) from a STEP file.
# Used 2026-10-28 for the MX17 mounting geometry (REPORT.md section 1). Run with the nTof_x17 venv python.
import re, sys, numpy as np
from collections import defaultdict
txt = open(sys.argv[1], encoding='latin-1').read()
ents = {}
for m in re.finditer(r'#(\d+)\s*=\s*(.*?);\s*(?=#\d+\s*=|ENDSEC)', txt, re.S):
    ents[int(m.group(1))] = m.group(2).replace('\n', '')
refs = lambda s: [int(x) for x in re.findall(r'#(\d+)', s)]
typ = lambda i: ents[i].split('(')[0].strip()
def nums(s): return [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', s)]
# product -> PDF -> PD -> PDS -> SDR -> rep
prod = {i: re.search(r"PRODUCT\('([^']*)'", e).group(1) for i, e in ents.items() if e.startswith('PRODUCT(')}
pdf = {i: refs(e)[-1] for i, e in ents.items() if typ(i) in ('PRODUCT_DEFINITION_FORMATION', 'PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE')}
pd = {}
for i, e in ents.items():
    if typ(i) == 'PRODUCT_DEFINITION':
        r = refs(e); pd[i] = prod.get(pdf.get(r[0]))
pds = {i: pd.get(refs(e)[0]) for i, e in ents.items() if typ(i) == 'PRODUCT_DEFINITION_SHAPE' and refs(e) and refs(e)[0] in pd}
rep_of = {}
for i, e in ents.items():
    if typ(i) == 'SHAPE_DEFINITION_REPRESENTATION':
        a, b = refs(e)[:2]
        if a in pds: rep_of.setdefault(pds[a], []).append(b)
# SRR links shape_rep <-> brep rep
srr = defaultdict(list)
for i, e in ents.items():
    if 'SHAPE_REPRESENTATION_RELATIONSHIP' in e and 'REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION' not in e:
        r = refs(e)
        if len(r) >= 2: srr[r[0]].append(r[1]); srr[r[1]].append(r[0])
def points(rep):
    seen, stack, pts = set(), [rep], []
    while stack:
        j = stack.pop()
        if j in seen or j not in ents: continue
        seen.add(j); t = typ(j)
        if t == 'CARTESIAN_POINT':
            pts.append(nums(ents[j].split(',', 1)[1])[:3]); continue
        if t in ('AXIS2_PLACEMENT_3D', 'DIRECTION'): continue
        for r in refs(ents[j]): stack.append(r)
    return np.array(pts)
want = sys.argv[2:]
for name in want:
    for rep in rep_of.get(name, []):
        reps = [rep] + srr.get(rep, [])
        L = [points(r) for r in reps]; L = [x for x in L if len(x)]
        if not L: continue
        P = np.vstack(L)
        P = P[np.isfinite(P).all(1)]
        print(name, len(P), 'min', P.min(0).round(2), 'max', P.max(0).round(2))
if '--ring' in sys.argv:
    for rep in rep_of.get('joint', []):
        L=[points(r) for r in [rep]+srr.get(rep,[])]; P=np.vstack([x for x in L if len(x)])
        c=(P.max(0)+P.min(0))/2; Q=P-c; r=np.maximum(np.abs(Q[:,0]),np.abs(Q[:,1]))
        sel=r>150
        print('centre',c.round(2),'ring radial extent',r[sel].min().round(2),r[sel].max().round(2),'z',Q[sel,2].min(),Q[sel,2].max())
        print(np.unique(np.round(r[sel],1))[:40])
