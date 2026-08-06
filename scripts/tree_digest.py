#!/usr/bin/env python3
"""
Order-sensitive digest of every TTree in a ROOT file.

Used to prove that a geometry refactor is bit-identical: run the sim twice with
the same fixed seed (`-s <n>`) before and after the change and compare digests.
A ROOT file cannot be compared byte-wise because it embeds timestamps and
compression state, so the branch values are hashed instead.

    python scripts/tree_digest.py out_before.root out_after.root
    python scripts/tree_digest.py out.root            # just print the digest

Trap this script exists to avoid: a char-array branch (TLeafC) comes back from
PyROOT as a cppyy LowLevelView whose repr() contains its memory ADDRESS, so
hashing str(value) reports spurious differences on every run. Those leaves must
be read with GetValueString().
"""

from __future__ import annotations

import argparse
import hashlib
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError


def digest_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {path}")
    out = {}
    for key in sorted({k.GetName() for k in f.GetListOfKeys()}):
        obj = f.Get(key)
        if not isinstance(obj, ROOT.TTree):
            continue
        h = hashlib.sha256()
        leaves = [obj.GetLeaf(b.GetName())
                  for b in obj.GetListOfBranches()]
        leaves = [lf for lf in leaves if lf]
        names = [lf.GetName() for lf in leaves]
        h.update(("|".join(names)).encode())
        for i in range(obj.GetEntries()):
            obj.GetEntry(i)
            for lf in leaves:
                if lf.IsA().GetName() == "TLeafC":
                    # char array: must NOT be str()'d, see module docstring
                    h.update(str(lf.GetValueString()).encode())
                else:
                    n = lf.GetLen()
                    for j in range(n):
                        h.update(repr(lf.GetValue(j)).encode())
        out[key] = (obj.GetEntries(), h.hexdigest())
    f.Close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    digs = [digest_file(p) for p in args.files]
    for p, d in zip(args.files, digs):
        print(p)
        for tree, (n, hx) in sorted(d.items()):
            print(f"   {tree:24s} entries {n:>8d}  {hx}")

    if len(digs) < 2:
        return 0
    ref = digs[0]
    ok = True
    for p, d in zip(args.files[1:], digs[1:]):
        if d == ref:
            print(f"\nIDENTICAL: {args.files[0]} == {p}")
        else:
            ok = False
            print(f"\nDIFFERS: {args.files[0]} != {p}")
            for tree in sorted(set(ref) | set(d)):
                if ref.get(tree) != d.get(tree):
                    print(f"   tree {tree}: {ref.get(tree)} vs {d.get(tree)}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
