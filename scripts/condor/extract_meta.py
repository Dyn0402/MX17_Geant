#!/usr/bin/env python3
"""Print the `meta` JSON of a greens_comb_*.npz without reading its arrays.

An S1 product is ~2.3 GB but its provenance/acceptance block is a few kB. npz
is a zip, so the meta member can be pulled on its own — which is what makes it
cheap to check acceptance at CERN and keep the bulk transfer to the one or two
products actually needed for downstream work.

    python3 extract_meta.py <product.npz> [> meta.json]
"""
import io
import sys
import zipfile

import numpy as np


def main():
    z = zipfile.ZipFile(sys.argv[1])
    m = np.load(io.BytesIO(z.read("meta.npy")), allow_pickle=True)
    print(str(m))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
