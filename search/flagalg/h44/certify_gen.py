"""Exact certificate (../n7/certify7.py, unchanged logic) for the threshold problem of gen.py: certify7.main with
the class list from data/l{LAM}/ and the raw scan's extension filter replaced by the threshold version
(popcount(h & pset) >= LAM on the p-sets through vertex 6).  Bound semantics as in certify7 (tau = 0 runs: plain).
Usage: python certify_gen.py LAM RUN_PREFIX [certify7 options]
"""
import sys

import gen

LAM = int(sys.argv[1])
import numpy as np  # noqa: E402

P = int(np.load(sys.argv[2] + ".best.npz")["p"])
gen.configure(P, LAM)

import certify7  # noqa: E402

certify7.DATA = gen.data_dir(LAM)
certify7.scan_raw_exact = gen.make_raw_kernels(LAM)
sys.argv = [sys.argv[0]] + sys.argv[2:]
print(f"# certify_gen: predicate 'every {P}-set spans >= {LAM} edges'; classes {certify7.DATA}", flush=True)
certify7.main()
