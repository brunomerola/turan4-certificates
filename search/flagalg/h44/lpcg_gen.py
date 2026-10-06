"""LP-CUT-CG (../n7/lpcg7.py, unchanged) for the threshold problem of gen.py.  Imports highspy (via lpcg7); never
import ortools here.
Usage: python lpcg_gen.py LAM P [lpcg7 options ...]     e.g.  python lpcg_gen.py 2 5 --no-stationarity --out results/x
The universe is data/l{LAM}/classes_p{P}.npy (enum_gen.py); everything else is lpcg7.main.
"""
import sys

import gen

LAM = int(sys.argv[1])
P = int(sys.argv[2])
gen.configure(P, LAM)

import lpcg7  # noqa: E402

lpcg7.DATA = gen.data_dir(LAM)
sys.argv = [sys.argv[0]] + sys.argv[2:]
print(f"# lpcg_gen: predicate 'every {P}-set spans >= {LAM} edges'; universe dir {lpcg7.DATA}", flush=True)
lpcg7.main()
