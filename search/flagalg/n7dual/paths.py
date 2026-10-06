"""Paths shared by the n7dual scripts (search for the N = 7 dual point).  Producer code, nothing here is certified."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FLAGALG = os.path.dirname(HERE)
N7 = os.path.join(FLAGALG, "n7")
for p in (FLAGALG, N7):
    if p not in sys.path:
        sys.path.insert(0, p)

# read-only inputs from the n7 directory (large files, not included in this package)
MAIN_N7 = os.environ.get("N7DUAL_MAIN_N7", N7)
CLASSES_P5 = os.path.join(MAIN_N7, "data", "classes_p5.npy")
STATE_CONV = os.path.join(MAIN_N7, "results", "lpcg_p5_conv_full.state.npz")
RES = os.path.join(HERE, "results")
LOGS = os.path.join(HERE, "logs")
