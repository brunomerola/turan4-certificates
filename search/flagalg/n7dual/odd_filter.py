"""Keep only ODD graphs (every 5-subset spans 1, 3 or 5 edges) of a mask list; writes a .npy for prep_support.py.
An odd-supported dual point bounds both the plain N = 7 method for t(5,4) and its odd-restricted variant.
Usage: python odd_filter.py IN.npz[:key] OUT.npy
"""
import sys
from itertools import combinations

import numpy as np

from prep_support import load_masks

E7 = sorted(combinations(range(7), 4), key=lambda e: e[::-1])
IDX = {e: i for i, e in enumerate(E7)}
P5 = [sum(1 << IDX[e] for e in combinations(S, 4)) for S in combinations(range(7), 5)]


def is_odd(m):
    return all(bin(m & p).count("1") % 2 == 1 for p in P5)


def main():
    masks = load_masks(sys.argv[1]).astype(np.int64)
    odd = np.array([is_odd(int(m)) for m in masks])
    print(f"{int(odd.sum())} odd of {masks.size}")
    np.save(sys.argv[2], masks[odd])


if __name__ == "__main__":
    main()
