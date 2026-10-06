"""All ODD admissible classes (every 5-subset spans 1, 3 or 5 edges) of the p = 5 class list (vectorised parity)."""
import sys
from itertools import combinations

import numpy as np

import paths

E7 = sorted(combinations(range(7), 4), key=lambda e: e[::-1])
IDX = {e: i for i, e in enumerate(E7)}
P5 = [sum(1 << IDX[e] for e in combinations(S, 4)) for S in combinations(range(7), 5)]


def popcount(x):
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
    return (x * 0x0101010101010101) >> 56 & 0xFF


def main():
    cls = np.load(paths.CLASSES_P5).astype(np.int64)
    odd = np.ones(cls.size, dtype=bool)
    for p in P5:
        odd &= (popcount(cls & np.int64(p)) & 1) == 1
    print(f"{int(odd.sum())} odd classes of {cls.size}")
    np.save(sys.argv[1], cls[odd])


if __name__ == "__main__":
    main()
