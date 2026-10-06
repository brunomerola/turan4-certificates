"""Independent checker (own code): 6-vertex representatives for the raw scan, and (optional) a class list as raw uint64.

1. All 2^15 4-graphs on {0..5}; canonical form = min over all 720 vertex permutations (brute force, numpy);
   admissible = every 5-subset of {0..5} contains an edge.  reps = canonical forms of admissible graphs, written as
   OWN 7-vertex masks (lex order of 4-subsets of {0..6}; vertex 6 unused) to local/reps6_p5.bin (uint64).
   Completeness of the raw scan: any admissible 7-vertex H has H - 6 admissible and isomorphic (by a permutation of
   {0..5}) to a rep R; applying it to H gives R + (link of 6), which the raw scan enumerates; the slack is an
   isomorphism invariant (rs_eval.c sums over ALL ordered root tuples).
2. Labelled count of admissible 7-vertex graphs: L7 = sum_R (720 / |Aut R|) * #{links : R + link admissible}.
3. local/classes_p5.u64: payload of the producer class list (uint64 COLEX masks), and its count (only if a class
   list path is given; the package does not include the class list).
"""
import json
import sys
from itertools import combinations, permutations

import numpy as np

E7 = list(combinations(range(7), 4))
I7 = {e: i for i, e in enumerate(E7)}
E6 = list(combinations(range(6), 4))
I6 = {e: i for i, e in enumerate(E6)}
PC = np.array([bin(v).count("1") for v in range(1 << 16)], dtype=np.int8)


def popc(a):
    return PC[a & 0xFFFF] + PC[(a >> 16) & 0xFFFF] + PC[(a >> 32) & 0xFFFF]


def main():
    cls_path = sys.argv[1] if len(sys.argv) > 1 else None
    X = np.arange(1 << 15, dtype=np.int64)
    canon = X.copy()
    for g in permutations(range(6)):
        img = [I6[tuple(sorted(g[v] for v in e))] for e in E6]
        Y = np.zeros_like(X)
        for j in range(15):
            Y |= ((X >> j) & 1) << img[j]
        canon = np.minimum(canon, Y)
    adm6 = np.ones(X.size, dtype=bool)
    for P in combinations(range(6), 5):
        adm6 &= (X & sum(1 << I6[e] for e in combinations(P, 4))) != 0
    reps6 = np.unique(canon[adm6])
    orbit = {int(r): int((canon == r).sum()) for r in reps6}
    assert sum(orbit.values()) == int(adm6.sum())
    to7 = [I7[e] for e in E6]
    reps7 = np.array([sum(1 << to7[j] for j in range(15) if (int(r) >> j) & 1) for r in reps6], dtype=np.uint64)
    reps7.tofile("local/reps6_p5.bin")
    link = [I7[t + (6,)] for t in combinations(range(6), 3)]
    Lr = np.arange(1 << 20, dtype=np.int64)
    LH = np.zeros_like(Lr)
    for j in range(20):
        LH |= ((Lr >> j) & 1) << link[j]
    PS = [sum(1 << I7[e] for e in combinations(P, 4)) for P in combinations(range(7), 5)]
    raw = 0
    L7 = 0
    for r6, r7 in zip(reps6, reps7):
        Hm = LH | np.int64(r7)
        ok = np.ones(Hm.size, dtype=bool)
        for m in PS:
            ok &= (Hm & m) != 0
        c = int(ok.sum())
        raw += c
        L7 += orbit[int(r6)] * c
    out = {"admissible_6vertex_labelled": int(adm6.sum()), "reps6": int(reps6.size),
           "raw_admissible_extensions": raw, "labelled_admissible_7vertex": L7}
    if cls_path is not None:
        cls = np.load(cls_path)
        assert cls.dtype == np.uint64
        cls.tofile("local/classes_p5.u64")
        out.update({"class_list": cls_path, "class_list_count": int(cls.size),
                    "class_list_unique": int(np.unique(cls).size)})
    json.dump(out, open("local/reps6_p5.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
