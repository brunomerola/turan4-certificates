"""Independent checker: CF1-CF3 with own code (st_lib).  Labelled Giraud systems on m = 7, 8, 9 vertices, their isomorphism
classes (orbit peeling), partition multiplicities, the co-Fano and SQS(8) classes, and the match of the 7-vertex
classes with the certificate's tight classes.  Also self-tests of st_lib (f identities, the two Giraud tests).
Output: results_cf.json, local/G7.txt, local/G8.txt, local/G9.txt (labelled masks), local/reps{7,8,9}.json.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from itertools import combinations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.getcwd()   # outputs (results_*.json, local/) go to the working directory
import st_lib as L  # noqa: E402

TIGHT = [2290237506, 17418883036, 19736566032, 26383237633, 26894418433, 29337331976, 31813832752, 31837153552,
         34326183967, 34359738367]     # sharp_klp5.tight.json (exact_tight_classes)


def selftests(rng):
    out = {}
    # f via 5-sets == definitional 2d - P on random graphs, n = 6..9
    ok = True
    for n in range(6, 10):
        for _ in range(15):
            mk = rng.getrandbits(comb(n, 4))
            ok &= L.f5(mk, n) == L.f_def(mk, n)
    # on 7 vertices: (24 e - P)/420
    for _ in range(100):
        mk = rng.getrandbits(35)
        e = bin(mk).count("1")
        P = 0
        for T in combinations(range(7), 3):
            c = sum(L.bit(mk, T + (x,)) for x in range(7) if x not in T)
            P += c * (c - 1)
        ok &= L.f5(mk, 7) == L.Fraction(24 * e - P, 420)
    out["f5_equals_definition"] = ok
    # the two Giraud tests agree: random Giraud systems, their single/double flips, random graphs (n = 7..9)
    agree = True
    npos = nneg = 0
    for n in (7, 8, 9):
        for _ in range(40):
            a = rng.randrange(0, n + 1)
            A = list(range(n))
            rng.shuffle(A)
            A = A[:a]
            M = {(x, y): rng.randrange(2) for x in A for y in range(n) if y not in A}
            g = L.giraud(n, A, M)
            for h in (g, g ^ (1 << rng.randrange(comb(n, 4))),
                      g ^ (1 << rng.randrange(comb(n, 4))) ^ (1 << rng.randrange(comb(n, 4))),
                      rng.getrandbits(comb(n, 4))):
                r1 = sorted(map(sorted, L.is_giraud_forced(h, n)))
                r2 = sorted(map(sorted, L.is_giraud_struct(h, n)))
                agree &= r1 == r2
                if r1:
                    npos += 1
                else:
                    nneg += 1
            agree &= bool(L.is_giraud_forced(g, n))
            agree &= L.is_odd(g, n)
    out["giraud_tests_agree"] = agree
    out["giraud_tests_pos_neg"] = [npos, nneg]
    return out


def multiplicity(fam, n):
    res = {}
    for mk, parts in fam.items():
        if len(parts) > 1:
            res.setdefault(len(parts), []).append(mk)
    return res


def main():
    t0 = time.time()
    rng = random.Random(7_20261006)
    out = {"selftests": selftests(rng)}
    print(out, flush=True)
    # ------------------------------------------------------------------ m = 7
    g7n = L.enum_giraud(7, normalised=True)
    g7f = L.enum_giraud(7, normalised=False)
    o7 = L.orbits(set(g7n), 7)
    tight_in = [m in g7n for m in TIGHT]
    tab7 = L.perm_quad_table(7)
    # orbit of every tight mask; which orbit each falls in
    orbit_of_tight = []
    orb_sets = []
    for rep, _ in o7:
        row = rep
        orb = {L.relabel(row, 7, p) for p in __import__("itertools").permutations(range(7))}
        orb_sets.append(orb)
    for m in TIGHT:
        orbit_of_tight.append([i for i, s in enumerate(orb_sets) if m in s])
    multi7 = multiplicity(g7n, 7)
    r7 = {"labelled_normalised": len(g7n), "labelled_all_matrices": len(g7f), "same_set": set(g7n) == set(g7f),
          "partition_sets_equal": all(g7n[k] == g7f[k] for k in g7n),
          "classes": len(o7), "orbit_sizes": sorted(s for _, s in o7), "all_odd": all(L.is_odd(m, 7) for m in g7n),
          "tight_masks_are_giraud": all(tight_in),
          "each_tight_mask_in_exactly_one_orbit_and_bijective":
              sorted(x[0] for x in orbit_of_tight if len(x) == 1) == list(range(10)),
          "multi_partition": {str(k): len(v) for k, v in multi7.items()}}
    # the multi-partition class(es)
    mp = [m for v in multi7.values() for m in v]
    cls_mp = {i for i, s in enumerate(orb_sets) for m in mp if m in s}
    r7["multi_partition_orbits"] = sorted(cls_mp)
    info = []
    for i in sorted(cls_mp):
        rep = min(orb_sets[i])
        parts = g7n[rep]
        e = bin(rep).count("1")
        comps = [tuple(sorted(set(range(7)) - set(L.subsets(7, 4)[j]))) for j in range(35) if (rep >> j) & 1]
        fano = all(sum(1 for c in comps if set(pr) <= set(c)) == 1 for pr in combinations(range(7), 2))
        info.append({"rep": rep, "edges": e, "partitions": len(parts),
                     "part_sizes": sorted(sorted((len(p), 7 - len(p))) for p in parts),
                     "five_counts_hist": [L.five_counts(rep, 7).count(j) for j in range(6)],
                     "complements_of_edges_form_Fano_plane": fano,
                     "orbit_size": len(orb_sets[i]),
                     "orbit_contains_tight_mask": [m for m in TIGHT if m in orb_sets[i]]})
    r7["multi_partition_info"] = info
    out["m7"] = r7
    print(json.dumps(r7), flush=True)
    open(os.path.join(OUT, "local", "G7.txt"), "w").write("\n".join(map(str, sorted(g7n))) + "\n")
    json.dump([[r, s] for r, s in o7], open(os.path.join(OUT, "local", "reps7.json"), "w"))
    # ------------------------------------------------------------------ m = 8
    g8n = L.enum_giraud(8, normalised=True)
    t1 = time.time()
    g8f = L.enum_giraud(8, normalised=False)
    print("m8 full enumeration", round(time.time() - t1, 1), flush=True)
    tab8 = L.perm_quad_table(8)
    o8 = L.orbits(set(g8n), 8, tab8)
    multi8 = multiplicity(g8n, 8)
    mp8 = sorted({m for v in multi8.values() for m in v})
    # orbits of the multi-partition graphs
    mp_orbits = []
    seen = set()
    for m in mp8:
        if m in seen:
            continue
        orb = set()
        import numpy as np
        row = np.array([(m >> i) & 1 for i in range(70)], dtype=bool)
        imgs = np.zeros(tab8.shape, dtype=bool)
        np.put_along_axis(imgs, tab8.astype(np.intp), np.broadcast_to(row, tab8.shape), axis=1)
        orb = set(L._pack_rows(np.unique(imgs, axis=0)))
        seen |= orb
        rep = min(orb)
        parts = g8n[rep]
        blocks = [L.subsets(8, 4)[j] for j in range(70) if (rep >> j) & 1]
        steiner = all(sum(1 for b in blocks if set(t) <= set(b)) == 1 for t in combinations(range(8), 3))
        # AG(3,2): with some labelling of the 8 points by F_2^3 the blocks are the affine planes: equivalently every
        # block is a 4-set and the symmetric difference of two disjoint blocks is empty / blocks closed under
        # complement (parallel classes); we test: complement of every block is a block, and the 7 parallel classes
        comp_closed = all(tuple(sorted(set(range(8)) - set(b))) in set(blocks) for b in blocks)
        mp_orbits.append({"rep": rep, "orbit_size": len(orb), "edges": len(blocks), "partitions": len(parts),
                          "part_sizes": sorted(sorted((len(p), 8 - len(p))) for p in parts),
                          "partition_parts_are_blocks": all(tuple(sorted(p)) in set(blocks) for p in parts),
                          "five_counts_hist": [L.five_counts(rep, 8).count(j) for j in range(6)],
                          "steiner_quadruple_system_S(3,4,8)": steiner, "blocks_closed_under_complement": comp_closed,
                          "all_multi_partition_graphs_in_this_orbit": set(mp8) <= orb})
    r8 = {"labelled_normalised": len(g8n), "labelled_all_matrices": len(g8f), "same_set": set(g8n) == set(g8f),
          "partition_sets_equal": all(g8n[k] == g8f[k] for k in g8n),
          "classes": len(o8), "orbit_sizes": sorted(s for _, s in o8), "all_odd": all(L.is_odd(m, 8) for m in g8n),
          "multi_partition": {str(k): len(v) for k, v in multi8.items()}, "multi_partition_orbits": mp_orbits}
    out["m8"] = r8
    print(json.dumps(r8), flush=True)
    del g8f
    open(os.path.join(OUT, "local", "G8.txt"), "w").write("\n".join(map(str, sorted(g8n))) + "\n")
    json.dump([[r, s] for r, s in o8], open(os.path.join(OUT, "local", "reps8.json"), "w"))
    # ------------------------------------------------------------------ m = 9
    t1 = time.time()
    g9n = L.enum_giraud(9, normalised=True)
    print("m9 enumeration", len(g9n), round(time.time() - t1, 1), flush=True)
    tab9 = L.perm_quad_table(9)
    print("tab9", tab9.shape, round(time.time() - t1, 1), flush=True)
    o9 = L.orbits(set(g9n), 9, tab9)
    print("orbits9", len(o9), round(time.time() - t1, 1), flush=True)
    multi9 = multiplicity(g9n, 9)
    # odd: sample (all 604,426 would be slow in Python); oddness of every Giraud system is Lemma odd (hand)
    smp = random.Random(9).sample(sorted(g9n), 2000)
    r9 = {"labelled_normalised": len(g9n), "classes": len(o9), "orbit_sizes": sorted(s for _, s in o9),
          "orbit_sizes_sum": sum(s for _, s in o9), "multi_partition": {str(k): len(v) for k, v in multi9.items()},
          "odd_on_2000_sample": all(L.is_odd(m, 9) for m in smp)}
    out["m9"] = r9
    print(json.dumps(r9), flush=True)
    open(os.path.join(OUT, "local", "G9.txt"), "w").write("\n".join(map(str, sorted(g9n))) + "\n")
    json.dump([[r, s] for r, s in o9], open(os.path.join(OUT, "local", "reps9.json"), "w"))
    out["time_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(OUT, "results_cf.json"), "w"), indent=1)
    print("DONE", out["time_s"])


if __name__ == "__main__":
    main()
