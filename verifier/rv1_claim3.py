"""Independent checker, CLAIM 3 (own code): the certificate does not force oddness.

 part A: the certificate's argmin class 27303369217 (colex): e(H), 5-set profile (does it contain an even 5-set?),
         exact d - c_Q with the own evaluator (must be exactly b).
 part B: all ODD labelled 4-graphs on {0..6}: an odd graph is the complement of a 3-cocycle z of the simplex over F_2
         (every 5-set contains an even number of z-sets); a cocycle is determined by its 20 values on the 4-sets
         containing vertex 0: z(abcd) = z(0bcd) + z(0acd) + z(0abd) + z(0abc).  Generates all 2^20, checks they are
         distinct and odd, writes their colex masks (uint64) for the C evaluator; class count by Burnside over S_7
         (fixed cocycles of each permutation = an F_2 null space).
 part C (after the C scan): exact own evaluation of the lowest odd graphs reported by rv_eval.
Usage: python rv1_claim3.py CERT_PREFIX OUTDIR {AB | C RVEVAL_OUT}
"""
import json
import os
import sys
import time
from fractions import Fraction
from itertools import combinations, permutations

import numpy as np

import rv1_lib as L
from rv1_dual import eval_N, load_cert


def f2_rank(rows):
    rows = [r for r in rows if r]
    rank = 0
    while rows:
        p = max(rows)
        hb = p.bit_length() - 1
        rows = [r ^ p if (r >> hb) & 1 else r for r in rows if r != p]
        rows = [r for r in rows if r]
        rank += 1
    return rank


def main():
    prefix, outdir, mode = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(outdir, f"claim3_{mode}.log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    cert, keys = load_cert(prefix)
    Mfac = int(cert["M"])
    DEN = 5040 * Mfac * Mfac
    b = Fraction(cert["bound"])
    if mode == "AB":
        t0 = time.time()
        h = L.colex_to_own(27303369217)
        prof = L.five_profile(h)
        evens = [P for P, c in zip(combinations(range(7), 5), prof) if c % 2 == 0]
        log("argmin 27303369217: e =", L.popcount(h), " admissible:", L.admissible5(h), " 5-set edge counts:",
            sorted(prof), " even 5-sets:", len(evens), evens)
        val = Fraction(eval_N(h, keys, Mfac), DEN)
        log("own evaluator: d - c =", val, " == b:", val == b)
        # part B
        S0 = [e for e in L.L7 if 0 in e]                 # the 20 4-sets containing 0
        others = [e for e in L.L7 if 0 not in e]
        pos = {e: i for i, e in enumerate(L.L7)}
        # bit tables: cocycle bit of each 4-set as an F_2-linear function (mask over the 20 free bits)
        lin = {}
        for i, e in enumerate(S0):
            lin[e] = 1 << i
        for e in others:
            x = 0
            for drop in range(4):
                f = (0,) + tuple(v for j, v in enumerate(e) if j != drop)
                x ^= lin[f]
            lin[e] = x
        free = np.arange(1 << 20, dtype=np.int64)
        z = np.zeros(1 << 20, dtype=np.int64)
        for e in L.L7:
            par = np.zeros(1 << 20, dtype=np.int64)
            x = lin[e]
            for i in range(20):
                if (x >> i) & 1:
                    par ^= (free >> i) & 1
            z |= par << pos[e]
        hs = (~z) & ((1 << 35) - 1)                    # odd graph = complement of the cocycle (own lex masks)
        log("generated:", len(hs), " distinct:", len(np.unique(hs)) == len(hs))
        oddok = np.ones(len(hs), dtype=bool)
        for F in L.FIVE:
            fm = sum(1 << i for i in F)
            cnt = np.zeros(len(hs), dtype=np.int64)
            x = hs & fm
            for i in F:
                cnt += (x >> i) & 1
            oddok &= (cnt % 2 == 1)
        log("all odd:", bool(oddok.all()))
        # converse sanity: random graphs that are odd must be in the set (checked on all graphs of the support of
        # the dual point and of the Giraud law, which are odd)
        hset = set(int(v) for v in hs)
        D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates",
                                     "n7_dual_point", "dual_Wodd_best.json")))
        log("dual-point support graphs contained:", all(L.colex_to_own(int(x["mask"])) in hset for x in D["support"]))
        G = json.load(open("giraud/giraud_values.json"))
        log("Giraud-law graphs contained:", all(L.colex_to_own(m) in hset for m in G["colex"]))
        # colex export
        cx = np.zeros(len(hs), dtype=np.uint64)
        for i, e in enumerate(L.L7):
            cx |= (((hs >> i) & 1).astype(np.uint64) << np.uint64(L.colex_rank(e)))
        cx.tofile(os.path.join(outdir, "odd_all.u64"))
        # Burnside: cocycle basis = delta(1_t) for the 20 triangles t avoiding 0; z_t(T) = [t subset T]
        basis = []
        for t in combinations(range(1, 7), 3):
            basis.append(sum(1 << pos[e] for e in L.L7 if set(t) <= set(e)))
        assert f2_rank(basis) == 20
        tot = 0
        for g in permutations(range(7)):
            pm = [pos[tuple(sorted(g[v] for v in e))] for e in L.L7]
            rows = []
            for zb in basis:
                img = 0
                for i in range(35):
                    if (zb >> i) & 1:
                        img |= 1 << pm[i]
                rows.append(img ^ zb)
            tot += 2 ** (20 - f2_rank(rows))
        assert tot % 5040 == 0
        log("odd labelled graphs 2^20 =", 1 << 20, "; isomorphism classes by Burnside:", tot // 5040)
        log(f"done in {time.time() - t0:.1f}s")
    else:
        rv = sys.argv[4]
        low = []
        for line in open(rv):
            if line.startswith("LOW "):
                _, Nn, m = line.split()
                low.append((int(Nn), int(m)))
            if line.startswith("NMIN"):
                parts = line.split()
                nmin, argc = int(parts[1]), int(parts[7])
            if line.startswith("TOTAL"):
                log(line.strip())
        low.sort()
        log("rv_eval odd scan: min (d - c) - b =", Fraction(nmin, DEN) - b, "=", float(Fraction(nmin, DEN) - b),
            " claimed 376193/28862180229120:", Fraction(nmin, DEN) - b == Fraction(376193, 28862180229120))
        pm = L.perm_edge_maps()
        seen = {}
        for Nn, m in low[:16]:
            h = L.colex_to_own(m)
            c = L.canonical_graph(h, pm)
            own = eval_N(h, keys, Mfac)
            seen.setdefault(c, []).append((m, Nn, own))
        for c, lst in seen.items():
            log("class", c, "e =", L.popcount(c), " labelled among 16 lowest:", len(lst),
                " own == C numerator:", all(o == n for _, n, o in lst),
                " (d - c) - b =", Fraction(lst[0][2], DEN) - b)
        log("argmin colex", argc, "own value - b =", Fraction(eval_N(L.colex_to_own(argc), keys, Mfac), DEN) - b)


if __name__ == "__main__":
    main()
