"""Independent checker: exhaustive checks of Lemma 1 and Lemma 2 (form F) for m <= 7 points, own code.

Form F (REPORT Sec. 2a) for a 3-graph L on U = range(m) and A' subset U, B = U - A':
 (F1) every triple inside A' is in L; (F2) no triple with exactly 2 points in A' is in L; (F3) no triple inside B
 is in L; (F4) for a in A' and b1 < b2 < b3 in B an odd number of a b1 b2, a b1 b3, a b2 b3 are in L.
Checks, for m = 3..7:
 (1) the set of 3-graphs generated from the definition of form F (all A', all cut choices P_a) equals the set of
     links of the last vertex of all labelled Giraud systems on m+1 vertices (Lemma 1 (i) and (ii));
 (2) the number of such links equals |G_{m+1}| (an odd 4-graph is determined by one link);
 (3) for every link of g = Gamma(A,B,M) in G_{m+1}: the part of v = m minus v satisfies the predicate (Lemma 1 (i));
 (4) for every form-F 3-graph, the set of ALL A' satisfying the predicate (all 2^m subsets tested) either has one
     element or all its elements have size <= 3 (Lemma 2);  m = 3..7 (m = 7: 32,981 graphs x 128 subsets).
 (5) m = 3, 4, 5: every one of the 2^C(m,3) 3-graphs is tested with the predicate against all subsets, and the
     form-F ones are exactly the generated ones (completeness of the generation).
Output: results_formF.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.getcwd()   # outputs (results_*.json, local/) go to the working directory
import st_lib as L  # noqa: E402


def pred(Lm, m, Ap):
    Ap = set(Ap)
    for t in combinations(range(m), 3):
        k = sum(1 for x in t if x in Ap)
        e = (Lm >> L.cidx(t)) & 1
        if (k == 3 and not e) or (k in (0, 2) and e):
            return False
    B = [x for x in range(m) if x not in Ap]
    for a in Ap:
        for b1, b2, b3 in combinations(B, 3):
            s = ((Lm >> L.cidx(tuple(sorted((a, b1, b2))))) & 1) + ((Lm >> L.cidx(tuple(sorted((a, b1, b3))))) & 1) \
                + ((Lm >> L.cidx(tuple(sorted((a, b2, b3))))) & 1)
            if s % 2 == 0:
                return False
    return True


def generate(m):
    out = {}
    for am in range(1 << m):
        Ap = [x for x in range(m) if (am >> x) & 1]
        B = [x for x in range(m) if not (am >> x) & 1]
        base = 0
        for t in combinations(Ap, 3):
            base |= 1 << L.cidx(t)
        opts = []
        for a in Ap:
            s = set()
            for pm in range(1 << len(B)):
                P = {B[j] for j in range(len(B)) if (pm >> j) & 1}
                v = 0
                for b, b2 in combinations(B, 2):
                    if (b in P) == (b2 in P):
                        v |= 1 << L.cidx(tuple(sorted((a, b, b2))))
                s.add(v)
            opts.append(s)
        cur = {base}
        for s in opts:
            cur = {x | y for x in cur for y in s}
        for v in cur:
            out.setdefault(v, set()).add(frozenset(Ap))
    return out


def main():
    t0 = time.time()
    res = {}
    for m in range(3, 8):
        n = m + 1
        if n == 7:
            G = [int(x) for x in open(os.path.join(OUT, "local", "G7.txt")).read().split()]
            fam = None
        elif n == 8:
            G = [int(x) for x in open(os.path.join(OUT, "local", "G8.txt")).read().split()]
            fam = None
        else:
            fam = L.enum_giraud(n)
            G = sorted(fam)
        n4 = comb(m, 4)
        links = {g >> n4 for g in G}
        gen = generate(m)
        r = {"G_{m+1}": len(G), "distinct_links": len(links), "generated_formF": len(gen),
             "generated_equals_links": set(gen) == links}
        # (3) Lemma 1 (i) on every labelled Giraud system (partitions from a fresh enumeration)
        famn = fam if fam is not None else L.enum_giraud(n)
        ok3 = True
        for g, parts in famn.items():
            Lm = g >> n4
            for A in parts:
                P = A if m in A else frozenset(range(n)) - A
                ok3 &= pred(Lm, m, P - {m})
        r["lemma1_i_ok"] = ok3
        # (4) Lemma 2 with ALL subsets
        multi = 0
        ok4 = True
        gen_reps_ok = True
        for Lm in sorted(links):
            reps = [frozenset(x for x in range(m) if (am >> x) & 1) for am in range(1 << m)
                    if pred(Lm, m, [x for x in range(m) if (am >> x) & 1])]
            gen_reps_ok &= set(reps) == gen.get(Lm, set())
            if len(reps) > 1:
                multi += 1
                ok4 &= all(len(A) <= 3 for A in reps)
        r.update({"with_several_reps": multi, "lemma2_ok": ok4, "generated_reps_equal_predicate_reps": gen_reps_ok})
        # (5) completeness for small m
        if m <= 5:
            allf = set()
            for Lm in range(1 << comb(m, 3)):
                if any(pred(Lm, m, [x for x in range(m) if (am >> x) & 1]) for am in range(1 << m)):
                    allf.add(Lm)
            r["all_3graphs_formF_equals_generated"] = allf == set(gen)
        res[str(m)] = r
        print(m, r, round(time.time() - t0, 1), flush=True)
    res["time_s"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(OUT, "results_formF.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
