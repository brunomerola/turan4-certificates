"""Independent checker: merge rs_eval.c scan logs (possibly several LO-HI parts), recombine the limbs EXACTLY, and report:
min slack, number of graphs at the bound, next value, per-F groups; compare the per-F maxima of Zsum with the
producer's stored groups (sharp_klp5.json, informational); classify the tight graphs (S == 0) against the orbits of
Giraud's 10 classes (all 5040 relabellings, own code).
Usage: python rs_analyze.py BLOB_JSON TARGET PRODUCER_JSON|- TIGHT_FILES(comma)|- LOG [LOG ...]
"""
import json
import sys
from fractions import Fraction
from itertools import combinations, permutations

import numpy as np

from rs_lib import colex_rank

LB = 45
GIRAUD = [2290237506, 17418883036, 19736566032, 26383237633, 26894418433, 29337331976, 31813832752, 31837153552,
          34326183967, 34359738367]


def val(limbs):
    return sum(int(v) << (LB * j) for j, v in enumerate(limbs))


def parse(log):
    out = {"G": {}, "LOW": []}
    for line in open(log):
        t = line.split()
        if not t:
            continue
        if t[0] == "TOTAL":
            out["TOTAL"] = {kv.split("=")[0]: kv.split("=")[1] for kv in t[1:]}
        elif t[0] in ("MIN", "NEXT", "LOW"):
            i = t.index("|")
            rec = (val(t[1:i]), int(t[i + 4]))
            if t[0] == "LOW":
                out["LOW"].append(rec)
            else:
                out[t[0]] = rec
        elif t[0] == "G":
            i = t.index("|")
            out["G"][int(t[1])] = (int(t[2]), int(t[3]), val(t[4:i]), int(t[i + 4]))
        elif t[0] == "RANGE":
            out["RANGE"] = (int(t[1]), int(t[2]))
    return out


def own_to_colex_table():
    E = list(combinations(range(7), 4))
    return [colex_rank(e) for e in E]


def orbit_colex(h):
    E = list(combinations(range(7), 4))
    es = [e for e in E if (h >> colex_rank(e)) & 1]
    return {sum(1 << colex_rank(tuple(sorted(g[v] for v in e))) for e in es) for g in permutations(range(7))}


def main():
    bj = json.load(open(sys.argv[1]))
    b = Fraction(sys.argv[2])
    pj = None if sys.argv[3] == "-" else json.load(open(sys.argv[3]))
    tfiles = [] if sys.argv[4] == "-" else sys.argv[4].split(",")
    logs = sys.argv[5:]
    D = int(bj["D"])
    T = 5040 * D * b
    assert T.denominator == 1
    T = int(T)
    parts = [parse(l) for l in logs]
    tot = {k: sum(int(p["TOTAL"][k].rstrip("s")) for p in parts) for k in ("admissible", "rejected", "zero", "negative")}
    mn = min((p["MIN"] for p in parts), key=lambda r: r[0])
    nx = [p["NEXT"] for p in parts if "NEXT" in p]
    nx = min(nx, key=lambda r: r[0]) if nx else None
    G = {}
    for p in parts:
        for F, (c, z, s, h) in p["G"].items():
            if F not in G:
                G[F] = [c, z, s, h]
            else:
                G[F][0] += c
                G[F][1] += z
                if s < G[F][2]:
                    G[F][2], G[F][3] = s, h
    slack = lambda S: b + Fraction(S, 5040 * D)
    rep = {"logs": logs, "ranges": [p.get("RANGE") for p in parts], "totals": tot,
           "min_S": str(mn[0]), "min_slack": str(slack(mn[0])), "min_slack_equals_b": slack(mn[0]) == b,
           "argmin_colex": mn[1]}
    if nx:
        rep.update({"next_S": str(nx[0]), "next_slack_minus_b": str(Fraction(nx[0], 5040 * D)),
                    "next_slack_minus_b_float": float(Fraction(nx[0], 5040 * D)), "next_argmin_colex": nx[1]})
    rep["groups"] = len(G)
    rep["groups_with_zero"] = sorted(F for F in G if G[F][1])
    rep["zero_counts_by_F"] = {F: G[F][1] for F in sorted(G) if G[F][1]}
    # per-F max Zsum = 12 D F - T - minS(F)
    if pj is not None:
        prod = {int(g[0]): (int(g[1]), int(g[2]), int(g[3])) for g in pj["groups"]}
        same_keys = set(prod) == set(G)
        same_Z = same_keys and all(12 * D * F - T - G[F][2] == prod[F][0] for F in G)
        same_cnt = same_keys and all(G[F][0] == prod[F][2] for F in G)
        rep["producer_groups_same_F"] = same_keys
        rep["producer_groups_same_maxZ"] = same_Z
        rep["producer_groups_same_counts"] = same_cnt
    # lowest values (class scan)
    lows = sorted(set(r for p in parts for r in p["LOW"]))
    if lows:
        rep["lowest"] = [(str(Fraction(S, 5040 * D)), float(Fraction(S, 5040 * D)), h) for S, h in lows[:30]]
    # tight graphs
    if tfiles:
        tm = np.concatenate([np.fromfile(f, dtype=np.uint64) for f in tfiles])
        assert tm.size == tot["zero"], (tm.size, tot["zero"])
        E = list(combinations(range(7), 4))
        cr = [colex_rank(e) for e in E]
        def o2c(h):
            return sum(1 << cr[i] for i in range(35) if (h >> i) & 1)
        orb = {g: orbit_colex(g) for g in GIRAUD}
        cls = {}
        uniq = set(int(x) for x in tm)
        for h in uniq:
            c = o2c(h)
            hit = [g for g in GIRAUD if c in orb[g]]
            assert len(hit) <= 1
            cls[hit[0] if hit else None] = cls.get(hit[0] if hit else None, 0) + 1
        rep["tight_graphs"] = int(tm.size)
        rep["tight_distinct_labelled"] = len(uniq)
        rep["tight_by_giraud_class"] = {str(k): v for k, v in cls.items()}
        rep["tight_all_giraud"] = None not in cls
        rep["tight_classes_hit"] = sum(1 for k in cls if k is not None)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
