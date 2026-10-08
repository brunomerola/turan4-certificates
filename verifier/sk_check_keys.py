#!/usr/bin/env python3
"""Self-contained checker of the sidecar key lists <prefix>.keys.json (producer code; no import from the
flag-algebra engine -- Python + numpy integer arrays only; < 100 MB).

For every keys.json given on the command line (default: */*.keys.json next to this script):
  1. the copied certificate files next to it have the sha256 recorded in keys.json, and so do the originals
     (source_dir, if present; the originals are not part of the data package);
  2. the admissibility predicate is re-stated HERE from the definitions (forbidden 3-graphs of BCL Table 2, or the r = 4
     threshold "every p-set spans >= lam edges") and must equal the one recorded in keys.json;
  3. for every block (s, m) of the certificate, the key list is REBUILT from the definitions: types = admissible
     r-graphs on [s] in canonical form (minimum mask over all permutations of [s]), increasing, kept iff they have a
     flag; flags of a type sigma = admissible r-graphs on [m] whose restriction to [s] is sigma, in canonical form
     (minimum mask over the permutations of [m] fixing 0..s-1), increasing; |Aut(sigma)| = number of permutations of
     [s] fixing sigma.  Masks: bit i = the i-th r-subset in colex order (so the r-subsets of [s] come first).
     Asserted equal to keys.json (order of keys, s, m, sigma, |Aut|, flags) and to cert.json's keys (s, m, sigma,
     n_flags);
  4. cert.npz: one integer array A<k> per key with A<k>.shape[1] = n_flags(k) and shape[0] = keys.json factor_rows;
     cert.json M = 2^M_bits.
Usage: python check_keys.py [KEYS.json ...]          exit code 0 iff every file passes
"""
import glob
import hashlib
import json
import math
import os
import sys
from itertools import combinations, permutations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIGCAT = os.path.dirname(HERE)

# BCL (arXiv:2108.10406v3) Table 2 edge lists (1-based), as forbidden (not necessarily induced) subgraphs
FORB = {
    "K4m": "123 124 134", "K5": "123 124 125 134 135 145 234 235 245 345",
    "K5m": "124 125 134 135 145 234 235 245 345", "K5eq": "124 125 134 135 145 234 235 245",
    "K5lt": "123 124 125 134 135 145 234 235", "F32": "123 145 245 345", "F5": "123 124 345",
    "J4": "123 124 125 134 135 145", "J5": "123 124 125 126 134 135 136 145 146 156", "C5": "123 234 345 145 125",
    "C5m": "123 234 345 145", "K4": "123 124 134 234",
    "K6": " ".join("".join(str(v + 1) for v in e) for e in combinations(range(6), 3)),
}
FAMILY = {"s3_K5lt": ["K5lt"], "s3_K5eq": ["K5eq"], "s3_C5": ["C5"], "s3_K5m": ["K5m"], "s3_K6": ["K6"],
          "s3_J4": ["J4"], "s3_J5": ["J5"]}
THRESH = {"c4_K5m": (5, 2), "c4_K6m": (6, 2), "c4_K6": (6, 1), "c4_K7": (7, 1)}


def colex(n, r):
    return sorted(combinations(range(n), r), key=lambda e: e[::-1])


def bits(masks, nb):
    return ((np.asarray(masks, dtype=np.int64)[:, None] >> np.arange(nb)[None, :]) & 1).astype(np.int64)


def permute(masks, n, r, g):
    """Masks after relabelling vertex v -> g[v] (vectorised)."""
    E = colex(n, r)
    pos = {e: i for i, e in enumerate(E)}
    tgt = np.array([pos[tuple(sorted(g[v] for v in e))] for e in E], dtype=np.int64)
    B = bits(masks, len(E))
    return (B << tgt[None, :]).sum(axis=1)


def canon(masks, n, r, perms):
    out = None
    for g in perms:
        p = permute(masks, n, r, g)
        out = p if out is None else np.minimum(out, p)
    return out


def admissible(masks, n, r, prob):
    masks = np.asarray(masks, dtype=np.int64)
    ok = np.ones(masks.size, dtype=bool)
    E = colex(n, r)
    pos = {e: i for i, e in enumerate(E)}
    if prob in THRESH:
        p, lam = THRESH[prob]
        if n < p:
            return ok
        for P in combinations(range(n), p):
            pm = sum(1 << pos[e] for e in combinations(P, r))
            ok &= np.bitwise_count(masks & pm) >= lam
        return ok
    for g in FAMILY[prob]:
        es = [tuple(int(ch) - 1 for ch in w) for w in FORB[g].split()]
        nv = 1 + max(max(e) for e in es)
        if nv > n:
            continue
        copies = set()
        for inj in permutations(range(n), nv):
            copies.add(sum(1 << pos[tuple(sorted(inj[v] for v in e))] for e in es))
        for c in copies:
            ok &= (masks & c) != c
    return ok


def rebuild(prob, r, s, m):
    """[(sigma, aut_size, flags)] for block (s, m), from the definitions."""
    ns, nm = math.comb(s, r), math.comb(m, r)
    lab = np.arange(1 << ns, dtype=np.int64)
    lab = lab[admissible(lab, s, r, prob)]
    perms_s = list(permutations(range(s)))
    types = np.unique(canon(lab, s, r, perms_s)) if ns else np.array([0], dtype=np.int64)
    fix = [tuple(range(s)) + p for p in permutations(range(s, m))]
    free = np.arange(1 << (nm - ns), dtype=np.int64)
    out = []
    for sig in types.tolist():
        cand = sig | (free << ns)
        cand = cand[admissible(cand, m, r, prob)]
        if cand.size == 0:
            continue
        fl = np.unique(canon(cand, m, r, fix))
        aut = sum(1 for g in perms_s if int(permute(np.array([sig]), s, r, g)[0]) == sig) if ns else math.factorial(s)
        out.append((int(sig), aut, [int(x) for x in fl]))
    return out


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def check(path):
    K = json.load(open(path))
    d = os.path.dirname(path)
    res = {"file": os.path.relpath(path, HERE), "problem": K["problem"], "objective": K["objective"]}
    errs = []
    for name, h in K["copied_files_sha256"].items():
        if sha(os.path.join(d, name)) != h:
            errs.append(f"copy {name} sha256 differs")
        orig = os.path.join(SIGCAT, K["source_dir"], name)
        if os.path.exists(orig) and sha(orig) != h:
            errs.append(f"original {name} sha256 differs")
    obj = K["objective"]
    cj = json.load(open(os.path.join(d, f"{K['certificate']}.{obj}.cert.json")))
    cz = np.load(os.path.join(d, f"{K['certificate']}.{obj}.cert.npz"))
    prob, r = K["problem"], K["r"]
    pr = K["predicate"]
    if prob in THRESH:
        if (pr.get("p"), pr.get("lam")) != THRESH[prob]:
            errs.append("predicate differs")
    elif pr.get("forbidden") != [FORB[g] for g in FAMILY[prob]]:
        errs.append("forbidden list differs")
    if cj["M"] != 2 ** cj["M_bits"] or K["M_bits"] != cj["M_bits"]:
        errs.append("M inconsistent")
    built = []
    for (s, m) in [tuple(b) for b in K["blocks"]]:
        for sig, aut, fl in rebuild(prob, r, s, m):
            built.append((s, m, sig, aut, fl))
    kk = K["keys"]
    if len(built) != len(kk) or len(cj["keys"]) != len(kk):
        errs.append(f"number of keys: rebuilt {len(built)}, keys.json {len(kk)}, cert.json {len(cj['keys'])}")
    nA = len([f for f in cz.files if f.startswith("A")])
    if nA != len(kk):
        errs.append(f"cert.npz has {nA} factor arrays for {len(kk)} keys")
    for k, (b, e, c) in enumerate(zip(built, kk, cj["keys"])):
        s, m, sig, aut, fl = b
        if (e["k"], e["s"], e["m"], e["sigma"], e["aut_size"], e["flags"], e["n_flags"]) != (k, s, m, sig, aut, fl, len(fl)):
            errs.append(f"key {k}: rebuilt list differs from keys.json")
        if (c["s"], c["m"], c["sigma"], c["n_flags"]) != (s, m, sig, len(fl)):
            errs.append(f"key {k}: cert.json key differs")
        A = cz[f"A{k}"]
        if not np.issubdtype(A.dtype, np.integer) or A.ndim != 2 or A.shape[1] != len(fl) or A.shape[0] != e["factor_rows"]:
            errs.append(f"key {k}: A{k} shape {A.shape} vs n_flags {len(fl)}, rows {e['factor_rows']}")
    res.update({"keys": len(kk), "flags": sum(len(b[4]) for b in built), "bound": K["bound"],
                "sigma_upper": K["sigma_upper"], "PASS": not errs, "errors": errs[:10]})
    return res


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, "*", "*.keys.json")))
    ok = True
    for f in files:
        r = check(f)
        ok &= r["PASS"]
        print(json.dumps(r), flush=True)
    print(f"{'ALL PASS' if ok else 'FAIL'}: {len(files)} key files checked", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
