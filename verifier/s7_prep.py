"""Independent prep (own code): rebuild types and flags, decode a certificate, write the evaluator blob and the 6-vertex
representatives.

usage: s7_prep.py PROBLEM CERT_PREFIX OUT_PREFIX
  CERT_PREFIX.cert.json / .cert.npz = the producer's certificate (data only: key list, M, integer factor rows A_k);
  writes OUT_PREFIX.blob (int64 stream for s7_eval), OUT_PREFIX.prep.json (summary + checks), and
  local/reps6_<PROBLEM>.npy (+ _aut.npy) in the working directory if absent.
Everything is computed here from the definitions in s7_lib.py; no producer module is imported.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from itertools import combinations, permutations
from math import comb, factorial

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s7_lib as L  # noqa: E402
OUT = os.getcwd()  # generated files go to the working directory

N = 7
DEN = 5040


def own_types_flags(prob: str, s: int, m: int):
    """Types (canonical admissible s-vertex graphs, increasing, only those with >= 1 flag) and their sorted flag
    lists (canonical under relabellings fixing 0..s-1) for block (s, m)."""
    r = L.PROBLEMS[prob][0]
    nsig = comb(s, r)
    nfree = comb(m, r) - nsig
    codes = np.arange(1 << nsig, dtype=np.int64)
    ok_s = L.admissible_np(prob, codes, s) if nsig else np.ones(1, dtype=bool)
    can_s = L.canon_all_np(codes, s, r) if nsig else np.zeros(1, dtype=np.int64)
    cand = sorted(set(int(x) for x in can_s[ok_s]))
    types, flags = [], []
    fc = np.arange(1 << nfree, dtype=np.int64)
    for sig in cand:
        masks = np.int64(sig) | (fc << nsig)
        ok = L.admissible_np(prob, masks, m)
        if not ok.any():
            continue
        can = L.canon_fix_np(masks[ok], m, r, s)
        fl = np.unique(can)
        types.append(sig)
        flags.append(fl)
    return types, flags


def permcode(c: int, S_order, s: int, r: int) -> int:
    """Labelled code of the ordering pi (tuple of positions into the sorted root set) given the sorted-order code c:
    position-edge e present iff the sorted-position set {pi[i] : i in e} is an edge of c."""
    out = 0
    for e in L.edges(s, r):
        j = L.cidx(tuple(S_order[i] for i in e))
        if (c >> j) & 1:
            out |= 1 << L.cidx(e)
    return out


def build_block(prob: str, s: int, m: int, key_list):
    r = L.PROBLEMS[prob][0]
    nsig = comb(s, r)
    nfree = comb(m, r) - nsig
    types, flags = own_types_flags(prob, s, m)
    tindex = {t: i for i, t in enumerate(types)}
    # compare with the certificate's key list for this block
    mine = [(s, m, i, t, int(fl.size)) for i, (t, fl) in enumerate(zip(types, flags))]
    theirs = [(k["s"], k["m"], k["type_index"], k["sigma"], k["n_flags"]) for k in key_list]
    assert mine == theirs, (s, m, mine, theirs)
    perms = list(permutations(range(s)))
    ncodes = 1 << nsig
    type_of = np.full(ncodes, -1, dtype=np.int64)
    match = [[] for _ in range(ncodes)]
    for c in range(ncodes):
        if nsig and not L.admissible(prob, c, s):
            continue
        can = L.canon_all(c, s, r) if nsig else 0
        if can not in tindex:
            continue
        type_of[c] = tindex[can]
        for pi_i, pi in enumerate(perms):
            if permcode(c, pi, s, r) == can:
                match[c].append(pi_i)
        assert len(match[c]) >= 1
    # flag index tables
    flagidx = []
    for t, (sig, fl) in enumerate(zip(types, flags)):
        fcs = np.arange(1 << nfree, dtype=np.int64)
        masks = np.int64(sig) | (fcs << nsig)
        ok = L.admissible_np(prob, masks, m)
        can = L.canon_fix_np(masks, m, r, s)
        j = np.searchsorted(fl, can)
        jj = np.minimum(j, fl.size - 1)
        idx = np.where(ok & (fl[jj] == can), jj, -1)
        assert np.array_equal(idx >= 0, ok)
        flagidx.append(idx.astype(np.int64))
    # geometry
    rootsets = list(combinations(range(N), s))
    sigpos = []
    freepos = []
    Em = L.edges(m, r)
    free_sub = Em[nsig:]
    nU = comb(N - s, m - s)
    pairs = None
    for S in rootsets:
        sigpos.append([L.cidx(tuple(S[x] for x in e)) for e in L.edges(s, r)])
        O = [v for v in range(N) if v not in S]
        Us = list(combinations(O, m - s))
        assert len(Us) == nU
        if pairs is None:
            pairs = [(i, j) for i in range(nU) for j in range(i + 1, nU) if not set(Us[i]) & set(Us[j])]
        for pi in perms:
            theta = tuple(S[pi[i]] for i in range(s))
            for U in Us:
                V = theta + tuple(U)
                freepos.append([L.cidx(tuple(V[x] for x in e)) for e in free_sub])
    denom = factorial(N) // factorial(N - s) * comb(N - s, m - s) * comb(N - m, m - s)
    assert denom == len(rootsets) * len(perms) * 2 * len(pairs)
    assert DEN % denom == 0
    return dict(s=s, m=m, r=r, nsig=nsig, nfree=nfree, types=types, flags=flags, type_of=type_of, match=match,
                flagidx=flagidx, sigpos=np.array(sigpos, dtype=np.int64).reshape(len(rootsets), nsig),
                freepos=np.array(freepos, dtype=np.int64).reshape(-1), nU=nU, pairs=pairs, nrs=len(rootsets),
                nperm=len(perms), denom=denom, wb=DEN // denom)


def reps6(prob: str):
    """Own canonical representatives (min over all 720 relabellings) of the admissible 6-vertex graphs, and |Aut|."""
    r = L.PROBLEMS[prob][0]
    path = os.path.join(OUT, "local", f"reps6_{prob}.npy")
    if os.path.exists(path):
        return np.load(path), np.load(path.replace(".npy", "_aut.npy"))
    nb = comb(6, r)
    allm = np.arange(1 << nb, dtype=np.int64)
    adm = allm[L.admissible_np(prob, allm, 6)]
    best = adm.copy()
    for g in permutations(range(6)):
        best = np.minimum(best, L.apply_perm_np(adm, 6, r, g))
    reps = np.unique(best)
    aut = np.zeros(reps.size, dtype=np.int64)
    for g in permutations(range(6)):
        aut += (L.apply_perm_np(reps, 6, r, g) == reps)
    np.save(path, reps)
    np.save(path.replace(".npy", "_aut.npy"), aut)
    # orbit-stabiliser check: sum 720/|Aut| = number of admissible labelled graphs
    assert int(sum(720 // int(a) for a in aut)) == adm.size and all(720 % int(a) == 0 for a in aut)
    return reps, aut


def objtab(obj: str):
    if obj == "cosig3":   # for each 4-set Q of [7]: mask of its four triples (r = 3)
        return [sum(1 << L.cidx(t) for t in combinations(Q, 3)) for Q in combinations(range(N), 4)]
    # sigma4: for each 3-set T of [7]: mask of the four 4-sets T + x (r = 4)
    return [sum(1 << L.cidx(tuple(sorted(T + (x,)))) for x in range(N) if x not in T) for T in combinations(range(N), 3)]


def main():
    prob, cert, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    t0 = time.time()
    r, mode, keys, obj = L.PROBLEMS[prob]
    J = json.load(open(cert + ".cert.json"))
    assert J["problem"] == prob and J["objective"] == obj, (J["problem"], J["objective"])
    assert J["tau"] == [0, 0]
    with open(cert + ".cert.npz", "rb") as f:
        npz_sha = hashlib.sha256(f.read()).hexdigest()
    assert npz_sha == J["cert_npz_sha256"], "npz hash differs from the json"
    Z = np.load(cert + ".cert.npz", allow_pickle=False)
    M = int(J["M"])
    assert M == 1 << int(J["M_bits"])
    M2 = M * M
    blocks = [tuple(b) for b in J["blocks"]]
    keylist = J["keys"]
    assert len(Z.files) == len(keylist) and set(Z.files) == {f"A{k}" for k in range(len(keylist))}
    built = []
    kpos = 0
    for (s, m) in blocks:
        kl = [k for k in keylist if (k["s"], k["m"]) == (s, m)]
        assert keylist[kpos:kpos + len(kl)] == kl, "keys not grouped by block in order"
        B = build_block(prob, s, m, kl)
        B["keys"] = list(range(kpos, kpos + len(kl)))
        kpos += len(kl)
        built.append(B)
        print(f"# block ({s},{m}): {len(B['types'])} types, flags {[int(f.size) for f in B['flags']]}, "
              f"denom {B['denom']}, wb {B['wb']}", flush=True)
    assert kpos == len(keylist)
    # Gram matrices, exactly
    Gs = []
    gstat = []
    for k in range(len(keylist)):
        A = Z[f"A{k}"]
        assert A.dtype == np.int64 and A.ndim == 2 and A.shape[1] == keylist[k]["n_flags"]
        amax = int(np.abs(A).max()) if A.size else 0
        assert amax * amax * max(1, A.shape[0]) < (1 << 62), "possible int64 overflow in A^T A"
        G = A.T @ A
        assert np.array_equal(G, G.T)
        if A.size:  # spot check against Python ints
            rng = np.random.default_rng(k)
            for _ in range(20):
                a, b = (int(x) for x in rng.integers(0, A.shape[1], 2))
                assert int(G[a, b]) == sum(int(A[i, a]) * int(A[i, b]) for i in range(A.shape[0]))
        Gs.append(G)
        gstat.append({"key": k, "rows": int(A.shape[0]), "n": int(A.shape[1]), "max_abs_A": amax,
                      "max_abs_G": int(np.abs(G).max()) if G.size else 0,
                      "max_diag_G": int(np.diag(G).max()) if G.size else 0})
    # a-priori |Z| bound with these weights: per block, sum of weights = DEN, |G| <= max diag (PSD Gram)
    zb = sum(DEN * max(g["max_diag_G"] for g in gstat if g["key"] in B["keys"]) for B in built)
    # blob
    out = []

    def put(*xs):
        out.extend(int(x) for x in xs)

    def arr(a):
        a = [int(x) for x in np.asarray(a).reshape(-1)]
        out.append(len(a))
        out.extend(a)

    cops = L.problem_copies(prob, N)
    put(0x53494737, r, 0 if mode == "sub" else 1, 0 if obj == "cosig3" else 1, L.DEN7[obj], M2)
    arr(cops)
    arr(objtab(obj))
    put(comb(6, r), comb(6, r - 1), len(built))
    goff = 0
    gflat = []
    for B in built:
        put(B["s"], B["m"], B["wb"], B["nsig"], B["nfree"], B["nU"], len(B["pairs"]), B["nrs"], B["nperm"],
            len(B["types"]))
        arr(np.array(B["pairs"], dtype=np.int64).reshape(-1))
        arr(B["sigpos"])
        arr(B["type_of"])
        off = [0]
        flat = []
        for lst in B["match"]:
            flat.extend(lst)
            off.append(len(flat))
        arr(off)
        arr(flat)
        arr(B["freepos"])
        for t, k in enumerate(B["keys"]):
            put(int(B["flags"][t].size))
            arr(B["flagidx"][t])
            put(goff)
            g = Gs[k]
            gflat.append(g.reshape(-1))
            goff += g.size
    allg = np.concatenate(gflat) if gflat else np.zeros(0, np.int64)
    put(allg.size)
    blob = np.concatenate([np.array(out, dtype=np.int64), allg.astype(np.int64)])
    blob.tofile(outp + ".blob")
    reps, aut = reps6(prob)
    reps.astype(np.int64).tofile(outp + ".reps.bin")
    summ = {"problem": prob, "cert_prefix": cert, "cert_npz_sha256": npz_sha, "M": M, "M2": M2,
            "bound_claimed": J["bound"], "blocks": blocks,
            "blocks_built": [{"s": B["s"], "m": B["m"], "denom": B["denom"], "wb": B["wb"], "types": B["types"],
                              "n_flags": [int(f.size) for f in B["flags"]], "nU": B["nU"], "npairs": len(B["pairs"]),
                              "aut_sizes": [len(B["match"][int(t)]) if B["nsig"] else 1 for t in B["types"]]}
                             for B in built],
            "keys_match_cert": True, "gram": gstat, "z_bound_apriori_mine": str(zb),
            "z_bound_apriori_cert": J["z_bound_apriori"], "n_copies7": len(cops),
            "reps6": int(reps.size), "labelled6": int(sum(720 // int(a) for a in aut)),
            "blob_sha256": hashlib.sha256(blob.tobytes()).hexdigest(), "time_s": time.time() - t0}
    json.dump(summ, open(outp + ".prep.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("gram", "blocks_built")}), flush=True)


if __name__ == "__main__":
    main()
