"""Independent streaming prep for the run directory (own code; low memory: Gram blocks are written key by key).

usage: s8_prepvm.py PROBLEM CERT_PREFIX NAME VMRUN_DIR
  decodes CERT_PREFIX.cert.{json,npz} exactly as s8_prep.py / s8_prep128.py (reviewer's own types and flags, compared
  with the json key list; npz hash = json), computes G_k = A_k^T A_k exactly and writes
    VMRUN_DIR/blobs/NAME.blob   magic 'SIG7' (int64 G) if max|A|^2 * rows < 2^62 for every key, else
                                magic 'SIG8' (three int64 limb matrices, G = 2^32 P + 2^16 Q + R; evaluator s8_eval128)
    VMRUN_DIR/reps/PROBLEM.reps.bin   the reviewer's 6-vertex representatives (int64)
    VMRUN_DIR/prep/NAME.prep.json     summary (bound claimed, Gram stats, a-priori bounds, blob sha256, evaluator)
  Every Gram entry check: 20 random entries per key (and all diagonals in limb mode) recomputed in Python integers.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s8_lib as L  # noqa: E402
from s8_prep import N, DEN, build_block, objtab, reps6  # noqa: E402
from s8_prep128 import limbs  # noqa: E402


def main():
    prob, cert, name, vm = sys.argv[1:5]
    t0 = time.time()
    r, mode, keys, obj = L.PROBLEMS[prob]
    J = json.load(open(cert + ".cert.json"))
    assert J["problem"] == prob and J["objective"] == obj and J["tau"] == [0, 0]
    with open(cert + ".cert.npz", "rb") as f:
        npz_sha = hashlib.sha256(f.read()).hexdigest()
    assert npz_sha == J["cert_npz_sha256"]
    Z = np.load(cert + ".cert.npz", allow_pickle=False)
    M = int(J["M"])
    assert M == 1 << int(J["M_bits"])
    M2 = M * M
    keylist = J["keys"]
    assert set(Z.files) - {"M_bits"} == {f"A{k}" for k in range(len(keylist))}
    blocks = [tuple(b) for b in J["blocks"]]
    built, kpos = [], 0
    for (s, m) in blocks:
        kl = [k for k in keylist if (k["s"], k["m"]) == (s, m)]
        assert keylist[kpos:kpos + len(kl)] == kl
        B = build_block(prob, s, m, kl)
        B["keys"] = list(range(kpos, kpos + len(kl)))
        kpos += len(kl)
        built.append(B)
        print(f"# block ({s},{m}): {len(B['types'])} types, flags {[int(f.size) for f in B['flags']]}", flush=True)
    assert kpos == len(keylist)
    # int64 possible?
    use64 = True
    for k in range(len(keylist)):
        A = Z[f"A{k}"]
        assert A.dtype == np.int64 and A.ndim == 2 and A.shape[1] == keylist[k]["n_flags"]
        amax = int(np.abs(A).max()) if A.size else 0
        if amax * amax * max(1, A.shape[0]) >= (1 << 62):
            use64 = False
    os.makedirs(os.path.join(vm, "blobs"), exist_ok=True)
    os.makedirs(os.path.join(vm, "reps"), exist_ok=True)
    os.makedirs(os.path.join(vm, "prep"), exist_ok=True)
    out = []

    def put(*xs):
        out.extend(int(x) for x in xs)

    def arr(a):
        a = [int(x) for x in np.asarray(a).reshape(-1)]
        out.append(len(a))
        out.extend(a)

    cops = L.problem_copies(prob, N)
    put(0x53494737 if use64 else 0x53494738, r, 0 if mode == "sub" else 1, 0 if obj == "cosig3" else 1,
        L.DEN7[obj], M2)
    arr(cops)
    arr(objtab(obj))
    put(comb(6, r), comb(6, r - 1), len(built))
    goff = 0
    for B in built:
        put(B["s"], B["m"], B["wb"], B["nsig"], B["nfree"], B["nU"], len(B["pairs"]), B["nrs"], B["nperm"],
            len(B["types"]))
        arr(np.array(B["pairs"], dtype=np.int64).reshape(-1))
        arr(B["sigpos"])
        arr(B["type_of"])
        off, flat = [0], []
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
            goff += int(B["flags"][t].size) ** 2
    put(goff)
    bpath = os.path.join(vm, "blobs", name + ".blob")
    h = hashlib.sha256()
    rng = np.random.default_rng(8)
    gstat = []
    with open(bpath, "wb") as f:
        hb = np.array(out, dtype=np.int64).tobytes()
        f.write(hb)
        h.update(hb)
        passes = [None] if use64 else ["P", "Q", "R"]
        for part in passes:
            for k in range(len(keylist)):
                A = Z[f"A{k}"]
                n = A.shape[1]
                if use64:
                    G = A.T @ A
                    assert np.array_equal(G, G.T)
                    piece = G
                else:
                    P, Q, R = limbs(A)
                    piece = {"P": P, "Q": Q, "R": R}[part]
                if part in (None, "P") and A.size:
                    rows = [[int(x) for x in row] for row in A.tolist()]
                    for _ in range(20):
                        a, b = (int(x) for x in rng.integers(0, n, 2))
                        g = int(G[a, b]) if use64 else (int(P[a, b]) << 32) + (int(Q[a, b]) << 16) + int(R[a, b])
                        assert g == sum(row[a] * row[b] for row in rows)
                    diag = [int(G[i, i]) for i in range(n)] if use64 else \
                        [(int(P[i, i]) << 32) + (int(Q[i, i]) << 16) + int(R[i, i]) for i in range(n)]
                    if not use64:
                        assert all(d == sum(row[i] * row[i] for row in rows) for i, d in enumerate(diag))
                    gstat.append({"key": k, "rows": int(A.shape[0]), "n": n, "max_abs_A": int(np.abs(A).max()),
                                  "max_diag_G": max(diag)})
                elif part in (None, "P"):
                    gstat.append({"key": k, "rows": 0, "n": n, "max_abs_A": 0, "max_diag_G": 0})
                pb = np.ascontiguousarray(piece, dtype=np.int64).tobytes()
                f.write(pb)
                h.update(pb)
                del piece
    zb = sum(DEN * max(g["max_diag_G"] for g in gstat if g["key"] in B["keys"]) for B in built)
    fmax = 210 if obj == "cosig3" else 840
    sb = fmax * DEN * M2 + L.DEN7[obj] * zb
    assert sb < (1 << 126)
    if use64:
        assert max(g["max_diag_G"] for g in gstat) < (1 << 62)
    reps, aut = reps6(prob)
    rpath = os.path.join(vm, "reps", prob + ".reps.bin")
    reps.astype(np.int64).tofile(rpath)
    summ = {"name": name, "problem": prob, "cert_prefix": os.path.relpath(cert, HERE), "cert_npz_sha256": npz_sha,
            "M": M, "bound_claimed": J["bound"], "evaluator": "s8_eval" if use64 else "s8_eval128",
            "keys": len(keylist), "keys_match_cert": True, "max_abs_A": max(g["max_abs_A"] for g in gstat),
            "max_diag_G_bits": max(g["max_diag_G"] for g in gstat).bit_length(),
            "z_bound_apriori_mine": str(zb), "z_bound_apriori_cert": J["z_bound_apriori"],
            "z_bound_bits": zb.bit_length(), "S_bound_bits": sb.bit_length(),
            "reps6": int(reps.size), "labelled6": int(sum(720 // int(a) for a in aut)),
            "raw_reps6_cert": J.get("raw_reps6"), "raw_extensions_cert": J.get("raw_extensions_scanned"),
            "classes_cert": J["classes_scanned"],
            "blob_sha256": h.hexdigest(), "blob_bytes": os.path.getsize(bpath), "time_s": time.time() - t0}
    json.dump(summ, open(os.path.join(vm, "prep", name + ".prep.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)


if __name__ == "__main__":
    main()
