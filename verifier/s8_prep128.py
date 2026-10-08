"""Independent checker, prep for certificates whose Gram entries exceed int64 (J_4 re-rounded at M = 2^30).

usage: s8_prep128.py PROBLEM CERT_PREFIX OUT_PREFIX
Same decoding as s8_prep.py (its build_block / reps6 / objtab are reused unchanged), plus:
  * the npz may carry a scalar 'M_bits' (must equal the json's M_bits);
  * the json now lists the flags per key: they must equal the reviewer's own flag lists exactly;
  * G_k = A_k^T A_k is computed EXACTLY without int64 overflow: A = 2^16 A1 + A0 with 0 <= A0 < 2^16 (floor split),
    G = 2^32 (A1^T A1) + 2^16 (A1^T A0 + A0^T A1) + A0^T A0, each product in int64 under an asserted a-priori bound;
    the three int64 limb matrices are written to the blob (magic 'SIG8'), and s8_eval128 recombines them in __int128.
    50 random entries per key are recomputed from the definition in Python integers.
  * a-priori bound for the evaluator: |Z| <= 5040 * sum over blocks of max diag G (PSD Gram), and
    |S| <= 210 * 5040 * M^2 + den7 * |Z| must be < 2^126.
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


def limbs(A: np.ndarray):
    A0 = np.mod(A, 1 << 16)
    A1 = (A - A0) >> 16
    assert np.array_equal(A1 * (1 << 16) + A0, A)
    rows = max(1, A.shape[0])
    b1 = int(np.abs(A1).max()) if A1.size else 0
    assert b1 * b1 * rows < (1 << 62) and 2 * b1 * (1 << 16) * rows < (1 << 62) and (1 << 32) * rows < (1 << 62)
    P = A1.T @ A1
    Q = A1.T @ A0 + A0.T @ A1
    R = A0.T @ A0
    return P, Q, R


def main():
    prob, cert, outp = sys.argv[1], sys.argv[2], sys.argv[3]
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
    if "M_bits" in Z.files:
        assert int(Z["M_bits"]) == int(J["M_bits"])
    M2 = M * M
    assert M2 < (1 << 63)
    keylist = J["keys"]
    assert set(Z.files) - {"M_bits"} == {f"A{k}" for k in range(len(keylist))}
    blocks = [tuple(b) for b in J["blocks"]]
    built = []
    kpos = 0
    flags_equal = True
    for (s, m) in blocks:
        kl = [k for k in keylist if (k["s"], k["m"]) == (s, m)]
        assert keylist[kpos:kpos + len(kl)] == kl
        B = build_block(prob, s, m, kl)       # asserts (s, m, type index, sigma, n_flags) = the json's
        for t, k in enumerate(kl):
            if "flags" in k:
                flags_equal &= [int(x) for x in k["flags"]] == [int(x) for x in B["flags"][t]]
            if "aut_size" in k and B["nsig"]:
                assert int(k["aut_size"]) == len(B["match"][int(B["types"][t])])
        B["keys"] = list(range(kpos, kpos + len(kl)))
        kpos += len(kl)
        built.append(B)
    assert kpos == len(keylist)
    assert flags_equal, "json flag lists differ from the reviewer's"
    Ps, Qs, Rs, gstat = [], [], [], []
    rng = np.random.default_rng(30)
    for k in range(len(keylist)):
        A = Z[f"A{k}"]
        assert A.dtype == np.int64 and A.ndim == 2 and A.shape[1] == keylist[k]["n_flags"]
        P, Q, R = limbs(A)
        n = A.shape[1]
        rowsA = [[int(x) for x in row] for row in A.tolist()]
        for _ in range(50 if A.size else 0):
            a, b = (int(x) for x in rng.integers(0, n, 2))
            g = (int(P[a, b]) << 32) + (int(Q[a, b]) << 16) + int(R[a, b])
            assert g == sum(row[a] * row[b] for row in rowsA)
        diag = [(int(P[i, i]) << 32) + (int(Q[i, i]) << 16) + int(R[i, i]) for i in range(n)]
        if A.size:
            assert all(d == sum(row[i] * row[i] for row in rowsA) for i, d in enumerate(diag))
        Ps.append(P.reshape(-1)); Qs.append(Q.reshape(-1)); Rs.append(R.reshape(-1))
        gstat.append({"key": k, "rows": int(A.shape[0]), "n": n,
                      "max_abs_A": int(np.abs(A).max()) if A.size else 0, "max_diag_G": max(diag) if diag else 0})
    zb = sum(DEN * max(g["max_diag_G"] for g in gstat if g["key"] in B["keys"]) for B in built)
    sb = 210 * DEN * M2 + L.DEN7[obj] * zb if obj == "cosig3" else 1260 * DEN * M2 + L.DEN7[obj] * zb
    assert sb < (1 << 126), "possible int128 overflow"
    out = []

    def put(*xs):
        out.extend(int(x) for x in xs)

    def arr(a):
        a = [int(x) for x in np.asarray(a).reshape(-1)]
        out.append(len(a))
        out.extend(a)

    cops = L.problem_copies(prob, N)
    put(0x53494738, r, 0 if mode == "sub" else 1, 0 if obj == "cosig3" else 1, L.DEN7[obj], M2)
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
            goff += Ps[k].size
    put(goff)
    head = np.array(out, dtype=np.int64)
    blob = np.concatenate([head] + [np.concatenate(x) if x else np.zeros(0, np.int64) for x in (Ps, Qs, Rs)])
    blob.astype(np.int64).tofile(outp + ".blob")
    reps, aut = reps6(prob)
    reps.astype(np.int64).tofile(outp + ".reps.bin")
    summ = {"problem": prob, "cert_prefix": cert, "cert_npz_sha256": npz_sha, "M": M, "M2": M2,
            "bound_claimed": J["bound"], "blocks": blocks, "keys_match_cert": True, "flags_equal_json": flags_equal,
            "max_abs_A": max(g["max_abs_A"] for g in gstat), "max_diag_G": str(max(g["max_diag_G"] for g in gstat)),
            "max_diag_G_bits": max(g["max_diag_G"] for g in gstat).bit_length(),
            "z_bound_apriori_mine": str(zb), "z_bound_bits": zb.bit_length(), "S_bound_bits": sb.bit_length(),
            "reps6": int(reps.size), "blob_sha256": hashlib.sha256(blob.tobytes()).hexdigest(),
            "time_s": time.time() - t0}
    json.dump(summ, open(outp + ".prep.json", "w"), indent=1)
    print(json.dumps(summ), flush=True)


if __name__ == "__main__":
    main()
