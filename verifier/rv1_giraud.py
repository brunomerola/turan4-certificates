"""Independent checker, CLAIM 2: exact 7-vertex law of Giraud's construction (own derivation) and exact d - c_Q on it.

Law (limit n -> infinity of a uniform 7-subset of the Giraud system with equal parts and an i.i.d. uniform 0/1 matrix
M): each of the 7 vertices is in part A or B independently with probability 1/2; the |A||B| entries M_ab (a in A,
b in B) are i.i.d. uniform bits.  Edges: 4-sets inside A, 4-sets inside B, and cross sets {a,a',b,b'} with even minor
M_ab + M_ab' + M_a'b + M_a'b'.  A configuration (A, M) has probability 2^-(7 + |A||B|).
Evaluation: own moment counts (rv1_lib), certificate Q = A^T A / M^2 (lpcg_p5_conv_full), exact integer numerators
N(H) = 144 e(H) M^2 - sum_b (5040/conf_b) Z_b(H), value d - c = N / (5040 M^2).
Also: positive control of the own PSD test (moment matrices of the law must be PSD: a genuine limit object), class
count, export of colex masks (batches of <= 64) for the C evaluator rv_eval (rv_eval.c, list mode).
Usage: python rv1_giraud.py CERT_PREFIX OUTDIR
"""
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

import rv1_lib as L
from rv1_dual import BLOCKS, build_Y, check_Y, eval_N, load_cert


def giraud_law():
    law = {}
    nconf = 0
    for Abits in range(1 << 7):
        Aset = [v for v in range(7) if (Abits >> v) & 1]
        Bset = [v for v in range(7) if not (Abits >> v) & 1]
        k, l = len(Aset), len(Bset)
        pidx = {(a, b): i for i, (a, b) in enumerate((a, b) for a in Aset for b in Bset)}
        w = Fraction(1, 2 ** (7 + k * l))
        # classify the 35 4-sets once per part assignment
        internal, cross = [], []
        for j, e in enumerate(L.L7):
            ina = [v for v in e if (Abits >> v) & 1]
            inb = [v for v in e if not (Abits >> v) & 1]
            if len(ina) in (0, 4):
                internal.append(j)
            elif len(ina) == 2:
                a, a2 = ina
                b, b2 = inb
                cross.append((j, pidx[(a, b)], pidx[(a, b2)], pidx[(a2, b)], pidx[(a2, b2)]))
        base = sum(1 << j for j in internal)
        for Mb in range(1 << (k * l)):
            nconf += 1
            h = base
            for (j, p, q, r, t) in cross:
                if not (((Mb >> p) ^ (Mb >> q) ^ (Mb >> r) ^ (Mb >> t)) & 1):
                    h |= 1 << j
            law[h] = law.get(h, 0) + w
    return law, nconf


def canon_numpy(hs):
    pm = np.array(L.perm_edge_maps(), dtype=np.int64)          # (5040, 35)
    pw = (np.int64(1) << pm)                                     # image bit of each edge under each permutation
    out = []
    for h in hs:
        bits = [i for i in range(35) if (h >> i) & 1]
        img = pw[:, bits].sum(axis=1) if bits else np.zeros(5040, dtype=np.int64)
        out.append(int(img.min()))
    return out


def main():
    t0 = time.time()
    prefix, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(outdir, "giraud.log"), "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    law, nconf = giraud_law()
    log("configurations:", nconf, " (sum_k C(7,k) 2^(k(7-k)) =",
        sum(__import__("math").comb(7, k) * 2 ** (k * (7 - k)) for k in range(8)), ")")
    log("total probability == 1:", sum(law.values()) == 1, "; distinct labelled graphs:", len(law))
    log("all odd:", all(L.is_odd(h) for h in law), "; all admissible:", all(L.admissible5(h) for h in law))
    Ed = sum(w * Fraction(L.popcount(h), 35) for h, w in law.items())
    log("E[d] =", Ed, "== 5/16:", Ed == Fraction(5, 16))
    prof = {}
    for h, w in law.items():
        for c in L.five_profile(h):
            prof[c] = prof.get(c, 0) + w / 21
    log("5-set profile (1,3,5 edges):", {k: str(v) for k, v in sorted(prof.items())})
    hs = sorted(law)
    cls = canon_numpy(hs)
    classes = {}
    for h, c in zip(hs, cls):
        classes.setdefault(c, []).append(h)
    log("isomorphism classes:", len(classes), " sizes:", sorted(len(v) for v in classes.values()))

    cert, keys = load_cert(prefix)
    Mfac = int(cert["M"])
    b = Fraction(cert["bound"])
    DEN = 5040 * Mfac * Mfac
    Wint = {h: int(w * 2 ** 19) for h, w in law.items()}
    assert all(Fraction(Wint[h], 2 ** 19) == law[h] for h in law)
    vals = {}
    Y = {bb: {} for bb in BLOCKS}
    t1 = time.time()
    for i, h in enumerate(hs):
        mcs = {}
        for (s, m) in BLOCKS:
            mc = L.moment_counts(h, s, m)
            mcs[(s, m)] = mc
            for tau, d in mc.items():
                Dd = Y[(s, m)].setdefault(tau, {})
                for key, c in d.items():
                    Dd[key] = Dd.get(key, 0) + Wint[h] * c
        vals[h] = eval_N(h, keys, Mfac, mc_cache=mcs)
        if i % 400 == 0:
            log(f"  evaluated {i}/{len(hs)} ({time.time() - t1:.0f}s)")
    v = {h: Fraction(n, DEN) for h, n in vals.items()}
    # isomorphism invariance of the own evaluator
    inv = all(len({vals[h] for h in g}) == 1 for g in classes.values())
    log("value constant on every isomorphism class (evaluator self-consistency):", inv)
    mn, mx = min(v.values()), max(v.values())
    Edc = sum(law[h] * v[h] for h in hs)
    Ec = Ed - Edc
    log("b =", b)
    log("min (d - c) - b =", mn - b, "=", float(mn - b), " claimed 1091855/15393162788864:",
        mn - b == Fraction(1091855, 15393162788864))
    log("max (d - c) - b =", mx - b, "=", float(mx - b), " claimed 52186103/115448720916480:",
        mx - b == Fraction(52186103, 115448720916480))
    log("E[d - c] =", Edc, "=", float(Edc), " claimed 173614396936021/2^49:",
        Edc == Fraction(173614396936021, 2 ** 49))
    log("E[c] =", Ec, "=", float(Ec), " claimed 2307463508139/2^49:", Ec == Fraction(2307463508139, 2 ** 49))
    log("E[d - c] - b =", float(Edc - b), "; 5/16 - b =", float(Fraction(5, 16) - b),
        "; E[c] + (E[d-c] - b) == 5/16 - b:", Ec + (Edc - b) == Fraction(5, 16) - b)
    log("all support graphs >= b:", mn >= b)
    log("per class: (canonical own mask, labelled count, P(class), e, (d-c)-b)")
    for c, g in sorted(classes.items(), key=lambda kv: v[kv[1][0]]):
        log("  ", c, len(g), str(sum(law[h] for h in g)), L.popcount(g[0]), float(v[g[0]] - b))

    log("# positive control: moment matrices of the Giraud 7-vertex law (exact PSD test, all labelled types)")
    universes = {}
    for (s, m) in BLOCKS:
        for tau in L.all_labelled_types(s):
            universes[(s, m, tau)] = L.flag_universe(s, m, tau)
    okg, _ = check_Y(Y, universes, log)
    log("Giraud law: all moment matrices PSD:", okg)

    # exports for rv_eval (colex masks, uint64, batches of <= 64)
    cm = [L.own_to_colex(h) for h in hs]
    nb = 0
    for i in range(0, len(cm), 64):
        np.array(cm[i:i + 64], dtype=np.uint64).tofile(os.path.join(outdir, f"giraud_batch{nb:02d}.u64"))
        nb += 1
    json.dump({"colex": cm, "N": [vals[h] for h in hs], "w_times_2^19": [Wint[h] for h in hs],
               "DEN": DEN, "b": str(b), "E_d_minus_c": str(Edc), "E_c": str(Ec), "min_minus_b": str(mn - b),
               "max_minus_b": str(mx - b), "classes": len(classes), "batches": nb},
              open(os.path.join(outdir, "giraud_values.json"), "w"))
    log(f"exported {nb} batches; done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
