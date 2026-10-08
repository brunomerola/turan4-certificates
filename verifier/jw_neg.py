#!/usr/bin/env python3
"""Independent checker: the reviewer's own negative / positive controls for the verifier jw_verify.py (reviewer code).
Writes corrupted (and one alternative valid) J_4 witnesses to out/neg/, runs jw_verify.py on each, on the producer's
three negative files, and with two deliberately WRONG moment definitions (--mutate), and compares exit codes with the
expected ones.  Exit 0 iff every control behaves as expected."""
from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jw_lib as L  # noqa: E402

PY = sys.executable
S1R = os.path.normpath(os.path.join(HERE, "..", "certificates", "sigma_catalogue", "J4_limit"))
OUT = os.path.join(os.getcwd(), "out", "neg")
os.makedirs(OUT, exist_ok=True)


def load(p):
    return json.load(open(p))


def write(name, W, note):
    W = dict(W)
    W["note"] = note
    W.pop("maker", None)
    p = os.path.join(OUT, name + ".json")
    json.dump(W, open(p, "w"), indent=0)
    return p


def run(path, extra=()):
    js = path.replace(".json", ".rv.json") if path.startswith(OUT) else os.path.join(OUT, os.path.basename(path).replace(".json", ".rv.json"))
    if extra:
        js = js.replace(".rv.json", f".rv_{extra[-1]}.json")
    r = subprocess.run([PY, "-B", os.path.join(HERE, "jw_verify.py"), path, "--json", js, *extra],
                       capture_output=True, text=True)
    rep = json.load(open(js)) if os.path.exists(js) else {}
    return r.returncode, rep.get("fails"), rep.get("delta")


def main():
    W = load(os.path.join(S1R, "witness_Fanoc.json"))
    cls = [int(c) for c in W["classes"]]
    mu = [Fraction(x) for x in W["mu"]]
    delta = Fraction(W["delta"])
    rnd = random.Random(20261007)
    cases = []
    # (a) move 1e-6 from the largest-mu class to the smallest positive one
    i_max = max(range(len(mu)), key=lambda i: mu[i])
    i_min = min((i for i in range(len(mu)) if mu[i] > 0), key=lambda i: mu[i])
    m2 = mu[:]
    eps = Fraction(1, 10 ** 6)
    m2[i_max] -= eps
    m2[i_min] += eps
    cases.append(("rv_moved", write("rv_moved", {**W, "mu": [str(x) for x in m2]},
                                    f"moved 1e-6 from class {cls[i_max]} to {cls[i_min]}"), 1))
    # (b) 1e-6 onto a J_4-free class NOT in supp(phi) (random search, own test)
    phi_can = set(L.canon_many([L.colex_to_lex(c) for c in cls]))
    while True:
        h = rnd.getrandbits(35) & rnd.getrandbits(35)      # sparse-ish random 3-graph
        if not L.contains_J4(h) and L.canon_many([h])[0] not in phi_can:
            break
    m3 = mu[:] + [eps]
    m3[i_max] -= eps
    cases.append(("rv_outside", write("rv_outside", {**W, "classes": cls + [L.lex_to_colex(h)],
                                                     "mu": [str(x) for x in m3]},
                                      f"moved 1e-6 onto a J4-free class outside supp(phi) (lex {h})"), 1))
    # (c) a negative entry: take 1e-9 from a class with mu = 0
    i_zero = next(i for i in range(len(mu)) if mu[i] == 0)
    m4 = mu[:]
    m4[i_zero] -= Fraction(1, 10 ** 9)
    m4[i_max] += Fraction(1, 10 ** 9)
    cases.append(("rv_negative", write("rv_negative", {**W, "mu": [str(x) for x in m4]},
                                       "mu = -1e-9 on a zero class"), 1))
    # (d) swap mu between two support classes with different mu (heaviest and lightest positive)
    i1, i2 = i_max, i_min
    assert mu[i1] != mu[i2]
    m5 = mu[:]
    m5[i1], m5[i2] = m5[i2], m5[i1]
    cases.append(("rv_swapped", write("rv_swapped", {**W, "mu": [str(x) for x in m5]},
                                      "mu values of the heaviest and lightest positive classes swapped"), 1))
    # (e) correct mu, wrong claimed delta
    cases.append(("rv_wrongdelta", write("rv_wrongdelta", {**W, "delta": str(delta + Fraction(1, 10 ** 12))},
                                         "file delta + 1e-12"), 1))
    # (f) POSITIVE: an alternative valid witness (mu + phi)/2 with delta/2 (conditions are linear in mu)
    import pickle
    phi = pickle.load(open(os.path.join(os.getcwd(), "local", "phi_Fanoc.pkl"), "rb"))["phi"]
    can = L.canon_many([L.colex_to_lex(c) for c in cls])
    m6 = [(mu[i] + phi[can[i]]) / 2 for i in range(len(mu))]
    cases.append(("rv_half", write("rv_half", {**W, "mu": [str(x) for x in m6], "delta": str(delta / 2)},
                                   "(mu + phi)/2, delta/2: valid"), 0))
    # (g) mu = phi
    m7 = [phi[can[i]] for i in range(len(mu))]
    cases.append(("rv_phi", write("rv_phi", {**W, "mu": [str(x) for x in m7], "delta": "0"}, "mu = phi"), 2))
    # producer's negative files
    for nm, ex in (("neg_moved", 1), ("neg_outside", 1), ("neg_phi", 2)):
        cases.append(("producer_" + nm, os.path.join(S1R, nm + ".json"), ex))
    allok = True
    lines = []
    for nm, path, ex in cases:
        rc, fails, d = run(path)
        ok = rc == ex
        allok &= ok
        lines.append(f"{'CONTROL OK ' if ok else 'CONTROL BAD'} {nm:22s} exit {rc} (expected {ex}); fails {fails}; delta {d}")
        print(lines[-1], flush=True)
    # mutations of the moment definition, TRUE witness.  edgeflip (moments of one class taken from a wrong graph) must
    # FAIL.  typeblind (all labelled types of a block merged into one matrix, flags identified across types) gives
    # conditions IMPLIED by the labelled ones (ker of a sum of PSD matrices = intersection of kernels), so it must PASS;
    # overlap (U1, U2 may overlap) makes M_OO(phi) positive definite on O (diagonal terms), no conditions: PASS, vacuous.
    for mut, ex in (("edgeflip", 1), ("typeblind", 0), ("overlap", 0)):
        rc, fails, d = run(os.path.join(S1R, "witness_Fanoc.json"), ("--mutate", mut))
        ok = rc == ex
        allok &= ok
        lines.append(f"{'CONTROL OK ' if ok else 'CONTROL BAD'} mutate_{mut:14s} exit {rc} (expected {ex}); fails {fails}")
        print(lines[-1], flush=True)
    print("ALL CONTROLS AS EXPECTED" if allok else "SOME CONTROL NOT AS EXPECTED")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
