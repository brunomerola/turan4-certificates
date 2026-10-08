"""Independent checker core (no producer imports): parse type-pattern certificates, Turan checks, exact cover check.

Conventions (Lemma pattern of the paper):
  type e = r-multiset over parts (count vector, |e| = r); block B = k-multiset (|B| = k)
  mass_e(x) = r! prod_i x_i^{e_i} / e_i!   (fraction of all r-sets of [n] that have type e, n -> oo)
  a(e,B)   = prod_i C(B_i, e_i)            (number of type-e r-subsets of a type-B k-set)
  (LP)     sum_B a(e,B) y_B >= mass_e(x) for e in T   =>   l(k,r,p) <= sum y,  c = C(k,r) l <= C(k,r) sum y.
Recursive patterns: parts listed in `inner` carry a (k,r,q) system of cost c_inner * x_i^r (c-units); Turan
domination is then "P_i >= q for some inner part i, or P dominates a type of T".
"""
import hashlib
import itertools
import math
from fractions import Fraction as Fr


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def parse(path):
    d = {"T": [], "y": [], "inner": None, "header": []}
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            d["header"].append(s)
            continue
        tok = s.split()
        key = tok[0]
        if key in ("k", "r", "p", "m"):
            assert key not in d, "duplicate " + key
            d[key] = int(tok[1])
        elif key == "x":
            d["x"] = [Fr(t) for t in tok[1:]]
        elif key == "T":
            d["T"].append(tuple(int(t) for t in tok[1:]))
        elif key == "y":
            d["y"].append((tuple(int(t) for t in tok[1:-1]), Fr(tok[-1])))
        elif key == "c":
            d["c"] = Fr(tok[1])
        elif key == "inner":
            # format: inner i j ... q Q c745 VALUE
            qi = tok.index("q")
            parts = [int(t) for t in tok[1:qi]]
            q = int(tok[qi + 1])
            val = Fr(tok[qi + 3])
            d["inner"] = {"parts": parts, "q": q, "cname": tok[qi + 2], "c": val}
        else:
            raise ValueError("unknown line: " + s)
    return d


def multisets(m, size):
    """All count vectors of length m summing to size."""
    for c in itertools.combinations_with_replacement(range(m), size):
        v = [0] * m
        for i in c:
            v[i] += 1
        yield tuple(v)


def mass(e, x):
    r = sum(e)
    v = Fr(math.factorial(r))
    for ei, xi in zip(e, x):
        v *= xi ** ei / math.factorial(ei)
    return v


def a_coef(e, B):
    return math.prod(math.comb(b, ei) for b, ei in zip(B, e))


def handled_by_inner(P, inner):
    return inner is not None and any(P[i] >= inner["q"] for i in inner["parts"])


def turan_multiset(d):
    """Every p-multiset is handled by an inner part or dominates a type of T. Returns list of failures."""
    m, p, T, inner = d["m"], d["p"], d["T"], d["inner"]
    bad = []
    for P in multisets(m, p):
        if handled_by_inner(P, inner):
            continue
        if not any(all(ei <= pi for ei, pi in zip(e, P)) for e in T):
            bad.append(P)
    return bad


def turan_brute(d, per=None):
    """Brute force on a blow-up with `per` (default p) labelled points per part: every p-subset of points has
    >= q points in an inner part, or contains an r-subset whose type is in T."""
    m, p, r, T, inner = d["m"], d["p"], d["r"], set(d["T"]), d["inner"]
    per = per or p
    lab = [i for i in range(m) for _ in range(per)]
    npts = len(lab)
    checked = 0
    for S in itertools.combinations(range(npts), p):
        checked += 1
        cnt = [0] * m
        for z in S:
            cnt[lab[z]] += 1
        if handled_by_inner(cnt, inner):
            continue
        ok = False
        for Q in itertools.combinations(S, r):
            v = [0] * m
            for z in Q:
                v[lab[z]] += 1
            if tuple(v) in T:
                ok = True
                break
        if not ok:
            return False, S, checked
    return True, None, checked


def check_cert(path, claimed_decimal=None, verbose=True):
    """Full exact check. Returns dict of results; raises nothing (failures recorded)."""
    d = parse(path)
    k, r, p, m = d["k"], d["r"], d["p"], d["m"]
    x, T, Y = d["x"], d["T"], d["y"]
    res = {"file": path, "k": k, "r": r, "p": p, "m": m, "fail": []}
    # basic sanity
    if len(x) != m:
        res["fail"].append("len(x) != m")
    if any(xi <= 0 for xi in x):
        res["fail"].append("x not > 0")
    if sum(x) != 1:
        res["fail"].append("sum x = %s != 1" % sum(x))
    if len(set(T)) != len(T):
        res["fail"].append("duplicate types")
    for e in T:
        if len(e) != m or sum(e) != r or min(e) < 0:
            res["fail"].append("bad type %s" % (e,))
    blocks = {}
    for B, v in Y:
        if len(B) != m or sum(B) != k or min(B) < 0:
            res["fail"].append("bad block %s" % (B,))
        if v < 0:
            res["fail"].append("negative y %s" % (B,))
        if B in blocks:
            res["fail"].append("duplicate block %s" % (B,))
        blocks[B] = blocks.get(B, Fr(0)) + v
    # total mass over all r-multisets must be 1 (normalisation sanity)
    tot = sum((mass(e, x) for e in multisets(m, r)), Fr(0))
    if tot != 1:
        res["fail"].append("sum of masses over all types = %s" % tot)
    # Turan
    bad = turan_multiset(d)
    res["turan_multiset_failures"] = len(bad)
    if bad:
        res["fail"].append("Turan multiset failures: %s" % bad[:5])
    # cover rows
    slack_min = None
    tight = 0
    for e in T:
        lhs = sum((a_coef(e, B) * v for B, v in blocks.items()), Fr(0))
        s = lhs - mass(e, x)
        if s < 0:
            res["fail"].append("row %s violated by %s" % (e, s))
        if s == 0:
            tight += 1
        slack_min = s if slack_min is None or s < slack_min else slack_min
    res["rows"] = len(T)
    res["tight"] = tight
    res["min_slack"] = slack_min
    sy = sum(blocks.values(), Fr(0))
    c_val = math.comb(k, r) * sy
    if d["inner"] is not None:
        inn = d["inner"]
        c_val += inn["c"] * sum((x[i] ** r for i in inn["parts"]), Fr(0))
    res["sum_y"] = sy
    res["c"] = c_val
    res["c_float"] = float(c_val)
    res["nblocks"] = len(blocks)
    if "c" in d and d["c"] != c_val:
        res["fail"].append("stated c %s != recomputed %s" % (d["c"], c_val))
    res["c_stated_matches"] = ("c" in d and d["c"] == c_val)
    if claimed_decimal is not None:
        cd = Fr(claimed_decimal)
        res["claimed_decimal"] = claimed_decimal
        res["decimal_ok"] = cd >= c_val
        if cd < c_val:
            res["fail"].append("claimed decimal %s is BELOW exact value %s (rounded down)" % (claimed_decimal,
                                                                                            float(c_val)))
    # density of the shadow at x (informative)
    res["density"] = float(sum((mass(e, x) for e in T), Fr(0)))
    return d, res


def decimal_up(fr, digits=9):
    """Smallest decimal with `digits` places that is >= fr."""
    s = 10 ** digits
    q = -((-fr.numerator * s) // fr.denominator)
    return "%d.%0*d" % (q // s, digits, q % s)
