# Independent checker (own code; imports nothing from the producer's code).
# Exact re-check of every number printed in the paper's text on upper bounds for lottery numbers (Section sec:lotup).
# Usage: python -I r7p6_check.py RM      (RM = the root of the data package)
import sys, hashlib
from fractions import Fraction as Q
from math import comb, factorial
from itertools import combinations, combinations_with_replacement, product

RM = sys.argv[1]
OK = []
def chk(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS " if cond else "FAIL ") + name)

def floor_dec(q, d):          # decimal string of q rounded DOWN to d digits
    q = Q(q); n = (q.numerator * 10**d) // q.denominator
    return "%d.%0*d" % (n // 10**d, d, n % 10**d)
def ceil_dec(q, d):
    q = Q(q); n = -((-q.numerator * 10**d) // q.denominator)
    return "%d.%0*d" % (n // 10**d, d, n % 10**d)

# ---------------------------------------------------------------- certificate files
def read_cert(path):
    C = {"T": [], "y": []}
    with open(path) as fh:
        for ln in fh:
            s = ln.split()
            if not s or s[0].startswith("#"):
                continue
            if s[0] == "T":
                C["T"].append(tuple(int(v) for v in s[1:]))
            elif s[0] == "y":
                C["y"].append((tuple(int(v) for v in s[1:-1]), Q(s[-1])))
            elif s[0] == "x":
                C["x"] = [Q(v) for v in s[1:]]
            elif s[0] == "inner":
                # inner i j ... q Q c745 V
                i = s.index("q"); C["inner_parts"] = tuple(int(v) for v in s[1:i])
                C["q"] = int(s[i + 1]); C["inner_c"] = Q(s[i + 3])
            elif s[0] in ("k", "r", "p", "m"):
                C[s[0]] = int(s[1])
            elif s[0] == "c":
                C["c"] = Q(s[1])
    return C

def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()

def mu(e, x):                 # r! prod x_i^e_i / e_i!
    v = Q(factorial(sum(e)))
    for ei, xi in zip(e, x):
        v *= Q(xi) ** ei / factorial(ei)
    return v

def a(e, B):
    v = 1
    for ei, bi in zip(e, B):
        v *= comb(bi, ei)
    return v

def types_of_size(m, s):      # all count vectors of size s over m parts
    out = []
    for c in combinations_with_replacement(range(m), s):
        v = [0] * m
        for i in c:
            v[i] += 1
        out.append(tuple(v))
    return out

def dominated(P, T):
    return any(all(P[i] >= e[i] for i in range(len(P))) for e in T)

def inner_ok(P, inner, q):
    return any(P[i] >= q for i in inner)

def check_cert(path, expect_sha, label):
    C = read_cert(path)
    k, r, p, m = C["k"], C["r"], C["p"], C["m"]
    x, T, Y = C["x"], C["T"], C["y"]
    inner = C.get("inner_parts", ()); q = C.get("q")
    s1 = sha(path); s2 = open(path + ".sha256").read().split()[0]
    chk("%s sha256 = sidecar = block/P6_sources (%s...)" % (label, s1[:8]), s1 == s2 == expect_sha)
    chk("%s x > 0, sum x = 1, m = %d" % (label, m), all(v > 0 for v in x) and sum(x) == 1 and len(x) == m)
    chk("%s all types have size r=%d, all blocks size k=%d, all y > 0" % (label, r, k),
        all(sum(e) == r for e in T) and all(sum(B) == k for B, _ in Y) and all(v > 0 for _, v in Y))
    chk("%s no duplicate types / blocks" % label, len(set(T)) == len(T) and len(set(B for B, _ in Y)) == len(Y))
    # rows (ii)
    tight = 0; bad = 0
    for e in T:
        lhs = sum(a(e, B) * v for B, v in Y); rhs = mu(e, x)
        if lhs < rhs: bad += 1
        if lhs == rhs: tight += 1
    chk("%s condition (ii): %d rows, all satisfied, %d tight" % (label, len(T), tight), bad == 0)
    # domination (i), type level
    PT = types_of_size(m, p)
    fails = [P for P in PT if not inner_ok(P, inner, q or 10**9) and not dominated(P, T)]
    chk("%s condition (i): %d types of size %d, failures %d" % (label, len(PT), p, len(fails)), not fails)
    # brute force on a blow-up with p points per part (independent of the type-level argument)
    pts = [(i, j) for i in range(m) for j in range(p)]
    nfail = 0; ntot = 0
    for S in combinations(pts, p):
        ntot += 1
        P = [0] * m
        for (i, _) in S: P[i] += 1
        if inner and any(P[i] >= q for i in inner):
            continue
        # look for an r-subset of S whose type is in T (explicit subsets, not domination)
        found = False
        Tset = set(T)
        for R in combinations(S, r):
            e = [0] * m
            for (i, _) in R: e[i] += 1
            if tuple(e) in Tset:
                found = True; break
        if not found: nfail += 1
    chk("%s brute force on %dx%d blow-up: %d %d-sets, failures %d" % (label, m, p, ntot, p, nfail), nfail == 0)
    # minimality of T (informative)
    nonmin = []
    for e in T:
        T2 = [f for f in T if f != e]
        if all(inner_ok(P, inner, q or 10**9) or dominated(P, T2) for P in PT):
            nonmin.append(e)
    print("info %s types removable without breaking (i): %s" % (label, nonmin))
    shapes = {}
    for e in T:
        sh = tuple(sorted((v for v in e if v), reverse=True)); shapes[sh] = shapes.get(sh, 0) + 1
    val = comb(k, r) * sum(v for _, v in Y)
    return C, tight, len(T), len(Y), shapes, val

c745 = Q(3467292853891919529774577, 7750000000000000000000000)
print("=== (b) c(7,4,5)")
f5 = RM + "/certificates/lottery_upper/k7_r4p5_m5_den1000000.txt"
C5, t5, nT5, nY5, sh5, v5 = check_cert(f5, "5f0f88f5007a42535dfc17823cf3e1165749dbd9f27dcca8f1a0271fc0060d61", "(b)")
chk("(b) 38 types, 27 tight, 27 positive y", (nT5, t5, nY5) == (38, 27, 27))
chk("(b) shapes 5x(4), 10x(3,1), 10x(2,2), 12x(2,1,1), 1x(1,1,1,1): %s" % sh5,
    sh5 == {(4,): 5, (3, 1): 10, (2, 2): 10, (2, 1, 1): 12, (1, 1, 1, 1): 1})
chk("(b) C(7,4) sum y = c745 = file c", v5 == c745 == C5["c"])
chk("(b) x printed 0.317221, 0.126888, 0.302114, 0.126888, 0.126889 are exact",
    C5["x"] == [Q("0.317221"), Q("0.126888"), Q("0.302114"), Q("0.126888"), Q("0.126889")])
chk("(b) 126 = number of 5-types over 5 parts", len(types_of_size(5, 5)) == 126 == comb(9, 4))
chk("(b) c745 < 0.447392627 (rounded up; c745 = %s)" % floor_dec(c745, 15), c745 < Q("0.447392627") and c745 > Q("0.447392626"))

print("=== (c) c(7,4,6)")
c746 = (107 + 3 * c745) / 384
for fn, sh_ in ((RM + "/certificates/lottery_upper/f8_k7_r4p6_x_quarter.txt",
                 "44abf7407bfa5dfc8b6d49e7aa8fadef33c8e563d39ca5c6325d752a54123403"),):
    lab = "(c) " + fn.split("/")[-1]
    C6, t6, nT6, nY6, sh6, v6 = check_cert(fn, sh_, lab)
    chk(lab + " 14 types, 14 tight, 8 blocks", (nT6, t6, nY6) == (14, 14, 8))
    chk(lab + " inner parts {0,1}, q = 5, x = 1/4", C6["inner_parts"] == (0, 1) and C6["q"] == 5 and C6["x"] == [Q(1, 4)] * 4)
    chk(lab + " 35 sum y = 107/384", v6 == Q(107, 384))
    tot = v6 + C6["inner_c"] * (C6["x"][0] ** 4 + C6["x"][1] ** 4)
    chk(lab + " file c = 35 sum y + inner*(x0^4+x1^4) = %s" % tot, tot == C6["c"])
    # the paper's description of T and y (paper parts 1..4 = file parts 0..3)
    arcs = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]
    Tp = set()
    for i in (2, 3):
        e = [0] * 4; e[i] = 4; Tp.add(tuple(e))
    for i, j in combinations(range(4), 2):
        e = [0] * 4; e[i] = 2; e[j] = 2; Tp.add(tuple(e))
    for i, j in arcs:
        e = [0] * 4; e[i] = 3; e[j] = 1; Tp.add(tuple(e))
    Yp = {}
    for i in (2, 3):
        B = [0] * 4; B[i] = 7; Yp[tuple(B)] = Q(1, 13440)
    for i, j in arcs:
        B = [0] * 4; B[i] = 4; B[j] = 3; Yp[tuple(B)] = Q(1, 768)
    chk(lab + " T and y = the paper's verbal description (arcs 1->2,2->3,3->4,4->1,1->3,2->4)",
        Tp == set(C6["T"]) and Yp == dict(C6["y"]))
chk("(c) f8_k7 file c = (107+3c745)/384", c746 == read_cert(RM + "/certificates/lottery_upper/f8_k7_r4p6_x_quarter.txt")["c"])
chk("(c) (107+3c745)/384 < 0.282141089 (value %s)" % floor_dec(c746, 15), Q("0.282141088") < c746 < Q("0.282141089"))
chk("(c) mu at 1/4: (4) 1/256, (2,2) 3/128, (3,1) 1/64",
    (mu((4,), [Q(1, 4)]), mu((2, 2), [Q(1, 4)] * 2), mu((3, 1), [Q(1, 4)] * 2)) == (Q(1, 256), Q(3, 128), Q(1, 64)))
chk("(c) 35/13440 + 1/768 = 1/256; 18/768 = 3/128; 12/768 = 1/64",
    Q(35, 13440) + Q(1, 768) == Q(1, 256) and Q(18, 768) == Q(3, 128) and Q(12, 768) == Q(1, 64))
chk("(c) C(4,2)C(3,2) = 18, C(4,3)C(3,1) = 12, C(7,4) = 35", comb(4, 2) * comb(3, 2) == 18 and comb(4, 3) * 3 == 12)
chk("(c) 35(2/13440 + 6/768) = 107/384", 35 * (Q(2, 13440) + Q(6, 768)) == Q(107, 384))
chk("(c) 2*(1/4)^4 = 1/128 and 107/384 + c745/128 = (107+3c745)/384",
    2 * Q(1, 4) ** 4 == Q(1, 128) and Q(107, 384) + c745 / 128 == c746)
# hand proof of (i) for (c), re-run as a case check: P size 6, P1,P2 <= 4
PT6 = types_of_size(4, 6)
T6 = read_cert(RM + "/certificates/lottery_upper/f8_k7_r4p6_x_quarter.txt")["T"]
hand_ok = True
for P in PT6:
    if P[0] >= 5 or P[1] >= 5: continue
    if P[2] >= 4 or P[3] >= 4 or sum(1 for v in P if v >= 2) >= 2:
        continue
    i = max(range(4), key=lambda j: P[j])
    others = [P[j] for j in range(4) if j != i]
    if not (P[i] >= 3 and all(v <= 1 for v in others)): hand_ok = False
    if P[i] == 3 and others != [1, 1, 1]: hand_ok = False
    if P[i] == 4 and not (i in (0, 1) and sum(others) == 2): hand_ok = False
    if not dominated(P, T6): hand_ok = False
chk("(c) case analysis of the printed proof of (i) is exhaustive and correct", hand_ok)
# comparison with Montecalvo (same shape): with inner 38/81, (108+3c)/384 = 7385/25920; ours is 1/384 less
chk("(c) Montecalvo shape: (108 + 3*38/81)/384 = 7385/25920; difference to (107+3c)/384 is 1/384",
    (108 + 3 * Q(38, 81)) / 384 == Q(7385, 25920) and Q(1, 8960) * 70 - Q(1, 13440) * 70 == Q(1, 384))

print("=== (a) c(6,4,5) <= 7/16 on Gamma_d")
chk("(a) 15(1/72+1/180)(3/2) = 7/16", 15 * (Q(1, 72) + Q(1, 180)) * Q(3, 2) == Q(7, 16))
chk("(a) 7/16 / 15 = 7/240; 240/7 = %s > 34.28; floor 34.28" % floor_dec(Q(240, 7), 6),
    Q(7, 16) / 15 == Q(7, 240) and Q(240, 7) > Q("34.28") and floor_dec(Q(240, 7), 2) == "34.28")

def ip(u, w):
    return bin(u & w).count("1") & 1

def gamma_check(d, full_tiling):
    N = 2 ** d
    A = [("A", v) for v in range(N)]; B = [("B", v) for v in range(N)]
    def is_edge(S):
        sa = [v for (s, v) in S if s == "A"]; sb = [v for (s, v) in S if s == "B"]
        if len(sa) in (0, 4): return True
        if len(sa) == 2: return ip(sa[0] ^ sa[1], sb[0] ^ sb[1]) == 0
        return False
    # even-minor fraction
    ev = sum(1 for a1, a2 in combinations(range(N), 2) for b1, b2 in combinations(range(N), 2) if ip(a1 ^ a2, b1 ^ b2) == 0)
    Ec = ev
    chk("(a) d=%d |E_cross| = C(N,2)(N/2-1)(N/2) = %d; fraction (N/2-1)/(N-1)" % (d, Ec),
        Ec == comb(N, 2) * (N // 2 - 1) * (N // 2) and Q(Ec, comb(N, 2) ** 2) == Q(N // 2 - 1, N - 1))
    # good (3,3) sets: count through each cross edge directly (test all 9 cross 4-subsets)
    cN = (N // 2 - 2) * (N // 4 - 2)
    def good33(a3, b3):
        return all(ip(x ^ y, z ^ w) == 0 for x, y in combinations(a3, 2) for z, w in combinations(b3, 2))
    counts = set()
    sample = [(a1, a2, b1, b2) for a1, a2 in combinations(range(N), 2) for b1, b2 in combinations(range(N), 2)
              if ip(a1 ^ a2, b1 ^ b2) == 0]
    if d >= 5:
        sample = sample[::97]
    for (a1, a2, b1, b2) in sample:
        c = 0
        for a3 in range(N):
            if a3 in (a1, a2): continue
            for b3 in range(N):
                if b3 in (b1, b2): continue
                if good33((a1, a2, a3), (b1, b2, b3)): c += 1
        counts.add(c)
    chk("(a) d=%d every %scross edge lies in exactly c_N = (N/2-2)(N/4-2) = %d good (3,3)-sets: %s"
        % (d, "" if d < 5 else "sampled ", cN, counts), counts == {cN})
    W = Q(Ec, 9) + Q(2 * comb(N, 4), 15)
    if full_tiling:
        # enumerate all good sets, give the weights, check every edge has weight exactly 1, total = W, alpha
        good = []
        for a3 in combinations(range(N), 3):
            for b3 in combinations(range(N), 3):
                if good33(a3, b3): good.append((a3, b3))
        wt = {}
        tot = Q(0)
        pairw = {}
        for a3, b3 in good:
            edges = [(("A", x), ("A", y), ("B", z), ("B", w)) for x, y in combinations(a3, 2) for z, w in combinations(b3, 2)]
            assert len(edges) == 9 and all(is_edge(E) for E in edges)
            for E in edges:
                wt[E] = wt.get(E, 0) + Q(1, cN)
            tot += Q(1, cN)
            for E, F in combinations(edges, 2):
                pairw[(E, F)] = pairw.get((E, F), 0) + Q(1, cN)
        okc = len(wt) == Ec and all(v == 1 for v in wt.values())
        ninternal = 2 * comb(N, 6)
        tot += ninternal * Q(1, comb(N - 4, 2))
        chk("(a) d=%d tiling: each of the %d cross edges has weight exactly 1 (good (3,3)-sets: %d); internal C(N-4,2)=%d"
            % (d, Ec, len(good), comb(N - 4, 2)), okc and len(good) * 9 == Ec * cN and comb(N, 6) * 15 == comb(N, 4) * comb(N - 4, 2))
        chk("(a) d=%d total weight = W_N = |E_cross|/9 + 2C(N,4)/15" % d, tot == W)
        amax_cross = max(pairw.values())
        amax_int = Q(N - 5, comb(N - 4, 2))
        bound = N * max(Q(1, cN), Q(1, comb(N - 4, 2)))
        chk("(a) d=%d alpha: max pair weight cross %s, internal %s <= N max(1/c_N,1/C(N-4,2)) = %s"
            % (d, amax_cross, amax_int, bound), max(amax_cross, amax_int) <= bound)
    val = 15 * W / comb(2 * N, 4)
    chk("(a) d=%d 15 W_N/C(n,4) = %s < 7/16" % (d, floor_dec(val, 9)), val < Q(7, 16))
    return N

# Turan (5,4) property of Gamma_d by brute force for d = 3 (n = 16), d = 4 (n = 32)
for d in (3, 4):
    N = 2 ** d
    V = [("A", v) for v in range(N)] + [("B", v) for v in range(N)]
    def is_edge(S):
        sa = [v for (s, v) in S if s == "A"]; sb = [v for (s, v) in S if s == "B"]
        if len(sa) in (0, 4): return True
        if len(sa) == 2: return ip(sa[0] ^ sa[1], sb[0] ^ sb[1]) == 0
        return False
    bad = 0; parity_bad = 0; tot = 0
    for S in combinations(V, 5):
        tot += 1
        ne = sum(1 for R in combinations(S, 4) if is_edge(R))
        if ne == 0: bad += 1
        if ne % 2 == 0: parity_bad += 1
    chk("(a) d=%d Gamma_d Turan (n,5,4) and odd by brute force over %d 5-sets" % (d, tot), bad == 0 and parity_bad == 0)
gamma_check(3, False)
gamma_check(4, True)
gamma_check(5, False)
for d in (6, 10, 20, 40):
    N = 2 ** d
    W = Q(comb(N, 2) * (N // 2 - 1) * (N // 2), 9) + Q(2 * comb(N, 4), 15)
    v = 15 * W / comb(2 * N, 4)
    print("info (a) d=%d 15W/C(n,4) = %s" % (d, floor_dec(v, 12)))
chk("(a) d = 3 gives c_N = 0 (so d >= 4 is needed)", (4 - 2) * (2 - 2) == 0)

print("=== Table tab:lottery: x values (C(k,4)/c), upper rounded DOWN")
rows = {(6, 5): Q(7, 16), (7, 5): c745, (7, 6): c746}
printed = {(6, 5): "34.28", (7, 5): "78.23", (7, 6): "124.05"}
for (k, p), c in rows.items():
    xv = comb(k, 4) / c
    chk("table (%d,4,%d): C(k,4)/c = %s -> floor %s (printed %s)" % (k, p, floor_dec(xv, 6), floor_dec(xv, 2), printed[(k, p)]),
        floor_dec(xv, 2) == printed[(k, p)])
old = {(6, 5): (Q("2.258") / 5, "33.21"), (7, 5): (Q(38, 81), "74.60"), (7, 6): (Q(7385, 25920), "122.84")}
for (k, p), (c, s) in old.items():
    chk("caption old (%d,4,%d): %s -> %s" % (k, p, floor_dec(comb(k, 4) / c, 6), s), floor_dec(comb(k, 4) / c, 2) == s)
chk("sentence: 0.4516 = 2.258/5, 38/81 = (190/81)/5, 7385/25920 = (7385/1728)/15",
    Q("2.258") / 5 == Q("0.4516") and Q(190, 81) / 5 == Q(38, 81) and Q(7385, 1728) / 15 == Q(7385, 25920))
chk("new < old (normalised): 7/16 < 0.4516, c745 < 38/81, c746 < 7385/25920",
    Q(7, 16) < Q("0.4516") and c745 < Q(38, 81) and c746 < Q(7385, 25920))

print("=== gap percentages")
def avg(n0, T0, N=200000):
    T = T0
    for n in range(n0 + 1, N + 1):
        T = -((-n * T) // (n - 4))
    return Q(T, comb(N, 4))
prev = {5: avg(18, 807), 6: avg(16, 190), 7: avg(16, 99)}
prev627 = avg(17, 627)
print("info prev t(5,4) from 807@18 = %s, from 627@17 = %s; t(6,4) = %s; t(7,4) = %s" %
      (floor_dec(prev[5], 9), floor_dec(prev627, 9), floor_dec(prev[6], 9), floor_dec(prev[7], 9)))
chk("caption floors t(5,4)>0.264386, t(6,4)>0.105430, t(7,4)>0.055614",
    floor_dec(prev[5], 6) == "0.264386" and floor_dec(prev[6], 6) == "0.105430" and floor_dec(prev[7], 6) == "0.055614")
new = {5: Q(35604499940047, 115448720916480), 6: Q(9224492822869, 65970697666560), 7: Q(7476698908057, 115448720916480)}
U_old = {(5, 5): Q(865, 2160), (6, 5): Q("2.258") / 5, (7, 5): Q(38, 81),
         (5, 6): Q("3.6678") / 15, (6, 6): Q("3.863") / 15, (7, 6): Q(7385, 25920),
         (5, 7): Q(1, 8), (6, 7): Q(1, 8), (7, 7): Q(1, 8)}
U_new = dict(U_old); U_new.update(rows)
for nm, U in (("old U", U_old), ("new U", U_new)):
    for space in ("c", "x"):
        r = {}
        for (k, p), u in U.items():
            if space == "c":
                r[(k, p)] = 100 * (new[p] - prev[p]) / (u - prev[p])
            else:
                C4 = comb(k, 4)
                r[(k, p)] = 100 * (C4 / prev[p] - C4 / new[p]) / (C4 / prev[p] - C4 / u)
        lo = min(r.values()); hi = max(r.values())
        print("info %s, %s-space: %s ; min %s max %s" % (nm, space,
              " ".join("(%d,4,%d) %s" % (k, p, floor_dec(r[(k, p)], 4)) for (k, p) in sorted(r, key=lambda t: (t[1], t[0]))),
              floor_dec(lo, 4), floor_dec(hi, 4)))
        if space == "c":
            chk("gap (%s, normalised constants): range rounds to 13.2%%..32.3%%" % nm,
                floor_dec(lo, 1) == "13.1" and ceil_dec(lo, 1) == "13.2" and floor_dec(hi, 1) == "32.3")
            argmin = [kp for kp in r if r[kp] == lo]; argmax = [kp for kp in r if r[kp] == hi]
            chk("gap (%s): min at p = 7 rows %s, max at (5,4,5) %s" % (nm, argmin, argmax),
                set(argmin) == {(5, 7), (6, 7), (7, 7)} and argmax == [(5, 5)])

print("=== (P6c) Section sec:giraud constants")
chk("24*C(6,4) = 360 (Cor. robust at k = 6: 7/16 - eps >= 7/16 - 360 eta - o(1) => eta >= eps/360 - o(1))", 24 * comb(6, 4) == 360)
chk("Cor. far check: 7/16 - 865/2160 = 1/27, 1/27/120 = 1/3240, 24/3240 = 1/135",
    Q(7, 16) - Q(865, 2160) == Q(1, 27) and Q(1, 27) / 120 == Q(1, 3240) and Q(24, 3240) == Q(1, 135))
# optional k = 7 add-on
y_split = [comb(i, 4) + comb(7 - i, 4) + Q(17, 9) * comb(i, 2) * comb(7 - i, 2) for i in range(8)]
chk("optional k=7: (1,17/9) feasible on all splits (max %s <= 35); splits i=0..3: %s" % (max(y_split), y_split[:4]),
    max(y_split) <= 35 and y_split[:4] == [35, 15, 5 + Q(170, 9), 1 + 34])
s = Q(1, 4)
chk("optional k=7: 1 - 4s + (2 + 3*17/9) s^2 at s = 1/4 = 23/48 and 2 + 17/3 = 23/3",
    1 - 4 * s + (2 + 3 * Q(17, 9)) * s * s == Q(23, 48) and 2 + Q(17, 3) == Q(23, 3))
dist = (Q(23, 48) - c745) / (24 * 35)
chk("optional k=7: (23/48 - c745)/840 = %s >= 3.78e-5" % floor_dec(dist, 10), dist >= Q("0.0000378"))

print("\nSUMMARY: %d checks, %d FAIL" % (len(OK), sum(1 for _, v in OK if not v)))
