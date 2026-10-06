"""Independent re-check (2026-10-06) of the producer value t(4,3) >= 0.42692378 from plain flag
algebras on N = 5 vertices (validation table of the producer code, not included). Own code from the definitions; imports only
numpy, scipy and clarabel (no producer module).

t(4,3) = lim min density of 3-graphs in which every 4-set contains an edge (= 1 - pi(K_4^(3))).
Certificate value b = min over ALL labelled admissible 3-graphs H on [5] of d(H) - sum_blocks c_Q(H), with
c_Q(H) = E_{(theta, U1, U2)} Q[flag(theta,U1), flag(theta,U2)] * [H[theta] = sigma], (theta, U1, U2) uniform among
injective s-tuples and ordered pairs of disjoint (m-s)-sets, for every block with 2m - s <= 5:
  (1,3): one type, flags = edge / no edge on the 3 vertices (2 flags);
  (2,3): one type, flags = edge / no edge (2 flags);
  (3,4): types = 3 labelled vertices without / with an edge; flags = link of the 4th vertex in the 3 root pairs (8).
Soundness (as in the paper, Proposition prop:fa): t(4,3) >= b for every PSD Q. The float SDP only proposes Q; the bound is
recomputed exactly (Fractions) for Q = A^T A / M^2 with integer A, PSD by construction.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial, sqrt

import numpy as np
import scipy.sparse as sp
import clarabel

V = range(5)
T5 = list(combinations(V, 3))                       # 10 triples
TI = {t: i for i, t in enumerate(T5)}


def has(H, tri):
    return (H >> TI[tuple(sorted(tri))]) & 1


def admissible(H):
    return all(any(has(H, t) for t in combinations(Q, 3)) for Q in combinations(V, 4))


ADM = [H for H in range(1 << 10) if admissible(H)]


def conf(s, m, n=5):
    return factorial(n) // factorial(n - s) * comb(n - s, m - s) * comb(n - m, m - s)


# blocks: (s, m, list of types); a type is a function H, theta -> bool; flag index function
def type_ok(s, sig, H, th):
    if s == 3:
        return has(H, th) == sig
    return True


def flag(s, m, H, th, U):
    if (s, m) == (1, 3):
        return has(H, (th[0],) + tuple(U))
    if (s, m) == (2, 3):
        return has(H, th + tuple(U))
    u = U[0]                                        # (3,4): link of u in the root pairs (th0th1, th0th2, th1th2)
    return has(H, (th[0], th[1], u)) | (has(H, (th[0], th[2], u)) << 1) | (has(H, (th[1], th[2], u)) << 2)


BLOCKS = [(1, 3, 0, 2), (2, 3, 0, 2), (3, 4, 0, 8), (3, 4, 1, 8)]   # (s, m, sigma, nflags)


def moment(H, blk):
    s, m, sig, nf = blk
    Mx = [[0] * nf for _ in range(nf)]
    for th in permutations(V, s):
        if not type_ok(s, sig, H, th):
            continue
        rest = [v for v in V if v not in th]
        for U1 in combinations(rest, m - s):
            for U2 in combinations([v for v in rest if v not in U1], m - s):
                Mx[flag(s, m, H, th, U1)][flag(s, m, H, th, U2)] += 1
    return Mx, conf(s, m)


MOM = {H: [moment(H, b) for b in BLOCKS] for H in ADM}
DENS = {H: Fraction(bin(H).count("1"), 10) for H in ADM}

# ---- float SDP with clarabel: x = (b, svec(Q_1), ..., svec(Q_4)), maximise b
idx, off = [], 1
for _, _, _, nf in BLOCKS:
    pos = {}
    for j in range(nf):                              # clarabel PSDTriangle: upper triangle, column-wise
        for i in range(j + 1):
            pos[(i, j)] = off
            off += 1
    idx.append(pos)
nv = off
rows, rhs = [], []
for H in ADM:
    r = np.zeros(nv)
    r[0] = 1.0
    for k, (blk, (Mx, cf)) in enumerate(zip(BLOCKS, MOM[H])):
        nf = blk[3]
        for i in range(nf):
            for j in range(nf):
                a, b_ = min(i, j), max(i, j)
                c = Mx[i][j] / cf
                r[idx[k][(a, b_)]] += c if a == b_ else c / sqrt(2)
    rows.append(r)
    rhs.append(float(DENS[H]))
A1 = np.array(rows)
nps = nv - 1
A2 = np.hstack([np.zeros((nps, 1)), -np.eye(nps)])
A = sp.csc_matrix(np.vstack([A1, A2]))
bvec = np.concatenate([rhs, np.zeros(nps)])
cones = [clarabel.NonnegativeConeT(len(ADM))] + [clarabel.PSDTriangleConeT(blk[3]) for blk in BLOCKS]
P = sp.csc_matrix((nv, nv))
q = np.zeros(nv)
q[0] = -1.0
st = clarabel.DefaultSettings()
st.verbose = False
sol = clarabel.DefaultSolver(P, q, A, bvec, cones, st).solve()
x = np.array(sol.x)
print("clarabel status", sol.status, "float b", x[0])


def unsvec(k, nf):
    Q = np.zeros((nf, nf))
    for (i, j), p in idx[k].items():
        v = x[p] if i == j else x[p] / sqrt(2)
        Q[i, j] = Q[j, i] = v
    return Q


# ---- exact rounding and exact evaluation
MSC = 1 << 24
Qex = []
for k, blk in enumerate(BLOCKS):
    Q = unsvec(k, blk[3])
    w, U = np.linalg.eigh(Q)
    keep = w > 1e-9 * max(w.max(), 1e-300)
    Aint = np.rint(MSC * (np.sqrt(w[keep])[:, None] * U[:, keep].T)).astype(np.int64)
    G = [[sum(int(Aint[r][i]) * int(Aint[r][j]) for r in range(Aint.shape[0])) for j in range(blk[3])]
         for i in range(blk[3])]
    Qex.append([[Fraction(G[i][j], MSC * MSC) for j in range(blk[3])] for i in range(blk[3])])
best, arg = None, None
for H in ADM:
    val = DENS[H]
    for k, (blk, (Mx, cf)) in enumerate(zip(BLOCKS, MOM[H])):
        nf = blk[3]
        val -= sum(Qex[k][i][j] * Mx[i][j] for i in range(nf) for j in range(nf)) / cf
    if best is None or val < best:
        best, arg = val, H
print("labelled admissible graphs", len(ADM))
print("exact b =", best, "=", float(best), "argmin edges", bin(arg).count("1"))
cl = (3 + sqrt(17)) / 12
print("Chung-Lu: pi(K4^3) <= %.9f, t(4,3) >= %.9f ; exact b > Chung-Lu: %s" % (cl, 1 - cl, float(best) > 1 - cl))
print("producer engine value 0.42692378221580 ; difference %.3e" % (float(best) - 0.42692378221580))
