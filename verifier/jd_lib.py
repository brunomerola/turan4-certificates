"""Independent checker library for the dual-point check of Remark j4limit (reviewer code).  Builds on jw_lib.py (the
reviewer's own encodings, J_4 test, objective, moment definition) and adds:
  * canonical labelled types (one labelled representative per isomorphism class = min LEX mask over S_s -- a
    different representative choice than the producer's min colex code; any choice is valid, see relabel lemma);
  * accumulation of Y_sigma = sum_H num_H * M_sigma(H) (integer counts) restricted to canonical types;
  * an exact positive-definiteness test by Sylvester's criterion with the leading principal minors computed EXACTLY
    by multimodular elimination + CRT (Hadamard bound), independent of any floating-point hint, and the exact
    fraction-free (Bareiss) PSD test of jw_lib for small or singular matrices;
  * the relabelling map of flags under a root permutation (for the spot check of the relabel lemma).
Imports nothing from the producer.
"""
from __future__ import annotations

import functools
import itertools
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jw_lib as L  # noqa: E402

N = 7
BLOCKS = L.BLOCKS
_PERMS_S = {s: list(itertools.permutations(range(s))) for s in range(0, 6)}


def relabel_mask(mask, n, perm):
    """mask on [n] -> mask of the graph with triple {perm[a], perm[b], perm[c]} for each triple {a, b, c}."""
    out = 0
    for (a, b, c) in L.edges_of(mask, n):
        out |= 1 << L.TIDX[n][tuple(sorted((perm[a], perm[b], perm[c])))]
    return out


CANON_TYPE = {}
for _s in range(0, 6):
    _tab = {}
    for _m in range(1 << len(L.TRI[_s])):
        _tab[_m] = min(relabel_mask(_m, _s, p) for p in _PERMS_S[_s])
    CANON_TYPE[_s] = _tab
IS_CANON = {s: {m for m, c in CANON_TYPE[s].items() if m == c} for s in CANON_TYPE}


def moment_counts_canon(mask, blocks=BLOCKS):
    """As jw_lib.moment_counts, restricted to theta whose type (lex mask of H[theta] on [s]) is the canonical
    representative of its isomorphism class (so each key is (s, m, canonical labelled type))."""
    E = L.adj3(mask)
    out = {}
    for (s, m) in blocks:
        r = m - s
        ft = L._FTRI[(s, m)]
        st = L.TRI[s]
        canon_ok = IS_CANON[s]
        for th in itertools.permutations(range(N), s):
            tm = 0
            for i, (a, b, c) in enumerate(st):
                if E[th[a]][th[b]][th[c]]:
                    tm |= 1 << i
            if tm not in canon_ok:
                continue
            rest = [v for v in range(N) if v not in th]
            codes = {}
            for U in itertools.combinations(rest, r):
                best = None
                for oU in itertools.permutations(U):
                    w = th + oU
                    code = 0
                    for i, (a, b, c) in enumerate(ft):
                        if E[w[a]][w[b]][w[c]]:
                            code |= 1 << i
                    if best is None or code < best:
                        best = code
                codes[U] = best
            cnt = out.setdefault((s, m, tm), Counter())
            Us = list(codes)
            for U1 in Us:
                s1 = set(U1)
                for U2 in Us:
                    if s1.isdisjoint(U2):
                        cnt[(codes[U1], codes[U2])] += 1
    return out


def flag_graph(sig, code, s, m):
    """3-adjacency on [m] of the flag with type sig (lex on [s]) and free-part code (lex bits of the triples meeting
    {s..m-1}, in the minimising order)."""
    E = [[[False] * m for _ in range(m)] for _ in range(m)]
    tris = [L.TRI[s][i] for i in range(len(L.TRI[s])) if (sig >> i) & 1]
    tris += [t for i, t in enumerate(L._FTRI[(s, m)]) if (code >> i) & 1]
    for t in tris:
        for a, b, c in itertools.permutations(t):
            E[a][b][c] = True
    return E


def flag_code(E, roots, U, s, m):
    best = None
    for oU in itertools.permutations(U):
        w = tuple(roots) + oU
        code = 0
        for i, (a, b, c) in enumerate(L._FTRI[(s, m)]):
            if E[w[a]][w[b]][w[c]]:
                code |= 1 << i
        if best is None or code < best:
            best = code
    return best


@functools.lru_cache(maxsize=None)
def relabel_flag(sig, code, s, m, pi):
    """sigma-flag (sig, code) -> the tau-flag code, tau = sig^pi (tau has {i,j,k} iff sig has {pi i, pi j, pi k}):
    root i of the tau-flag is root pi(i) of the sigma-flag; free vertices unchanged."""
    E = flag_graph(sig, code, s, m)
    roots = tuple(pi[i] for i in range(s))
    return flag_code(E, roots, tuple(range(s, m)), s, m)


# ------------------------------------------------------------------ primes, multimodular leading minors
def is_prime32(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7):            # deterministic for n < 3,215,031,751
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


_PRIMES = []


def primes(k):
    q = (1 << 31) - 1 if not _PRIMES else _PRIMES[-1] - 2
    while len(_PRIMES) < k:
        if is_prime32(q):
            _PRIMES.append(q)
        q -= 2
    return _PRIMES[:k]


class Limbs:
    """Integer matrix as 30-bit limbs (numpy int64) for fast reduction mod many primes."""

    def __init__(self, A):
        self.n = len(A)
        mx = max((abs(x) for row in A for x in row), default=0)
        self.nl = max(1, (mx.bit_length() + 29) // 30)
        self.sign = np.array([[1 if x < 0 else 0 for x in row] for row in A], dtype=bool)
        absA = [[abs(x) for x in row] for row in A]
        self.limbs = np.zeros((self.nl, self.n, self.n), dtype=np.int64)
        for j in range(self.nl):
            sh = 30 * j
            self.limbs[j] = np.array([[(x >> sh) & ((1 << 30) - 1) for x in row] for row in absA], dtype=np.int64)

    def mod(self, p):
        acc = np.zeros((self.n, self.n), dtype=np.int64)
        for j in range(self.nl):
            w = pow(2, 30 * j, p)
            acc = (acc + (self.limbs[j] % p) * w) % p
        acc[self.sign] = (p - acc[self.sign]) % p
        return acc


def leading_minors_mod_p(Amod, p):
    """All leading principal minors of A mod p (A given reduced mod p as int64), elimination without pivoting;
    None if a pivot vanishes mod p."""
    A = Amod.copy()
    n = A.shape[0]
    out = []
    d = 1
    for k in range(n):
        piv = int(A[k, k])
        if piv == 0:
            return None
        d = d * piv % p
        out.append(d)
        if k + 1 < n:
            inv = pow(piv, p - 2, p)
            f = (A[k + 1:, k] * inv) % p
            A[k + 1:, k + 1:] = (A[k + 1:, k + 1:] - (f[:, None] * A[k, k + 1:][None, :]) % p) % p
    return out


def leading_minors_exact(A):
    """Exact leading principal minors det(A_k), k = 1..n, of an integer matrix: multimodular + CRT, number of primes
    from the Hadamard bound of the full rows (an upper bound for every leading minor)."""
    n = len(A)
    log2H = sum(0.5 * math.log2(max(1, sum(x * x for x in row))) for row in A)
    need_bits = log2H + 2 + 64              # margin
    res, mods = [], []
    k = 0
    bits = 0.0
    cand = primes(int(need_bits / 30) + 64)
    used = 0
    LA = Limbs(A)
    for p in cand:
        if bits >= need_bits:
            break
        r = leading_minors_mod_p(LA.mod(p), p)
        if r is None:
            continue
        res.append(r)
        mods.append(p)
        bits += math.log2(p)
        used += 1
    if bits < need_bits:
        raise RuntimeError("not enough primes")
    # CRT (Garner, incremental) per minor
    M = 1
    xs = [0] * n
    for r, p in zip(res, mods):
        Minv = pow(M % p, p - 2, p)
        for i in range(n):
            t = ((r[i] - xs[i]) % p) * Minv % p
            xs[i] += M * t
        M *= p
    half = M // 2
    xs = [x - M if x > half else x for x in xs]
    return xs, {"primes": used, "bits_needed": round(need_bits), "bits_used": round(bits)}


def bareiss_leading_minors(A):
    """Exact leading principal minors by fraction-free elimination without pivoting (for cross-checks; stops at the
    first zero minor)."""
    n = len(A)
    M = [list(r) for r in A]
    prev = 1
    out = []
    for k in range(n):
        d = M[k][k]
        out.append(d)
        if d == 0:
            return out
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                q, rr = divmod(M[i][j] * d - M[i][k] * M[k][j], prev)
                assert rr == 0
                M[i][j] = q
        prev = d
    return out
