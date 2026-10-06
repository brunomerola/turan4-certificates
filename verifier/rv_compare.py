"""Independent checker: compare an rv_eval output (list or raw scan) with the producer certificate json.
Checks: exact bound = NMIN / DEN (reduced) equals the certificate's "bound"; group table (e, cdd): counts (list scan of
the producer's class list only) and maxima of Zsum equal the certificate's groups; argmin graph.
Usage: python rv_compare.py CERT_JSON RVLOG [--raw]
"""
import json
import sys
from fractions import Fraction

cert = json.load(open(sys.argv[1]))
raw = "--raw" in sys.argv
G, low = {}, []
nmin = den = None
for line in open(sys.argv[2]):
    t = line.split()
    if not t:
        continue
    if t[0] == "NMIN":
        nmin, den = int(t[1]), int(t[3])
        argc = int(t[7])
    elif t[0] == "G":
        G[(int(t[1]), int(t[2]))] = (int(t[3]), int(t[4]), int(t[5]))
    elif t[0] == "TOTAL":
        total = t
    elif t[0] == "LOW":
        low.append((int(t[1]), int(t[2])))
low.sort()
# rv_eval2 binary (used for the p=5/6 list and raw runs) printed the SUMS numerator in the NMIN field (buffer reuse,
# print-only bug, fixed in rv_eval.c); the exact minimum is the smallest entry of the LOW list (the KLOW lowest
# numerators over all scanned graphs), and its graph must be the reported argmin.
if low[0][0] != nmin:
    print(f"note: NMIN field {nmin} != min LOW {low[0][0]} (rv_eval2 print bug); using min LOW")
    nmin = low[0][0]
_a = [n for n, h in low if h == argc]
assert not _a or _a[0] == nmin, (low[0], argc)   # raw scans: many labelled copies tie at the minimum
b = Fraction(nmin, den)
cb = Fraction(cert["bound"])
print("TOTAL", " ".join(total[1:]))
print(f"rv bound {b} = {float(b):.15f}; cert bound {cb} = {float(cb):.15f}; EQUAL: {b == cb}")
print(f"rv argmin (colex) {argc}; cert argmin {cert['argmin_graph']}")
CG = {(g[0], g[1]): (int(g[4]), int(g[2]), int(g[3])) for g in cert["groups"]}
same_keys = set(G) == set(CG)
mz_ok = all(G[k][1] == CG[k][1] for k in CG if k in G)
cnt_ok = all(G[k][0] == CG[k][0] for k in CG if k in G)
print(f"groups: rv {len(G)} cert {len(CG)} same set: {same_keys}; maxZ equal on all: {mz_ok}; counts equal: {cnt_ok}"
      + (" (counts not comparable for a raw scan)" if raw else ""))
bad = [k for k in CG if k in G and G[k][1] != CG[k][1]][:5]
if bad:
    print("first maxZ mismatches:", [(k, G[k][1], CG[k][1]) for k in bad])
low.sort()
D = den
print("lowest values N/DEN:", [f"{nmin_ / D:.12f}" for nmin_, _ in low[:10]])
margin = [(n - nmin) / D for n, _ in low]
print("distinct H among lowest-64 list with value < bound + 1e-6:", sum(1 for n, _ in low if (n - nmin) / D < 1e-6))
