"""Independent checker: compare the pure-Python exact S (rs_pycheck.py) with rs_eval.c limbs ("each" mode), graph by graph."""
import sys
C = {}
for l in open(sys.argv[1]):
    t = l.split()
    if t and t[0] == "EACH":
        C[int(t[2])] = sum(int(v) << (45 * j) for j, v in enumerate(t[6:]))
n = ok = 0
for l in open(sys.argv[2]):
    t = l.split()
    if t and t[0] == "PY":
        h, S = int(t[1]), int(t[-1])
        n += 1
        ok += C.get(h) == S
        print(h, "PY S == C S:", C.get(h) == S)
print(f"{ok}/{n} identical")
