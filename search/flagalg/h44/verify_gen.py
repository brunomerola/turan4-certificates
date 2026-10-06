"""Independent sample verifier (../n7/verify7.py, pure Python) for the threshold problem: verify7.main with its
admissibility test replaced by 'every p-set spans >= LAM edges' (used for flag re-enumeration and random samples).
Usage: python verify_gen.py LAM RUN_PREFIX [--random K]
"""
import os
import sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "n7"))
import verify7  # noqa: E402

LAM = int(sys.argv[1])


def admissible(es, n, p):
    return all(sum(frozenset(f) in es for f in combinations(P, 4)) >= LAM for P in combinations(range(n), p))


verify7.admissible = admissible
sys.argv = [sys.argv[0]] + sys.argv[2:]
print(f"# verify_gen: every p-set spans >= {LAM} edges", flush=True)
verify7.main()
