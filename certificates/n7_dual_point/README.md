# n7_dual_point: the seven-vertex optimum for t(5,4)

**Proves** Theorem n7opt of the paper: no plain flag-algebra certificate on seven vertices proves

    t(5,4) > V = 1336237682928914994292138923/(7*2^89) = 0.30840120119463679...,

so the best plain seven-vertex bound lies in [0.308400990997, 0.308401201195]. The same holds for the problem
restricted to odd 4-graphs, because every support graph is odd.

**File** `dual_Wodd_best.json`: `den` = 2^110; `value` = V (stored as 1336237682928914994292138923/
4332790137498830962146934784); `support`, 131 entries `{mask, num}` with `mask` a 4-graph on {0, ..., 6} (colex) and
weight y_H = num/den. The field `source` names an intermediate file of the producer's search, not included and not
needed.

**Checked by** part (d) of `verifier/verify_all.sh`:

* `rv1_dual.py`: den = 2^110, sum of weights 1, all positive; 131 pairwise non-isomorphic, admissible, odd graphs;
  V recomputed exactly; for every labelled type of the six pairs (1,4), (2,4), (3,4), (3,5), (4,5), (5,6) the moment
  matrix (own flag enumeration) is symmetric and PSD by exact elimination; a negative control (a point mass on the
  certificate's argmin) must fail; weak-duality consistency with the certificate of `K5_4/`.
* `rv1_lift.py`: the lifting identity that makes moment matrices of smaller or mixed flag sizes sums of the checked
  ones (17,045 exact identities on the support and on 30 random graphs).
* `rv1_selftest.py`: self-tests of the library (encodings; flag counts of the paper's block table).
