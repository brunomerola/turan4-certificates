# n7_dual_points: what seven vertices can give for the other certificates

**Proves** the rows of Table n7limits of the paper other than t(5,4) (the t(5,4) row is Theorem n7opt, with its dual
point in `../n7_dual_point/`): no plain flag-algebra certificate on seven vertices (Definition cert: any blocks with
2m - s <= 7, any flag lists, mixed flag sizes allowed) proves more than V, where

| file | problem (p, lambda) | blocks | support | V (exact) | V rounded up | certified b (Theorem main) |
|---|---|---|---|---|---|---|
| `dual_p6_f3.json` | t(6,4) (6, 1) | all six | 884 | 37867960095539974857512229/(7*2^85) | 0.139837689411 | 9224492822869/(15*2^42) |
| `dual_p7_f1.json` | t(7,4) (7, 1) | all six | 1,058 | 6413961863433425283788748011/(5*2^94) | 0.064764464150 | 7476698908057/115448720916480 |
| `dual_h44_nt3.json` | t_2(5,4) (5, 2) | all six | 1,583 | 384021679403583125317285924719/(35*2^94) | 0.553946861054 | 48641756816487/(5*2^44) |
| `dual_k6m_W.json` | t_2(6,4) (6, 2) | all six | 2,126 | 153033054540662097928475789999/(35*2^94) | 0.220748423194 | 202684381843379/(105*2^43) |
| `dual_p5no56_W.json` | t(5,4) (5, 1) | no (5,6) | 1,569 | 390498160148386371476706381487/(35*2^95) | 0.281644555065 | --- |

So the best plain seven-vertex bound lies in [b, V]: V - b = 1.0567e-5 for t(6,4), 2.3875e-6 for t(7,4), 9.5411e-4 for
t_2(5,4) and 1.2956e-3 for t_2(6,4). In Turan form: no plain seven-vertex certificate proves pi(K_6^(4)) <
0.860162310589, pi(K_7^(4)) < 0.935235535850, pi(K_5^(4)-) < 0.446053138946 or pi(K_6^(4)-) < 0.779251576806.

The last row is the statement "without the (5,6) block": no plain seven-vertex certificate for t(5,4) whose types of
size five enter only through 1 x 1 terms (types of size at most 4 with any flag sizes, plus the 1 x 1 terms with
m <= 3 or m = s) proves t(5,4) > 0.281644555065, i.e. pi(K_5^(4)) < 0.718355444935. With the certificate of `../K5_4/`
(b = 0.308400990997...) the (5,6) block is therefore worth at least 0.026756435932.

**Supplementary file** (not used in the paper): `dual_p5no56_Wodd.json`, 55 graphs, all odd, for the same blocks as
`dual_p5no56_W.json`, with V' = 24406135009720244409729340413/(35*2^91) < 0.281644555070 (V' - V = 5.15e-12); it
bounds the odd-restricted problem without (5,6). It is checked like the others and is one of the inputs of
`rv2_lift.py`.

**Format** (as `../n7_dual_point/dual_Wodd_best.json`): `p`, `lam`, `blocks` (the blocks the dual point is for),
`den` = 2^110, `value` = V as a fraction, `support`, a list of `{mask, num}` with `mask` a 4-graph on {0, ..., 6}
(colex) and weight y_H = num/den. The field `source` names an intermediate file of the producer's search, not
included and not needed.

**Checked by** part (i) of `verifier/verify_all.sh` (fast), with `verifier/rv2_dual.py`, one run per file:

* den = 2^110, sum of the weights 1, all positive; the support graphs are distinct, pairwise non-isomorphic
  (canonical forms over all 5040 permutations) and admissible for "every p-set spans at least lambda edges"; the
  file's p, lambda and blocks equal the claimed ones (fixed in the program, not read from the file);
* V recomputed exactly and compared with the file and with the claimed fraction; decimals ceil12(V), floor12(1 - V);
* for every labelled type of every listed block the moment matrix Y_tau = sum_H num_H count_H,tau, built from the
  definition (every ordered injective theta, every ordered pair of disjoint flag sets) with the library `rv1_lib.py`
  of `rv1_dual.py`: all occurring flags in the own admissible flag universe, symmetric, rows with zero diagonal
  entirely zero, and the remaining part positive semidefinite by an exact test: symmetric elimination over the
  rationals (dimension <= 64), Bareiss fraction-free elimination (<= 160), and, for every size, a floating-point
  Cholesky hint followed by an exact integer residual that must be diagonally dominant (Gershgorin). Each test alone
  is a proof; where several apply they must agree. The largest live matrices have dimension 768;
* controls: on the largest matrix of each point the residual test must reject A - 2 lambda_min x x^T (proved
  indefinite exactly) and accept A - (lambda_min/2) x x^T; for the two no-(5,6) points the (5,6) block is built
  as well and must fail;
* weak-duality consistency with the certificate of Theorem main (`../K6_4/`, `../K7_4/`, `../K5_4minus/`,
  `../K6_4minus/`): b <= sum_H y_H (d(H) - c_Q(H)) <= V, every <Q_k, Y_k> >= 0, and b equal to the certificate's bound.

Also: `rv2_selftest.py` (flag universes against `rv1_lib.py`, canonical forms, admissibility against brute force, the
three PSD tests on random exact matrices), `rv2_lift.py` (the lifting identity that reduces mixed or smaller flag
sizes to the checked blocks, 27,254 exact identities on the 40 heaviest support graphs of each file and 20 random
graphs) and `rv2_decimals.py` (the decimals and differences quoted above). The residual test uses numpy's long
double for its hint; where long double is plain double (e.g. some ARM platforms) the hint may be too weak and the
test may fail to prove a true statement, but it cannot accept a false one.

**Provenance.** Found by the producer's dual-point search (facial reduction and a primal-dual interior-point solve
followed by rounding to a dyadic point; a second version of `../../search/flagalg/n7dual/`, not included), recorded
in commit d9a6a0c8 (branch research/turan4-n7) of the private working repository. The files are byte-identical to
those of the independent review R7_DUAL2 (2026-10-06, all six PASS; its notes concern the wording of the no-(5,6)
statement, used above), whose program `rv2_dual.py` and expected outputs `verify_all.sh` reproduces.
The dual points do not depend on the certificates. When the certificates of Theorem main (b), (d), (e) were
replaced by larger ones (review R7_B2), the claimed bounds b in `rv2_dual.py` and `rv2_decimals.py` and the paths
of the weak-duality partners were updated to the new certificates (`verifier/CHANGES.txt`, item 11); R7_DUAL2 had
run the weak-duality checks with the earlier certificates, and R7_B2 checked b <= V for the new ones.
