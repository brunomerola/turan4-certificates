# sigma_catalogue/J4_limit: what plain seven-vertex certificates can give for sigma(J_4)

**Proves** Remark j4limit of the paper (Section "Further codegree-squared densities"). Setting: plain flag-algebra
certificates on N = 7 vertices for J_4-free 3-graphs with the objective f = 1 - gamma (`cosig3`, see `../README.md`):
for every block (s, m) with m >= 3 and 2m - s <= 7 and every labelled type sigma a real PSD matrix Q on any set of
sigma-flags of size m (mixed flag sizes of one type and 1 x 1 terms included), and a real b with
f(H) - sum <Q, M(H)> >= b for every J_4-free 3-graph H on [7]. Such a certificate proves sigma(J_4) <= 1 - b. Not
covered: stationarity or other extra linear constraints, more than seven vertices, other objectives.

**1. The exact dual point** (`dual_W2.json`; the strict statement of the remark). A probability vector
y = num / 2^110 on 757 pairwise non-isomorphic J_4-free 3-graphs on [7] whose moment matrices Y = sum_H y_H M(H) are
positive semidefinite for all nine blocks and all types (exactly). By weak duality (as in Lemma dual of the paper)
every plain seven-vertex certificate has

    b <= V = sum_H y_H f(H) = 598379856892795294218114109955/831895706399775544732211478528 = 0.71929672468490...,

so it proves at best sigma(J_4) <= 1 - V = 0.28070327531...; no such certificate proves any bound below
0.280703275315 (1 - V rounded down), and 16/57 is missed by 41/57 - V = 1.5209e-6. With the J_4 certificate of
`../J4/` (b = 23220181141640969669/32281802128991715328 <= V; V - b = 4.15e-7), the plain seven-vertex optimum lies in
(0.280703275315, 0.280703690307].

**2. The structural obstruction** (`witness_Fanoc.json`). phi is the exact 7-vertex law of the iterated blow-up of the
complement of the Fano plane: 93 isomorphism classes, E_phi f = 41/57 (sigma = 16/57). The witness is a rational
probability vector mu on 87 of these classes with X^T M(mu) X = 0 for every block, every type and every basis X of the
kernel of M(phi) on its occurring flags. Hence every plain seven-vertex certificate with <Q, M(phi)> = 0 for every
block (which a certificate of value 41/57 would need) has

    b <= E_mu f = 41/57 - delta,   delta = 316891/14980607589 = 2.1153e-5,

i.e. it proves at best sigma(J_4) <= 16/57 + delta; delta is the exact optimum of this argument (a linear programme
over the null space of the conditions, solved exactly). In particular no plain seven-vertex certificate has b = 41/57.

**Files**

| file | role | format |
|---|---|---|
| `dual_W2.json` | the dual point of part 1 | `den` = 2^110, `value` = V, `support`: list of {`mask` (colex 3-subsets of [7]), `num`} |
| `witness_Fanoc.json` | the witness of part 2 | `construction` = Fanoc, `target` = 41/57, `classes` (colex masks of the 93 classes of phi), `mu` (fractions), `delta`; `maker`: a record of the producer's run |
| `dual_W1.json` | a second, weaker dual point (649 graphs, 41/57 - V = 1.5167e-6) | as `dual_W2.json`; used as a positive control |
| `ctrl_pos_Sn.json`, `ctrl_pos_Fanoc.json` | exact 7-vertex laws of S_n (V = 25/27) and of phi (V = 41/57) | as `dual_W2.json`; positive controls |
| `ctrl_neg_lp.json`, `ctrl_neg_dead.json`, `ctrl_neg_j4.json` | a rounded LP point (not PSD), a point with 2^-40 moved onto a graph outside the face of the dual point (a zero diagonal entry with a nonzero row), a point with 2^-40 on K_7^(3) (not J_4-free) | as `dual_W2.json`; negative controls (must fail) |
| `neg_moved.json`, `neg_outside.json`, `neg_phi.json` | the witness with 1e-6 moved between two classes, with 1e-6 moved outside the support of phi (onto K_7^(3)), and mu = phi (delta = 0) | as `witness_Fanoc.json`; controls (exit 1, 1, 2) |

All files are byte-identical to the producer's files that the independent reviews read (their sha256 are in
`SHA256SUMS`; e.g. dual_W2.json 3d024afe..., witness_Fanoc.json 6398d47e...). The fields `note`, `source` and `maker`
are records of the producer's run.

**Checked by** `verifier/verify_all.sh` (`verifier/EXPECTED.txt`, part (n)), with the programs of the two independent
reviews (Python with integers and fractions; numpy only for canonical forms and rank estimates mod p; scipy only to
find the active set of the linear programme, which is then solved exactly):

* fast: `j4_s2` (`jd_verify_dual.py`): the 757 graphs pairwise non-isomorphic and J_4-free, sum y = 1, V exactly, the
  34 occurring (block, type) moment matrices built from the definition and PSD exactly (fraction-free elimination for
  n <= 40; for n > 40 all leading principal minors positive, computed exactly by multimodular elimination and Chinese
  remaindering with a Hadamard bound; one cross-check against Bareiss minors), the relabelling lemma on 12 graphs;
  `j4_s1` (`jw_verify.py`): phi recomputed from the construction (24,727 labelled graphs, 93 classes, E_phi f = 41/57),
  mu on the support of phi, the 191 (block, labelled type) matrices M(phi) PSD, their kernels, the conditions on mu,
  delta, and the exact linear programme (the null space of the conditions has dimension 7, the optimum is the
  witness); `j4_s1_controls` (`jw_neg.py`: corrupted witnesses fail, a second valid witness and mu = phi behave as
  expected, deliberately wrong moment definitions); `j4_s1_lift` (`jw_lift.py`: the lifting identity for mixed flag
  sizes on 33 graphs); `j4_window` (b <= V and the two decimals above);
* `--raw j4_controls` (`jd_controls.py`, 4 to 11 minutes): the positive and negative controls of the dual-point
  checker, on the files above and on points it builds itself;
* optional `--producer-j4`: the producer's own self-contained checkers `search/sigma_catalogue/verify_witness.py` and
  `verify_dual3.py` on the two files (not independent; same verdicts).

**Provenance.** Producer steps of 2026-10-06/07 (commits 0f92052c, the witness, and fa77dde5, the dual point; branch
research/turan4-n7 of the private working repository). Independent reviews (by one reviewer who did not produce
them): R7_J4S1 of the witness, PASS-with-notes (wording: part 2 alone shows only that 41/57 is not attained), and
R7_J4S2 of the dual point, PASS-with-notes (wording; its window used the J_4 certificate, reviewed separately in
R7_SIG7 addendum A1, PASS).
