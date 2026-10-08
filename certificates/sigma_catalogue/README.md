# sigma_catalogue: further codegree-squared densities at N = 7

**Proves** Theorem sigmacat (the upper bounds of Table sigmacat; exact values in Appendix A) and Remark j4limit of the
paper (Section "Further codegree-squared densities"). One subdirectory per certificate, each with its own README
(claim, files, checks, provenance), and one for the two files of Remark j4limit:

| directory | statement | forbidden graph | objective | M | exact bound b | sigma(F) <= 1 - b | < (12 digits) | Table sigmacat | review |
|---|---|---|---|---|---|---|---|---|---|
| `J4/` | Theorem sigmacat, line 1 | J_4 (3-graph) | cosig3 = F/210 | 2^30 | 23220181141640969669/32281802128991715328 | 9061620987350745659/32281802128991715328 | 0.280703690307 | 0.2807037 | R7_SIG7 (A1) |
| `K5lt/` | Theorem sigmacat, line 2 | K_5^< (3-graph) | cosig3 | 2^19 | 2747492282707/4123168604160 | 1375676321453/4123168604160 | 0.333645420191 | 0.3336455 | R7_SIG7 |
| `K5eq/` | Table sigmacat | K_5^= (3-graph) | cosig3 | 2^19 | 106600906627513/173173081374720 | 66572174747207/173173081374720 | 0.384425652179 | 0.3844257 | R7_SIG8 |
| `K5_3minus/` | Table sigmacat | K_5^{3-} (3-graph) | cosig3 | 2^19 | 1143076695549/1924145348608 | 781068653059/1924145348608 | 0.405930172388 | 0.4059302 | R7_SIG8 |
| `C5/` | Table sigmacat | C_5, tight 5-cycle (3-graph) | cosig3 | 2^18 | 6422651703/8589934592 | 2167282889/8589934592 | 0.252304935013 | 0.2523050 | R7_SIG8 |
| `K6_3/` | Table sigmacat | K_6^(3) (3-graph) | cosig3 | 2^19 | 29749027005683/115448720916480 | 85699693910797/115448720916480 | 0.742318262433 | 0.7423183 | R7_SIG8 |
| `K5_4minus/` | Theorem sigmacat, line 3 | K_5^(4)- (4-graph) | sigma4 = F/420 | 2^21 | 24504905122391/30786325577728 | 6281420455337/30786325577728 | 0.204032807991 | 0.2040329 | R7_SIG7 |
| `K6_4minus/` | Table sigmacat | K_6^(4)- (4-graph) | sigma4 | 2^21 | 358317739722203/923589767331840 | 565272027609637/923589767331840 | 0.612037993062 | 0.6120380 | R7_SIG8 |
| `K6_4/` | Table sigmacat | K_6^(4) (4-graph) | sigma4 | 2^20 | 140873242475/549755813888 | 408882571413/549755813888 | 0.743753064695 | 0.7437531 | R7_SIG8 |
| `K7_4/` | Table sigmacat | K_7^(4) (4-graph) | sigma4 | 2^20 | 9454851490593/76965813944320 | 67510962453727/76965813944320 | 0.877155180904 | 0.8771552 | R7_SIG8 |
| `J4_limit/` | Remark j4limit | J_4 | cosig3 | | dual point: every plain N = 7 certificate has b <= V; witness: 41/57 - delta | 1 - V > 0.280703275315 | | | R7_J4S1, R7_J4S2 |

Denominators: 32281802128991715328 = 7 * 2^62, 4123168604160 = 15 * 2^38, 173173081374720 = 315 * 2^39,
1924145348608 = 7 * 2^38, 8589934592 = 2^33, 115448720916480 = 105 * 2^40, 30786325577728 = 7 * 2^42,
923589767331840 = 105 * 2^43, 549755813888 = 2^39, 76965813944320 = 35 * 2^41. The decimals are upper bounds rounded
up (the lower end 0.280703275315 of Remark j4limit is rounded down). The rows of Table sigmacat in the order of the
paper are J_4, K_5^<, K_5^=, K_5^{3-}, C_5, K_6^(3) (3-graphs; J_5 has no seven-vertex bound) and K_5^(4)-, K_6^(4)-,
K_6^(4), K_7^(4) (4-graphs). The cut generation stopped at its time limit for K_5^<, K_5^= and K_6^(3), and for the
other certificates a "converged" cut LP need not reach the seven-vertex optimum either (findings L3 of R7_SIG7 and
R7_SIG8): the values are valid bounds, not limits of the method (for J_4 see `J4_limit/`).

## Objectives

The codegree-squared density sigma(F) is defined as in the paper (for 3-graphs with pairs in place of triples and the
normalisation C(n,2)(n-2)^2). Both objectives are densities of configurations on at most five vertices, so
Remark objective of the paper applies, and on a 7-vertex graph H they are F(H)/210 or F(H)/420 for an integer F(H):

* `cosig3` (3-graphs, F-free graphs themselves): f(H) = 1 - gamma(H) = F(H)/210 with
  F(H) = 210 - sum over the 35 4-subsets Q of [7] of C(e_H(Q), 2); gamma is the probability that P + x and P + y are
  both edges (P a uniform pair, (x, y) a uniform ordered pair of distinct vertices outside P). A certificate of value b
  gives gamma(G) <= 1 - b + O(1/n) for every F-free G, hence sigma(F) <= 1 - b by the identity
  co_2(G)/(C(n,2)(n-2)^2) = gamma(G) - (gamma(G) - d(G))/(n-2).
* `sigma4` (4-graphs, complement form): H is the complement of an F-free G; kappa(H) = 2 d(H) - gamma(H) = F(H)/420
  with F(H) = 24 e(H) - sum over the 35 3-subsets T of [7] of c_T (c_T - 1) (Lemma co2 of the paper). For
  F = K_5^(4)-, H spans at least two edges in every 5-set. A certificate of value b gives sigma(F) <= 1 - b.

## Certificate format (`<prefix>.<objective>.cert.json` + `.cert.npz`, and `<prefix>.keys.json`)

* `.cert.npz`: one int64 array `A<k>` per key k (shape (rows, n_flags)); the certificate matrix is
  Q_k = A_k^T A_k / M^2, PSD by construction. `M` = 2^`M_bits` is in the json (and, for J4, also in the npz).
* `.cert.json`: `problem` (s3_J4, s3_K5lt, c4_K5m), `objective`, `blocks`, `M`, `M_bits`, `tau` = [0, 0] (no
  stationarity term), `keys` (one entry per key in matrix order: `s`, `m`, `type_index`, `sigma`, `n_flags`; the J4
  file also has `aut_size` and `flags`), the exact `bound` b, `argmin_graph`, the producer's table `groups`
  ([F, max Z, argmax, count] for every value F of the objective numerator, over its class list),
  `raw_extensions_scanned` and `cert_npz_sha256`; the other fields are records of the producer's run.
* `.keys.json` (written for the package from the producer's flag model, and compared with the independent
  reviewer's own types, flags and |Aut| for all ten certificates: R7_SIG8, addenda A1 and A2): the full key list in
  matrix order, `k`, `s`, `m`, `type_index`, `sigma`
  (canonical labelled type mask on [s]), `aut_size` = |Aut(sigma)|, `n_flags`, `factor_rows` (rows of A_k), `flags`
  (canonical flag masks on [m], increasing = the row and column order of Q_k); the predicate (forbidden 3-graphs, 1-based
  triples, or p and lambda), `M`, `bound`, `sigma_upper` = 1 - b, a `format` string with the formulas, and the sha256 of
  the three certificate files next to it.
* `.verify.json`: the producer's own sampling check; a record, no role in the proof.

Masks: bit i of the mask of an r-graph on {0, ..., n-1} is the i-th r-subset in colex order (as in `../README.md`).
A type mask is the minimum over all relabellings of [s]; a flag mask the minimum over the relabellings of [m] that fix
0, ..., s-1; types of a block in increasing mask order (only those with a flag), flags in increasing mask order.
Writing Z(H) = 5040 M^2 sum_k c_{Q_k}(H) (an integer), the slack of a 7-vertex graph H is
(F(H) 5040 M^2 - den7 Z(H)) / (den7 5040 M^2) with den7 = 210 (cosig3) or 420 (sigma4), and b is its minimum over all
admissible H.

Validity needs only: every Q_k PSD (by construction), the flags of a key pairwise non-isomorphic sigma-flags, the types
of a block pairwise non-isomorphic, and the exact minimum over all admissible 7-vertex graphs. Completeness of the
flag lists is not needed for validity (the checkers verify it anyway).

## Checks

`verifier/verify_all.sh`, parts (m) and (n) (`verifier/EXPECTED.txt`); `--fast` runs them with everything else,
`--sigcat` runs only them (6 to 8 minutes). The programs are those of the independent reviews R7_SIG7 (J4, K5lt,
K5_4minus: `s7_*.py`, `s7_eval.c`, `s7_eval128.c`), R7_SIG8 (the other seven: `s8_*.py`, `s8_eval.c`; and the key
lists of all ten, `s8_keys_cmp.py` and `s8_keys_cmp_stage1.py`), R7_J4S1 (the witness: `jw_*.py`) and R7_J4S2 (the
dual point: `jd_*.py`), copied with label and path edits only (`verifier/CHANGES.txt`, items 12 and 14); the key lists
are also checked by the producer's self-contained `sk_check_keys.py` (not part of a review). The exhaustive raw scans
are `--raw sigcat_<directory>` (or `--raw sigcat` for all ten): for J4, K5lt and K5_4minus about 27, 20 and 8 minutes
on 2 threads in their review, with up to 1.5 GB of memory (the decoding of the J4 certificate); for the seven of
R7_SIG8 20 s to 13.6 min on 15 threads of a cloud machine in that review (the 3-graphs about 3 to 35 minutes, the
4-graphs about 2.2 to 2.3 hours on 2 threads by its measured rates; less than 0.3 GB of memory, and the evaluation
tables written by the decoding take 0.05 to 0.29 GB of disk each). The class-list scans of the reviews need the
producer's class lists (not included) and are not repeated; the raw scan is the decisive one.
