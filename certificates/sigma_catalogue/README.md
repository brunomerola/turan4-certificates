# sigma_catalogue: further codegree-squared densities at N = 7

**Proves** Theorem sigmacat and Remark j4limit of the paper (Section "Further codegree-squared densities"). One
subdirectory per certificate, each with its own README (claim, files, checks, provenance), and one for the two files of
Remark j4limit:

| directory | statement | forbidden graph | objective | exact bound b | sigma(F) <= 1 - b | < |
|---|---|---|---|---|---|---|
| `J4/` | Theorem sigmacat, line 1 | J_4 (3-graph) | cosig3 = F/210 | 23220181141640969669/32281802128991715328 | 9061620987350745659/32281802128991715328 | 0.280703690307 |
| `K5lt/` | Theorem sigmacat, line 2 | K_5^< (3-graph) | cosig3 = F/210 | 2747492282707/4123168604160 | 1375676321453/4123168604160 | 0.333645420191 |
| `K5_4minus/` | Theorem sigmacat, line 3 | K_5^(4)- (4-graph) | sigma4 = F/420 | 24504905122391/30786325577728 | 6281420455337/30786325577728 | 0.204032807991 |
| `J4_limit/` | Remark j4limit | J_4 | cosig3 | dual point: every plain N = 7 certificate has b <= V; witness: 41/57 - delta | 1 - V > 0.280703275315 | |

Denominators: 32281802128991715328 = 7 * 2^62, 4123168604160 = 15 * 2^38, 30786325577728 = 7 * 2^42. The decimals are
upper bounds rounded up (the lower end 0.280703275315 of Remark j4limit is rounded down).

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
* `.keys.json` (written for the package): the full key list in matrix order, `k`, `s`, `m`, `type_index`, `sigma`
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
`--sigcat` runs only them (about 3 minutes). The programs are those of the independent reviews R7_SIG7 (the three
certificates: `s7_*.py`, `s7_eval.c`, `s7_eval128.c`), R7_J4S1 (the witness: `jw_*.py`) and R7_J4S2 (the dual point:
`jd_*.py`), copied with label and path edits only (`verifier/CHANGES.txt`, item 12); the key lists are checked by the
producer's self-contained `sk_check_keys.py` (not part of a review). The exhaustive raw scans are `--raw sigcat_J4`,
`sigcat_K5lt`, `sigcat_K5_4minus` (or `--raw sigcat`): on 2 threads about 27, 20 and 8 minutes in the review, with up
to 1.5 GB of memory (the decoding of the J4 certificate). The class-list scans of the review need the producer's class
lists (not included) and are not repeated; the raw scan is the decisive one.
