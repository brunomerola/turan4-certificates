# Verification guide

[Overview](README.md) · [Results and certificates](RESULTS.md)

Run both `--fast` and `--raw all` to complete the independent checks supplied in this package.
The first omits the long raw scans; the second runs those scans without rerunning the fast suite.
Optional producer checks provide additional comparisons. They are not independent checks.

## How to verify

Requirements: Python 3.11 or newer with numpy (2.0 or newer for part (m)); gcc with OpenMP (on macOS install gcc and
set `CC`, e.g. `CC=gcc-14`); bash; sha256sum or shasum. Check (f) also needs scipy and clarabel, and part (n) scipy
(`pip install numpy scipy clarabel`). The search code in `search/` has further dependencies but is not needed.

From the package root:

    verifier/verify_all.sh --fast                     # every check except the long raw scans
    verifier/verify_all.sh --raw cat_p6_lam11         # one full exhaustive raw scan (seconds)
    verifier/verify_all.sh --raw K5_4minus            # one full exhaustive raw scan (5 to 30 minutes)
    verifier/verify_all.sh --raw main --threads 8     # the five raw scans of Theorem main (hours)
    verifier/verify_all.sh --raw catalogue --threads 8   # the ten raw scans of Theorem catalogue (hours)
    verifier/verify_all.sh --raw sigma,stab_c10       # Theorem sigma (all raw extensions) and fact C10 (about 35 min)
    verifier/verify_all.sh --sigcat                   # only the fast checks of Theorem sigmacat and Remark j4limit (6-8 min)
    verifier/verify_all.sh --raw sigcat --threads 16  # the ten raw scans of Theorem sigmacat (about 9 h on 2 threads)
    verifier/verify_all.sh --lotup                    # only the checks of Proposition lotup (1-2 min)
    verifier/verify_all.sh --producer-n6              # optional: first implementation's check of the sharp N = 6 certificates
    verifier/verify_all.sh --producer-j4              # optional: the producer's own checks of the files of Remark j4limit

`--raw` takes a comma-separated list of names (K5_4, K6_4, K7_4, K5_4minus, K6_4minus, cat_p5_lam3, ..., cat_p7_lam4,
sigma, stab_c10, sigcat_J4, sigcat_K5lt, sigcat_K5_4minus, sigcat_K5eq, sigcat_K5_3minus, sigcat_C5, sigcat_K6_3,
sigcat_K6_4minus, sigcat_K6_4, sigcat_K7_4, j4_controls; see `--help`), or `main`, `catalogue`,
`sigcat` or `all`. Options can be combined. The driver writes
everything to `work/` (option `--work`), prints one PASS/FAIL line per step, keeps the full output of every step in
`work/logs/`, and exits with status 0 only if every check passed. Set `PYTHON` to choose the interpreter.
`verifier/EXPECTED.txt` lists the exact expected outputs and running times.

| part | what it checks | time |
|---|---|---|
| integrity | `SHA256SUMS`; each `.npz` against the `cert_npz_sha256` field of its `.json` (26 files: 15 certificates, the sigma key list and the ten certificates of Theorem sigmacat) | seconds |
| (a), fast | for each 4-graph certificate of Theorem main: own 6-vertex representatives and Burnside class counts (`rv_reps.py`), flag lists decoded, complete, non-isomorphic and admissible (`rv_prep.py`), exact value at the certificate's minimiser by a pure-Python evaluator (`rv_pycheck.py`) | about 15 min |
| (a), raw | the decisive step: exact minimum over **all** admissible one-vertex extensions with the C evaluator (`rv_eval.c`), comparison with the certificate (`rv_compare.py`), pure-Python value at the scan's minimiser | K5_4minus 7 to 11 min on 2 threads; each of the others 1 to 3.5 h on 2 threads, about 5 min on 28 to 32 threads |
| (b) | the three 5-graph certificates: exact PSD test and exact minimum over all 2^21 labelled 5-graphs (`rv_r5.py`) | about 5 min |
| (c) | the four six-vertex dual points (`rv_dual6.py`) | seconds |
| (d) | the seven-vertex dual point of Theorem n7opt (`rv1_dual.py`, `rv1_lift.py`, `rv1_selftest.py`) | about 2 min |
| (e) | Giraud's construction against the t(5,4) certificate, and the parity remark (`rv1_giraud.py`, `rv1_claim3.py`, a C scan of all 2^20 odd graphs) | about 5 min |
| (f) | the r = 3 comparison on five vertices (`indep_k43_n5.py`) | seconds |
| (g) | optional, `--producer-n6`: the first implementation's `search/flagalg/verify_cert.py` on the sharp six-vertex certificates | seconds |
| (h), fast | the ten catalogue certificates, with the same programs as (a) and lambda = 2 to 11; their exact bounds and the decimals of the paper's table (`cat_decimals`) | about 21 min |
| (h), raw | as (a), raw, for `cat_*` | cat_p5_lam3, cat_p6_lam11 seconds; cat_p6_lam9 13 min on 2 threads; the other seven several hours on 2 threads |
| (i) | the five dual points of Table n7limits and a supplementary one (`rv2_dual.py`, `rv2_selftest.py`, `rv2_lift.py`, `rv2_decimals.py`) | about 17 min |
| (j) | the three sharp six-vertex certificates, independently of the first implementation (`rv_sharp6.py`), and a negative control | 30 s |
| (k), fast | the sigma certificate: key list, PSD of X exactly, flag lists, the reduction identities and Giraud's law (`rs_prep.py`, `rs_reduction.py`, `rs_giraud7.py`), exact values at 26 graphs in Python and in C (`rs_pycheck.py`, `rs_eval.c`) | about 5 min |
| (k), raw | `--raw sigma`: exact minimum over all 86,952,880 admissible one-vertex extensions (`rs_eval.c`, `rs_analyze.py`): 33/64, tight only on Giraud's classes, next value 33/64 + mu_0 | about 10 min with `--threads 2` (15 to 20 CPU-min) |
| (l) | stability: Giraud systems on 7, 8, 9 vertices (`r_cf.py`), (C8), (C9), (C10) on the 33 class representatives (`r_lg.py`), the constants of Proposition margin and Theorem stab (`r_thm1.py`), Lemmas formF and rep on small cases (`r_formF.py`), the corollary numerics (`r_cor.py`) | about 6 min |
| (l), long | `--raw stab_c10`: (C10) over all 604,426 labelled bases (`r_lg10_full.py`) | about 23 min |
| (m), fast | the ten certificates of Theorem sigmacat: the key lists compared with an independent reviewer's own types and flags (`s8_keys_cmp.py`, `s8_keys_cmp_stage1.py`) and rebuilt by the producer's `sk_check_keys.py` with its negative controls, the reductions behind the objectives and the copy-list predicates (`s7_reduction.py`, `s8_predicates.py`), exact pure-Python values at 155 spot-check graphs equal to the independent reviews' (`s7_pycheck.py`, `s8_pycheck.py`), the decimals of the theorem and of Table sigmacat | about 1.5 min |
| (m), raw | `--raw sigcat_<name>` for the ten names: own decoding of the certificate (`s7_prep.py`, `s7_prep128.py`; `s8_prepvm.py` for the seven of R7_SIG8), C and Python values equal at the spot-check graphs, exact minimum over **all** admissible one-vertex extensions (`s7_eval.c`, `s7_eval128.c` with 128-bit Gram entries for J4, `s8_eval.c`), comparison with the certificate (`s7_analyze.py`, `s8_analyze.py`) | J4, K5lt, K5_4minus 27, 20 and 8 min on 2 threads (their review); the seven others 20 s to 13.6 min on 15 threads (their review), about 3 min to 2.3 h on 2 threads (estimate); up to 1.5 GB of memory |
| (n) | Remark j4limit: the exact dual point (`jd_verify_dual.py`: 757 graphs, all moment matrices PSD exactly) and the obstruction witness (`jw_verify.py`, its controls `jw_neg.py`, the lifting identity `jw_lift.py`), and the window b <= V | about 3 min |
| (n), long | `--raw j4_controls`: positive and negative controls of the dual-point checker (`jd_controls.py`) | 4 to 11 min |
| (o) | optional, `--producer-j4`: the producer's own checkers `search/sigma_catalogue/verify_witness.py` and `verify_dual3.py` on the two files of (n) | under 1 min |
| (p) | Proposition lotup (`--lotup`, also in `--fast`): the two type-pattern certificates (`rv3b_certs.py` with `rv3b_core.py`: conditions (i) and (ii) exactly, brute force on a blow-up, value, negative controls), the composition of (c) (`rv3b_derived.py`), the numerical companion of (a) on the Giraud host (`rv3b_giraud.py`), and 69 exact checks of the numbers of the section (`r7p6_check.py`) | 1 to 2 min |

Times were measured on one core of an AMD Ryzen 7 5800H laptop unless stated otherwise; memory stays below
600 MB, except for `--raw sigcat_J4` (1.5 GB) and `--raw sigcat_K5lt` (0.85 GB). The current `--fast` mode has 161
steps. The table gives timings per part; historical whole-run timings and the history of the package's development
versions are in `verifier/CHANGES.txt`, item 19. All fast checks and raw scans were run on clean cloud machines
before the initial public releases; Bruno also ran `--fast` (161/161 checks passed).

## Trust model

* **Exact certificates.** A 4-graph certificate is a list of integer matrices A_k with Q_k = A_k^T A_k / M^2, so
  every Q_k is positive semidefinite by construction. The 5-graph and sharp six-vertex certificates store rational
  matrices Q, and their positive semidefiniteness is tested exactly (LDL^T over the rationals). Dual points are
  rational vectors whose moment matrices are tested exactly: by elimination over the rationals, and, for the large
  matrices of `certificates/n7_dual_points/` (dimension up to 768), also by fraction-free (Bareiss) elimination up
  to dimension 160 and by a floating-point Cholesky hint whose exact integer residual must be diagonally dominant.
* **Exact value over all admissible graphs.** The bound b of a certificate is the minimum, over all admissible
  7-vertex graphs H, of d(H) - sum c_Q(H). The checker recomputes it in exact integer arithmetic (128-bit integers
  in C, fractions in Python). It enumerates every admissible graph as a labelled one-vertex extension R + L of its
  own 6-vertex representatives R (found by brute force over all 720 permutations of all 2^15 labelled 6-vertex
  graphs), with all 2^20 links L. Completeness follows from heredity alone (Lemma ext of the paper): no canonical
  forms on seven vertices are used. For each graph it enumerates every ordered injective root map with an exact
  labelled type test and every ordered pair of disjoint flag sets, without automorphism shortcuts.
* **Independent implementations.** The checkers in `verifier/` were written separately from the search code in
  `search/` and share no code with it. They read only the certificate files and use their own edge order, flag
  encoding, representatives and evaluators. Their results were compared with the producer's (same minimum, same
  maxima of Z on every group of graphs with equal edge count and equal second invariant), and the minimisers are
  re-evaluated by a second, pure-Python evaluator with explicit isomorphism search.
* **Not formally verified.** No proof assistant has checked the certificates or the programs. The proofs rely on
  the paper's Proposition fa and Lemma ext, on the correctness of the checkers, and on the compiler and hardware
  that ran them.
* **The sigma certificate** stores Q_k = W_k^T W_k / 2^72 + B_k^T X_k B_k: the first term is PSD by construction, the
  rational X_k are tested exactly (LDL^T); its objective is kappa = (24e - P)/420 instead of the edge density, and its
  checkers (part (k)) evaluate it in exact multi-limb integer arithmetic in C and with fractions in Python.
* **The certificates of Theorem sigmacat** have the format Q_k = A_k^T A_k / M^2 of the 4-graph certificates (so every
  Q_k is PSD by construction), for F-free 3-graphs (nine blocks; F = J_4, K_5^<, K_5^=, K_5^{3-}, C_5, K_6^(3)) and,
  in complement form, for F = K_5^(4)-, K_6^(4)-, K_6^(4), K_7^(4); the objectives are 1 - gamma (3-graphs) and
  kappa = 2d - gamma (4-graphs), see `certificates/sigma_catalogue/README.md`. Their checkers (part (m)) evaluate in
  exact 128-bit integer arithmetic in C (with 128-bit Gram entries for the J_4 certificate, M = 2^30) and with
  fractions in Python. The key lists (`*.keys.json`, written for the package because most certificate files do not
  list the flags) were compared with an independent reviewer's own types and flags for all ten certificates; the
  independent checkers of the certificates rebuild types and flags from the definitions and do not need them.
  The dual point of Remark j4limit is tested exactly (fraction-free elimination, or all leading principal minors by
  multimodular arithmetic), the witness by exact elimination and exact kernels.
* **Proposition lotup** (upper bounds for lottery numbers) rests on Kahn's fractional form of the
  Frankl-Rodl-Pippenger theorem (Kahn 1996, Theorem 1.2(b)) and on the existence of the limits l(k,r,p) (Sidorenko
  2023), which are cited, not checked here. Its certificates (`certificates/lottery_upper/`) are finite lists of exact
  rationals; part (p) checks the two conditions of Lemma pattern in exact arithmetic, by type-level domination and by
  brute force, with the programs of two independent reviews. Part (a) is a construction proved by hand; the package
  checks its counts numerically.
* **Limits.** The exact fraction printed by check (f) depends on the floating-point SDP solution and hence on the
  clarabel version; the check itself only requires the bound to beat Chung and Lu's. Of the ten catalogue raw
  scans, the independent review repeated four on its own machine; the other six were run on the cloud machines
  with the same checker programs, and the review checked their outputs against the certificates (see
  `certificates/catalogue/`). The raw scan of the t(6,4) certificate (Theorem main (b)) was run by its independent
  review with the review's own build and inputs on a temporary cloud machine, not on the review's own machine (see
  `certificates/K6_4/`). The raw scans of Theorem sigmacat were run by their independent reviews with the programs of
  part (m): J4, K5lt and K5_4minus on the reviewer's own machine, the other seven on a temporary cloud machine
  (16 CPUs). The initial packaging test covered decoding and C spot values for K5_4minus, C5 and K7_4, reproducing
  the reviews' evaluation tables byte for byte. The later cloud validation ran all raw scans through the package
  driver (see `verifier/CHANGES.txt`, item 18, for two driver fixes found during that validation). The reviews'
  second scans over the producer's class lists are not repeated here because those lists are not included.
  Statements that the paper labels as floating-point observations are not checked here.
