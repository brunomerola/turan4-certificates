# Turán densities of complete 4-graphs: certificates and checkers

Data package for the paper

> Bruno Mérola Corrêa, *Turán densities of complete 4-graphs via flag algebras on seven vertices*,
> preprint, 2026 (arXiv identifier to be added after submission).

It contains the exact rational certificates behind the paper's computer-assisted bounds, independent programs that
re-check them from the certificate files alone, a single driver that runs every check, and, for transparency, the
code that found the certificates.

## Results and files

| statement in the paper | bound | files | checks (`verifier/verify_all.sh`) |
|---|---|---|---|
| Theorem main (a) | t(5,4) >= 35604499940047/115448720916480 > 0.308400990997, so pi(K_5^(4)) < 0.691599009003 | `certificates/K5_4/` | (a) |
| Theorem main (b) | t(6,4) >= 9224492822869/65970697666560 > 0.139827122482, so pi(K_6^(4)) < 0.860172877518 | `certificates/K6_4/` | (a) |
| Theorem main (c) | t(7,4) >= 7476698908057/115448720916480 > 0.064762076605, so pi(K_7^(4)) < 0.935237923395 | `certificates/K7_4/` | (a) |
| Theorem main (d) | t_2(5,4) >= 48641756816487/87960930222080 > 0.552992751369, so pi(K_5^(4)-) < 0.447007248631 | `certificates/K5_4minus/` | (a) |
| Theorem main (e) | t_2(6,4) >= 202684381843379/923589767331840 > 0.219452823117, so pi(K_6^(4)-) < 0.780547176883 | `certificates/K6_4minus/` | (a) |
| Theorem catalogue (Table catalogue) | t_lambda(p,4) for (p, lambda) = (5,3), (6,3), (6,4), (6,5), (6,6), (6,9), (6,11), (7,2), (7,3), (7,4); exact fractions in `certificates/catalogue/README.md` | `certificates/catalogue/` | (h) |
| Theorem five | pi(K_6^(5)) < 0.767400802744, pi(H^5_4) < 0.389708097035, pi(K_7^(5)) < 0.927917451132 (exact fractions in Appendix A) | `certificates/fivegraphs/` | (b) |
| Theorem n6 | on six vertices the plain method gives at most 1/4, 1/10, 0 and 1/2 for t(5,4), t(6,4), t(7,4), t_2(5,4) | `certificates/n6_dual_points/` | (c) |
| Remark n6sharp | six-vertex certificates attaining 1/4, 1/10 and 1/2 | `certificates/n6_dual_points/cert_*_sharp.json` | (j); also (g), first implementation |
| Theorem n7opt | no plain seven-vertex certificate proves t(5,4) > 1336237682928914994292138923/(7*2^89) | `certificates/n7_dual_point/` | (d) |
| Table n7limits | no plain seven-vertex certificate proves t(6,4) > 0.139837689411, t(7,4) > 0.064764464150, t_2(5,4) > 0.553946861054, t_2(6,4) > 0.220748423194, or, without the (5,6) block, t(5,4) > 0.281644555065 (exact fractions in `certificates/n7_dual_points/README.md`) | `certificates/n7_dual_points/` | (i) |
| Theorem sigma | sigma(K_5^(4)) = 31/64: a sharp certificate of value exactly 33/64 for kappa = 2d - gamma, tight exactly on Giraud's ten classes | `certificates/sigma_K5_4/` | (k) |
| Section "Stability for the codegree-squared density": (G7), (G8), (C8), (C9), (C10) and the margin mu_0 of Proposition margin | 2,404 / 10 and 32,981 / 22 labelled / classes; locally Giraud => Giraud for n = 8, 9, 10 (15,636,107 on ten vertices); mu_0 = 176165518826891630107/6280747422216628134215680 | `certificates/sigma_K5_4/` (for mu_0 and the constants) | (l); C10 over all bases and the minimality of mu_0: `--raw stab_c10`, `--raw sigma` |
| Theorem sigmacat (Table sigmacat) | the ten upper bounds of the table: sigma(J_4) <= 9061620987350745659/32281802128991715328 < 0.280703690307, sigma(K_5^<) <= 1375676321453/4123168604160 < 0.333645420191, sigma(K_5^(4)-) <= 6281420455337/30786325577728 < 0.204032807991, and sigma(K_5^=) < 0.384425652179, sigma(K_5^{3-}) < 0.405930172388, sigma(C_5) < 0.252304935013, sigma(K_6^(3)) < 0.742318262433, sigma(K_6^(4)-) < 0.612037993062, sigma(K_6^(4)) < 0.743753064695, sigma(K_7^(4)) < 0.877155180904 (exact fractions in `certificates/sigma_catalogue/README.md`) | `certificates/sigma_catalogue/` (one subdirectory per certificate) | (m); the exhaustive scans: `--raw sigcat` |
| Remark j4limit | every plain seven-vertex certificate for J_4 proves at best sigma(J_4) <= 1 - V, 1 - V > 0.280703275315 (exact dual point); certificates vanishing on the moments of the extremal construction prove at best 16/57 + 316891/14980607589 (exact witness) | `certificates/sigma_catalogue/J4_limit/` | (n); controls `--raw j4_controls` |
| Proposition lotup (Section "Upper bounds for lottery numbers") | (b) C(7,4) l(7,4,5) <= c745 = 3467292853891919529774577/7750000000000000000000000 < 0.447392627; (c) C(7,4) l(7,4,6) <= (107 + 3 c745)/384 < 0.282141089 (type-pattern certificates); (a) C(6,4) l(6,4,5) <= 7/16 is an explicit construction (numerical companion checks only) | `certificates/lottery_upper/` | (p) |
| Section "Remarks on the seven-vertex optimum": Giraud's construction against the t(5,4) certificate; the certificate does not force parity | exact values, see `verifier/EXPECTED.txt` | `certificates/K5_4/` | (e) |
| Section "What six and seven vertices can give": the r = 3 comparison t(4,3) > 0.426923 on five vertices | 60084175574225/2^47 (see the note in EXPECTED.txt) | `verifier/indep_k43_n5.py` | (f) |

Here t_lambda(p,r) is the least limit density of an r-graph in which every p vertices span at least lambda edges,
t = t_1, and pi(K_p^(r)) = 1 - t(p,r) (complement form; see `certificates/README.md` for the formats). In the
catalogue, F_{p,lambda} is the family of p-vertex 4-graphs with at least C(p,4) - lambda + 1 edges and
pi(F_{p,lambda}) = 1 - t_lambda(p,4). Theorem and table names are the labels of the paper's LaTeX source. Each
subdirectory of `certificates/` has a README that states what its files prove and which checker reads them.

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
600 MB, except for `--raw sigcat_J4` (1.5 GB) and `--raw sigcat_K5lt` (0.85 GB). The whole `--fast` run (124 steps)
took 81 minutes with `--threads 2` on that laptop while other heavy jobs were running (62 minutes in the test of the
fourth version, which replaced the certificates of Theorem main (b), (d), (e)), and `--raw sigma,stab_c10` 37 minutes.
The fifth version added parts (m) and (n) to `--fast` (16 steps, 3 minutes with `--sigcat` in its test, peak memory
128 MB), the raw scans of part (m) (run by the independent review, not in the test of this package) and the optional
part (o); the sixth version the seven further certificates of Theorem sigmacat (part (m) now 28 steps, 33 with part
(n); `--sigcat` 8.3 minutes in its test while other jobs were running; 157 steps in `--fast`); the seventh version
part (p), Proposition lotup (4 steps, 1.7 minutes in its test; 161 steps in `--fast`).

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
  part (m) -- for J4, K5lt and K5_4minus on the reviewer's own machine, for the other seven on a temporary cloud
  machine (16 CPUs) -- not in the test of this package, which ran the decoding and the C spot values for three of
  them (K5_4minus, C5, K7_4) and reproduced the reviews' evaluation tables byte for byte; the reviews' second scans
  over the producer's class lists are not repeated (the class lists are not included). Statements that the paper labels as floating-point
  observations are not checked here.

## Contents

    README.md               this file
    LICENSE-CODE            MIT license (code)
    LICENSE-DATA            CC BY 4.0 (data and documentation)
    SHA256SUMS              SHA-256 of every file of the package
    certificates/           the certificates and dual points, one subdirectory per result, each with a README
    verifier/               the independent checkers, verify_all.sh, EXPECTED.txt, CHANGES.txt
    search/                 the producer code (not needed to verify), with its own README and CHANGES.txt

`verifier/CHANGES.txt` and `search/CHANGES.txt` record every change made to the programs when they were copied into
this package (all are comment edits or path adaptations; the SHA-256 of each original is given), and the programs
written for the package.

## AI assistance

Code in this package was written with AI coding agents built on several models (including Claude and GPT); the verifiers were written by separate agent
sessions that did not share code with the search code.

## License

The code (`verifier/`, `search/`) is released under the MIT license, see `LICENSE-CODE`. The data
(`certificates/`) and the documentation are released under the Creative Commons Attribution 4.0 International
license, see `LICENSE-DATA`, which lists the files it covers.

## How to cite

Please cite the paper:

    @misc{MerolaCorrea-turan4,
      author = {Bruno M{\'e}rola Corr{\^e}a},
      title  = {Tur{\'a}n densities of complete 4-graphs via flag algebras on seven vertices},
      year   = {2026},
      note   = {Preprint; arXiv identifier to be added}
    }

and, for the data, this package: https://github.com/brunomerola/turan4-certificates (release v1.0.0; see also
`CITATION.cff`).

## Provenance fields

Some certificate files (`*.cert.json`, `*.verify.json` under `certificates/sigma_catalogue/`) keep the `source` and
path fields written by the program that produced them, which name directories of the author's private research
repository. They are kept unchanged because these files are byte-identical to the files the independent reviews
checked (their SHA-256 digests are recorded in the reviews and in the key lists); the paths play no role in any check.
