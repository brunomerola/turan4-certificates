# sigma_K5_4: the codegree-squared density of K_5^(4)

**Proves** Theorem sigma of the paper (Section "The codegree-squared density"): sigma(K_5^(4)) = 31/64, the upper
bound. The lower bound is Giraud's construction (by hand in the paper; its exact seven-vertex law is checked by
`rs_giraud7.py`). Section "Stability for the codegree-squared density" uses this certificate through its margin mu_0
(Proposition margin; see `verifier/EXPECTED.txt`, part (l)).

**The claim.** In complement form (H is the complement of a K_5^(4)-free 4-graph, so every 5 vertices of H span an
edge: p = 5, lambda = 1) the objective on a 7-vertex 4-graph H is

    kappa(H) = (24 e(H) - P(H)) / 420,   P(H) = sum over 3-sets T of c_T (c_T - 1),   c_T = #{x not in T : T + x in H},

that is kappa = 2d - gamma (Lemma co2 of the paper; 420 = C(7,3) * 4 * 3 ordered triples (T, x, y)). The certificate
consists of PSD matrices Q_k, one per key k (a block (s, m) and a type sigma), and

    min over all admissible 7-vertex 4-graphs H of  kappa(H) - sum_k c_{Q_k}(H)  =  33/64   exactly,

with c_Q as in Definition cert of the paper. The minimum is attained exactly on the ten isomorphism classes of
7-vertex subgraphs of Giraud's construction (colex masks 2290237506, 17418883036, 19736566032, 26383237633,
26894418433, 29337331976, 31813832752, 31837153552, 34326183967, 34359738367) and nowhere else; every other admissible
graph exceeds it by at least mu_0 = 176165518826891630107/6280747422216628134215680 > 2.8048e-5 (attained at the
11-edge class 26887858440). With Lemma co2 this gives sigma(K_5^(4)) <= 1 - 33/64 = 31/64.

**The matrices.** For every key,

    Q_k = W_k^T W_k / 2^72  +  B_k^T X_k B_k,      X_k = Xq_k / q,

with q = 57 * 2^80, so that D Q_k = q W_k^T W_k + 2^72 B_k^T Xq_k B_k is an integer matrix for D = 2^72 q = 57 * 2^152.
W_k are integer matrices (rows 1, 0, 0, 19, 4, 7, 272, 193, 160, 44, 106 for the keys 0 to 10). The second term is
present only for the keys 6, 8 and 10 (the (5,6) types with 1, 3 and 5 edges), with integer B_k of 20, 5 and 1 rows
and symmetric rational X_k (dimensions 20, 5, 1), which must be checked to be PSD; the first term is PSD by
construction.

**Files**

* `sharp_klp5.cert.npz`: the certificate, byte-identical to the producer's file (sha256 caefdc1c...). Fields: `kb` =
  36, `q` and `D` (decimal strings), `W0` ... `W10` (decimal strings, one row per entry of W_k; an empty array means
  no rows), `B6`, `B8`, `B10` (int64), `corr` (JSON: for keys 6, 8, 10 the integer matrix `X_times_q` = q X_k and
  records of the rounding), `eps` = 2^-20 (the target diagonal of X_k).
* `sharp_klp5.cert.json`: the **key list** that the .npz does not carry: for each key k, in the index order of the
  matrices, `s`, `m`, `type_index`, `sigma` (colex mask on [s]) and `flags` (colex masks on [m], roots 0, ..., s-1),
  `n_flags`; plus the formula above, `kb`, `q`, `D`, the claimed bound and the sha256 of the .npz
  (`cert_npz_sha256`, checked by the step `npz_hashes`). The keys are copied unchanged from
  `../K5_4/lpcg_p5_conv_full.cert.json`: both certificates use the same flag model (blocks (1,4), (2,4), (3,4),
  (3,5), (4,5), (5,6); 2, 2, 2, 23, 15, 16, 809, 854, 904, 960, 1024 flags), exactly as the independent review
  evaluated this certificate. This file was written for the package (review finding F1: the .npz alone does not
  determine the flags).
* `sharp_klp5.json`: the producer's record of the run (byte-identical): its per-F group table over the 113 values of
  F = 24e - P (max of Z and counts), the tight classes, scan times. It is used only for an informational comparison
  (`rs_analyze.py`) and plays no role in the proof.

Validity needs only: every Q_k PSD, the flags of a key pairwise non-isomorphic sigma-flags, and the exact minimum
over all admissible 7-vertex graphs. Tightness on Giraud's classes is not needed for validity.

**Checked by** `verifier/verify_all.sh`, part (k) (the checkers of the independent review of this certificate,
copied with path and label edits only, see `verifier/CHANGES.txt`):

* fast: `sigma_keys` (the key list equals the K5_4 key list; `kb`, `q`, `D` agree with the .npz);
  `sigma_reduction` (`rs_reduction.py`: the identities of Lemma co2 and the exact 7-vertex averaging, by brute force);
  `sigma_giraud7` (`rs_giraud7.py`: Giraud's exact 7-vertex law, 2,404 labelled graphs in the ten classes, E kappa =
  33/64); `sigma_reps` (`rs_reps.py`: own 6-vertex representatives, 122 classes, 86,952,880 raw extensions);
  `sigma_prep` (`rs_prep.py`: the integer Gram matrices rebuilt exactly, X_6, X_8, X_10 symmetric and PSD by exact
  LDL^T, every flag list decoded, root part = sigma, pairwise non-isomorphic, admissible, complete, the types of a block
  pairwise non-isomorphic; writes the evaluation tables); `sigma_pycheck` (`rs_pycheck.py`: pure-Python exact value at
  26 graphs: Giraud's ten classes (exactly 33/64), the margin graph (33/64 + mu_0), five further low classes and ten
  random or relabelled graphs); `sigma_each` (the C evaluator `rs_eval.c` at the same 26 graphs, equal to the Python
  values);
* `--raw sigma`: the exhaustive scan with the symmetrised C evaluator over all 86,952,880 admissible one-vertex
  extensions of the representatives (completeness by heredity alone), split over `--threads` processes, then
  `rs_analyze.py`: minimum exactly 33/64, attained by 105 labelled extensions, all in Giraud's ten classes and
  covering all ten; next value exactly 33/64 + mu_0; the per-F maxima equal the producer's table. About 15 to 20 CPU minutes (10 minutes with `--threads 2`).

**Provenance.** Found by the producer's custom-objective pipeline (LP cut and column generation for the objective
kappa, then an exact "sharp" rounding with a forced kernel and a correction term; not included in `search/`), recorded
in commit eaab195e (branch research/turan4-n7) of the private working repository. Independent review R7_SIGMA
(2026-10-06): PASS-with-notes (notes on presentation and archiving only; F1 is addressed by `sharp_klp5.cert.json`).
The weaker, non-sharp certificate from the cloud run (b = 2267741359833/2^42) is not included.
