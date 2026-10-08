# sigma_catalogue/K6_4minus: the codegree-squared density of K_6^(4)-

**Proves** the row K_6^(4)- of Table sigmacat (Theorem sigmacat) of the paper (Section "Further codegree-squared
densities"; the exact value is in Appendix A):

    sigma(K_6^(4)-) <= 565272027609637/923589767331840 < 0.612037993062      (923589767331840 = 105 * 2^43)

i.e. sigma(K_6^(4)-) < 0.6120380 in the table (rounded up). K_6^(4)- is K_6^(4) minus one edge (14 edges on six
vertices). We found no published bound beyond sigma <= pi (0.852398 with the previous bound for pi, Table sigmacat);
the lower bound > 0.5012201 is a three-part pattern construction (Remark sigma4other of the paper).

**The claim.** Complement form, as for Theorem sigma: a 4-graph G is K_6^(4)--free iff its complement H spans at least
two edges in every 6-set (p = 6, lambda = 2, the admissibility of `../../K6_4minus/`). Objective `sigma4`, kappa(H) =
2 d(H) - gamma(H) = F(H)/420 with F(H) = 24 e(H) - sum over the 35 3-subsets T of [7] of c_T (c_T - 1) (Lemma co2 of
the paper; see `../K5_4minus/README.md`). The certificate consists of PSD matrices Q_k and

    min over all 7-vertex 4-graphs H with >= 2 edges in every 6-set of  kappa(H) - sum_k c_{Q_k}(H)  =  b
                                                                         =  358317739722203/923589767331840

exactly (b = 0.387962006938...). With Lemma co2 this gives sigma(K_6^(4)-) <= 1 - b.

**Blocks and matrices.** Blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6), 12 keys, 6,194 flags (largest Gram matrix
1024); Q_k = A_k^T A_k / M^2 with M = 2^21 and integer A_k (max |A| 6,206,619; 64-bit Gram entries suffice).

The minimum is attained by 36 labelled one-vertex extensions (the producer's argmin is 2522585087); the next value is
larger by 1.79e-8.

**Files**

* `sigma4_c4_K6m_f.sigma4.cert.npz`: integer arrays `A0` ... `A11` (int64).
* `sigma4_c4_K6m_f.sigma4.cert.json`: blocks, `M`, `M_bits`, `tau` = [0, 0] (no stationarity term), the key list
  without the flag lists (`s`, `m`, `type_index`, `sigma`, `n_flags`), the exact `bound`, `argmin_graph`, the a-priori
  bound `z_bound_apriori`, the producer's table `groups` ([F, max Z, argmax, count] for every value F of the objective
  numerator over its class list), `raw_extensions_scanned` = 160,859,017 and `cert_npz_sha256`; the other fields are
  records of the producer's run.
* `sigma4_c4_K6m_f.keys.json`: the full key list (types, |Aut|, factor rows and every flag list, colex masks), the
  formula, the predicate and the sha256 of the three certificate files, in the format of `../README.md`. Written for
  the package because the certificate json does not list the flags (finding M1 of the independent review); compared
  with the reviewer's own types and flags in the review's addendum A1 (all equal), which resolved M1.
* `sigma4_c4_K6m_f.sigma4.verify.json`: the producer's own sampling check; a record, no role in the proof.

All four files are byte-identical to the files the independent review read (sha256 of the certificate files in
`verifier/s8_review_hashes.txt`, of the key list in the review's addendum A1).

**Checked by** `verifier/verify_all.sh` (name K6_4minus; `verifier/EXPECTED.txt`, part (m)), with the programs of the
independent review R7_SIG8:

* fast: `sigcat_keys_indep7` (`s8_keys_cmp.py`, the review's comparison of the key list with its own types, flags and
  |Aut|, of the predicate, the format string and the hashes) and `sigcat_keys` (the producer's `sk_check_keys.py`);
  `sigcat8_predicates` (`s8_predicates.py`: the copy-list predicate against an independent definitional test);
  `sigcat_pycheck_K6_4minus` and `sigcat_pyvals_K6_4minus` (`s8_pycheck.py`, pure-Python values from the definition at
  16 graphs, equal to the review's, `../../../verifier/sigcat_spot_values.txt`; 3 at b); `sigcat_decimals`;
* `--raw sigcat_K6_4minus`: `s8_prepvm.py` (own types and flags equal to the json's key list, exact Gram matrices
  written key by key, own 154 six-vertex representatives; the evaluation table has the review's sha256), `s8_eval` at
  the spot-check graphs (equal to the Python values, `s8_cmp.py`), the exhaustive scan over all 161,480,704 one-vertex
  extensions of the representatives (160,859,017 admissible; completeness by heredity), and `s8_analyze.py`: minimum
  exactly b, attained 36 times, all 149 group maxima equal the producer's table, 34,244,802,014 labelled admissible
  7-vertex graphs (the counts of `../../K6_4minus/`, which has the same admissibility). In the review the scan took
  797 s on 15 threads of a 16-CPU cloud machine (peak memory 54 MB); on 2 threads it should take about 2.2 hours
  (estimate from the review's rates, not measured). The decoding takes a few seconds and less than 0.12 GB.

**Provenance.** Phase F (frozen cuts, suffix `_f`) certificate of the cloud run paperb-a5 (2026-10-07) (results in
commit 8a0657f3, branch research/turan4-n7 of the private working repository); key list from commit da3e33b0.
Independent review R7_SIG8 (2026-10-08): PASS-with-notes (M1, resolved by A1; L1, decimals of the brief only) (the
review's own decoding, its own evaluator on a cloud machine for the raw scan; no class list exists for this problem,
so the raw scan is the only one, pure-Python spot checks and the exact law of an extremal construction where one
exists); addendum A1: the key list checked, M1 resolved. Note (finding L3): the frozen-cut LP 'converged', which need
not be the seven-vertex SDP optimum.
