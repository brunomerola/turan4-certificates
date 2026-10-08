# sigma_catalogue/K5_3minus: the codegree-squared density of K_5^{3-}

**Proves** the row K_5^{3-} of Table sigmacat (Theorem sigmacat) of the paper (Section "Further codegree-squared
densities"; the exact value is in Appendix A):

    sigma(K_5^{3-}) <= 781068653059/1924145348608 < 0.405930172388      (1924145348608 = 7 * 2^38)

i.e. sigma(K_5^{3-}) < 0.4059302 in the table (rounded up). K_5^{3-} is K_5^(3) minus one triple (1-based: 124 125 134
135 145 234 235 245 345; missing 123). The previous bound was 0.41961743 (Balogh, Clemen and Lidicky); the lower bound
2909/8127 = 0.3579427 is the construction for K_5^= (also K_5^{3-}-free); the gap shrinks from 6.17e-2 to 4.80e-2.

**The claim.** Admissible graphs: K_5^{3-}-free 3-graphs on 7 vertices (no copy, not necessarily induced). Objective
`cosig3`, f(H) = 1 - gamma(H) = F(H)/210 with F(H) = 210 - sum over the 35 4-subsets Q of [7] of C(e_H(Q), 2), as for
J_4 (see `../README.md`). The certificate consists of PSD matrices Q_k and

    min over all K_5^{3-}-free 7-vertex 3-graphs H of  f(H) - sum_k c_{Q_k}(H)  =  b  =  1143076695549/1924145348608

exactly (b = 0.594069827612...). Hence sigma(K_5^{3-}) <= 1 - b, as in the proof of Theorem sigmacat.

**Blocks and matrices.** The nine blocks (0,3), (1,3), (2,3), (1,4), (2,4), (3,4), (3,5), (4,5), (5,6), 46 keys,
32,089 flags (largest Gram matrix 1024); Q_k = A_k^T A_k / M^2 with M = 2^19 and integer A_k (max |A| 6,672,660; the
Gram entries fit in 64 bits, max diagonal below 2^47, as the checker asserts).

The minimum is attained by 2 labelled one-vertex extensions (the producer's argmin is 242467); the next value is
larger by 4.92e-8.

**Files**

* `cosig3_s3_K5m_f.cosig3.cert.npz`: integer arrays `A0` ... `A45` (int64).
* `cosig3_s3_K5m_f.cosig3.cert.json`: blocks, `M`, `M_bits`, `tau` = [0, 0] (no stationarity term), the key list
  without the flag lists (`s`, `m`, `type_index`, `sigma`, `n_flags`), the exact `bound`, `argmin_graph`, the a-priori
  bound `z_bound_apriori`, the producer's table `groups` ([F, max Z, argmax, count] for every value F of the objective
  numerator over its class list), `raw_extensions_scanned` = 57,357,848 and `cert_npz_sha256`; the other fields are
  records of the producer's run.
* `cosig3_s3_K5m_f.keys.json`: the full key list (types, |Aut|, factor rows and every flag list, colex masks), the
  formula, the predicate and the sha256 of the three certificate files, in the format of `../README.md`. Written for
  the package because the certificate json does not list the flags (finding M1 of the independent review); compared
  with the reviewer's own types and flags in the review's addendum A1 (all equal), which resolved M1.
* `cosig3_s3_K5m_f.cosig3.verify.json`: the producer's own sampling check; a record, no role in the proof.

All four files are byte-identical to the files the independent review read (sha256 of the certificate files in
`verifier/s8_review_hashes.txt`, of the key list in the review's addendum A1).

**Checked by** `verifier/verify_all.sh` (name K5_3minus; `verifier/EXPECTED.txt`, part (m)), with the programs of the
independent review R7_SIG8:

* fast: `sigcat_keys_indep7` (`s8_keys_cmp.py`, the review's comparison of the key list with its own types, flags and
  |Aut|, of the predicate, the format string and the hashes) and `sigcat_keys` (the producer's `sk_check_keys.py`);
  `sigcat8_predicates` (`s8_predicates.py`: the copy-list predicate against an independent definitional test);
  `sigcat_pycheck_K5_3minus` and `sigcat_pyvals_K5_3minus` (`s8_pycheck.py`, pure-Python values from the definition at
  16 graphs, equal to the review's, `../../../verifier/sigcat_spot_values.txt`; 3 at b); `sigcat_decimals`;
* `--raw sigcat_K5_3minus`: `s8_prepvm.py` (own types and flags equal to the json's key list, exact Gram matrices
  written key by key, own 1,970 six-vertex representatives; the evaluation table has the review's sha256), `s8_eval`
  at the spot-check graphs (equal to the Python values, `s8_cmp.py`), the exhaustive scan over all 64,552,960
  one-vertex extensions of the representatives (57,357,848 admissible; completeness by heredity), and `s8_analyze.py`:
  minimum exactly b, attained 2 times, all 118 group maxima equal the producer's table, 28,854,605,654 labelled
  admissible 7-vertex graphs. In the review the scan took 194 s on 15 threads of a 16-CPU cloud machine (peak memory
  247 MB); on 2 threads it should take about 30 minutes (estimate from the review's rates, not measured). The decoding
  takes about 50 seconds and less than 0.12 GB.

**Provenance.** Phase F (frozen cuts, suffix `_f`) certificate of the cloud run paperb-a6 (2026-10-07/08) (results in
commit 505cc10d, branch research/turan4-n7 of the private working repository); key list from commit da3e33b0.
Independent review R7_SIG8 (2026-10-08): PASS-with-notes (M1, resolved by A1; L1, decimals of the brief only) (the
review's own decoding, its own evaluator on a cloud machine for the raw scan, and a second scan over the producer's
class list of 5,868,713 classes with the same groups and class counts, pure-Python spot checks and the exact law of an
extremal construction where one exists); addendum A1: the key list checked, M1 resolved. Note (finding L3): the
frozen-cut LP 'converged', which need not be the seven-vertex SDP optimum; the value is a valid bound, not a limit of
the method.
