# sigma_catalogue/C5: the codegree-squared density of C_5

**Proves** the row C_5 of Table sigmacat (Theorem sigmacat) of the paper (Section "Further codegree-squared
densities"; the exact value is in Appendix A):

    sigma(C_5) <= 2167282889/8589934592 < 0.252304935013      (8589934592 = 2^33)

i.e. sigma(C_5) < 0.2523050 in the table (rounded up). C_5 is the tight 5-cycle (1-based: 123 234 345 145 125). The
previous bound was 0.25310457 (Balogh, Clemen and Lidicky); the lower bound 0.2519408 comes from their iterated
construction G_it (see the note below); the gap shrinks from 1.16e-3 to 3.64e-4.

**The claim.** Admissible graphs: C_5-free 3-graphs on 7 vertices (no copy, not necessarily induced). Objective
`cosig3`, f(H) = 1 - gamma(H) = F(H)/210 with F(H) = 210 - sum over the 35 4-subsets Q of [7] of C(e_H(Q), 2), as for
J_4 (see `../README.md`). The certificate consists of PSD matrices Q_k and

    min over all C_5-free 7-vertex 3-graphs H of  f(H) - sum_k c_{Q_k}(H)  =  b  =  6422651703/8589934592

exactly (b = 0.747695064987...). Hence sigma(C_5) <= 1 - b, as in the proof of Theorem sigmacat.

**Blocks and matrices.** The nine blocks (0,3), (1,3), (2,3), (1,4), (2,4), (3,4), (3,5), (4,5), (5,6), 40 keys,
13,879 flags (largest Gram matrix 1024); Q_k = A_k^T A_k / M^2 with M = 2^18 and integer A_k (max |A| 5,235,778; the
Gram entries fit in 64 bits, max diagonal below 2^47, as the checker asserts).

The minimum is attained by a single labelled one-vertex extension (the empty graph, mask 0, also the producer's
argmin); the next value is larger by 2.79e-7.

**Files**

* `cosig3_s3_C5_f.cosig3.cert.npz`: integer arrays `A0` ... `A39` (int64).
* `cosig3_s3_C5_f.cosig3.cert.json`: blocks, `M`, `M_bits`, `tau` = [0, 0] (no stationarity term), the key list
  without the flag lists (`s`, `m`, `type_index`, `sigma`, `n_flags`), the exact `bound`, `argmin_graph`, the a-priori
  bound `z_bound_apriori`, the producer's table `groups` ([F, max Z, argmax, count] for every value F of the objective
  numerator over its class list), `raw_extensions_scanned` = 5,635,379 and `cert_npz_sha256`; the other fields are
  records of the producer's run.
* `cosig3_s3_C5_f.keys.json`: the full key list (types, |Aut|, factor rows and every flag list, colex masks), the
  formula, the predicate and the sha256 of the three certificate files, in the format of `../README.md`. Written for
  the package because the certificate json does not list the flags (finding M1 of the independent review); compared
  with the reviewer's own types and flags in the review's addendum A1 (all equal), which resolved M1.
* `cosig3_s3_C5_f.cosig3.verify.json`: the producer's own sampling check; a record, no role in the proof.

All four files are byte-identical to the files the independent review read (sha256 of the certificate files in
`verifier/s8_review_hashes.txt`, of the key list in the review's addendum A1).

**Checked by** `verifier/verify_all.sh` (name C5; `verifier/EXPECTED.txt`, part (m)), with the programs of the
independent review R7_SIG8:

* fast: `sigcat_keys_indep7` (`s8_keys_cmp.py`, the review's comparison of the key list with its own types, flags and
  |Aut|, of the predicate, the format string and the hashes) and `sigcat_keys` (the producer's `sk_check_keys.py`);
  `sigcat8_predicates` (`s8_predicates.py`: the copy-list predicate against an independent definitional test);
  `sigcat_pycheck_C5` and `sigcat_pyvals_C5` (`s8_pycheck.py`, pure-Python values from the definition at 14 graphs,
  equal to the review's, `../../../verifier/sigcat_spot_values.txt`; 1 at b); `sigcat_decimals`;
* `--raw sigcat_C5`: `s8_prepvm.py` (own types and flags equal to the json's key list, exact Gram matrices written key
  by key, own 835 six-vertex representatives; the evaluation table has the review's sha256), `s8_eval` at the
  spot-check graphs (equal to the Python values, `s8_cmp.py`), the exhaustive scan over all 27,361,280 one-vertex
  extensions of the representatives (5,635,379 admissible; completeness by heredity), and `s8_analyze.py`: minimum
  exactly b, attained once, all 76 group maxima equal the producer's table, 2,342,442,087 labelled admissible 7-vertex
  graphs. In the review the scan took 20 s on 15 threads of a 16-CPU cloud machine (peak memory 68 MB); on 2 threads
  it should take about 3 minutes (estimate from the review's rates, not measured). The decoding takes about 10 to 50
  seconds and less than 0.12 GB.

**Provenance.** Phase F (frozen cuts, suffix `_f`) certificate of the cloud run paperb-a5 (2026-10-07) (results in
commit 8a0657f3, branch research/turan4-n7 of the private working repository); key list from commit da3e33b0.
Independent review R7_SIG8 (2026-10-08): PASS-with-notes (M1, resolved by A1) (the review's own decoding, its own
evaluator on a cloud machine for the raw scan, and a second scan over the producer's class list of 497,783 classes
with the same groups and class counts, pure-Python spot checks and the exact law of an extremal construction where one
exists); addendum A1: the key list checked, M1 resolved. Note (review finding L2): in Balogh, Clemen and Lidicky's
description of G_it (Section 3.17 of arXiv v3) the part sizes and the part of the iteration are interchanged; their
displayed formula and the value 0.25194 fit the reading |B| = bn with the edges xyz, x, y in B, z in A, iterated
inside A, which is C_5-free (the paper states this in the caption of Table sigmacat). The value is not known to be the
seven-vertex optimum: the frozen-cut LP 'converged', which does not show that seven vertices cannot do better (finding
L3).
