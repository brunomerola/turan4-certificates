# K6_4minus: t_2(6,4), K_6^(4) minus an edge

**Proves** Theorem main (e) of the paper:

    t_2(6,4) >= 202684381843379/923589767331840 > 0.219452823117,  hence  pi(K_6^(4)-) < 0.780547176883.

(923589767331840 = 105 * 2^43.) A 4-graph contains K_6^(4)- (K_6^(4) minus one edge) iff some six vertices span at
least 14 edges. Admissible (complement form): every 6 vertices span at least **two** edges (p = 6, lambda = 2). The
certificate file does not record lambda; the checkers must be run with `RV_LAM=2` and `rv_reps.py 6 2`
(verifier/verify_all.sh does this).

**Files**

* `c6l2_f.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table (361 groups).
* `c6l2_f.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites; it is larger than the certificate
of the preceding phase A, which is not included.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 0 to 5 edges (1013, 1023, 1024, 1024, 1024, 1024 admissible flags).
Ranks (rows of A): (1,4) 1; (3,5) 16; (5,6) 77, 214, 280, 133, 77, 51 for the types with 0 to 5 edges; the other
blocks are zero. The minimum is attained at a 12-edge graph with degrees 8, 8, 7, 7, 6, 6, 6, a K_5^(4) on five
vertices plus 7 edges through the remaining pair of vertices (`argmin_graph` 17444171808), and only on its
isomorphism class.

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 6 2` (154 six-vertex representatives; 6,986,573 classes on 7 vertices by Burnside;
  34,244,802,014 labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic,
  admissible), `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions);
* `--raw K6_4minus`: `rv_eval.c` over all 160,859,017 admissible one-vertex extensions (621,687 links rejected), then
  `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. This raw scan was run by the independent review
  R7_B2 on its own machine (3 h 24 min on 2 threads), not in the test of this package. Its output (exact minimum
  4864425164241096/22166154415964160 = the bound, the SUMS checksum, all 361 group maxima equal to the certificate's)
  is the expected output in `verifier/verify_all.sh`.
* The certificate is also the weak-duality partner of the dual point `../n7_dual_points/dual_k6m_W.json` (part (i)).

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-b2 (c3d-highcpu-30)
on 2026-10-06 (cloud run paperb-b2), which resumed the saved LP state of the earlier t_2(6,4) run; the results are
recorded in commit ece523f8 (branch research/turan4-n7) of the private working repository. Independent review R7_B2
(2026-10-06/07): PASS. The review re-derived the counts, re-ran the flag-list checks with lambda = 2, ran the full raw
scan locally with its own inputs, compared it and the cloud machine's raw and class-list scans with the certificate
on every group maximum, and re-evaluated the certificate exactly, with two pure-Python evaluators, at eight graphs.
All expected values above are taken from its outputs. It supersedes the certificate `lpcg_k6m_f.cert.*` of the
earlier versions of this package (bound 964431985683/2^42; see `verifier/CHANGES.txt`, item 11).
