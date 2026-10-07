# K6_4: t(6,4), complete 4-graph on six vertices

**Proves** Theorem main (b) of the paper:

    t(6,4) >= 9224492822869/65970697666560 > 0.139827122482,  hence  pi(K_6^(4)) < 0.860172877518.

(65970697666560 = 15 * 2^42.) Admissible: every 6 vertices span at least one edge (p = 6, lambda = 1).

**Files**

* `c6l1.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table (372 groups).
* `c6l1.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

This is the certificate of the time-limited cut and column generation (phase A, no suffix), the one the paper cites:
it is larger, by 2.5e-8, than the certificate of the following frozen-cut phase (phase F, `c6l1_f`, bound
4304762543451/30786325577728 > 0.139827097344), which is not included.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 0 to 5 edges (1023, 1024, 1024, 1024, 1024, 1024 admissible flags).
Ranks (rows of A): (3,5) 12; (5,6) 64, 182, 194, 116, 54, 21 for the types with 0 to 5 edges; the other blocks are
zero. The minimum is attained at the 7-edge 4-graph whose edges are the complements of the seven lines of a Fano
plane (`argmin_graph` 2290237506), and only on its isomorphism class.

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 6 1` (155 six-vertex representatives; 7,011,184 classes on 7 vertices by Burnside;
  34,352,419,335 labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic,
  admissible), `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions);
* `--raw K6_4`: `rv_eval.c` over all 162,481,445 admissible one-vertex extensions (47,835 links rejected), then
  `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. This raw scan was run by the independent review
  R7_B2, not in the test of this package (on 2 threads it takes hours): with its own build of `rv_eval.c` and its own
  evaluation tables and representatives, on a temporary 32-thread cloud machine (323 s), because a crash of the
  reviewing laptop had interrupted its local raw scans. Its output (exact minimum 3099429588483984/22166154415964160
  = the bound, all 372 group maxima equal to the certificate's) is the expected output in `verifier/verify_all.sh`;
  the review checked it against the certificate and against its own local scan of all 7,011,184 classes.
* The certificate is also the weak-duality partner of the dual point `../n7_dual_points/dual_p6_f3.json` (part (i)).

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", here lambda = 1; exact rounding; exact class and raw scans) on the cloud VM paperb-b2
(c3d-highcpu-30) on 2026-10-06 (cloud run paperb-b2), which resumed the saved LP state of the earlier t(6,4) run;
the results are recorded in commit ece523f8 (branch research/turan4-n7) of the private working repository.
Independent review R7_B2 (2026-10-06/07): PASS. The review re-derived the counts, re-ran the flag-list checks,
ran the full raw scan and a class-list scan with its own inputs, and re-evaluated the certificate exactly, with two
pure-Python evaluators, at six graphs (the certificate's argmin, the raw scan's argmin and four others). All expected
values above are taken from its outputs. It supersedes the certificate `lpcg_p6_full.cert.*` of the earlier versions
of this package (bound 614350111275/2^42; see `verifier/CHANGES.txt`, item 11).
