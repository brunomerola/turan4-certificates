# K5_4minus: t_2(5,4), K_5^(4) minus an edge

**Proves** Theorem main (d) of the paper:

    t_2(5,4) >= 48641756816487/87960930222080 > 0.552992751369,  hence  pi(K_5^(4)-) < 0.447007248631.

(87960930222080 = 5 * 2^44.) K_5^(4)- = H^4_4 is the 4-graph on five vertices with four edges. Admissible: every 5
vertices span at least **two** edges (p = 5, lambda = 2). The certificate file does not record lambda; the checkers
must be run with `RV_LAM=2` and `rv_reps.py 5 2` (verifier/verify_all.sh does this).

**Files**

* `c5l2_f.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table (170 groups).
* `c5l2_f.cert.npz`: integer factor matrices A0..A9 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites; it is larger than the certificate
of the preceding phase A, which is not included.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 2 to 5 edges (373, 463, 588, 768 admissible flags). Ranks (rows of A):
(3,5) 9; (5,6) 189, 225, 201, 124 for the types with 2 to 5 edges; the other blocks are zero. The minimum is attained
at a 25-edge graph with degrees 19, 14, 14, 14, 13, 13, 13 (`argmin_graph` 30440134655), and only on its isomorphism
class.

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 5 2` (62 six-vertex representatives; 360,127 classes on 7 vertices by Burnside; 1,741,143,596
  labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions);
* `--raw K5_4minus`: `rv_eval.c` over all 10,058,620 admissible one-vertex extensions (54,953,092 links rejected;
  7 to 11 min on 2 threads), then `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. Run in the test of
  this package (428 s on 2 threads); the expected values (exact minimum 12257722717754724/22166154415964160 = the
  bound, the SUMS checksum, all 170 group maxima equal to the certificate's) are those of the independent review's
  own raw scan, and the test reproduced them.
* The certificate is also the weak-duality partner of the dual point `../n7_dual_points/dual_h44_nt3.json` (part (i)).

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-b2 (c3d-highcpu-30)
on 2026-10-06 (cloud run paperb-b2), which resumed the saved LP state of the earlier t_2(5,4) run; the results are
recorded in commit ece523f8 (branch research/turan4-n7) of the private working repository. Independent review R7_B2
(2026-10-06/07): PASS. The review re-derived the counts, re-ran the flag-list checks with lambda = 2, ran the full raw
scan locally with its own inputs (634 s on 2 threads), compared it and the cloud machine's raw and class-list scans
with the certificate on every group maximum, and re-evaluated the certificate exactly, with two pure-Python
evaluators, at eight graphs. All expected values above are taken from its outputs. It supersedes the certificate
`lpcg_h44e.cert.*` of the earlier versions of this package (bound 2430536617277/2^42; see `verifier/CHANGES.txt`,
item 11).
