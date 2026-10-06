# catalogue/p6_lam5: t_5(6,4), the family F_{6,5}

**Proves** the row (p, lambda) = (6, 5) of Theorem catalogue (Table catalogue) of the paper:

    t_5(6,4) >= 64873027783677/153931627888640 > 0.421440536123,  hence  pi(F_{6,5}) < 0.578559463877.

F_{6,5} is the family of 6-vertex 4-graphs with at least C(6,4) - 5 + 1 = 11 edges; its minimal members (9 up to
isomorphism) are the copies of K_6^(4) minus 4 edges. Admissible (complement form): every 6 vertices span at least
**5** edges (p = 6, lambda = 5). The certificate file does not record lambda; the checkers must be run with
`RV_LAM=5` and `rv_reps.py 6 5` (verifier/verify_all.sh does this).

**Files**

* `c6l5_f.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table (272 groups).
* `c6l5_f.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 0 to 5 edges (638, 848, 968, 1013, 1023, 1024 admissible flags). Ranks
(rows of A): (3,5) 16; (4,5) 2, 1 (empty and one-edge type); (5,6) 1, 149, 348, 185, 115, 113 for the types with 0 to
5 edges; the other blocks ((1,4), (2,4), (3,4)) are zero. The minimum is attained at a 25-edge graph (`argmin_graph`
31803041519).

**Checked by** `verifier/verify_all.sh` (name `cat_p6_lam5`):

* fast: `rv_reps.py 6 5` (138 six-vertex representatives; 5,059,432 classes on 7 vertices by Burnside; 24,868,266,169
  labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions), and `cat_decimals` (the exact bound
  and the decimals of the paper's table);
* `--raw cat_p6_lam5`: `rv_eval.c` over all 112,218,914 admissible one-vertex extensions (32,484,574 links rejected),
  then `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. The full raw scan was run on the cloud VM with
  the same checker programs (`rv_eval.c`, `rv_compare.py`, 28 threads) and checked by the independent review, not
  repeated locally (on 2 threads it takes several hours).

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-a2 (c3d-highcpu-30)
on 2026-10-06; the results are recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): PASS with notes (wording only). The review re-derived the counts,
re-ran the flag-list checks with lambda = 5, compared the raw and class-list scans with the certificate on every
group maximum, and re-evaluated the certificate exactly at four graphs with two further evaluators. All expected
values above are taken from its outputs.
