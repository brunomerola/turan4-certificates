# catalogue/p7_lam3: t_3(7,4), the family F_{7,3}

**Proves** the row (p, lambda) = (7, 3) of Theorem catalogue (Table catalogue) of the paper:

    t_3(7,4) >= 4043518541347/30786325577728 > 0.131341381781,  hence  pi(F_{7,3}) < 0.868658618219.

F_{7,3} is the family of 7-vertex 4-graphs with at least C(7,4) - 3 + 1 = 33 edges; its minimal members (3 up to
isomorphism) are the copies of K_7^(4) minus 2 edges. Admissible (complement form): every 7 vertices span at least
**3** edges (p = 7, lambda = 3). The certificate file does not record lambda; the checkers must be run with
`RV_LAM=3` and `rv_reps.py 7 3` (verifier/verify_all.sh does this).

For p = 7 the only 7-subset of a 7-vertex graph is its whole vertex set, so a 7-vertex graph is admissible iff it has
at least 3 edges, and N = 7 is the first level that sees the constraint. The six-vertex condition is vacuous, so the
representatives are all 156 4-graphs on six vertices.

**Files**

* `c7l3_f.cert.json`: keys, flag lists, M = 2^20, the exact bound, the producer's group table (374 groups).
* `c7l3_f.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 0 to 5 edges (1024, 1024, 1024, 1024, 1024, 1024 admissible flags).
Ranks (rows of A): (3,5) 9; (5,6) 68, 142, 134, 88, 41, 23 for the types with 0 to 5 edges; the other blocks ((1,4),
(2,4), (3,4), (4,5)) are zero. The minimum is attained at a 3-edge graph (`argmin_graph` 30064771072).

**Checked by** `verifier/verify_all.sh` (name `cat_p7_lam3`):

* fast: `rv_reps.py 7 3` (156 six-vertex representatives; 7,013,315 classes on 7 vertices by Burnside; 34,359,737,737
  labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions), and `cat_decimals` (the exact bound
  and the decimals of the paper's table);
* `--raw cat_p7_lam3`: `rv_eval.c` over all 163,577,622 admissible one-vertex extensions (234 links rejected), then
  `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. The full raw scan was run on the cloud VM with the
  same checker programs (`rv_eval.c`, `rv_compare.py`, 28 threads) and checked by the independent review, not
  repeated locally (on 2 threads it takes several hours).

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-a1 (c3d-highcpu-30)
on 2026-10-06; the results are recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): PASS with notes (wording only). The review re-derived the counts,
re-ran the flag-list checks with lambda = 3, compared the raw and class-list scans with the certificate on every
group maximum, and re-evaluated the certificate exactly at four graphs with two further evaluators. All expected
values above are taken from its outputs.
