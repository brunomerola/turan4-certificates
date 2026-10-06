# catalogue/p6_lam11: t_11(6,4), the family F_{6,11}

**Proves** the row (p, lambda) = (6, 11) of Theorem catalogue (Table catalogue) of the paper:

    t_11(6,4) >= 24267758577011/28862180229120 > 0.840815156178,  hence  pi(F_{6,11}) < 0.159184843822.

F_{6,11} is the family of 6-vertex 4-graphs with at least C(6,4) - 11 + 1 = 5 edges; its minimal members (15 up to
isomorphism) are the copies of K_6^(4) minus 10 edges, i.e. the 6-vertex 4-graphs with 5 edges. Admissible
(complement form): every 6 vertices span at least **11** edges (p = 6, lambda = 11). The certificate file does not
record lambda; the checkers must be run with `RV_LAM=11` and `rv_reps.py 6 11` (verifier/verify_all.sh does this).

**Files**

* `c6l11_f.cert.json`: keys, flag lists, M = 2^19, the exact bound, the producer's group table (47 groups).
* `c6l11_f.cert.npz`: integer factor matrices A0..A10 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 1 to 5 edges (1, 11, 56, 176, 386 admissible flags). Ranks (rows of A):
(3,5) 10; (5,6) 0, 2, 41, 72, 34 for the types with 1 to 5 edges; the other blocks ((1,4), (2,4), (3,4), (4,5)) are
zero. The minimum is attained at a 27-edge graph (`argmin_graph` 2146992127).

**Checked by** `verifier/verify_all.sh` (name `cat_p6_lam11`):

* fast: `rv_reps.py 6 11` (18 six-vertex representatives; 1,939 classes on 7 vertices by Burnside; 7,009,743 labelled
  admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions), and `cat_decimals` (the exact bound
  and the decimals of the paper's table);
* `--raw cat_p6_lam11`: `rv_eval.c` over all 155,364 admissible one-vertex extensions (18,719,004 links rejected),
  then `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. The full raw scan was reproduced by the
  independent review on 2 threads (16 s), equal to the VM's scan on every deterministic field.

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-a3 (n2-highcpu-16) on
2026-10-06; the results are recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): PASS with notes (wording only). The review re-derived the counts,
re-ran the flag-list checks with lambda = 11, compared the raw and class-list scans with the certificate on every
group maximum, and re-evaluated the certificate exactly at four graphs with two further evaluators. All expected
values above are taken from its outputs.
