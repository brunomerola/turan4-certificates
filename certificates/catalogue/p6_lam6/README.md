# catalogue/p6_lam6: t_6(6,4), the family F_{6,6}

**Proves** the row (p, lambda) = (6, 6) of Theorem catalogue (Table catalogue) of the paper:

    t_6(6,4) >= 612132591153137/1231453023109120 > 0.497081561103,  hence  pi(F_{6,6}) < 0.502918438897.

F_{6,6} is the family of 6-vertex 4-graphs with at least C(6,4) - 6 + 1 = 10 edges; its minimal members (15 up to
isomorphism) are the copies of K_6^(4) minus 5 edges. Admissible (complement form): every 6 vertices span at least
**6** edges (p = 6, lambda = 6). The certificate file does not record lambda; the checkers must be run with
`RV_LAM=6` and `rv_reps.py 6 6` (verifier/verify_all.sh does this).

**Files**

* `c6l6_f.cert.json`: keys, flag lists, M = 2^22, the exact bound, the producer's group table (230 groups).
* `c6l6_f.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 0 to 5 edges (386, 638, 848, 968, 1013, 1023 admissible flags). Ranks
(rows of A): (3,5) 16; (4,5) 4, 4 (empty and one-edge type); (5,6) 0, 108, 356, 224, 140, 95 for the types with 0 to
5 edges; the other blocks ((1,4), (2,4), (3,4)) are zero. The minimum is attained at a 17-edge graph (`argmin_graph`
34342977567).

**Checked by** `verifier/verify_all.sh` (name `cat_p6_lam6`):

* fast: `rv_reps.py 6 6` (123 six-vertex representatives; 3,263,333 classes on 7 vertices by Burnside; 16,021,163,016
  labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions), and `cat_decimals` (the exact bound
  and the decimals of the paper's table);
* `--raw cat_p6_lam6`: `rv_eval.c` over all 74,015,441 admissible one-vertex extensions (54,959,407 links rejected),
  then `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. The full raw scan was reproduced by the
  independent review on 2 shared threads (4 h 32 min wall clock), equal to the VM's scan on every deterministic
  field.

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-a2 (c3d-highcpu-30)
on 2026-10-06; the results are recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): PASS with notes (wording only). The review re-derived the counts,
re-ran the flag-list checks with lambda = 6, compared the raw and class-list scans with the certificate on every
group maximum, and re-evaluated the certificate exactly at four graphs with two further evaluators. All expected
values above are taken from its outputs.
