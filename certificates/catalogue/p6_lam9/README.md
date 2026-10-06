# catalogue/p6_lam9: t_9(6,4), the family F_{6,9}

**Proves** the row (p, lambda) = (6, 9) of Theorem catalogue (Table catalogue) of the paper:

    t_9(6,4) >= 217134104143267/307863255777280 > 0.705293990330,  hence  pi(F_{6,9}) < 0.294706009670.

F_{6,9} is the family of 6-vertex 4-graphs with at least C(6,4) - 9 + 1 = 7 edges; its minimal members (24 up to
isomorphism) are the copies of K_6^(4) minus 8 edges, i.e. the 6-vertex 4-graphs with 7 edges. Admissible (complement
form): every 6 vertices span at least **9** edges (p = 6, lambda = 9). The certificate file does not record lambda;
the checkers must be run with `RV_LAM=9` and `rv_reps.py 6 9` (verifier/verify_all.sh does this).

**Files**

* `c6l9_f.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table (112 groups).
* `c6l9_f.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

This is the frozen-cut (phase F, suffix `_f`) certificate, the one the paper cites.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 0 to 5 edges (11, 56, 176, 386, 638, 848 admissible flags). Ranks (rows
of A): (3,5) 16; (4,5) 2, 0 (empty and one-edge type); (5,6) 0, 0, 101, 251, 258, 100 for the types with 0 to 5
edges; the other blocks ((1,4), (2,4), (3,4)) are zero. The minimum is attained at a 26-edge graph (`argmin_graph`
4294711295).

**Checked by** `verifier/verify_all.sh` (name `cat_p6_lam9`):

* fast: `rv_reps.py 6 9` (54 six-vertex representatives; 122,149 classes on 7 vertices by Burnside; 574,200,468
  labelled admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions), and `cat_decimals` (the exact bound
  and the decimals of the paper's table);
* `--raw cat_p6_lam9`: `rv_eval.c` over all 4,517,328 admissible one-vertex extensions (52,105,776 links rejected),
  then `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. The full raw scan was reproduced by the
  independent review on 2 threads (13 min), equal to the VM's scan on every deterministic field.

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-a3 (n2-highcpu-16) on
2026-10-06; the results are recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): PASS with notes (wording only). The review re-derived the counts,
re-ran the flag-list checks with lambda = 9, compared the raw and class-list scans with the certificate on every
group maximum, and re-evaluated the certificate exactly at four graphs with two further evaluators. All expected
values above are taken from its outputs.
