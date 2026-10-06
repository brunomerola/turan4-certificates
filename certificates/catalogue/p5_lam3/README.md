# catalogue/p5_lam3: t_3(5,4), the family F_{5,3}

**Proves** the row (p, lambda) = (5, 3) of Theorem catalogue (Table catalogue) of the paper:

    t_3(5,4) >= 3298534811415/4398046511104 > 0.749999983648,  hence  pi(F_{5,3}) < 0.250000016352.

F_{5,3} is the family of 5-vertex 4-graphs with at least C(5,4) - 3 + 1 = 3 edges; its minimal members (1 up to
isomorphism) are the copies of H^4_3 = K_5^(4) minus 2 edges (5 vertices, 3 edges). Admissible (complement form):
every 5 vertices span at least **3** edges (p = 5, lambda = 3). The certificate file does not record lambda; the
checkers must be run with `RV_LAM=3` and `rv_reps.py 5 3` (verifier/verify_all.sh does this).

This case validates the method for lambda >= 3: the exact value is known, t_3(5,4) = 3/4, that is pi(H^4_3) = 1/4
(see the paper). The certified b is below 3/4, as it must be: 3/4 - b = 1.635e-8 (the paper's "up to 1.7e-8").

**Files**

* `c5l3.cert.json`: keys, flag lists, M = 2^22, the exact bound, the producer's group table (69 groups).
* `c5l3.cert.npz`: integer factor matrices A0..A8 (sha256 = `cert_npz_sha256` in the json).

This is the phase-A certificate (no `_f` suffix): for this case the frozen-cut phase F found no better point and
wrote no certificate, so the phase-A certificate is the one the paper cites.

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 3 to 5 edges (88, 145, 253 admissible flags). Ranks (rows of A): (2,4)
1; (4,5) 1, 3 (empty and one-edge type); (5,6) 0, 2, 2 for the types with 3 to 5 edges; the other blocks ((1,4),
(3,4), (3,5)) are zero. The minimum is attained at K_7^(4) (`argmin_graph` 34359738367).

**Checked by** `verifier/verify_all.sh` (name `cat_p5_lam3`):

* fast: `rv_reps.py 5 3` (19 six-vertex representatives; 3,065 classes on 7 vertices by Burnside; 12,420,940 labelled
  admissible 7-vertex graphs), `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible),
  `rv_pycheck.py` at the certificate's argmin (value = bound, exact fractions), and `cat_decimals` (the exact bound
  and the decimals of the paper's table);
* `--raw cat_p5_lam3`: `rv_eval.c` over all 172,932 admissible one-vertex extensions (19,750,012 links rejected),
  then `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin. The full raw scan was reproduced by the
  independent review on 2 threads (17.6 s).

**Provenance.** Found by the search code of `search/flagalg/h44/` (LP cut and column generation for "every p-set
spans at least lambda edges", exact rounding, exact class and raw scans) on the cloud VM paperb-a1 (c3d-highcpu-30)
on 2026-10-06; the results are recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): PASS with notes (the notes concern wording only; one is that this
case has no frozen-cut certificate, see Files). The review re-derived the counts, re-ran the flag-list checks with
lambda = 3, compared the raw and class-list scans with the certificate on every group maximum, and re-evaluated the
certificate exactly at four graphs with two further evaluators. All expected values above are taken from its outputs.
