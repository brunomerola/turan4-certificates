# K5_4: t(5,4), complete 4-graph on five vertices

**Proves** Theorem main (a) of the paper:

    t(5,4) >= 35604499940047/115448720916480 > 0.308400990997,  hence  pi(K_5^(4)) < 0.691599009003.

Admissible: every 5 vertices span at least one edge (p = 5, lambda = 1).

**Files**

* `lpcg_p5_conv_full.cert.json`: keys, flag lists, M = 2^20, the exact bound, the producer's group table.
* `lpcg_p5_conv_full.cert.npz`: integer factor matrices A0..A10 (sha256 = `cert_npz_sha256` in the json).

**Certificate.** Plain flag-algebra certificate on N = 7 vertices with the blocks (s, m) = (1,4), (2,4), (3,4),
(3,5), (4,5), (5,6). Ranks (rows of A): (3,5) 11; (4,5) 1, 1; (5,6) 103, 119, 91, 33, 35 for the types with 1 to 5
edges; the other blocks are zero. The minimum is attained at a 15-edge graph (`argmin_graph` 27303369217).

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 5 1` (122 six-vertex representatives; 3,908,438 classes on 7 vertices by Burnside),
  `rv_prep.py` (flag lists decoded, complete, non-isomorphic, admissible), `rv_pycheck.py` at the certificate's
  argmin (value = bound, exact fractions);
* `--raw K5_4`: `rv_eval.c` over all 86,952,880 admissible one-vertex extensions (the exact minimum), then
  `rv_compare.py --raw` and `rv_pycheck.py` at the scan's argmin.

The same certificate is also the input of `rv1_dual.py` (weak-duality consistency with the dual point of
`n7_dual_point/`), `rv1_giraud.py` (Giraud's construction against the certificate) and `rv1_claim3.py` (the
certificate does not force parity), which check remarks of the paper.
