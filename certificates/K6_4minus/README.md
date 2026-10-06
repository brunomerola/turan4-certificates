# K6_4minus: t_2(6,4), K_6^(4) minus an edge

**Proves** Theorem main (e) of the paper:

    t_2(6,4) >= 964431985683/2^42 > 0.219286445299,  hence  pi(K_6^(4)-) < 0.780713554701.

A 4-graph contains K_6^(4)- (K_6^(4) minus one edge) iff some six vertices span at least 14 edges. Admissible
(complement form): every 6 vertices span at least **two** edges (p = 6, lambda = 2). The certificate file does not
record lambda; the checkers must be run with `RV_LAM=2` and `rv_reps.py 6 2` (verifier/verify_all.sh does this).

**Files**

* `lpcg_k6m_f.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table.
* `lpcg_k6m_f.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6).
Ranks: (1,4) 1; (3,5) 15; (5,6) 80, 226, 230, 121, 80, 43 for the types with 0 to 5 edges; the other blocks are zero.
The minimum is attained at K_6^(4) plus an isolated vertex (`argmin_graph` 31837153552).

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 6 2` (154 representatives; 6,986,573 classes), `rv_prep.py`, `rv_pycheck.py` at the argmin;
* `--raw K6_4minus`: `rv_eval.c` over all 160,859,017 admissible one-vertex extensions, `rv_compare.py --raw`,
  `rv_pycheck.py` at the scan's argmin.
