# K6_4: t(6,4), complete 4-graph on six vertices

**Proves** Theorem main (b) of the paper:

    t(6,4) >= 614350111275/2^42 > 0.139687042809,  hence  pi(K_6^(4)) < 0.860312957191.

Admissible: every 6 vertices span at least one edge (p = 6, lambda = 1).

**Files**

* `lpcg_p6_full.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table.
* `lpcg_p6_full.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6).
Ranks: (3,5) 12; (5,6) 54, 122, 116, 56, 26, 5 for the types with 0 to 5 edges; the other blocks are zero. The
minimum is attained at K_7^(4) (`argmin_graph` 34359738367 = all 35 edges).

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 6 1` (155 representatives; 7,011,184 classes), `rv_prep.py`, `rv_pycheck.py` at the argmin;
* `--raw K6_4`: `rv_eval.c` over all 162,481,445 admissible one-vertex extensions, `rv_compare.py --raw`,
  `rv_pycheck.py` at the scan's argmin.
