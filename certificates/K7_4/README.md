# K7_4: t(7,4), complete 4-graph on seven vertices

**Proves** Theorem main (c) of the paper:

    t(7,4) >= 7476698908057/115448720916480 > 0.064762076605,  hence  pi(K_7^(4)) < 0.935237923395.

Admissible: every 7 vertices span at least one edge (p = 7, lambda = 1); on 7 vertices this excludes only the empty
graph.

**Files**

* `lpcg_p7_full.cert.json`: keys, flag lists, M = 2^20, the exact bound, the producer's group table.
* `lpcg_p7_full.cert.npz`: integer factor matrices A0..A11 (sha256 = `cert_npz_sha256` in the json).

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6).
Ranks: (3,5) 11; (5,6) 35, 94, 73, 58, 41, 27 for the types with 0 to 5 edges; the other blocks are zero. The
minimum is attained at a 13-edge graph (`argmin_graph` 33940849153).

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 7 1` (156 representatives; 7,013,319 classes = 7,013,320 - 1), `rv_prep.py`, `rv_pycheck.py`
  at the argmin;
* `--raw K7_4`: `rv_eval.c` over all 163,577,855 admissible one-vertex extensions, `rv_compare.py --raw`,
  `rv_pycheck.py` at the scan's argmin.
