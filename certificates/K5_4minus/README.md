# K5_4minus: t_2(5,4), K_5^(4) minus an edge

**Proves** Theorem main (d) of the paper:

    t_2(5,4) >= 2430536617277/2^42 > 0.552640043969,  hence  pi(K_5^(4)-) < 0.447359956031.

K_5^(4)- = H^4_4 is the 4-graph on five vertices with four edges. Admissible: every 5 vertices span at least
**two** edges (p = 5, lambda = 2). The certificate file does not record lambda; the checkers must be run with
`RV_LAM=2` and `rv_reps.py 5 2` (verifier/verify_all.sh does this).

**Files**

* `lpcg_h44e.cert.json`: keys, flag lists, M = 2^21, the exact bound, the producer's group table.
* `lpcg_h44e.cert.npz`: integer factor matrices A0..A9 (sha256 = `cert_npz_sha256` in the json).

**Certificate.** Plain flag-algebra certificate on N = 7 vertices, blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6);
the (5,6) types are the 5-vertex graphs with 2 to 5 edges. Ranks: (3,5) 9; (5,6) 175, 199, 176, 113; the other
blocks are zero. The minimum is attained at K_7^(4) (`argmin_graph` 34359738367).

**Checked by** `verifier/verify_all.sh`:

* fast: `rv_reps.py 5 2` (62 representatives; 360,127 classes), `rv_prep.py`, `rv_pycheck.py` at the argmin;
* `--raw K5_4minus`: `rv_eval.c` over all 10,058,620 admissible one-vertex extensions (5 to 30 min), `rv_compare.py
  --raw`, `rv_pycheck.py` at the scan's argmin.
