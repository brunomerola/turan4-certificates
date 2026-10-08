# Certificates

One subdirectory per result. Every data file here is byte-identical to the file produced by the search code and
re-checked by the independent checkers, except the key lists written for the package from the producer's data
(`sigma_K5_4/sharp_klp5.cert.json` and `sigma_catalogue/*/*.keys.json`; each records the SHA-256 of the certificate
files it describes); `SHA256SUMS` at the package root lists all digests, and each certificate `.json` records the
SHA-256 of its `.npz` in the field `cert_npz_sha256`.

| directory | result in the paper | exact value |
|---|---|---|
| `K5_4/` | Theorem main (a): t(5,4), so pi(K_5^(4)) | 35604499940047/115448720916480 |
| `K6_4/` | Theorem main (b): t(6,4), so pi(K_6^(4)) | 9224492822869/65970697666560 |
| `K7_4/` | Theorem main (c): t(7,4), so pi(K_7^(4)) | 7476698908057/115448720916480 |
| `K5_4minus/` | Theorem main (d): t_2(5,4), so pi(K_5^(4)-) | 48641756816487/87960930222080 |
| `K6_4minus/` | Theorem main (e): t_2(6,4), so pi(K_6^(4)-) | 202684381843379/923589767331840 |
| `fivegraphs/` | Theorem five: t(6,5), t_3(6,5), t(7,5) | see `fivegraphs/README.md` |
| `catalogue/` | Theorem catalogue (Table catalogue): t_lambda(p,4) for ten (p, lambda), one subdirectory each | see `catalogue/README.md` |
| `n6_dual_points/` | Theorem n6 (six-vertex barrier) and Remark n6sharp | 1/4, 1/10, 0, 1/2 |
| `n7_dual_point/` | Theorem n7opt (seven-vertex optimum for t(5,4)) | 1336237682928914994292138923/(7*2^89) |
| `n7_dual_points/` | Table n7limits: seven-vertex limits for t(6,4), t(7,4), t_2(5,4), t_2(6,4) and for t(5,4) without the (5,6) block | see `n7_dual_points/README.md` |
| `sigma_K5_4/` | Theorem sigma: sigma(K_5^(4)) = 31/64 (sharp certificate, value 33/64 for kappa = 2d - gamma); also the margin mu_0 of Proposition margin (stability section) | 33/64 |
| `sigma_catalogue/` | Theorem sigmacat: sigma(J_4), sigma(K_5^<), sigma(K_5^(4)-) (one subdirectory each); Remark j4limit: the dual point and the witness for J_4 (`J4_limit/`) | see `sigma_catalogue/README.md` |

Theorem, table and remark names are the LaTeX labels of the paper (`thm:main`, `thm:five`, `thm:catalogue`,
`tab:catalogue`, `thm:n6`, `thm:n7opt`, `tab:n7limits`, `rem:n6sharp`, `thm:sigma`, `prop:margin`, `thm:sigmacat`,
`rem:j4limit`);
`115448720916480 = 105 * 2^40`, `65970697666560 = 15 * 2^42`, `87960930222080 = 5 * 2^44`,
`923589767331840 = 105 * 2^43`. The computer facts of the stability section ((G7), (G8), (C8), (C9), (C10)) need no
data file: they are recomputed from the definitions by part (l) of `verifier/verify_all.sh`.

## Conventions

* `t_lambda(p, r)` is the least limit density of an r-graph in which every p-set of vertices spans at least
  `lambda` edges; `t(p, r) = t_1(p, r)`. Then `pi(K_p^(r)) = 1 - t(p, r)`, `pi(K_5^(4)-) = 1 - t_2(5,4)` and
  `pi(K_6^(4)-) = 1 - t_2(6,4)` (complement form). A graph satisfying the condition is called admissible.
* Graphs are bit masks. Bit i of the mask of an r-graph on {0, ..., n-1} is the i-th r-subset in **colex** order;
  the colex rank of {a < b < c < d} is C(a,1) + C(b,2) + C(c,3) + C(d,4) (and likewise for r = 5).
* A sigma-flag on m vertices has its roots on 0, ..., s-1; its mask is a graph on {0, ..., m-1} whose restriction to
  {0, ..., s-1} is the type sigma.

## 4-graph certificates (`*.cert.json` + `*.cert.npz`)

* `.json`: `p`; `keys`, one entry per (block (s, m), type): `s`, `m`, `type_index`, `sigma` (colex mask on [s]),
  `n_flags`, `flags` (colex masks on [m]); `M` (scale, 2^20 or 2^21; 2^19 to 2^22 in `catalogue/`);
  `tau1_num = tau2_num = 0` (no stationarity term); `bound` (the exact value b); `argmin_graph` (a minimiser, colex
  mask on [7]); `groups`, the producer's table `[e, cdd, max Z, argmax, count]` over its class list, grouped by edge
  count e and by cdd, the number of pairs of edges meeting in exactly one vertex. The other fields (`float_b_*`, scan times, ...) are records of the
  producer's run and play no role in the proof.
* `.npz`: one int64 array `A{k}` per key k, of shape (rank, n_flags); the certificate matrix is
  `Q_k = A_k^T A_k / M^2`, positive semidefinite by construction.
* For a 7-vertex 4-graph H, `Z(H) = 5040 M^2 sum_k c_{Q_k}(H)` is an integer and
  `d(H) - sum_k c_{Q_k}(H) = (144 M^2 e(H) - Z(H)) / (5040 M^2)`; the value b is its minimum over all admissible H.
* The files do not record lambda: it is 1 for `K5_4`, `K6_4`, `K7_4`, 2 for `K5_4minus`, `K6_4minus`, and L for
  `catalogue/pP_lamL`. The checkers receive it through `RV_LAM` and as an argument of `rv_reps.py`
  (verifier/verify_all.sh does this).

## The sigma certificate (`sigma_K5_4/`)

It has a different format: rational matrices Q_k = W_k^T W_k / 2^72 + B_k^T X_k B_k in `sharp_klp5.cert.npz`, with the
key list (types and flags, same conventions as above) in `sharp_klp5.cert.json`; see `sigma_K5_4/README.md`.

## The certificates of Theorem sigmacat (`sigma_catalogue/`)

They have the matrix format of the 4-graph certificates (Q_k = A_k^T A_k / M^2, one int64 array `A<k>` per key in
`<prefix>.<objective>.cert.npz`), but for 3-graphs (J_4, K_5^<; nine blocks) or for the 4-graph K_5^(4)- in complement
form, and with the objectives 1 - gamma (`cosig3`, F/210) or kappa = 2d - gamma (`sigma4`, F/420) instead of the edge
density. The json records the forbidden graph through `problem` (s3_J4, s3_K5lt, c4_K5m); the full key lists with the
flag masks are in `<prefix>.keys.json`. The dual point and the witness of Remark j4limit are rational vectors on
7-vertex 3-graphs (colex masks). See `sigma_catalogue/README.md`.
