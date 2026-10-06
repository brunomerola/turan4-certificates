# n6_dual_points: what six vertices can give

**Proves** Theorem n6 of the paper: no plain flag-algebra certificate on six vertices proves t(5,4) > 1/4,
t(6,4) > 1/10, t(7,4) > 0 or t_2(5,4) > 1/2.

**Dual points** (field `dual` of each file: `value`, and `support`, a list of `y`, `edges`, `graph_edges` with the
edges of a 4-graph on {0, ..., 5}):

| file | problem | sum y_H d(H) |
|---|---|---|
| `n6_control_p5.json` | t(5,4) (p = 5, lambda = 1) | 1/4 |
| `n6_control_p6.json` | t(6,4) (p = 6, lambda = 1) | 1/10 |
| `n6_control_p7.json` | t(7,4) (p = 7, lambda = 1); support = the empty graph | 0 |
| `n6_control_p5_l2.json` | t_2(5,4) (p = 5, lambda = 2) | 1/2 |

The other fields of these files (`sdp`, `sharp`, `blocks`, ...) are records of the producer's run.

**Checked by** `verifier/rv_dual6.py` (fast; part (c) of `verify_all.sh`): y > 0, sum y = 1, every support graph
admissible, the value recomputed exactly, and all ten moment matrices of the pairs (s, m) with 2m - s <= 6
(types and flags enumerated from the definitions) symmetric and PSD by an exact LDL^T test.

**Sharp six-vertex certificates** (Remark n6sharp): `cert_turan_r4_p5_N6_all_sharp.json` (value 1/4),
`cert_turan_r4_p6_N6_all_sharp.json` (1/10) and `cert_turan_r4_p5_l2_N6_all_sharp.json` (1/2), in the format of
`../fivegraphs/` with r = 4, N = 6 (`allowed` = null means at least one edge, [2, 3, 4, 5] means at least two). With
the dual points they show that these values are the exact six-vertex optima.

They are checked by `verifier/rv_sharp6.py` (fast; part (j) of `verify_all.sh`), an independent checker written from
the paper's definitions that shares no code with the search code or with the other checkers: every Q is symmetric
and PSD (an exact LDL^T decomposition over the rationals, multiplied back and compared with Q); every listed flag is
a sigma-flag and the listed flags are pairwise non-isomorphic; all 2^15 labelled 4-graphs on six vertices are
enumerated, the admissible ones (27,449 / 32,767 / 12,068) split into isomorphism classes by orbits over all 720
permutations (122 / 155 / 62 classes); and the value e(H)/15 - sum c_Q(H) is computed exactly from the definition for
every class (and again after a random relabelling). Its minimum equals the claimed bound, and it is attained by 3, 4
and 4 classes, the sizes of the supports of the corresponding dual points. A negative control (`sharp6_negative`)
checks that a changed bound, a non-PSD Q and a doubled Q are rejected. When the paper was written they had been
checked by the first implementation only; `verify_all.sh --producer-n6` still runs that checker
(`search/flagalg/verify_cert.py`).
