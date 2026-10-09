# Results and certificates

[Overview](README.md) · [Verification guide](VERIFICATION.md)

This index maps the paper's results to certificate files and checks in `verifier/verify_all.sh`.
Theorem, proposition, remark and table names below are the labels in the paper's LaTeX source.
Each certificate directory gives the exact values, file formats and provenance.

Here t_lambda(p,r) is the least limit density of an r-graph in which every p vertices span at least lambda edges,
t = t_1, and pi(K_p^(r)) = 1 - t(p,r). In the catalogue, F_{p,lambda} is the family of p-vertex 4-graphs with
at least C(p,4) - lambda + 1 edges, so pi(F_{p,lambda}) = 1 - t_lambda(p,4).
See [conventions](certificates/README.md#conventions).

Decimal summaries use six places, rounded up for upper bounds and down for lower bounds.
Twelve places are retained for the seven-vertex method limits and the J_4 comparison, where small gaps matter.
Exact fractions and the more detailed certificate tables remain the reference values.

## Turán densities of 4-graphs

| statement in the paper | bound | files | checks |
|---|---|---|---|
| Theorem main (a) | t(5,4) >= 35604499940047/115448720916480 > 0.308400, so pi(K_5^(4)) < 0.691600 | [certificates/K5_4/](certificates/K5_4/) | (a) |
| Theorem main (b) | t(6,4) >= 9224492822869/65970697666560 > 0.139827, so pi(K_6^(4)) < 0.860173 | [certificates/K6_4/](certificates/K6_4/) | (a) |
| Theorem main (c) | t(7,4) >= 7476698908057/115448720916480 > 0.064762, so pi(K_7^(4)) < 0.935238 | [certificates/K7_4/](certificates/K7_4/) | (a) |
| Theorem main (d) | t_2(5,4) >= 48641756816487/87960930222080 > 0.552992, so pi(K_5^(4)-) < 0.447008 | [certificates/K5_4minus/](certificates/K5_4minus/) | (a) |
| Theorem main (e) | t_2(6,4) >= 202684381843379/923589767331840 > 0.219452, so pi(K_6^(4)-) < 0.780548 | [certificates/K6_4minus/](certificates/K6_4minus/) | (a) |
| Theorem catalogue (Table catalogue) | t_lambda(p,4) for (p, lambda) = (5,3), (6,3), (6,4), (6,5), (6,6), (6,9), (6,11), (7,2), (7,3), (7,4); exact fractions in `certificates/catalogue/README.md` | [certificates/catalogue/](certificates/catalogue/) | (h) |

## Exact codegree-squared density and stability

| statement in the paper | bound | files | checks |
|---|---|---|---|
| Theorem sigma | sigma(K_5^(4)) = 31/64: a sharp certificate of value exactly 33/64 for kappa = 2d - gamma, tight exactly on Giraud's ten classes | [certificates/sigma_K5_4/](certificates/sigma_K5_4/) | (k) |
| Section "Stability for the codegree-squared density": (G7), (G8), (C8), (C9), (C10) and the margin mu_0 of Proposition margin | 2,404 / 10 and 32,981 / 22 labelled / classes; locally Giraud => Giraud for n = 8, 9, 10 (15,636,107 on ten vertices); mu_0 = 176165518826891630107/6280747422216628134215680 | [certificates/sigma_K5_4/](certificates/sigma_K5_4/) (for mu_0 and the constants) | (l); C10 over all bases and the minimality of mu_0: `--raw stab_c10`, `--raw sigma` |

## Further codegree-squared bounds

Theorem sigmacat: each entry bounds sigma(F). See the [exact fractions](certificates/sigma_catalogue/README.md).
Check (m) runs the fast checks; `--raw sigcat` runs the exhaustive scans.

| F | upper bound | certificate |
|---|---|---|
| J_4 | < 0.280703690307 | [J4](certificates/sigma_catalogue/J4/) |
| K_5^< | < 0.333646 | [K5lt](certificates/sigma_catalogue/K5lt/) |
| K_5^= | < 0.384426 | [K5eq](certificates/sigma_catalogue/K5eq/) |
| K_5^{3-} | < 0.405931 | [K5_3minus](certificates/sigma_catalogue/K5_3minus/) |
| C_5, tight 5-cycle | < 0.252305 | [C5](certificates/sigma_catalogue/C5/) |
| K_6^(3) | < 0.742319 | [K6_3](certificates/sigma_catalogue/K6_3/) |
| K_5^(4)- | < 0.204033 | [K5_4minus](certificates/sigma_catalogue/K5_4minus/) |
| K_6^(4)- | < 0.612038 | [K6_4minus](certificates/sigma_catalogue/K6_4minus/) |
| K_6^(4) | < 0.743754 | [K6_4](certificates/sigma_catalogue/K6_4/) |
| K_7^(4) | < 0.877156 | [K7_4](certificates/sigma_catalogue/K7_4/) |

## Limits of the plain method

| statement in the paper | bound | files | checks |
|---|---|---|---|
| Theorem n6 | on six vertices the plain method gives at most 1/4, 1/10, 0 and 1/2 for t(5,4), t(6,4), t(7,4), t_2(5,4) | [certificates/n6_dual_points/](certificates/n6_dual_points/) | (c) |
| Remark n6sharp | six-vertex certificates attaining 1/4, 1/10 and 1/2 | [certificates/n6_dual_points/cert_*_sharp.json](certificates/n6_dual_points/) | (j); also (g), first implementation |
| Theorem n7opt | no plain seven-vertex certificate proves t(5,4) > 1336237682928914994292138923/(7*2^89) | [certificates/n7_dual_point/](certificates/n7_dual_point/) | (d) |
| Table n7limits | no plain seven-vertex certificate proves t(6,4) > 0.139837689411, t(7,4) > 0.064764464150, t_2(5,4) > 0.553946861054, t_2(6,4) > 0.220748423194, or, without the (5,6) block, t(5,4) > 0.281644555065 (exact fractions in `certificates/n7_dual_points/README.md`) | [certificates/n7_dual_points/](certificates/n7_dual_points/) | (i) |
| Remark j4limit | every plain seven-vertex certificate for J_4 proves at best sigma(J_4) <= 1 - V, 1 - V > 0.280703275315 (exact dual point); certificates vanishing on the moments of the extremal construction prove at best 16/57 + 316891/14980607589 (exact witness) | [certificates/sigma_catalogue/J4_limit/](certificates/sigma_catalogue/J4_limit/) | (n); controls `--raw j4_controls` |

## Other bounds and construction checks

| statement in the paper | bound | files | checks |
|---|---|---|---|
| Theorem five | pi(K_6^(5)) < 0.767401, pi(H^5_4) < 0.389709, pi(K_7^(5)) < 0.927918 (exact fractions in Appendix A) | [certificates/fivegraphs/](certificates/fivegraphs/) | (b) |
| Proposition lotup (Section "Upper bounds for lottery numbers") | (b) C(7,4) l(7,4,5) <= c745 = 3467292853891919529774577/7750000000000000000000000 < 0.447393; (c) C(7,4) l(7,4,6) <= (107 + 3 c745)/384 < 0.282142 (type-pattern certificates); (a) C(6,4) l(6,4,5) <= 7/16 is an explicit construction (numerical companion checks only) | [certificates/lottery_upper/](certificates/lottery_upper/) | (p) |
| Giraud's construction against the t(5,4) certificate; the certificate does not force parity | exact values, see `verifier/EXPECTED.txt` | [certificates/K5_4/](certificates/K5_4/) | (e) |
| the r = 3 comparison t(4,3) > 0.426923 on five vertices | 60084175574225/2^47 (see the note in EXPECTED.txt) | [verifier/indep_k43_n5.py](verifier/indep_k43_n5.py) | (f) |

## Citation

Please cite the paper:

    @misc{MerolaCorrea-turan4,
      author = {Bruno M{\'e}rola Corr{\^e}a},
      title  = {Tur{\'a}n densities of complete 4-graphs via flag algebras on seven vertices},
      year   = {2026},
      note   = {Preprint; arXiv identifier to be added}
    }

and, for the data, this package: https://github.com/brunomerola/turan4-certificates (releases v1.0.0 and later; see also
`CITATION.cff`).
