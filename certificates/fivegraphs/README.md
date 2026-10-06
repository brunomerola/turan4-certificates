# fivegraphs: three 5-graph bounds on seven vertices

**Proves** Theorem five of the paper (exact values in its Appendix A):

| file | problem (complement form) | exact bound | consequence |
|---|---|---|---|
| `cert_k65_N7_all.json` | every 6-set spans >= 1 edge | t(6,5) >= 5410603883828909241569858123831800227/23261489926236027775816623554906030080 > 0.232599197256 | pi(K_6^(5)) < 0.767400802744 |
| `cert_h45_N7_all.json` | every 6-set spans >= 3 edges | t_3(6,5) >= 28392597905796184739036021237749789573/46522979852472055551633247109812060160 > 0.610291902965 | pi(H^5_4) < 0.389708097035 |
| `cert_k75_N7_all.json` | every 7-set spans >= 1 edge | t(7,5) >= 1006048490611906127582298374928111023/13956893955741616665489974132943618048 > 0.072082548868 | pi(K_7^(5)) < 0.927917451132 |

**Format.** `r = 5`, `p`, `N = 7`, `bound`, `allowed` (the admissible edge counts of a p-set; `null` means at least
one), and `blocks`, each with `s`, `m`, `sigma`, `flags` (colex masks of 5-subsets) and the rational matrix
`Q = Qnum / den`. These certificates store Q itself, not a factor, so positive semidefiniteness is a property to be
checked, not a property by construction.

**Checked by** `verifier/rv_r5.py` (fast; part (b) of `verify_all.sh`): every Q is tested symmetric and PSD exactly
(LDL^T over the rationals); all 2^21 labelled 5-graphs on seven vertices are scanned, the admissible ones are split
into classes by orbit marking over all 5040 permutations (888, 150 and 1043 classes), and the exact minimum of
e(H)/21 - sum c_Q(H) is compared with `bound`.
