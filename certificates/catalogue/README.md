# catalogue: threshold families F_{p,lambda} at N = 7

**Proves** Theorem catalogue of the paper (Table catalogue). Let F_{p,lambda} be the family of p-vertex 4-graphs with
at least C(p,4) - lambda + 1 edges. A 4-graph is F_{p,lambda}-free iff its complement spans at least lambda edges in
every p-set, so pi(F_{p,lambda}) = 1 - t_lambda(p,4). One subdirectory per case, each with its own README (family,
exact bound, certificate, checks, provenance):

| directory | (p, lambda) | minimal members of F_{p,lambda} | exact bound b <= t_lambda(p,4) | t > | pi(F) < |
|---|---|---|---|---|---|
| `p5_lam3/` | (5, 3) | H^4_3 (validation: t_3(5,4) = 3/4 is known) | 3298534811415/4398046511104 | 0.749999983648 | 0.250000016352 |
| `p6_lam3/` | (6, 3) | K_6^(4) minus 2 edges | 11166973764957/38482906972160 | 0.290180099259 | 0.709819900741 |
| `p6_lam4/` | (6, 4) | K_6^(4) minus 3 edges | 391547931101/1099511627776 | 0.356110768826 | 0.643889231174 |
| `p6_lam5/` | (6, 5) | K_6^(4) minus 4 edges | 64873027783677/153931627888640 | 0.421440536123 | 0.578559463877 |
| `p6_lam6/` | (6, 6) | K_6^(4) minus 5 edges | 612132591153137/1231453023109120 | 0.497081561103 | 0.502918438897 |
| `p6_lam9/` | (6, 9) | 6-vertex 4-graphs with 7 edges | 217134104143267/307863255777280 | 0.705293990330 | 0.294706009670 |
| `p6_lam11/` | (6, 11) | 6-vertex 4-graphs with 5 edges | 24267758577011/28862180229120 | 0.840815156178 | 0.159184843822 |
| `p7_lam2/` | (7, 2) | K_7^(4)- | 111867473425/1099511627776 | 0.101742874380 | 0.898257125620 |
| `p7_lam3/` | (7, 3) | K_7^(4) minus 2 edges | 4043518541347/30786325577728 | 0.131341381781 | 0.868658618219 |
| `p7_lam4/` | (7, 4) | K_7^(4) minus 3 edges | 11199502214035/69269232549888 | 0.161680760732 | 0.838319239268 |

The decimals are those of the paper's table (t floored, pi ceiled to 12 digits). The cases (5,1), (5,2), (6,1), (6,2)
and (7,1) are K_5^(4), K_5^(4)-, K_6^(4), K_6^(4)- and K_7^(4), in `../K5_4/` ... `../K6_4minus/`.

**Format and checks.** The certificates have the format of the 4-graph certificates of Theorem main (see
`../README.md`) and are checked by the same programs, with lambda passed through `RV_LAM` and as the second argument
of `rv_reps.py`. `verifier/verify_all.sh --fast` runs, for every case, the representatives and Burnside counts, the
flag-list checks, the pure-Python value at the certificate's argmin, and the step `cat_decimals` (exact bounds and
table decimals); `--raw cat_pP_lamL` (or `--raw catalogue` for all ten) runs the exhaustive raw scan. For
p5_lam3 and p6_lam11 the raw scan takes seconds and for p6_lam9 about 15 minutes; the other seven take hours on
2 threads. See `verifier/EXPECTED.txt`, part (h).

**Which certificate.** Every case was run as a time-limited LP cut and column generation (phase A) and then
continued over a frozen cut set (phase F). The paper cites the phase-F certificates (`*_f.cert.*`), except for (5,3),
where phase F found no better point and wrote no certificate, so the phase-A certificate `c5l3.cert.*` is cited.
The certificates not cited are not included.

**Provenance.** Cloud run of 2026-10-06 on three VMs (paperb-a1: (5,3), (6,3), (7,2), (7,3); paperb-a2: (6,4),
(6,5), (6,6); paperb-a3: (7,4), (6,9), (6,11)) with the search code of `../../search/flagalg/h44/` (identical to the
code in this package); results recorded in commit 51154f47 (branch research/turan4-n7) of the private working
repository. Independent review R7_CAT (2026-10-06): all ten PASS (nine of them "with notes", on wording outside the
certified numbers); the expected values in `verifier/verify_all.sh` and `verifier/EXPECTED.txt` are taken from its
outputs.
