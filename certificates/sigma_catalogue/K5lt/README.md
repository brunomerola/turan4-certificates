# sigma_catalogue/K5lt: the codegree-squared density of K_5^<

**Proves** the second line of Theorem sigmacat of the paper (Section "Further codegree-squared densities"):

    sigma(K_5^<) <= 1375676321453/4123168604160 < 0.333645420191      (4123168604160 = 15 * 2^38)

K_5^< is the 3-graph on five vertices with eight edges whose two missing edges meet in exactly two vertices (1-based:
123 124 125 134 135 145 234 235; missing 245 and 345). Balogh, Clemen and Lidicky conjecture sigma(K_5^<) = 1/3 (Turán's
construction for the tetrahedron); the gap is 3.12e-4.

**The claim.** Admissible graphs: K_5^<-free 3-graphs on 7 vertices (no copy of K_5^<, not necessarily induced).
Objective `cosig3`, f(H) = 1 - gamma(H) = F(H)/210 with F(H) = 210 - sum over the 35 4-subsets Q of [7] of
C(e_H(Q), 2), exactly as for J_4 (see `../J4/README.md` and `../README.md`). The certificate consists of PSD matrices
Q_k and

    min over all K_5^<-free 7-vertex 3-graphs H of  f(H) - sum_k c_{Q_k}(H)  =  b  =  2747492282707/4123168604160

exactly (b = 0.666354579809...). Hence sigma(K_5^<) <= 1 - b. The minimum is attained by 7 labelled one-vertex
extensions (one isomorphism class; the producer's argmin is 1073740992); the next value is larger by 5.30e-7.

**Blocks and matrices.** The nine blocks of `../J4/`, 45 keys, 28,108 flags (largest Gram matrix 1024);
Q_k = A_k^T A_k / M^2 with M = 2^19 and integer A_k (64-bit Gram entries suffice; the checker asserts this).

**Files**

* `cosig3_s3_K5lt_f.cosig3.cert.npz`: integer arrays `A0` ... `A44` (int64).
* `cosig3_s3_K5lt_f.cosig3.cert.json`: blocks, `M`, `M_bits`, `tau` = [0, 0], the key list **without the flag lists**
  (`s`, `m`, `type_index`, `sigma`, `n_flags`), the exact `bound`, `argmin_graph`, the a-priori bound `z_bound_apriori`,
  the producer's table `groups`, `raw_extensions_scanned` = 40,527,893 and `cert_npz_sha256`; the other fields are
  records of the run.
* `cosig3_s3_K5lt_f.keys.json`: the full key list (types, |Aut|, factor rows and every flag list, colex masks), the
  formula, the predicate and the sha256 of the three certificate files; written for the package, because the
  certificate json does not list the flags (finding M1 of the independent review).
* `cosig3_s3_K5lt_f.cosig3.verify.json`: the producer's sampling check (308 graphs); a record, no role in the proof.

The .cert.* and .verify.json files are byte-identical to the files reviewed (sha256 d7f20322..., 5c5a3d1b...,
5b8025c6...).

**Checked by** `verifier/verify_all.sh` (name K5lt; `verifier/EXPECTED.txt`, part (m)):

* fast: `sigcat_keys` (the key list rebuilt from the definitions equals `keys.json` and the json's keys), the
  reductions `sigcat_reduction` (including the K_5^<-freeness test against an independent definitional test),
  `sigcat_pycheck_K5lt` and `sigcat_pyvals_K5lt` (pure-Python values at 17 graphs, equal to the independent review's;
  4 at b), `sigcat_decimals`;
* `--raw sigcat_K5lt`: `s7_prep.py` (own types and flags, equal to the json's key list; exact Gram matrices; the
  a-priori |Z| bound equal to the json's), `s7_eval` at the 17 spot-check graphs (equal to the Python values), the
  exhaustive scan over all 56,623,104 one-vertex extensions of the checker's own 1,728 six-vertex representatives
  (40,527,893 admissible), and `s7_analyze.py`: minimum exactly b, attained 7 times, all 108 group maxima equal the
  producer's table, 20,294,982,806 labelled admissible 7-vertex graphs. 1,175 s on 2 threads in the review; peak
  memory about 0.85 GB (the prep; the scan 0.2 GB).

**Provenance.** Phase F (frozen cuts, suffix `_f`) certificate of the cloud run paperb-a4 (2026-10-06; results in
commit bf6de3a7, branch research/turan4-n7 of the private working repository); key list from commit da3e33b0.
Independent review R7_SIG7 (2026-10-06/07): PASS-with-notes (notes on archiving -- M1, the flag lists, addressed by
`keys.json` -- and on decimals and interpretation; none affects validity). The weaker phase-A certificate is not
included. The run stopped before convergence, so it neither shows nor rules out that seven vertices reach 1/3.
