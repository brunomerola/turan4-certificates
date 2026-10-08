# sigma_catalogue/K5_4minus: the codegree-squared density of K_5^(4)-

**Proves** the third line of Theorem sigmacat of the paper (Section "Further codegree-squared densities"):

    sigma(K_5^(4)-) <= 6281420455337/30786325577728 < 0.204032807991      (30786325577728 = 7 * 2^42)

K_5^(4)- is the 4-graph on five vertices with four edges. The lower bound 1/5 is Sidorenko's doubling construction;
the gap is 4.0e-3 (the previous bound was 1/4).

**The claim.** Complement form, as for Theorem sigma (`../../sigma_K5_4/`): a 4-graph G is K_5^(4)--free iff its
complement H spans at least two edges in every 5-set (p = 5, lambda = 2, the admissibility of `../../K5_4minus/`).
The objective (`sigma4` in the files) is

    kappa(H) = 2 d(H) - gamma(H) = F(H)/420,   F(H) = 24 e(H) - sum over the 35 3-subsets T of [7] of c_T (c_T - 1),

c_T = number of edges of H containing T (Lemma co2 of the paper; 420 = C(7,3) * 4 * 3 ordered triples (T, x, y)).
The certificate consists of PSD matrices Q_k and

    min over all 7-vertex 4-graphs H with >= 2 edges in every 5-set of  kappa(H) - sum_k c_{Q_k}(H)  =  b
                                                                         =  24504905122391/30786325577728

exactly (b = 0.795967192009...). With Lemma co2 and gamma(G) - d(G) >= -1 this gives sigma(K_5^(4)-) <= 1 - b, as in
the proof of Theorem sigma. The minimum is attained by 7 labelled one-vertex extensions (one isomorphism class:
mask 33554431, the complement of the ten 4-sets through a fixed pair); the next value is larger by 2.57e-8.

**Blocks and matrices.** Blocks (1,4), (2,4), (3,4), (3,5), (4,5), (5,6), 10 keys, 2,243 flags (largest Gram matrix
768); Q_k = A_k^T A_k / M^2 with M = 2^21 and integer A_k.

**Files**

* `sigma4_c4_K5m_f.sigma4.cert.npz`: integer arrays `A0` ... `A9` (int64).
* `sigma4_c4_K5m_f.sigma4.cert.json`: blocks, `M`, `M_bits`, `tau` = [0, 0], the key list without the flag lists, the
  exact `bound`, `argmin_graph`, `z_bound_apriori`, the producer's table `groups`, `raw_extensions_scanned` =
  10,058,620 and `cert_npz_sha256`; the other fields are records of the run.
* `sigma4_c4_K5m_f.keys.json`: the full key list with every flag list, the formula, the predicate and the sha256 of
  the three certificate files; written for the package (finding M1 of the independent review R7_SIG7), and compared
  with an independent reviewer's own types, flags and |Aut| (review R7_SIG8, addendum A2: all equal; M1 resolved).
* `sigma4_c4_K5m_f.sigma4.verify.json`: the producer's sampling check (256 graphs); a record, no role in the proof.

The .cert.* and .verify.json files are byte-identical to the files reviewed (sha256 bf8ab264..., b0e8a528...,
c205d6ce...).

**Checked by** `verifier/verify_all.sh` (name K5_4minus; `verifier/EXPECTED.txt`, part (m)):

* fast: `sigcat_keys_indep3` (review R7_SIG8 A2) and `sigcat_keys` (producer), `sigcat_reduction` (including the
  identities of Lemma co2 on random 4-graphs and the equivalence "K_5^(4)--free <=> every 5-set of the complement
  spans >= 2 edges"), `sigcat_pycheck_K5_4minus` and `sigcat_pyvals_K5_4minus` (pure-Python values at 16 graphs,
  equal to the independent review's; 3 at b), `sigcat_decimals`;
* `--raw sigcat_K5_4minus`: `s7_prep.py`, `s7_eval` at the 16 spot-check graphs (equal to the Python values), the
  exhaustive scan over all 65,011,712 one-vertex extensions of the checker's own 62 six-vertex representatives
  (10,058,620 admissible), and `s7_analyze.py`: minimum exactly b, attained 7 times, all 56 group maxima equal the
  producer's table, 1,741,143,596 labelled admissible 7-vertex graphs (the same counts as `../../K5_4minus/`, which has
  the same admissibility). 487 s on 2 threads in the review; less than 0.1 GB of memory.

**Provenance.** Phase F certificate of the cloud run paperb-a4 (2026-10-06; commit bf6de3a7, branch
research/turan4-n7 of the private working repository); key list from commit da3e33b0. Independent review R7_SIG7
(2026-10-06/07): PASS-with-notes (archiving, decimals, interpretation; none affects validity; the archiving note M1
was resolved by the independent comparison of the key list in R7_SIG8 addendum A2, 2026-10-08). The phase-A
certificate is not included.
