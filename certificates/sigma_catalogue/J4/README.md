# sigma_catalogue/J4: the codegree-squared density of J_4

**Proves** the first line of Theorem sigmacat of the paper (Section "Further codegree-squared densities"):

    sigma(J_4) <= 9061620987350745659/32281802128991715328 < 0.280703690307      (32281802128991715328 = 7 * 2^62)

J_4 is the 3-graph on five vertices with the six triples through one vertex u (1-based: 123 124 125 134 135 145); a
3-graph is J_4-free iff no vertex link contains K_4. The lower bound 16/57 = 0.2807017... (iterated blow-up of the
complement of the Fano plane, Proposition jt) is by hand; the gap is 1.936e-6. Remark j4limit (directory
`../J4_limit/`) shows that plain certificates on seven vertices cannot go below 0.280703275315.

**The claim.** Admissible graphs: J_4-free 3-graphs on 7 vertices. Objective (`cosig3` in the files): on a 3-graph H
on [7],

    f(H) = 1 - gamma(H) = F(H)/210,      F(H) = 210 - sum over the 35 4-subsets Q of [7] of C(e_H(Q), 2),

where gamma is the probability that P + x and P + y are both edges for a uniform pair P and a uniform ordered pair
(x, y) of distinct vertices outside P (two edges inside a 4-set Q meet in a pair; 210 = 21 pairs x 10 unordered
{x, y}). By the identity after Lemma co2 of the paper, co_2(G)/(C(n,2)(n-2)^2) = gamma(G) - (gamma(G) - d(G))/(n-2),
and f is a density on four vertices, so Remark objective applies. The certificate consists of PSD matrices Q_k, one per
key k (a block (s, m) and a type sigma), and

    min over all J_4-free 7-vertex 3-graphs H of  f(H) - sum_k c_{Q_k}(H)  =  b  =  23220181141640969669/32281802128991715328

exactly (b = 0.719296309693...), with c_Q as in Definition cert of the paper. Then gamma(G) <= 1 - b + O(1/n) for every
J_4-free G, hence sigma(J_4) <= 1 - b. The minimum is attained by 5 labelled one-vertex extensions (one isomorphism
class; the producer's argmin is 9134677119); the next value is larger by 1.79e-11.

**Blocks and matrices.** The nine blocks (0,3), (1,3), (2,3), (1,4), (2,4), (3,4), (3,5), (4,5), (5,6) (all (s, m) with
m >= 3 and 2m - s <= 7), 43 keys, 24,402 flags in all (largest Gram matrix 958). For every key, Q_k = A_k^T A_k / M^2
with M = 2^30, A_k an integer matrix (rows = `factor_rows` in the json), so Q_k is PSD by construction. The entries of
A_k reach 1.5e10 (about 2^34), so the Gram entries exceed 64 bits; the checker splits A_k into 16-bit limbs and
evaluates in 128-bit integers with asserted a-priori bounds (`s7_prep128.py`, `s7_eval128.c`).

**Files**

* `cosig3_s3_J4_r30.cosig3.cert.npz`: integer arrays `A0` ... `A42` (int64) and the scalar `M_bits` = 30.
* `cosig3_s3_J4_r30.cosig3.cert.json`: blocks, `M`, `M_bits`, `tau` = [0, 0] (no stationarity term), the key list
  (`s`, `m`, `type_index`, `sigma`, `n_flags`, `aut_size` and `flags`, colex masks), a `format` string (the formula
  above), the exact `bound` and `sigma_upper` = 1 - b, `argmin_graph`, the producer's per-objective-value table `groups`
  ([F, max Z, argmax, count] over its class list), its raw-scan count `raw_extensions_scanned` = 29,092,719, and
  `cert_npz_sha256` (checked by the step `npz_hashes`). The other fields (float values, times, the private path of the
  certificate it was re-rounded from) are records of the producer's run.
* `cosig3_s3_J4_r30.keys.json`: the key list again, with the formula, the predicate and the sha256 of the three files
  above, in the format shared by all certificates of `../` (see `../README.md`); written for the package.
* `cosig3_s3_J4_r30.cosig3.verify.json`: the producer's own sampling check (306 graphs, all >= b); a record, no role in
  the proof.

The .cert.* and .verify.json files are byte-identical to the files reviewed (sha256 95277da6..., b5c4035c...,
8398e731...).

**Checked by** `verifier/verify_all.sh` (name J4; see `verifier/EXPECTED.txt`, part (m)):

* fast: `sigcat_keys` (`sk_check_keys.py`, the producer's self-contained key-list checker: the key list rebuilt from
  the definitions equals `keys.json` and the json, shapes of the `A_k`, sha256 of the copies); `sigcat_reduction`
  (`s7_reduction.py`: the identities behind the objective, the 7-vertex numerator F, exact averaging over 7-subsets, the
  J_4-freeness test against an independent definitional test); `sigcat_pycheck_J4` and `sigcat_pyvals_J4`
  (`s7_pycheck.py`, a pure-Python evaluation from the definition at 17 graphs, every value equal to the independent
  review's, `../../../verifier/sigcat_spot_values.txt`; 4 of them at b); `sigcat_decimals` (b, 1 - b and the ceiling
  of the paper);
* `--raw sigcat_J4`: `s7_prep128.py` (own types and flags, equal to the json's key list including every flag list;
  exact Gram matrices; the a-priori bounds 2^81 for |Z| and 2^89 for the numerator of the slack), `s7_eval128` at the 17
  spot-check graphs (equal to the Python values, `s7_cmp.py`), then the exhaustive scan over all 49,577,984 one-vertex
  extensions of the checker's own 1,513 six-vertex representatives (29,092,719 admissible; completeness by heredity, as
  in Lemma ext of the paper), and `s7_analyze.py`: minimum exactly b, attained 5 times, all 106 group maxima equal the
  producer's table, 14,704,790,947 labelled admissible 7-vertex graphs. 1,598 s on 2 threads in the review (one of them
  shared), peak memory about 1.5 GB (the prep; the scan 0.8 GB).

**Provenance.** The cloud run paperb-a4 (2026-10-06) produced a J_4 certificate at M = 2^19 with bound
810174058051/2886218022912 (sigma(J_4) <= 0.280704386024); this certificate is the same floating-point dual rounded
again at M = 2^30 (producer step recorded in commit ed36f398, branch research/turan4-n7 of the private working
repository), which removes most of the rounding loss. The key list `keys.json` was written from the producer's flag
model and checked against a rebuild from the definitions (commit da3e33b0). Independent review R7_SIG7, addendum A1
(2026-10-07): PASS (own decoding, own 128-bit evaluator, raw scan and class-list scan, pure-Python spot checks,
construction law). The M = 2^19 certificate is superseded and not included.
