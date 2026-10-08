# lottery_upper: upper bounds for lottery numbers (Proposition lotup)

**Proves** parts (b) and (c) of Proposition lotup of the paper (Section "Upper bounds for lottery numbers"), through
Lemma pattern (type patterns) and Lemma kahn (Kahn's fractional form of the Frankl-Rodl-Pippenger theorem):

    (b)  C(7,4) l(7,4,5) <= c745 = 3467292853891919529774577/7750000000000000000000000 < 0.447392627
    (c)  C(7,4) l(7,4,6) <= (107 + 3 c745)/384 = 839651878561675758589323731/2976000000000000000000000000
                          < 0.282141089

Here l(k,r,p) = lim L(n,k,r,p)/C(n,r) is the normalised lottery number (Sidorenko 2023); the paper's Table lottery
gives the corresponding values C(k,4)/c: 78.23 for (7,4,5) and 124.05 for (7,4,6) (rounded down). Part (a),
C(6,4) l(6,4,5) <= 7/16, needs no certificate: it is an explicit construction (the Giraud system Gamma_d on two
copies of F_2^d with weights on its "good" 6-sets), proved by hand in the paper; its counts are checked numerically
by `lotup_giraud` and `lotup_check` (see below).

**The claim of a certificate.** A file describes the data of Lemma pattern: k, r, p, the number m of parts, part
weights x in (0,1]^m with sum 1, a set T of types of size r (count vectors over the m parts), values y_B >= 0 for
types B of size k, and optionally inner parts I with an inner requirement q and a constant for l(k,r,q). It proves
C(k,r) l(k,r,p) <= c with c = C(k,r) (sum_B y_B + c_inner sum_{i in I} x_i^r) provided

    (i)  every type P of size p with P_i < q for all i in I satisfies P >= e (coordinatewise) for some e in T, and
    (ii) sum_B a(e,B) y_B >= mu_e(x) for every e in T,
         with a(e,B) = prod_i C(B_i, e_i) and mu_e(x) = r! prod_i x_i^e_i / e_i!.

**File format** (plain text; lines starting with `#` are comments, records of the producer's run):

    k K | r R | p P | m M                        one line each
    inner i j ... q Q NAME VALUE                 optional: inner parts, requirement q, name and value of the constant
    x x_1 ... x_m                                part weights (exact rationals)
    T e_1 ... e_m                                one line per type of T
    y B_1 ... B_m VALUE                          one line per type B of size k with y_B > 0
    c VALUE                                      the claimed value c (exact rational)

**File parts 0, ..., m-1 are the paper's parts 1, ..., m.** In (c), for example, the inner parts `inner 0 1` are the
paper's parts 1 and 2 (I = {1, 2}), and the file's types `3 1 0 0`, `0 3 1 0`, `0 0 3 1`, `1 0 0 3`, `3 0 1 0`,
`0 3 0 1` are the paper's arcs 1->2, 2->3, 3->4, 4->1, 1->3, 2->4.

**Files**

| file | claim | m, x | T, y | value c |
|---|---|---|---|---|
| `k7_r4p5_m5_den1000000.txt` | (b): (k,r,p) = (7,4,5), no inner parts | m = 5, x = 0.317221, 0.126888, 0.302114, 0.126888, 0.126889 | 38 types (5 with four points in one part, 10 of shape (3,1), 10 of shape (2,2), 12 of shape (2,1,1), 1 of shape (1,1,1,1)); 27 values y_B | c745 (27 of the 38 rows of (ii) tight) |
| `f8_k7_r4p6_x_quarter.txt` | (c): (k,r,p) = (7,4,6), inner parts 0, 1 with q = 5 and constant c745 of (b) | m = 4, x = 1/4 each | 14 types: (4) in parts 2, 3, the six (2,2), the six arcs (3,1); y = 1/13440 on the 7-sets in parts 2, 3 and 1/768 on the six arc blocks (4,3) | (107 + 3 c745)/384 (all 14 rows tight) |

Each `.sha256` file holds the SHA-256 of its certificate (5f0f88f5... and 44abf740...), as recorded by the producer and
by the reviews. The files are byte-identical to the producer's. Only this certificate is shipped for (c): the producer's
other file for (c), with the same T and y but an older inner constant (37/81, value 1463/5184), is not needed.

**Checked by** `verifier/verify_all.sh`, part (p) (`--lotup`, also in `--fast`; `verifier/EXPECTED.txt`):

* `lotup_certs` (`rv3b_certs.py` with `rv3b_core.py`, the independent review R7_T3b's checker): for each file the
  sha256 sidecar, x > 0 with sum 1, condition (i) by domination over all types of size p and by brute force on a
  blow-up with p points per part (53,130 5-sets for (b), 134,596 6-sets for (c)), every row of (ii) exactly, the value
  c, and four negative controls per file (y shrunk by 1e-6, a needed type dropped, x mass moved, c lowered by 1e-12),
  all rejected;
* `lotup_derived` (`rv3b_derived.py`, R7_T3b): the composition of (c), 35 sum y = 107/384 and the inner weight
  x_1^4 + x_2^4 = 1/128 (file parts 0, 1), so c = 107/384 + c745/128 = (107 + 3 c745)/384 < 0.282141089;
* `lotup_giraud` (`rv3b_giraud.py`, R7_T3b): the numerical companion of (a): Gamma_d is a Turán (n,5,4)-system
  (brute force, d = 2, 3, 4), every cross edge lies in exactly (N/2-2)(N/4-2) good (3,3)-sets (exhaustively, d = 3, 4,
  5), the internal counts, the pair weights, the exact value 15 W_N / C(2N,4) as a rational function of N with limit
  7/16, and an end-to-end check that the good 6-sets form a (32,6,4,5)-lottery system at d = 4;
* `lotup_check` (`r7p6_check.py`, the independent review R7_P6 of the paper's text): 69 exact checks of the numbers of
  Section "Upper bounds for lottery numbers" -- both certificates again with an own parser (sha256, rows, domination,
  brute force, value, T and y equal to the paper's verbal description of (c)), the case analysis of the printed proof
  of (c)(i), the counts and the tiling of (a) at d = 3, 4, 5, the table values, the gap percentages and the constants
  of Section "Giraud systems are expensive to cover".

**Trust base.** Lemma kahn rests on Kahn 1996, Theorem 1.2(b), and the existence of l(k,r,p) on Sidorenko 2023; both
are cited, not re-proved. The independent review R7_T3b (2026-10-04) checked the certificate of (b) (row 3a), the
composition of (c) with the inner constant c745 (row 4; on the producer's other file for (c), which has the same T and
y) and part (a) (row 6); R7_F2 (2026-10-05, claim 2) checked the file for (c) shipped here, with the constant c745;
R7_P6 (2026-10-08, PASS-with-notes; its notes are applied in the paper) checked the text of Section "Upper bounds for
lottery numbers" with all its numbers and both files. In this package `rv3b_certs.py` also runs on the file for (c).
The checks of (i) and (ii) for (b) were done in exact arithmetic by two independent programs (`rv3b_certs.py`,
`r7p6_check.py`), as the paper states. Conditions (i) and (ii) are all a certificate needs; the producer's own
sampling of near-optimal x and its LP play no role.
