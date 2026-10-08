"""Independent checker: F8 propagation c(7,4,6) <= (107 + 3 c745)/384 and split bounds c(7,4,8), c(6,4,8) -- exact / rigorous.

Split (Montecalvo): parts with requirements q = (4,5), p = 1 - 2 + 9 = 8; inner costs c(k,4,4) = 1 (Rodl) and c(k,4,5).
   c(k,4,8) <= f(a) := a^4 * 1 + (1-a)^4 * c45 for every real a in (0,1)   (part sizes floor(a n)).
   min_a f = c45 / (1 + c45^(1/3))^3 at a* = s/(1+s), s = c45^(1/3).
Rigour: (i) exhibit a RATIONAL a with f(a) <= claimed decimal (exact Fractions) -- this alone proves the claim;
        (ii) the exact minimum to 30 digits (Decimal) to show the claimed decimal is the min rounded UP.
"""
import os
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv3b_core import parse, decimal_up, mass, a_coef

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.normpath(os.path.join(HERE, "..", "certificates", "lottery_upper"))


def cert_value(fn):
    d = parse(os.path.join(CERTS, fn))
    return d


def cbrt_dec(x):
    x = Decimal(x.numerator) / Decimal(x.denominator)
    g = Decimal(float(x) ** (1 / 3))
    for _ in range(200):
        g = g - (g ** 3 - x) / (3 * g * g)
    return g


def split_check(name, c45, claimed):
    s = cbrt_dec(c45)
    cdec = Decimal(c45.numerator) / Decimal(c45.denominator)
    mn = cdec / (1 + s) ** 3
    astar = s / (1 + s)
    # rational a near a*
    best = None
    for den in (10 ** 6, 10 ** 9, 10 ** 12):
        a = Fr(int((astar * den).to_integral_value()), den)
        f = a ** 4 + (1 - a) ** 4 * c45
        if best is None or f < best[1]:
            best = (a, f)
    a, f = best
    cl = Fr(claimed)
    print("%s: c45 = %.12f ; exact min = %s ; a* = %.12f" % (name, float(c45), str(mn)[:22], float(astar)))
    print("    rational a = %s gives f(a) = %.15f <= claimed %s : %s ; claimed - min = %.3e ; min rounded up 9dp = %s"
          % (a, float(f), claimed, f <= cl, float(cl) - float(mn), decimal_up_dec(mn)))
    return f <= cl and Decimal(claimed) >= mn


def decimal_up_dec(v, digits=9):
    q = (v * Decimal(10) ** digits).to_integral_value(rounding="ROUND_CEILING")
    return str(q / Decimal(10) ** digits)


if __name__ == "__main__":
    # c745 values (from certificates, values re-derived by rv3b_certs.py)
    d5 = parse(os.path.join(CERTS, "k7_r4p5_m5_den1000000.txt"))
    c745_m5 = d5["c"]          # verified equal to recomputed value in rv3b_certs.log
    c745_37 = Fr(37, 81)
    # F8: LP part at x = 1/4 from the F8 certificate, recomputed here
    d8 = parse(os.path.join(CERTS, "f8_k7_r4p6_x_quarter.txt"))
    sy = sum((v for _, v in d8["y"]), Fr(0))
    lp_part = 35 * sy
    inner_w = sum((d8["x"][i] ** 4 for i in d8["inner"]["parts"]), Fr(0))
    print("F8: 35*sum(y) = %s (== 107/384: %s) ; inner weight x0^4 + x1^4 = %s (== 1/128: %s)"
          % (lp_part, lp_part == Fr(107, 384), inner_w, inner_w == Fr(1, 128)))
    for lab, c in (("37/81", c745_37), ("m5 cert", c745_m5)):
        v = lp_part + c * inner_w
        print("   c(7,4,6) <= 35 sum y + c745/128 with c745 = %s: %s = %.12f ; == (107+3c)/384: %s ; up 9dp %s"
              % (lab, v if c == c745_37 else "...", float(v), v == (107 + 3 * c) / 384, decimal_up(v)))
    v = lp_part + c745_m5 * inner_w
    print("   claimed 0.282141089 >= exact: %s" % (Fr("0.282141089") >= v))
    print()
    ok1 = split_check("c(7,4,8) [q=(4,5), c45 = c745 m5 cert]", c745_m5, "0.081392194")
    ok2 = split_check("c(6,4,8) [q=(4,5), c45 = 7/16 (F9)]", Fr(7, 16), "0.080365839")
    ok3 = split_check("c(7,4,8) [with 37/81, old row]", c745_37, "0.082355105")
    print("split claims valid:", ok1, ok2, ok3)
    # pigeonhole for split: p = 1 - m + sum q
    print("split domination: p = 1 - 2 + (4 + 5) =", 1 - 2 + 4 + 5)
