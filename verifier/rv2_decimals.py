"""Independent checker: exact decimal checks of the numbers quoted for the N = 7 dual points (V from the JSONs,
recomputed exactly; b from the certificates' "bound" fields, as claimed)."""
import json
import os
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction

import rv1_lib as L

getcontext().prec = 90
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates", "n7_dual_points")


def Vof(tag):
    D = json.load(open(f"{RES}/dual_{tag}.json"))
    den = int(D["den"])
    return Fraction(sum(int(x["num"]) * L.popcount(L.colex_to_own(int(x["mask"]))) for x in D["support"]), 35 * den)


def dd(x, nd, mode=ROUND_FLOOR):
    return (Decimal(x.numerator) / Decimal(x.denominator)).quantize(Decimal(1).scaleb(-nd), rounding=mode)


def trunc20(x):
    return str(dd(x, 20, ROUND_FLOOR))


chk = []


def claim(name, ok, got):
    chk.append(ok)
    print(f"{'OK ' if ok else 'BAD'} {name}: {got}")


V = {t: Vof(t) for t in ["p5no56_W", "p5no56_Wodd", "p6_f3", "p7_f1", "h44_nt3", "k6m_W"]}
long = {"p5no56_W": "0.28164455506449893547", "p5no56_Wodd": "0.28164455506964684454",
        "p6_f3": "0.13983768941054627141", "p7_f1": "0.06476446414900658646", "h44_nt3": "0.55394686105381260646",
        "k6m_W": "0.22074842319302066019"}
for t, s in long.items():
    claim(f"V({t}) 20-digit expansion {s}...", trunc20(V[t]) == s, trunc20(V[t]))
claim("V(Wodd) - V(W) ~ 5.1e-12 (> 0)", V["p5no56_Wodd"] > V["p5no56_W"] and
      abs(float(V["p5no56_Wodd"] - V["p5no56_W"]) - 5.1e-12) < 0.05e-12, float(V["p5no56_Wodd"] - V["p5no56_W"]))
b56 = Fraction(35604499940047, 115448720916480)
claim("floor12(b with (5,6)) = 0.308400990997", str(dd(b56, 12)) == "0.308400990997", dd(b56, 12))
claim("ceil12(1 - b with (5,6)) = 0.691599009003", str(dd(1 - b56, 12, ROUND_CEILING)) == "0.691599009003",
      dd(1 - b56, 12, ROUND_CEILING))
gap = Decimal("0.308400990997") - Decimal("0.281644555065")
claim("0.308400990997 - 0.281644555065 = 0.026756435932", str(gap) == "0.026756435932", gap)
claim("b56 - V(e) >= 0.026756435932 (exact)", b56 - V["p5no56_W"] >= Fraction("0.026756435932"),
      float(b56 - V["p5no56_W"]))
B = {"p6_f3": (Fraction(9224492822869, 15 * 2 ** 42), "0.139827122482362", "1.0567e-05", "1.056693e-5"),
     "p7_f1": (Fraction(7476698908057, 115448720916480), "0.064762076605993", "2.3875e-06", "2.387545e-6"),
     "h44_nt3": (Fraction(48641756816487, 5 * 2 ** 44), "0.552992751369026", "9.5411e-04", "9.541097e-4"),
     "k6m_W": (Fraction(202684381843379, 105 * 2 ** 43), "0.219452823117469", "1.2956e-03", "1.295600e-3")}
for t, (b, bs, vmb, wid) in B.items():
    claim(f"b({t}) 15-digit {bs}...", str(dd(b, 15)) == bs, dd(b, 15))
    claim(f"V - b ({t}) = {vmb}", f"{float(V[t] - b):.4e}" == vmb, f"{float(V[t] - b):.6e}")
    w = Fraction(str(dd(V[t], 12, ROUND_CEILING))) - Fraction(str(dd(b, 12)))
    claim(f"width ceil12(V) - floor12(b) ({t}) = {wid}", f"{float(w):.6e}".replace("e-0", "e-") == wid, float(w))
    claim(f"b <= V ({t})", b <= V[t], "")
print("ALL OK" if all(chk) else f"{chk.count(False)} BAD")
