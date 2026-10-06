"""Negative controls for verify_dual.py (both must FAIL):
 neg_lp.json   : the frozen-cut LP pseudo-density y_LP (rounded, support 398) -- PSD only on the 1,337 cuts;
 neg_dead.json : the final dual point with weight 2^-40 moved onto one graph of W that facial reduction killed
                 (it puts mass on a row whose diagonal is zero)."""
import json

import numpy as np

import paths

K = 100
den = 1 << K


def write(fn, masks, y, src):
    num = [max(1, int(float(v) * 2.0 ** 60) << (K - 60)) for v in y]
    i = int(np.argmax(y))
    num[i] += den - sum(num)
    assert num[i] > 0 and sum(num) == den
    json.dump({"den": den, "source": src, "support": [{"mask": int(m), "num": n} for m, n in zip(masks, num)]},
              open(fn, "w"))


def main():
    d = np.load(paths.RES + "/lp_y.npz")
    W, y = d["W"], d["y"]
    s = y > 0
    write(paths.RES + "/neg_lp.json", W[s], y[s] / y[s].sum(), "y_LP (frozen-cut LP)")
    fin = json.load(open(paths.RES + "/dual_W_final.json"))
    alive = np.load(paths.RES + "/facial_W.npz")["alive"]
    killed = W[~alive & s][0]
    dfin = int(fin["den"])
    eps = dfin >> 40                                   # weight 2^-40
    sup = fin["support"] + [{"mask": int(killed), "num": eps}]
    big = max(range(len(fin["support"])), key=lambda i: sup[i]["num"])
    sup[big] = {"mask": sup[big]["mask"], "num": sup[big]["num"] - eps}
    assert sum(r["num"] for r in sup) == dfin
    json.dump({"den": dfin, "source": "dual_W_final + 2^-40 on a killed graph", "support": sup},
              open(paths.RES + "/neg_dead.json", "w"))
    print("written neg_lp.json, neg_dead.json; killed graph", int(killed))


if __name__ == "__main__":
    main()
