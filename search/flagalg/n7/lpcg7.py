"""Matrix-free LP-CUT-CG for t(p, 4) at N = 7 (r = 4), p in {5, 6, 7}.  Imports highspy; never import ortools here.

Method (Jeong-Park-Im-Lee-Yang 2026, arXiv 2609.27495, Section 4, code github.com/taeyool/tetrahedron-turan,
Apache-2.0; Ahmadi-Dash-Hall 2017 for the LP inner approximation of the PSD cone by rank-one columns).  This is an
independent re-implementation of the idea; no code was copied.

Bound side (what is certified):  maximise b over lambda, tau >= 0 such that for every admissible 7-vertex H
    b + sum_v lambda_v v^T M_{sigma(v)}(H) v + tau1 s1(H) + tau2 s2(H) <= d(H),
    s1(H) = b0 d(H) - c_dd(H),   s2(H) = b0 (1 - d(H)) - (d(H) - c_dd(H))     (stationarity, see below).
Q_sigma = sum lambda_v v v^T is PSD by construction, so a complete scan of ALL admissible H (numba, matrix-free:
M_sigma(H) is never stored for H outside the working set) gives a valid bound for any lambda, tau.

Restricted LP actually solved (its LP dual, Jeong's orientation, so that the simplex basis has the size of the cut
set rather than of the graph set):  columns y_H >= 0 for H in a working set W (+ slacks z, w1, w2),
    min  sum_H y_H d(H) + cap z + tau_max (w1 + w2)
    s.t. sum_H y_H + z = 1                          (dual: b <= cap)
         sum_H y_H v^T M(H) v >= 0   for each cut v  (dual: lambda_v >= 0)
         sum_H y_H s_i(H) + w_i >= 0                 (dual: 0 <= tau_i <= tau_max).
Cuts = eigenvectors of the most negative eigenvalues of the moment matrices Y_sigma = sum_H y_H M_sigma(H); columns
= the H with the most negative reduced cost slack(H) - b found by the scan.  b, lambda, tau are the row duals.

Stationarity (degree regularity; Razborov's differential method, here in its elementary cloning form, A0_METHODS
§(a).2 and Sidorenko's Lemma 2): in an n-vertex graph G with the maximum number of edges and no K_p^(4), deleting a
vertex u and adding a non-adjacent twin of v loses deg(u) edges, gains >= deg(v) - C(n-2, 2), and creates no K_p
(a K_p through both twins would need an edge containing both; one through the twin only maps to one through v).
Hence all degrees differ by O(n^2) = o(n^3): extremal graphs (and their complements, the minimisers here) are
degree-regular, and for every limit t and every 1-flag F >= 0:  [[(d - t) F]]_1 = 0, d the rooted edge flag.
At N = 7 the product rho * F does not fit (rho needs 4 vertices, F >= 4), so we use the relaxation valid when
t <= b0:  [[(d - b0) F]]_1 = (t - b0) [[F]] <= 0, for F = d (tau1) and F = 1 - d (tau2), all expressible at N = 7:
    [[d . d]]_1 = c_dd(H) = 72 cdd(H)/5040,   [[d (1 - d)]]_1 = d - c_dd,   [[d]] = d(H).
For fixed Q and tau the bound is  t >= max_{b0} min(b0, min_H slack(H; b0)) = min_H a_H / (1 + beta_H)  with
a_H = slack(H; 0), beta_H = tau1 d + tau2 (1 - d) (exact fixed point; see certify7.py).  b0 follows a Dinkelbach
iteration (b0 <- current fixed-point bound); the plain problem (tau = 0) is solved first.

Usage: python lpcg7.py P [--seconds 3600] [--sample 0.1] [--blocks all|std] ...
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import time

import numpy as np
import scipy.linalg

from common import DATA, build_model
from kernels import coef_columns, moment_raw, qsum_flat, scan_float, tables, terms

import highspy  # noqa: E402  (after numba imports; no ortools in this process)

INF = highspy.kHighsInf
STD_BLOCKS = ((1, 4), (3, 5), (5, 6))
ALL_BLOCKS = ((1, 4), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6))
NFIX = 3          # fixed rows (norm, stat1, stat2) and fixed columns (z, w1, w2)


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def popcount_array(x):
    e = np.zeros(x.size, dtype=np.int64)
    u = x.copy()
    while u.any():
        e += u & 1
        u >>= 1
    return e


class LPCG7:
    def __init__(self, model, universe, b0, a):
        self.m, self.a = model, a
        self.U = universe
        self.eU = popcount_array(universe)
        self.tab = tables(model)
        self.nkeys = len(model.keys)
        self.keyid = np.full((len(model.blocks), 6), -1, dtype=np.int64)
        for k, (bi, t) in enumerate(model.keys):
            self.keyid[bi, t] = k
        self.P = [model.blocks[bi].aut[t] for bi, t in model.keys]
        self.dims = [int(model.qdim[bi, t]) for bi, t in model.keys]
        self.b0 = float(b0)
        self.W, self.Wset = [], set()
        self.e = np.zeros(0)
        self.cdd = np.zeros(0)
        self.T = [[np.zeros(0, np.int64)] * 4 for _ in range(self.nkeys)]   # graph-column terms per key
        self.cuts = []                                    # row NFIX + j -> (key, unit vector, birth iteration)
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("primal_feasibility_tolerance", a.lp_tol)
        h.setOptionValue("dual_feasibility_tolerance", a.lp_tol)
        h.setOptionValue("threads", a.lp_threads)
        if a.lp_threads > 1:
            h.setOptionValue("parallel", "on")
        h.setOptionValue("solver", a.lp_solver)
        if a.simplex_strategy >= 0:
            h.setOptionValue("simplex_strategy", a.simplex_strategy)
        h.addRows(NFIX, np.array([1.0, 0.0, 0.0]), np.array([1.0, INF, INF]), 0, np.zeros(NFIX, np.int32),
                  np.zeros(0, np.int32), np.zeros(0))
        h.addCols(NFIX, np.array([1.0, 0.0, 0.0]), np.zeros(NFIX), np.full(NFIX, INF), NFIX,
                  np.arange(NFIX, dtype=np.int32), np.arange(NFIX, dtype=np.int32), np.ones(NFIX))
        self.h = h
        self.it = 0
        self.stat_on = False
        self.best = {"b_star": -np.inf}
        self.hist = []
        self.hot = None

    # ------------------------------------------------------------------------------------------- graphs (columns)
    def s_coefs(self, e, cdd, b0):
        d = e / 35.0
        c = 72.0 * cdd / 5040.0
        return b0 * d - c, b0 * (1.0 - d) - (d - c)

    def add_graphs(self, masks):
        masks = [int(x) for x in masks if int(x) not in self.Wset]
        if not masks:
            return 0
        H = np.array(masks, dtype=np.int64)
        rows, keys, aa, bb, ww, ee, cc = terms(H, *self.tab[:12], self.keyid, self.m.cdd_pos)
        code = ((rows * self.nkeys + keys) * 1024 + aa) * 1024 + bb
        order = np.argsort(code, kind="stable")
        code, ww = code[order], ww[order]
        u, start = np.unique(code, return_index=True)
        wsum = np.add.reduceat(ww, start)
        bb = u % 1024
        aa = (u // 1024) % 1024
        keys = (u // (1024 * 1024)) % self.nkeys
        rows = u // (1024 * 1024 * self.nkeys)
        c0 = len(self.W)
        nnew = H.size
        dense = np.zeros((NFIX + len(self.cuts), nnew))
        dense[0] = 1.0
        dense[1], dense[2] = self.s_coefs(ee, cc, self.b0)
        bykey = {}
        for j, (k, v, _) in enumerate(self.cuts):
            bykey.setdefault(k, []).append(j)
        for k in range(self.nkeys):
            sel = keys == k
            if not sel.any():
                continue
            rk, ak, bk, wk = rows[sel], aa[sel], bb[sel], wsum[sel]
            if k in bykey:
                js = bykey[k]
                V = np.array([self.cuts[j][1] for j in js])
                dense[NFIX + np.array(js)] = coef_columns(rk, ak, bk, wk, self.P[k], V, nnew)
            T = self.T[k]
            self.T[k] = [np.concatenate([T[0], rk + c0]), np.concatenate([T[1], ak]), np.concatenate([T[2], bk]),
                         np.concatenate([T[3], wk])]
        D = dense.T                                        # (nnew, nrows): column-wise entries
        nz = np.abs(D) > 1e-15
        starts = np.concatenate([[0], np.cumsum(nz.sum(axis=1))[:-1]]).astype(np.int32)
        idx = np.nonzero(nz)
        self.h.addCols(nnew, ee / 35.0, np.zeros(nnew), np.full(nnew, INF), int(nz.sum()), starts,
                       idx[1].astype(np.int32), D[idx])
        self.W.extend(masks)
        self.Wset.update(masks)
        self.e = np.concatenate([self.e, ee])
        self.cdd = np.concatenate([self.cdd, cc])
        return nnew

    def set_b0(self, b0):
        self.b0 = float(b0)
        s1, s2 = self.s_coefs(self.e, self.cdd, self.b0)
        for i in range(len(self.W)):
            self.h.changeCoeff(1, NFIX + i, float(s1[i]))
            self.h.changeCoeff(2, NFIX + i, float(s2[i]))

    # ------------------------------------------------------------------------------------------- cuts (rows)
    def add_cuts(self, items):
        if not items:
            return 0
        nW = len(self.W)
        bykey = {}
        for k, v in items:
            bykey.setdefault(k, []).append(v)
        starts, idx, val, kk, new = [], [], [], 0, []
        for k, vs in bykey.items():
            T = self.T[k]
            C = coef_columns(T[0], T[1], T[2], T[3], self.P[k], np.array(vs), nW)
            for v, c in zip(vs, C):
                if self.a.row_scale:
                    # scale the cut row to max |coef| = 1 by storing v * sqrt(s) (Q = sum lambda v v^T unchanged)
                    sc = 1.0 / max(float(np.abs(c).max()), 1e-12)
                    c = c * sc
                    v = v * np.sqrt(sc)
                nz = np.nonzero(np.abs(c) > 1e-15)[0]
                starts.append(kk)
                idx.append((nz + NFIX).astype(np.int32))
                val.append(c[nz])
                kk += nz.size
                new.append((k, v, self.it))
        n = len(new)
        self.h.addRows(n, np.zeros(n), np.full(n, INF), kk, np.array(starts, dtype=np.int32),
                       np.concatenate(idx), np.concatenate(val))
        self.cuts.extend(new)
        return n

    def prune_cuts(self, lam):
        if len(self.cuts) <= self.a.keep_cuts:
            return 0
        old = np.array([self.it - c[2] > self.a.prune_age for c in self.cuts])
        dead = np.nonzero((lam <= 1e-14) & old)[0]
        if dead.size == 0:
            return 0
        self.h.deleteRows(dead.size, (dead + NFIX).astype(np.int32))
        keep = np.ones(len(self.cuts), dtype=bool)
        keep[dead] = False
        self.cuts = [c for c, kp in zip(self.cuts, keep) if kp]
        return int(dead.size)

    # ------------------------------------------------------------------------------------------- Grams, moments
    def moments(self, y):
        Ys = []
        for k in range(self.nkeys):
            T = self.T[k]
            Rm = moment_raw(T[0], T[1], T[2], T[3], y, self.dims[k])
            P = self.P[k]
            Y = np.zeros_like(Rm)
            for g in range(P.shape[0]):
                Y += Rm[np.ix_(P[g], P[g])]
            Ys.append(Y)
        return Ys

    def grams(self, lam):
        Qs = [np.zeros((n, n)) for n in self.dims]
        for (k, v, _), l in zip(self.cuts, lam):
            if l > 0:
                Qs[k] += l * np.outer(v, v)
        return Qs

    def factors(self, lam):
        F = [[] for _ in self.dims]
        for (k, v, _), l in zip(self.cuts, lam):
            if l > 0:
                F[k].append(np.sqrt(l) * v)
        return [np.array(f).reshape(-1, n) for f, n in zip(F, self.dims)]

    # ------------------------------------------------------------------------------------------- pricing
    def scan(self, lam, tau1, tau2, subset=None):
        qs = qsum_flat(self.m, self.grams(lam))
        U = self.U if subset is None else self.U[subset]
        eU = self.eU if subset is None else self.eU[subset]
        sl, bad = scan_float(U, *self.tab, qs, tau1, tau2, self.b0, self.a.nchunks)
        if bad.any():
            raise RuntimeError(f"{int(bad.sum())} inadmissible masks in the universe")
        d = eU / 35.0
        beta = tau1 * d + tau2 * (1.0 - d)
        bstar = float(np.min((sl + self.b0 * beta) / (1.0 + beta)))
        return sl, bstar

    def price(self, lam, tau1, tau2, b_lp, full):
        a = self.a
        if full or self.hot is None:
            sl, bstar = self.scan(lam, tau1, tau2)
            nh = min(a.hot_size, sl.size)
            self.hot = np.argpartition(sl, nh - 1)[:nh] if nh < sl.size else np.arange(sl.size)
            idx = np.arange(sl.size)
            info = {"scan": "full", "b_valid_at_b0": float(sl.min()), "b_star": bstar}
        else:
            sl, _ = self.scan(lam, tau1, tau2, self.hot)
            idx = self.hot
            info = {"scan": "hot", "hot_min_at_b0": float(sl.min())}
        viol = np.nonzero(sl < b_lp - a.price_tol)[0]
        added = 0
        if viol.size:
            pick = viol[np.argsort(sl[viol])[:a.rows_per_price * 3]]
            pick = [h for h in self.U[idx[pick]] if int(h) not in self.Wset][:a.rows_per_price]
            added = self.add_graphs(pick)
        return added, info

    # ------------------------------------------------------------------------------------------- control
    def seed_n6(self):
        """Cuts from the N = 6 flag SDP (engine, Clarabel): its blocks (s, m) with 2m - s = 6 are blocks of this
        N = 7 model with the same flags, so its solution is feasible here."""
        from flags import build_problem
        from hypergraphs import EveryPSetHasEdge
        from sdp import solve_flag_sdp
        prob = build_problem(6, 4, EveryPSetHasEdge(self.m.p), type_mode="standard")
        num, den = prob.density_num_den()
        sol = solve_flag_sdp(prob, num / den)
        print(f"# N=6 seed SDP: {sol['status']} b = {sol['b']:.12f}", flush=True)
        items = []
        for blk, Q in zip(prob.blocks, sol["Q"]):
            for k, (bi, t) in enumerate(self.m.keys):
                b = self.m.blocks[bi]
                if (b.s, b.m) == (blk.fam.s, blk.fam.m) and b.types[t] == int(blk.fam.sigma):
                    assert np.array_equal(b.flags[t], blk.fam.flags)
                    w, V = np.linalg.eigh(Q)
                    items += [(k, V[:, i]) for i in range(w.size) if w[i] > 1e-9 * max(w.max(), 1e-300)]
        return self.add_cuts(items)

    def set_stationarity(self, on):
        self.stat_on = bool(on)
        tm = self.a.tau_max if on else 0.0
        self.h.changeColCost(1, tm)
        self.h.changeColCost(2, tm)
        self.set_cap()

    def set_cap(self):
        self.h.changeColCost(0, self.b0 + self.a.b_cap if self.stat_on else 1.0)

    def run(self):
        a = self.a
        t0 = time.time()
        k0 = min(a.init_rows, self.U.size)
        self.add_graphs(self.U[np.argsort(self.eU, kind="stable")[:k0]])
        if a.resume:
            st = np.load(a.resume + ".state.npz")
            self.add_graphs(st["W"])
            ck = st["cutkeys"]
            items = []
            for k in range(self.nkeys):
                V = st[f"V{k}"]
                items += [(k, v) for v in V]
            self.add_cuts(items)
            if "b0" in st.files and a.b0 is None:
                self.set_b0(float(st["b0"]))
            print(json.dumps({"resumed_graphs": len(self.W), "resumed_cuts": len(self.cuts), "b0": self.b0}),
                  flush=True)
        elif a.seed_n6:
            print(json.dumps({"seed_n6_cuts": self.seed_n6()}), flush=True)
        self.set_stationarity(a.tau_max > 0 and (a.stat_after < 0 or (a.resume and a.resume_stat)))
        since_price, since_full = 10 ** 9, 10 ** 9
        status = "running"
        while True:
            self.it += 1
            t = time.time()
            self.h.run()
            t_lp = time.time() - t
            st = self.h.getModelStatus()
            if st != highspy.HighsModelStatus.kOptimal:
                status = f"lp_status_{self.h.modelStatusToString(st)}"
                break
            sol = self.h.getSolution()
            info_lp = self.h.getInfo()
            yall = np.array(sol.col_value)
            rd = np.array(sol.row_dual)
            b_lp = float(rd[0])
            tau1, tau2 = max(float(rd[1]), 0.0), max(float(rd[2]), 0.0)
            lam = np.maximum(rd[NFIX:], 0.0)
            y = np.maximum(yall[NFIX:], 0.0)
            t = time.time()
            Ys = self.moments(y)
            cand, mins = [], []
            for k, Y in enumerate(Ys):
                nev = min(a.per_key, self.dims[k])
                w, V = scipy.linalg.eigh(Y, subset_by_index=[0, nev - 1], driver="evr", check_finite=False)
                mins.append(float(w[0]))
                cand += [(float(w[i]), k, V[:, i]) for i in range(w.size) if w[i] < -a.psd_tol]
            cand.sort(key=lambda z: z[0])
            t_eig = time.time() - t
            row = {"it": self.it, "t": round(time.time() - t0, 1), "b_lp": b_lp,
                   "obj": float(info_lp.objective_function_value), "tau1": tau1, "tau2": tau2, "b0": self.b0,
                   "stat": self.stat_on, "z": float(yall[0]), "min_eig": min(mins), "graphs": len(self.W),
                   "cuts": len(self.cuts), "simplex_its": int(info_lp.simplex_iteration_count),
                   "lp_s": round(t_lp, 2), "eig_s": round(t_eig, 2)}
            added_g, b0_changed, full = 0, False, False
            if since_price >= a.price_every or not cand:
                full = since_full >= a.full_every or not cand
                t = time.time()
                added_g, info = self.price(lam, tau1, tau2, b_lp, full)
                if info["scan"] == "hot" and added_g == 0 and not cand:
                    added_g, info = self.price(lam, tau1, tau2, b_lp, True)
                full = info["scan"] == "full"
                row.update(info)
                row["price_s"] = round(time.time() - t, 2)
                since_price = 0
                since_full = 0 if full else since_full + 1
                if full:
                    bstar = info["b_star"]
                    if bstar > self.best["b_star"]:
                        self.best = {"b_star": bstar, "it": self.it, "t": round(time.time() - t0, 1), "b_lp": b_lp,
                                     "tau1": tau1, "tau2": tau2, "b0": self.b0}
                        self.save_best(lam, tau1, tau2, b_lp)
                    if self.stat_on and a.update_b0:
                        if bstar > self.b0 + a.b0_step or (added_g == 0 and not cand and bstar < self.b0 - a.b0_step):
                            self.set_b0(bstar)
                            self.set_cap()
                            b0_changed = True
                            row["b0_updated"] = bstar
            else:
                since_price += 1
            pruned = self.prune_cuts(lam) if self.it % a.prune_every == 0 else 0
            t = time.time()
            added_c = self.add_cuts([(k, v / np.linalg.norm(v)) for _, k, v in cand[:a.max_cuts]])
            row["addcut_s"] = round(time.time() - t, 2)
            row.update({"added_graphs": added_g, "added_cuts": added_c, "pruned": pruned, "rss_mb": round(rss_mb())})
            stalled = added_g == 0 and added_c == 0 and not b0_changed and full
            if not self.stat_on and a.tau_max > 0 and a.stat_after >= 0 and (
                    stalled or time.time() - t0 > a.stat_after):
                if np.isfinite(self.best["b_star"]):
                    self.set_b0(max(self.best["b_star"], self.b0))
                self.set_stationarity(True)
                row["stationarity_on"] = self.b0
                stalled = False
            self.hist.append(row)
            print(json.dumps(row), flush=True)
            if a.save_every and self.it % a.save_every == 0:
                self.save_state()
            if stalled:
                status = "converged"
                break
            if time.time() - t0 > a.seconds:
                status = "time_limit"
                break
            if self.it >= a.max_iter:
                status = "iteration_limit"
                break
        return status

    def save_best(self, lam, tau1, tau2, b_lp):
        if not self.a.out:
            return
        F = self.factors(lam)
        tmp = self.a.out + ".best.tmp.npz"                 # atomic replace (a Spot preemption may kill the writer)
        np.savez(tmp, tau1=tau1, tau2=tau2, b0=self.b0, b_star=self.best["b_star"], b_lp=b_lp,
                 keys=np.array(self.m.keys), **{f"F{k}": f for k, f in enumerate(F)},
                 blocks=np.array([(b.s, b.m) for b in self.m.blocks]), p=self.m.p)
        os.replace(tmp, self.a.out + ".best.npz")

    def save_state(self):
        if not self.a.out:
            return
        tmp = self.a.out + ".state.tmp.npz"
        np.savez(tmp, W=np.array(self.W, dtype=np.int64), b0=self.b0, stat=self.stat_on,
                 cutkeys=np.array([c[0] for c in self.cuts]),
                 **{f"V{k}": np.array([c[1] for c in self.cuts if c[0] == k]).reshape(-1, n)
                    for k, n in enumerate(self.dims)})
        os.replace(tmp, self.a.out + ".state.npz")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("p", type=int)
    ap.add_argument("--blocks", choices=["all", "std"], default="all")
    ap.add_argument("--blocks-list", default=None,
                    help="explicit block subset 's,m;s,m;...' (overrides --blocks; e.g. '5,6' or '3,5;5,6')")
    ap.add_argument("--sample", type=float, default=1.0, help="profiling: scan only a random fraction of classes")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--b0", type=float, default=None)
    ap.add_argument("--no-stationarity", action="store_true")
    ap.add_argument("--update-b0", type=int, default=1)
    ap.add_argument("--b0-step", type=float, default=1e-6)
    ap.add_argument("--tau-max", type=float, default=1000.0)
    ap.add_argument("--stat-after", type=float, default=float("inf"),
                    help="seconds before the stationarity rows are enabled (they are also enabled when the plain "
                         "LP-CUT-CG stalls); negative = from the start")
    ap.add_argument("--b-cap", type=float, default=2e-3, help="with stationarity: b <= b0 + cap")
    ap.add_argument("--hot-size", type=int, default=300000)
    ap.add_argument("--seed-n6", type=int, default=1)
    ap.add_argument("--full-every", type=int, default=2)
    ap.add_argument("--seconds", type=float, default=3600)
    ap.add_argument("--max-iter", type=int, default=100000)
    ap.add_argument("--init-rows", type=int, default=200)
    ap.add_argument("--rows-per-price", type=int, default=500)
    ap.add_argument("--price-every", type=int, default=1)
    ap.add_argument("--price-tol", type=float, default=1e-9)
    ap.add_argument("--per-key", type=int, default=12)
    ap.add_argument("--max-cuts", type=int, default=30)
    ap.add_argument("--psd-tol", type=float, default=1e-10)
    ap.add_argument("--keep-cuts", type=int, default=1500)
    ap.add_argument("--prune-every", type=int, default=5)
    ap.add_argument("--prune-age", type=int, default=10)
    ap.add_argument("--lp-tol", type=float, default=1e-9)
    ap.add_argument("--lp-solver", default="choose", choices=["choose", "simplex", "ipm", "pdlp"])
    ap.add_argument("--lp-threads", type=int, default=1, help="HiGHS threads (> 1 turns on parallel dual simplex)")
    ap.add_argument("--simplex-strategy", type=int, default=-1, help="HiGHS: 1 dual, 4 primal; -1 default")
    ap.add_argument("--nchunks", type=int, default=64)
    ap.add_argument("--out", default=None, help="prefix for .best.npz/.state.npz/.json")
    ap.add_argument("--save-every", type=int, default=10, help="write .state.npz every k iterations")
    ap.add_argument("--row-scale", type=int, default=1, help="normalise each cut row to max |coefficient| 1")
    ap.add_argument("--resume", default=None, help="prefix of a previous run: reload its working set and cuts")
    ap.add_argument("--resume-stat", type=int, default=0, help="with --resume: stationarity on from the start")
    ap.add_argument("--init-best", default=None, help="prefix whose .best.npz b_star is the initial best (resume)")
    a = ap.parse_args()
    t0 = time.time()
    if a.blocks_list:
        sm = tuple(tuple(int(x) for x in t.split(",")) for t in a.blocks_list.split(";"))
        assert all(blk in ALL_BLOCKS for blk in sm), sm
        sm = tuple(blk for blk in ALL_BLOCKS if blk in sm)          # canonical order
    else:
        sm = ALL_BLOCKS if a.blocks == "all" else STD_BLOCKS
    model = build_model(a.p, sm)
    cls = np.load(os.path.join(DATA, f"classes_p{a.p}.npy")).astype(np.int64)
    if a.sample < 1.0:
        rng = np.random.default_rng(a.seed)
        low = cls[np.argsort(popcount_array(cls), kind="stable")[:5000]]   # keep the sparsest graphs in every sample
        cls = np.unique(np.concatenate([low, rng.choice(cls, size=int(cls.size * a.sample), replace=False)]))
    b0 = a.b0 if a.b0 is not None else {5: 0.25, 6: 0.1, 7: 0.05}[a.p]
    if a.no_stationarity:
        a.tau_max = 0.0
    print(f"# lpcg7 p={a.p} blocks={a.blocks} universe {cls.size:,} Gram dims {model.gram_dims()} "
          f"b0={b0} setup {time.time() - t0:.1f}s", flush=True)
    cg = LPCG7(model, cls, b0, a)
    if a.init_best and os.path.exists(a.init_best + ".best.npz"):
        # keep a previous run's best (e.g. before a preemption): .best.npz is overwritten only by a better b_star
        prev = np.load(a.init_best + ".best.npz")
        cg.best = {"b_star": float(prev["b_star"]), "from": a.init_best}
        print(json.dumps({"init_best": cg.best}), flush=True)
    status = cg.run()
    cg.save_state()
    out = {"p": a.p, "args": vars(a), "status": status, "best": cg.best, "iterations": cg.it,
           "graphs": len(cg.W), "cuts": len(cg.cuts), "universe": int(cls.size), "time_s": time.time() - t0,
           "peak_rss_mb": rss_mb(), "history": cg.hist}
    if a.out:
        with open(a.out + ".json", "w") as f:
            json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "history"}), flush=True)


if __name__ == "__main__":
    main()
