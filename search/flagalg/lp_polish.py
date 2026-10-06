"""High-accuracy polish of a flag-algebra SDP by LP column generation with HiGHS (imports highspy; never import
ortools in this process).

Idea (Ahmadi-Dash-Hall 2017 "DSOS/SDSOS + column generation"; used for flag algebras by Jeong-Park-Im-Lee-Yang 2026,
arXiv 2609.27495, github.com/taeyool/tetrahedron-turan, Apache-2.0 -- this file is an independent re-implementation
of the same idea, no code copied):  restrict every PSD block to  Q = sum_v lambda_v v v^T,  lambda >= 0, over a
finite set of vectors v.  The problem becomes an LP whose every feasible point is SDP-feasible (PSD by construction)
and which HiGHS solves to ~1e-9 accuracy.  Pricing: with LP duals y_H (sum = 1) the reduced cost of a column v of
block sigma is v^T Y_sigma v, Y_sigma = sum_H y_H M_sigma(H); add the eigenvectors of the most negative eigenvalues.
At convergence (all Y_sigma, Z PSD up to tolerance) the LP value is the SDP value.  The initial vectors are all
eigenvectors of the Clarabel solution.

Usage:
  python lp_polish.py turan --r 3 --p 4 --N 6
  python lp_polish.py lottery --k 4 --r 3 --p 4 --N 6 --m 5 --tau-sos
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from fractions import Fraction
from math import comb

import highspy
import numpy as np

from sdp import solve_flag_sdp

INF = highspy.kHighsInf


class LPCG:
    def __init__(self, prob, d, Rrho=None, Crow=None, ta=None, tb=None, n1=0, verbose=True, link_perms=None):
        """Rows: H rows  b - Rrho c + sum lambda v^T M(H) v <= d(H);  packing rows  Crow c + sum mu w_a w_b <= 1."""
        self.prob, self.verbose = prob, verbose
        self.nH = prob.nH
        self.nc = 0 if Rrho is None else Rrho.shape[1]
        self.nT = 0 if Crow is None else Crow.shape[0]
        self.ta, self.tb, self.n1 = ta, tb, n1
        self.link_perms = link_perms      # symmetric mode: W columns are group averages of w w^T
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("primal_feasibility_tolerance", 1e-10)
        h.setOptionValue("dual_feasibility_tolerance", 1e-10)
        nrow = self.nH + self.nT
        lower = np.full(nrow, -INF)
        upper = np.concatenate([np.asarray(d, dtype=float), np.ones(self.nT)])
        h.addRows(nrow, lower, upper, 0, np.zeros(nrow, dtype=np.int32), np.zeros(0, dtype=np.int32),
                  np.zeros(0))
        self.h = h
        # b
        self._add_cols(np.array([-1.0]), np.array([-INF]), np.array([INF]),
                       [(np.arange(self.nH), np.ones(self.nH))])
        # c
        if self.nc:
            Rc = Rrho.tocsc()
            Cc = Crow.tocsc() if self.nT else None
            cols = []
            for j in range(self.nc):
                ri = Rc.indices[Rc.indptr[j]:Rc.indptr[j + 1]]
                rv = -Rc.data[Rc.indptr[j]:Rc.indptr[j + 1]]
                ci = Cc.indices[Cc.indptr[j]:Cc.indptr[j + 1]] + self.nH
                cv = Cc.data[Cc.indptr[j]:Cc.indptr[j + 1]]
                cols.append((np.concatenate([ri, ci]), np.concatenate([rv, cv])))
            self._add_cols(np.zeros(self.nc), np.zeros(self.nc), np.full(self.nc, INF), cols)
        self.vecs = [[] for _ in prob.blocks]     # list of vectors per block
        self.colidx = [[] for _ in prob.blocks]
        self.wvecs, self.wcol = [], []

    def _add_cols(self, cost, lo, hi, cols):
        n0 = self.h.getNumCol()
        starts, idx, val = [], [], []
        k = 0
        for ri, rv in cols:
            starts.append(k)
            keep = np.abs(rv) > 0
            idx.append(np.asarray(ri)[keep].astype(np.int32))
            val.append(np.asarray(rv, dtype=float)[keep])
            k += int(keep.sum())
        self.h.addCols(len(cols), cost, lo, hi, k, np.array(starts, dtype=np.int32),
                       np.concatenate(idx) if idx else np.zeros(0, np.int32),
                       np.concatenate(val) if val else np.zeros(0))
        return list(range(n0, n0 + len(cols)))

    def add_block_vectors(self, bi, V):
        blk = self.prob.blocks[bi]
        cols = []
        for v in V:
            v = v / np.linalg.norm(v)
            w = blk.cnt / blk.denom * v[blk.a] * v[blk.b]
            coef = np.bincount(blk.h, weights=w, minlength=self.nH)
            nz = np.nonzero(np.abs(coef) > 1e-15)[0]
            cols.append((nz, coef[nz]))
            self.vecs[bi].append(v)
        n = len(cols)
        self.colidx[bi] += self._add_cols(np.zeros(n), np.zeros(n), np.full(n, INF), cols)

    def wcol_matrix(self, w):
        if self.link_perms is None:
            return np.outer(w, w)
        acc = np.zeros((self.n1, self.n1))
        for img in self.link_perms:
            wp = w[np.argsort(img)]
            acc += np.outer(wp, wp)
        return acc / len(self.link_perms)

    def add_w_vectors(self, V):
        cols = []
        for w in V:
            w = w / np.linalg.norm(w)
            coef = self.wcol_matrix(w)[self.ta, self.tb]
            nz = np.nonzero(np.abs(coef) > 1e-15)[0]
            cols.append((nz + self.nH, coef[nz]))
            self.wvecs.append(w)
        n = len(cols)
        self.wcol += self._add_cols(np.zeros(n), np.zeros(n), np.full(n, INF), cols)

    def solve(self):
        self.h.run()
        sol = self.h.getSolution()
        x = np.array(sol.col_value)
        rd = np.array(sol.row_dual)
        return x, rd, str(self.h.getModelStatus())

    def dual_matrices(self, rd):
        y = -rd[:self.nH]
        Ys = []
        for blk in self.prob.blocks:
            Y = np.zeros((blk.size, blk.size))
            np.add.at(Y, (blk.a, blk.b), y[blk.h] * blk.cnt / blk.denom)
            Ys.append(Y)
        Z = None
        if self.n1:
            z = -rd[self.nH:]
            Z = np.zeros((self.n1, self.n1))
            np.add.at(Z, (self.ta, self.tb), z * np.where(self.ta == self.tb, 1.0, 0.5))
            np.add.at(Z, (self.tb, self.ta), z * np.where(self.ta == self.tb, 0.0, 0.5))
            if self.link_perms is not None:   # reduced cost of a symmetrised column = w^T Z_sym w
                Z = sum(Z[np.ix_(img, img)] for img in self.link_perms) / len(self.link_perms)
        return y, Ys, Z

    def run(self, Q0, W0=None, iters=40, per_block=6, tol=1e-10):
        for bi, Q in enumerate(Q0):
            _, V = np.linalg.eigh(Q)
            self.add_block_vectors(bi, list(V.T))
        if self.n1 and W0 is not None:
            _, V = np.linalg.eigh(W0)
            self.add_w_vectors(list(V.T))
        hist = []
        for it in range(iters):
            t = time.time()
            x, rd, st = self.solve()
            b = float(x[0])
            y, Ys, Z = self.dual_matrices(rd)
            mins = []
            added = 0
            for bi, Y in enumerate(Ys):
                w, V = np.linalg.eigh(Y)
                mins.append(float(w[0]))
                neg = [V[:, i] for i in range(min(per_block, w.size)) if w[i] < -tol]
                if neg:
                    self.add_block_vectors(bi, neg)
                    added += len(neg)
            zmin = None
            if Z is not None:
                w, V = np.linalg.eigh(Z)
                zmin = float(w[0])
                neg = [V[:, i] for i in range(min(per_block, w.size)) if w[i] < -tol]
                if neg:
                    self.add_w_vectors(neg)
                    added += len(neg)
            hist.append({"it": it, "status": st, "b": b, "sum_y": float(y.sum()), "min_eig_Y": min(mins),
                         "min_eig_Z": zmin, "added": added, "t": time.time() - t})
            if self.verbose:
                print(f"  LP-CG it {it}: {st} b = {b:.13f}  sum y = {y.sum():.6f}  min eig Y = {min(mins):.2e}"
                      f"  min eig Z = {zmin}  added {added}  ({time.time() - t:.1f}s)", flush=True)
            if added == 0:
                break
            if len(hist) >= 4 and hist[-1]["b"] - hist[-4]["b"] < 1e-12:   # stalled at LP precision
                break
        x, rd, st = self.solve()
        return self.extract(x), hist

    def extract(self, x):
        Qs = []
        for bi, blk in enumerate(self.prob.blocks):
            Q = np.zeros((blk.size, blk.size))
            for v, ci in zip(self.vecs[bi], self.colidx[bi]):
                if x[ci] > 0:
                    Q += x[ci] * np.outer(v, v)
            Qs.append(Q)
        W = None
        if self.n1:
            W = np.zeros((self.n1, self.n1))
            for w, ci in zip(self.wvecs, self.wcol):
                if x[ci] > 0:
                    W += x[ci] * self.wcol_matrix(w)
        return {"b": float(x[0]), "c": np.maximum(x[1:1 + self.nc], 0.0), "Q": Qs, "W": [W] if W is not None else []}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("problem", choices=["turan", "lottery"])
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--p", type=int, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--tau-sos", action="store_true")
    ap.add_argument("--tau-sym", action="store_true")
    ap.add_argument("--per-block", type=int, default=6)
    ap.add_argument("--iters", type=int, default=40)
    ap.add_argument("--bits", type=int, default=40)
    ap.add_argument("--out", default=None)
    ap.add_argument("--cert", default=None)
    a = ap.parse_args()
    t0 = time.time()
    out = {"args": vars(a)}
    if a.problem == "turan":
        from exact import certify_eig
        from flags import build_problem
        from hypergraphs import EveryPSetHasEdge
        prob = build_problem(a.N, a.r, EveryPSetHasEdge(a.p))
        num, den = prob.density_num_den()
        d = num / den
        sol0 = solve_flag_sdp(prob, d)
        print(f"# clarabel: {sol0['status']} b = {sol0['b']:.13f} dual = {-sol0['obj_dual']:.13f}", flush=True)
        cg = LPCG(prob, d)
        sol, hist = cg.run(sol0["Q"], iters=a.iters)
        ce = certify_eig(prob, [Fraction(int(x), den) for x in num], sol["Q"], bits=a.bits)
        out.update({"clarabel_b": sol0["b"], "clarabel_dual": -sol0["obj_dual"], "lp_b": sol["b"], "history": hist,
                    "b_exact": str(ce["b_exact"]), "b_exact_float": float(ce["b_exact"]),
                    "psd": [list(x) for x in ce["psd"]]})
        print(f"# LP-CG b = {sol['b']:.13f}; exact (eig rounding) = {float(ce['b_exact']):.13f}", flush=True)
        if a.cert:
            from turan import write_cert
            write_cert(a.cert, prob, ce["Qexact"], ce["b_exact"], a.r, a.p, a.N)
    else:
        from lottery import build_lottery, certify_lottery
        ctx = build_lottery(a.k, a.r, a.p, a.N, a.m, a.tau_sos or a.tau_sym, tau_sym=a.tau_sym)
        prob = ctx["prob"]
        Rrho = ctx["Rm"] / ctx["rden"]
        Crow = ctx["Cm"].astype(float) / ctx["cdiv"]
        sol0 = solve_flag_sdp(prob, np.zeros(prob.nH), R=Rrho, G=ctx["Gm"], h=np.ones(ctx["nT"]),
                              extra_psd=ctx["extra"], psd_maps=ctx.get("psd_maps", ()))
        print(f"# clarabel: {sol0['status']} b = {sol0['b']:.13f} dual = {-sol0['obj_dual']:.13f}", flush=True)
        num, den = prob.density_num_den()
        solT = solve_flag_sdp(prob, num / den)
        cgT = LPCG(prob, num / den, verbose=False)
        solTp, _ = cgT.run(solT["Q"], iters=a.iters)
        bate = solTp["b"] / comb(a.k, a.r)
        print(f"# Turan same N (LP-CG polished): {solTp['b']:.13f}; Bate = {bate:.13f}", flush=True)
        cg = LPCG(prob, np.zeros(prob.nH), Rrho=Rrho, Crow=Crow, ta=ctx["ta"], tb=ctx["tb"], n1=ctx["n1"],
                  link_perms=ctx.get("link_perms"))
        sol, hist = cg.run(sol0["Q"], sol0["W"][0] if sol0["W"] else None, iters=a.iters, per_block=a.per_block)
        cl = certify_lottery(ctx, sol, a.bits, a.cert)
        out.update({"clarabel_b": sol0["b"], "clarabel_dual": -sol0["obj_dual"], "lp_b": sol["b"],
                    "turan_lp_b": solTp["b"], "bate_same_N": bate, "gain_lp": sol["b"] - bate, "history": hist,
                    "b_exact": str(cl["b_exact"]), "b_exact_float": float(cl["b_exact"]),
                    "alpha": float(cl["alpha"]), "psd": [list(x) for x in cl["psd"]], "W_psd": cl["w_psd"]})
        print(f"# LP-CG b = {sol['b']:.13f} (gain over Bate {sol['b'] - bate:.3e}); exact = "
              f"{float(cl['b_exact']):.13f}", flush=True)
    out["time_s"] = time.time() - t0
    out["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    if a.out:
        with open(a.out, "w") as f:
            json.dump(out, f, indent=1, default=str)


if __name__ == "__main__":
    main()
