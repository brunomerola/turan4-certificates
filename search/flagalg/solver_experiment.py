"""Clarabel settings experiment on the K4^3 complement Turan SDP at N = 6 (accuracy vs settings)."""
import sys
import time

import numpy as np

from flags import build_problem
from hypergraphs import EveryPSetHasEdge
from sdp import solve_flag_sdp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
prob = build_problem(N, 3, EveryPSetHasEdge(4))
num, den = prob.density_num_den()
d = num / den
trials = [
    ("default-tol1e-8", 1e-8, {}),
    ("tol1e-9", 1e-9, {}),
    ("tol1e-10", 1e-10, {}),
    ("tol1e-10 no-equil", 1e-10, {"equilibrate_enable": False}),
    ("tol1e-10 maxiter1000 refine", 1e-10, {"iterative_refinement_reltol": 1e-14, "iterative_refinement_abstol": 1e-14,
                                           "iterative_refinement_max_iter": 50, "max_iter": 1000}),
    ("tol1e-10 static_reg 1e-10", 1e-10, {"static_regularization_constant": 1e-10}),
    ("tol1e-10 max_step 0.95", 1e-10, {"max_step_fraction": 0.95}),
]
for name, tol, st in trials:
    t = time.time()
    sol = solve_flag_sdp(prob, d, tol=tol, settings=st)
    print(f"{name:35s} {sol['status']:14s} b={sol['b']:.12f} dual={-sol['obj_dual']:.12f} "
          f"gap={-sol['obj_dual'] - sol['b']:.2e} it={sol['iterations']} rp={sol['r_prim']:.1e} "
          f"rd={sol['r_dual']:.1e} t={time.time() - t:.1f}s", flush=True)
