# search/ -- the code that found the certificates (not needed to verify them)

This directory holds the producer code (the "first implementation" of the paper), for transparency. **Nothing in
it is needed to verify the results**: a certificate is valid because its matrices are positive semidefinite and
its value is the exact minimum over all admissible graphs, and both facts are re-checked from the certificate files
alone by the independent checkers in `../verifier/`, which share no code with this directory. How a certificate
was found plays no role in its correctness.

The code is provided as it was run, with the sanitising edits listed in `CHANGES.txt`. It expects input files that
are not included (class lists of all admissible 7-vertex graphs, about 137 MB, and the saved states of the linear
programs), so most scripts will not run out of the box. Python dependencies, beyond numpy: numba, highspy
(HiGHS), clarabel, scipy and pynauty.

## Map

* `flagalg/` -- general flag-algebra engine for small N (types, flags, densities: `flags.py`, `hypergraphs.py`;
  SDP with Clarabel: `sdp.py`; exact rounding: `exact.py`; Turán problems in complement form: `turan.py`, which
  also wrote the 5-graph certificates of `../certificates/fivegraphs/` and the sharp six-vertex certificates;
  `verify_cert.py`, the first implementation's exact checker of such certificates, run by
  `verify_all.sh --producer-n6`; unit tests `test_flagalg.py`, validation ladder `run_validation.py`). `lottery.py`,
  `check_l1.py`, `lp_polish.py`, `estimate_n7.py` and `solver_experiment.py` are earlier experiments of the same
  project.
* `flagalg/n7/` -- the N = 7, r = 4 machinery for t(p,4), p = 5, 6, 7: class lists by one-vertex extension with
  nauty canonical forms (`enum7.py`), tables and numba kernels (`common.py`, `kernels.py`), the LP cutting-plane and
  column-generation solver following Jeong et al. (`lpcg7.py`), exact rounding and exact scans (`certify7.py`), a
  sampling checker (`verify7.py`), tests (`test_n7.py`) and a throughput profile (`profile7.py`).
* `flagalg/h44/` -- the same machinery for "every p-set spans at least lambda edges" (`gen.py`, `enum_gen.py`,
  `lpcg_gen.py`, `certify_gen.py`, `verify_gen.py`, `test_gen.py`), used for K6_4 (p = 6, lambda = 1), K5_4minus
  (p = 5, lambda = 2), K6_4minus (p = 6, lambda = 2) and the ten cases of `../certificates/catalogue/` (the files
  here are those that ran on the cloud machines for the catalogue and for these three certificates);
  `plugin_check.py` evaluates a certificate on Sidorenko's construction.
* `flagalg/n6control/` -- the six-vertex control: SDP, sharp rounding and the rational dual points
  (`n6_control.py`) and the producer's own exact check of the dual points (`dual_check.py`).
* `flagalg/n7dual/` -- the search for the seven-vertex dual point of Theorem n7opt: re-solve of the restricted LP
  (`lp_resolve.py`), facial reduction (`facial.py`), dual-SDP data (`prep_support.py`), barrier and primal-dual
  interior-point solvers (`barrier.py`, `pdip.py`), rounding to a dyadic point (`round_dual.py`), the producer's
  exact check (`verify_dual.py`, `crosscheck_engine.py`, `make_negative_tests.py`), and explorations
  (`analyse_support.py`, `lp_subset.py`, `near_tight.py`, `odd_classes.py`, `odd_filter.py`, `giraud_exact.py`).
  Floating-point steps here certify nothing; only the exact check of the final rational point matters, and it is
  repeated independently by `../verifier/rv1_dual.py`. The dual points of `../certificates/n7_dual_points/` were
  found by a second version of this search (generalised to "every p-set spans at least lambda edges" and to a
  given block list), which is not included; their exact check is `../verifier/rv2_dual.py`. Likewise the
  custom-objective pipeline and the sharp rounding that produced `../certificates/sigma_K5_4/` are not included;
  its exact checks are `../verifier/rs_*.py` and `rs_eval.c`.

All of this code is under the MIT license (`../LICENSE-CODE`).
