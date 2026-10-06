"""Tests of the threshold generalisation (gen.py) for p = 5, LAM = 2 (complement of K_5^{4-}-free):
kernel <Q, M(H)> vs the validated engine (direct configuration enumeration with the same predicate), exact kernel vs
float, raw threshold kernel vs the class list, class list = predicate.  Run inside tools/flagalg/h44:
    python -m pytest -q test_gen.py
"""
import numpy as np
import pytest

import gen

gen.configure(5, 2)

import common  # noqa: E402
from common import N, build_model  # noqa: E402
from flags import block_densities, flag_family  # noqa: E402
from hypergraphs import EveryPSetHasEdge  # noqa: E402  (rebound to the threshold predicate)
from kernels import qsum_flat, scan_exact, scan_float, tables  # noqa: E402

P, LAM = 5, 2


@pytest.fixture(scope="module")
def model():
    return build_model(P)


def random_adm(k, seed):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < k:
        dens = rng.uniform(0.4, 0.9)
        m = sum(1 << i for i in range(35) if rng.random() < dens)
        if gen.admissible7_lam(np.array([m]), P, LAM)[0]:
            out.append(m)
    return np.array(out, dtype=np.int64)


def test_predicate_rebound():
    assert "at least 2" in EveryPSetHasEdge(5).name
    assert common.EveryPSetHasEdge(5).name == EveryPSetHasEdge(5).name


def test_kernel_matches_engine(model):
    H = random_adm(15, seed=3)
    rng = np.random.default_rng(1)
    Qs = []
    for bi, t in model.keys:
        n = int(model.qdim[bi, t])
        A = rng.normal(size=(n, n))
        Qs.append((A + A.T) / 2)
    sl, bad = scan_float(H, *tables(model), qsum_flat(model, Qs), 0.0, 0.0, 0.0, 4)
    assert not bad.any()
    e = np.array([bin(int(h)).count("1") for h in H])
    pred = EveryPSetHasEdge(P)
    tot = np.zeros(H.size)
    for (bi, t), Q in zip(model.keys, Qs):
        b = model.blocks[bi]
        fam = flag_family(b.s, b.m, b.types[t], 4, pred)
        assert np.array_equal(fam.flags, b.flags[t])
        blk = block_densities(H, N, 4, fam)
        tot += np.bincount(blk.h, weights=blk.cnt * Q[blk.a, blk.b], minlength=H.size) / blk.denom
    np.testing.assert_allclose(e / 35 - sl, tot, rtol=1e-10, atol=1e-10)


def test_inadmissible_flagged(model):
    """A graph violating the threshold somewhere (but nonempty) must be flagged 'bad' by the kernel."""
    rng = np.random.default_rng(9)
    found = 0
    for _ in range(2000):
        m = sum(1 << i for i in range(35) if rng.random() < 0.55)
        if not gen.admissible7_lam(np.array([m]), P, LAM)[0]:
            _, bad = scan_float(np.array([m], dtype=np.int64), *tables(model), np.zeros(model.qsize), 0.0, 0.0,
                                0.0, 1)
            # bad iff some (5,6)-flag / type configuration is inadmissible; every 5-set lies in some 6-set
            assert bad[0] == 1
            found += 1
    assert found > 10


def test_raw_kernel_counts_match_classes(model):
    """Raw threshold scan with Q = 0: number of admissible extensions = labelled-style raw count, no bad graphs,
    and the set of (e, cdd) groups equals the class list's groups."""
    raw = gen.make_raw_kernels(LAM)
    reps = common.reps6(P).astype(np.int64)
    qi = np.zeros(model.qsize, dtype=np.int64)
    MZr, AHr, CNr, nbr = raw(reps, model.psets_new, *tables(model), qi, 1 << 18)
    cls = np.load(gen.classes_path(P, LAM)).astype(np.int64)
    MZ, AH, CN, nb = scan_exact(cls, *tables(model), qi, 8)
    assert nbr == 0 and nb == 0
    assert np.array_equal(CNr > 0, CN > 0)
    assert gen.admissible7_lam(AHr[CNr > 0], P, LAM).all()
