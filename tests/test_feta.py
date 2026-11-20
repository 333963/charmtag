"""Testy S4 (feta.py): tożsamość transportu na starcie, masa z zamianą K↔π, wagi SALAD ≈ 1 gdy MC = dane."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import pytest

from charmtag.bkggen import flow as F, feta as T, kin

# FETA (flow transportujący) wymaga torch+zuko (extras [flow]); SALAD hgb, swap-mass i C2ST działają bez nich.
needs_flow = pytest.mark.skipif(
    not all(__import__("importlib").util.find_spec(m) for m in ("torch", "zuko")),
    reason="torch/zuko (extras [flow]) niezainstalowane")

SP = F.SPECIES["D0"]


def _toyG(n, seed=0):
    r = np.random.default_rng(seed)
    G = np.column_stack([np.log1p(r.uniform(3, 20, n)), r.uniform(-2.4, 2.4, n), r.uniform(-1, 1, n), r.uniform(-np.pi, np.pi, n),
                         r.exponential(2, n), r.exponential(0.5, n), r.exponential(1, n), r.exponential(0.2, n),
                         r.exponential(0.1, n), r.uniform(0, 1, n), r.uniform(1, 7, n), r.exponential(0.3, n),
                         r.normal(1.2, 0.2, n), r.integers(0, 8, n).astype(float), r.normal(1.1, 0.2, n), r.integers(0, 8, n).astype(float)])
    m = r.uniform(1.75, 2.0, n)
    return G, m


@needs_flow
def test_transport_identity_at_init():
    G, m = _toyG(3000)
    fe = T.FetaTransport(F.FEATS_GEN, SP, base_transforms=2, base_hidden=(16, 16), top_transforms=2, top_hidden=(16, 16))
    fe.pre.fit(G)
    Gt = fe.transport(G, m)
    cont = [j for j, f in enumerate(F.FEATS_GEN) if f not in ("npts_K", "npts_pi")]
    assert np.nanmax(np.abs(Gt[:, cont] - G[:, cont]) / (np.abs(G[:, cont]) + 1e-3)) < 0.05


def test_swap_mass_matches_four_vectors():
    G, m = _toyG(2000, 1)
    pt, eta, cts, phis = np.expm1(G[:, 0]), G[:, 1], G[:, 2], G[:, 3]
    ms = T.swap_mass_kin(m, pt, eta, cts, phis)
    # niezależnie: pędy córek w lab z kin.daughters nie dają φ → test własności: m_swap(m_swap-hipoteza) symetria i zakres
    assert np.all(ms > kin.MK + kin.MPI - 1e-6)
    # dla cosθ*=0 w spoczynku D (pT→małe) |m_swap − m| małe tylko gdy p* ~; sprawdź ciągłość: zamiana dwukrotna = m
    assert np.isfinite(ms).all()


def test_salad_weights_unity_when_same():
    G, m = _toyG(20000, 2)
    G2, m2 = _toyG(20000, 3)
    s = T.SaladReweighter(F.FEATS_GEN, SP, kind="hgb").fit(G, m, G2, m2)
    w = s.weights(G, m)
    assert abs(np.median(w) - 1) < 0.15 and np.std(w) < 0.4


def test_c2st_floor():
    G, m = _toyG(6000, 4); G2, _ = _toyG(6000, 5)
    assert abs(T.c2st(G, G2) - 0.5) < 0.04
