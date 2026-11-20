"""Testy scalenia AA: jedno API bkggen, jeden c2st, dispatche CLI, placeholdery, hak V."""
import importlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from charmtag.bkggen import api
from charmtag.bkggen.c2st import c2st as _c2st
from charmtag.bkggen import flow as F2


def _synth_sb(n=6000, seed=0, species="D0"):
    """SB syntetyczne w cechach FEATS_GEN (+m) — tylko mechanika, nie fizyka."""
    rng = np.random.default_rng(seed)
    sp = dict(F2.SPECIES[species])
    lo, hi = sp["mass_range"]
    m = rng.uniform(lo, hi, 3 * n)
    m = m[np.abs(m - sp["m0"]) > sp["sb_min"]][:n]
    d = dict(log_pt=np.log1p(rng.uniform(3, 20, n)), eta=rng.uniform(-2, 2, n), cts=rng.uniform(-1, 1, n),
             phis=rng.uniform(-np.pi, np.pi, n), chi2=rng.exponential(2, n), Lxy=rng.exponential(0.5, n),
             impact_sig=rng.exponential(1, n), abs_impact_z0=rng.exponential(0.2, n), d_max_abs_d0=rng.exponential(0.05, n),
             r_d0=rng.uniform(0, 1, n), pointing=rng.exponential(1, n), d_max_abs_z0=rng.exponential(0.3, n),
             dedx_K=rng.normal(1.1, 0.2, n), npts_K=rng.integers(1, 9, n).astype(float),
             dedx_pi_max=rng.normal(1.0, 0.2, n), npts_pi=rng.integers(1, 9, n).astype(float))
    df = pd.DataFrame(d)
    df["m"] = m
    return df, sp


def test_registry_lists_all_adapters():
    assert {"copula", "flow", "s_copula", "s_flow", "salad", "candrot"} <= set(api.available())


def test_s2_copula_fit_sample_derive_H():
    SB, sp = _synth_sb()
    g = api.make_generator("copula", sp=sp).fit(SB)
    m = g.mass_model.sample(500, region="SR", rng=np.random.default_rng(1))
    X = g.sample(m)
    assert list(X.columns[:len(F2.FEATS_H)]) == F2.FEATS_H
    assert len(X) == 500 and np.all(np.abs(X["m"] - sp["m0"]) < sp["sr_half"])
    # więzy: abs_impact_sig == impact_sig (IPsig ≥ 0), d_min ≤ d_max
    assert np.all(X["d_min_abs_d0"] <= X["d_max_abs_d0"] + 1e-12)
    assert np.allclose(X["abs_impact_sig"], np.abs(X["impact_sig"]))


def test_s_copula_adapter_outputs_H_subset():
    SB, sp = _synth_sb()
    G = F2.derive_H(SB[F2.FEATS_GEN].to_numpy(float), SB["m"].to_numpy(float))
    SBH = pd.DataFrame(G, columns=F2.FEATS_H)
    SBH["m"] = SB["m"].to_numpy()
    g = api.make_generator("s_copula", feats=F2.FEATS_H, sp=sp).fit(SBH)
    X = g.sample(np.full(200, sp["m0"]))
    assert len(X) == 200 and set(F2.FEATS_H) <= set(X.columns)


def test_sample_before_fit_raises():
    _, sp = _synth_sb()
    with pytest.raises(RuntimeError):
        api.make_generator("copula", sp=sp).sample([sp["m0"]])


def test_c2st_single_implementation_and_forwarding():
    rng = np.random.default_rng(0)
    A = rng.normal(0, 1, (3000, 4))
    B = rng.normal(0, 1, (3000, 4))
    a1 = _c2st(A, B, seed=0)
    assert 0.45 < a1 < 0.56
    # validate.c2st_auc i flow.auc_two_samples → ta sama implementacja (identyczny wynik przy tym samym seed)
    from charmtag.bkggen import validate as V
    fa = [f"f{i}" for i in range(4)]
    a2 = V.c2st_auc(pd.DataFrame(A, columns=fa), pd.DataFrame(B, columns=fa), fa, seed=0, folds=3, n_max=150_000)
    a3 = F2.auc_two_samples(A, B, seed=0, n_folds=3, max_n=150_000)
    assert a1 == pytest.approx(a2, abs=1e-12) and a1 == pytest.approx(a3, abs=1e-12)
    # stary silnik S dostępny jawnie
    a4 = V.c2st_auc(pd.DataFrame(A, columns=fa), pd.DataFrame(B, columns=fa), fa, engine="xgb")
    assert 0.4 < a4 < 0.6


def test_c2st_detects_shift():
    rng = np.random.default_rng(1)
    A = rng.normal(0, 1, (3000, 3))
    B = rng.normal(0.8, 1, (3000, 3))
    assert _c2st(A, B) > 0.65


def test_cli_dispatch_has_no_duplicates():
    import inspect
    from charmtag import cli
    # REPO: dispatch tabelą _SUBCOMMANDS (bez łańcucha if); każda podkomenda raz, alias `candrot` usunięty (= bkggen --method candrot)
    assert set(cli._SUBCOMMANDS) == {"pairs", "fit2d", "syst", "s5compare", "bkggen", "s8salad"}
    assert "candrot" not in cli._SUBCOMMANDS
    assert inspect.getsource(cli.main).count("_SUBCOMMANDS.get") == 1


def test_cli_bkggen_help_and_candrot_choice():
    r = subprocess.run([sys.executable, "-m", "charmtag.cli", "bkggen", "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "candrot" in r.stdout


def test_hand_control_importable_and_no_variants_ext():
    from charmtag import hand
    assert callable(hand.score) and "lxy_sig_t" in hand.PREPARED   # ręczny dyskryminant (kontrola w validate)
    assert importlib.util.find_spec("charmtag.variants_ext") is None  # REPO: Hdec odrzucony (X-val) → moduł usunięty


def test_dstar_refit_hook_importable():
    DR = importlib.import_module("charmtag.dstar_refit")
    assert callable(DR.control_row) and callable(DR.fit_passfail)
    import inspect
    from charmtag import pipeline
    assert "dstar_refit" in inspect.getsource(pipeline.dstar_control)
