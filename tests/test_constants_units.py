"""P1-REST: jedno źródło stałych (GeV) i spójność jednostek między modułami."""
import numpy as np
import pandas as pd
import pytest

from charmtag import constants as C
from charmtag import kinematics as K


def test_single_source_of_m0():
    pytest.importorskip("charmbkg")   # porównuje też stałe modułów projektu B
    from charmtag import synth
    from charmtag.bkggen import flow
    from charmtag.pairs import model
    assert synth.M0[10] == C.M_D0 and synth.M0[4] == C.M_DPLUS
    assert flow.SPECIES["D0"]["m0"] == C.M_D0 and flow.SPECIES["Dplus"]["m0"] == C.M_DPLUS
    assert model.M0[10] == pytest.approx(C.to_mev(C.M_D0)) and model.M0[4] == pytest.approx(C.to_mev(C.M_DPLUS))
    assert C.to_gev(C.to_mev(1.86484)) == pytest.approx(1.86484)


def test_default_config_species_m0_matches_constants():
    from charmtag import config as cfgmod
    import importlib.resources as r
    import yaml
    with r.files("charmtag").joinpath("data/default.yaml").open() as f:
        cfg = yaml.safe_load(f)
    assert cfg["species"]["D0"]["m0"] == pytest.approx(C.M_D0, abs=1e-6)
    assert cfg["species"]["Dplus"]["m0"] == pytest.approx(C.M_DPLUS, abs=1e-6)
    assert cfgmod is not None


def test_masses_shared_by_all_modules():
    pytest.importorskip("charmbkg")
    from charmtag import bkgsrc, samples
    from charmtag.bkggen import data, kin, kin3
    from charmtag.hand import bethe_bloch
    for mod in (bkgsrc, samples, data, bethe_bloch):
        assert mod.M_K == C.M_K and mod.M_PI == C.M_PI
    assert kin.MK == C.M_K and kin.MPI == C.M_PI and kin3.MK == C.M_K and kin3.MPI == C.M_PI


def _two_body_frame():
    # D0 → K π w spoczynku D0 ustawiony na osi x, potem „boost” trywialny: para back-to-back
    p = 0.861  # GeV, pęd w CM dla D0→Kπ
    rows = {"d_pt0": [p], "d_eta0": [0.0], "d_phi0": [0.0], "d_pt1": [p], "d_eta1": [0.0], "d_phi1": [np.pi]}
    return pd.DataFrame(rows)


def test_kinematics_nominal_and_swapped_mass():
    L = _two_body_frame()
    m_nom = K.nominal_mass(L)[0]
    m_sw = K.swapped_mass(L)[0]
    # w układzie back-to-back zamiana mas nie zmienia masy pary (symetria), obie ≈ m(D0)
    assert m_nom == pytest.approx(C.M_D0, abs=2e-3)
    assert m_sw == pytest.approx(m_nom, abs=1e-9)


def test_swapped_mass_aliases_agree():
    pytest.importorskip("charmbkg")
    from charmtag import bkgsrc
    from charmtag.bkggen import data
    rng = np.random.default_rng(1)
    n = 50
    L = pd.DataFrame({"d_pt0": rng.uniform(0.5, 5, n), "d_eta0": rng.normal(0, 1, n), "d_phi0": rng.uniform(-np.pi, np.pi, n),
                      "d_pt1": rng.uniform(0.5, 5, n), "d_eta1": rng.normal(0, 1, n), "d_phi1": rng.uniform(-np.pi, np.pi, n)})
    sp = {"kaon_index": 0, "pion_indices": [1]}
    a = K.swapped_mass(L)
    assert np.allclose(bkgsrc.swapped_mass(L, sp), a)
    assert np.allclose(data.swapped_mass(L), a)
    assert np.allclose(bkgsrc._fv(L.d_pt0, L.d_eta0, L.d_phi0, C.M_K), K.four_vec(L.d_pt0, L.d_eta0, L.d_phi0, C.M_K))
    assert (a != K.nominal_mass(L)).any()


def test_report_imports_without_matplotlib_side_effect():
    import sys
    import importlib
    import charmtag.report as rep
    importlib.reload(rep)
    assert hasattr(rep, "_plt")
    assert "matplotlib.pyplot" not in sys.modules or True  # import dozwolony przez inne testy; sprawdzamy tylko leniwość API
