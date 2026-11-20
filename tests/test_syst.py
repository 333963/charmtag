"""Testy modułu systematyk (etap U): Asimov z globalnymi obserwablami, SF, normalizacja, pliki, raport."""
import json

import numpy as np
import pytest
import pandas as pd

from charmtag.pairs import CHANNELS
from charmtag.pairs.model import Model, region_mask
from charmtag.syst import core, files, norm, report
from charmtag.pairs.model import Inputs, M0


# toy_inputs: kopia z tests/test_pairs.py (bez importu pipeline/xgboost)
def toy_inputs(seed=1, nleg=40000, npair=400.0):
    rng = np.random.default_rng(seed)
    E = {10: np.linspace(1750, 2000, 26), 4: np.linspace(1700, 2200, 26)}
    C_ = {k: 0.5 * (v[1:] + v[:-1]) for k, v in E.items()}
    ncat = {10: 5, 4: 3}
    sig, dl, Y = {}, {}, {}
    for dt in (10, 4):
        for c in range(ncat[dt]):
            s = np.exp(-0.5 * ((C_[dt] - M0[dt]) / (12 + 3 * c)) ** 2) + (0.15 * np.exp(-0.5 * ((C_[dt] - M0[dt]) / 60) ** 2) if dt == 10 else 0)
            sig[(dt, c)] = s / s.sum()
            u = (C_[dt] - C_[dt][0]) / (C_[dt][-1] - C_[dt][0])
            b = np.exp(-1.2 * u) * nleg * (c + 1) / 6
            dl[(dt, c)] = rng.poisson(b).astype(float)
            Y[(dt, c)] = 0.2
    ncat_pairs = {}
    I0 = Inputs(edges=E, ncat=ncat, sig=sig, sig_comp={}, Y=Y, dataleg=dl, nev_legs=1e6, ndata=1e6, n_real={10: 3000.0, 4: 2000.0},
                n_pair={"D0D0": npair, "DpD0": npair, "DpDp": npair}, sbsb={}, meta={})
    # sbsb z niezależnych nóg × czynnik korelacji 1
    from charmtag.pairs.model import Shapes, cat_pairs
    S = Shapes(I0)
    for ch, (a, b) in CHANNELS.items():
        for i, j in cat_pairs(ncat, a, b):
            BB2 = np.outer(S.fshape(a, S.TH[(a, i)]), S.fshape(b, S.TH[(b, j)]))
            fr = BB2[np.ix_(I0.SBM[a], I0.SBM[b])].sum()
            ncat_pairs[(ch, i, j)] = int(rng.poisson(1e6 * S.NF[(a, i)] * S.NF[(b, j)] * (0.5 if (a == b and i == j) else 1) * fr * 3))
    I0.sbsb = ncat_pairs
    return I0


def _blind_hists(M, seed=3):
    rng = np.random.default_rng(seed)
    H = {}
    for ch in CHANNELS:
        cf = M.channel(ch)
        lam = cf.lam(cf.x0())
        n = rng.poisson(lam).astype(float)
        n[:, region_mask(M.I, ch, "centre")] = 0.0
        H[ch] = n
    return H


def test_asimov_nominal_mu1_and_sf():
    M = Model(toy_inputs())
    H = _blind_hists(M)
    cf, x, src = core.base_point(M, H, "DpDp")
    assert src.startswith("blind_fit")
    a = core.asimov_fit(cf, cf.lam(x), x)
    assert abs(a["mu"] - 1) < 2e-3, a["mu"]
    up = core.truth_yscale(M, "DpDp", x, {k: 1.1 for k in M.I.Y if k[0] == 4})
    b = core.asimov_fit(cf, up, x)
    assert abs(b["mu"] - 1.21) < 0.03, b["mu"]          # SS ∝ Y_a·Y_b → μ' = 1,1²


@pytest.mark.slow
def test_run_channel_dpdp_sources():
    M = Model(toy_inputs())
    H = _blind_hists(M)
    r = core.run_channel(M, H, "DpDp", sources=["width_Dp", "mscale_Dp", "refl_Ds", "sf_Dp", "kappa"])
    s = r["sources"]
    assert abs(r["nominal"]["mu"] - 1) < 2e-3
    assert s["width_Dp"]["up"]["dmu"] < 0 < s["width_Dp"]["down"]["dmu"]      # szerszy pik w prawdzie → μ w dół (okno)
    assert s["refl_Ds"]["up"]["dmu"] > 0 and "down" not in s["refl_Ds"]
    assert "kappa_channel" in s and "kappa_partner05" in s


def test_norm_collect_fallbacks():
    meta = {"norm_mode": "frozen", "dstar_fit_data": {"N": 3000.0, "N_err": 300.0}}
    r = norm.collect(meta, None)
    assert abs(r["mc_stat"]["delta"]["D0D0"] - 0.35) < 1e-9
    assert abs(r["ndstar_fit"]["delta"]["DpD0"] - 0.10) < 1e-9
    assert r["tune_MPI"]["pending"] and r["eff_ntrk"]["pending"]
    assert r["diffraction_x1p5"]["delta"]["DpD0"] == 0.019
    assert abs(norm.dmu_rel(0.25) + 0.2) < 1e-12


def test_files_assign_and_overdispersion():
    ft = pd.DataFrame({"file": ["a", "b", "c"], "job": ["1", "2", "3"], "rows": [10, 20, 5], "events": [1, 2, 1]})
    L = pd.DataFrame({"cand": [0, 9, 10, 29, 30, 34]})
    assert files.assign_files(L, ft).file_idx.tolist() == [0, 0, 1, 1, 2, 2]
    M = Model(toy_inputs())
    rng = np.random.default_rng(5)
    nf = 30
    w = np.full(nf, 1 / nf)
    Hf = {}
    for ch in CHANNELS:
        lam = M.channel(ch).lam(M.channel(ch).x0())
        Hf[ch] = np.array([rng.poisson(lam * w[k]) for k in range(nf)]).astype(float)
    od = files.overdispersion(Hf, M.I, w)
    assert all(0.4 < v["phi"] < 2.0 for v in od.values()), od


def test_report_build(tmp_path):
    nomin = {"mu": 1.0, "err": 1.0, "Z": 1.0}
    for ch, e in zip(CHANNELS, (1.0, 0.9, 3.9)):
        r = {"nominal": dict(nomin, err=e), "sources": {"width_D0": {"up": {"dmu": -0.1, "dZ": -0.1}, "down": {"dmu": 0.05, "dZ": 0.05},
                                                                    "dmu_sym": 0.1, "dZ_min": -0.1}}}
        json.dump(r, open(tmp_path / f"core_{ch}.json", "w"))
    json.dump(norm.collect({"norm_mode": "frozen"}, None), open(tmp_path / "norm.json", "w"))
    md, js = report.build(str(tmp_path))
    assert "szerokość piku D0" in md and "statystyka MC par" in md and (tmp_path / "syst_blind.md").exists()
    assert js["total_all"]["D0D0"] > 0.35
