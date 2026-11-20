"""Testy B1 (salad_v1.py): konfiguracja, regiony, closure gładki na syntetyku (C2ST < 0,52), closure adwersarialny
(skład zlokalizowany pod pikiem — liczba REFERENCYJNA jako regresja, nie próg), zapis/odczyt wag, rejestracja api."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")
from charmtag.bkggen import salad_v1 as SV  # noqa: E402

FEATS = ["log_pt", "abs_eta", "chi2", "Lxy", "impact_sig", "abs_impact_z0", "d_max_abs_d0"]
M0 = 1.86966


def _toy(n, seed, comp_shift=0.0, local=False, smooth=0.0):
    """Tło D⁺ syntetyczne: dwa składniki (A: Lxy mała, B: Lxy duża), masa ~ exp. `comp_shift` zwiększa udział B
    (globalnie lub tylko pod pikiem gdy local); `smooth` — gładki w m dryf chi2."""
    r = np.random.default_rng(seed)
    m = 1.70 + r.exponential(0.25, n); m = m[(m < 2.20)]; n = len(m)
    g = np.exp(-0.5 * ((m - M0) / 0.035) ** 2) if local else np.ones(n)
    pB = np.clip(0.3 + comp_shift * g, 0, 1)
    isB = r.random(n) < pB
    t = (m - M0) / 0.1
    df = pd.DataFrame({
        "log_pt": np.log1p(r.uniform(3, 20, n)), "abs_eta": r.uniform(0, 2.4, n),
        "chi2": r.exponential(2 * (1 + smooth * (0.5 + 0.3 * t)), n), "Lxy": np.where(isB, r.exponential(1.5, n), r.exponential(0.3, n)),
        "impact_sig": r.exponential(1, n), "abs_impact_z0": r.exponential(0.2, n), "d_max_abs_d0": r.exponential(0.1, n),
        "m": m, "truth": np.where(isB, 9, 12)})
    return df


def _cfg():
    cfg = SV.load_cfg()
    cfg["features"] = FEATS; cfg["max_n"] = 30_000; cfg["hgb"]["max_iter"] = 80; cfg["c2st"]["max_n"] = 15_000
    return cfg


def test_cfg_and_regions():
    cfg = SV.load_cfg()
    assert cfg["kind"] in SV.KINDS and cfg["kind_crosscheck"] in SV.KINDS and len(cfg["features"]) == 18
    sp = SV.species_params(cfg, "Dplus")
    r = SV.regions([M0, M0 + 0.05, M0 + 0.1, 2.5], sp)
    assert list(r) == ["SR", "band", "SB", "out"]


def test_closure_smooth_and_weights_roundtrip(tmp_path):
    cfg = _cfg()
    mc, data = _toy(70_000, 1), _toy(70_000, 2, smooth=0.8)
    rw = SV.fit_salad_v1(cfg, "Dplus", data, mc, c2st_band=True)
    assert 0.5 <= rw.diag["val_auc_SB"] < 0.8 and rw.diag["neff_over_N_SR"] > 0.5
    w = SV.salad_weights(rw, mc)
    assert abs(w.mean() - 1) < 1e-6 and w.min() >= cfg["clip"][0] / w.mean() * 0.999
    sp = SV.species_params(cfg, "Dplus")
    from charmtag.bkggen.feta import c2st
    rM, rD = SV.regions(mc["m"], sp), SV.regions(data["m"], sp)
    auc = c2st(data.loc[rD == "SR", FEATS].to_numpy(float), mc.loc[rM == "SR", FEATS].to_numpy(float), None, w[rM == "SR"], max_n=15_000)
    raw = c2st(data.loc[rD == "SR", FEATS].to_numpy(float), mc.loc[rM == "SR", FEATS].to_numpy(float), max_n=15_000)
    assert auc < 0.52 and auc < raw - 0.02, (auc, raw)
    assert rw.diag["c2st_band"] < 0.53 and rw.diag["ok"]
    p = SV.run(cfg, "Dplus", data, mc, str(tmp_path), sample="toy")
    w2, meta = SV.load_weights(p)
    assert np.allclose(w2, w, atol=1e-6) and meta["version"] == SV.VERSION == "1.1" and meta["sha256_mc"] == SV.input_hash(mc, FEATS + ["m"])


def test_closure_adversarial_local_reference():
    """Skład zlokalizowany pod pikiem (B: 30 % → 60 % w SR, SB bez zmian): SALAD z SB NIE może tego skorygować.
    Liczby referencyjne (regresja; zmiana algorytmu = zmiana liczb, nie „błąd"): ε_B(Lxy>1) pred/true ≈ 0,6–0,8,
    C2ST SR po SALAD > 0,53. Wariant globalny (ten sam skład w SB) musi być skorygowany (pred/true 0,95–1,05)."""
    cfg = _cfg()
    sp = SV.species_params(cfg, "Dplus")
    out = {}
    for name, local in (("local", True), ("global", False)):
        mc, data = _toy(80_000, 3), _toy(80_000, 4, comp_shift=0.3, local=local)
        rw = SV.fit_salad_v1(cfg, "Dplus", data, mc, c2st_band=False)
        w = SV.salad_weights(rw, mc)
        rM, rD = SV.regions(mc["m"], sp), SV.regions(data["m"], sp)
        cut = lambda d: (d["Lxy"].to_numpy(float) > 1.0)  # noqa: E731
        true = cut(data[rD == "SR"]).mean()
        pred = np.average(cut(mc[rM == "SR"]), weights=w[rM == "SR"])
        out[name] = pred / true
    print("B1 ref pred/true:", out)
    assert 0.9 < out["global"] < 1.1, out            # skład stały w m → korekta pełna
    assert 0.5 < out["local"] < 0.85, out            # skład lokalny → błąd ~ −30 % (referencja B1 syntetyk)


def test_api_registered():
    from charmtag.bkggen import api
    assert "salad_v1" in api.available()
    cfg = _cfg()
    mc, data = _toy(40_000, 5), _toy(40_000, 6)
    gen = api.make_generator("salad_v1", feats=FEATS, species="Dplus", seed=0)
    gen.fit(data[np.abs(data["m"] - M0) > 0.07], mc=mc, cfg=cfg)
    df = gen.sample(np.full(500, M0))
    assert len(df) == 500 and "w" in df and set(FEATS) <= set(df.columns)


def test_mc_truth_min_per_species_v11():
    """v1.1 (B1-D0): D0 bierze tło MC truth ≥ 2 (z częściowo dopasowanym charmem 2–7 i swapami 100/101), D⁺ ≥ 8."""
    cfg = SV.load_cfg()
    assert cfg["version"] == "1.1" and SV.mc_truth_min(cfg, "D0") == 2 and SV.mc_truth_min(cfg, "Dplus") == 8
    assert SV.mc_truth_min({"mc_truth_min": 8}, "D0") == 8          # stary yaml bez sekcji per gatunek
    mc = pd.DataFrame({"m": np.full(6, M0), "truth": [0, 1, 3, 9, 12, 100]})
    assert sorted(SV._mc_bkg(cfg, mc, "D0")["truth"]) == [3, 9, 12, 100]
    assert sorted(SV._mc_bkg(cfg, mc, "Dplus")["truth"]) == [9, 12, 100]
    assert sorted(SV._mc_bkg(cfg, mc)["truth"]) == [9, 12, 100]
