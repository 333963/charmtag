"""Test normalizacji --norm truth (etap T) na syntetycznym truth_norm.json."""
import json

import pytest

from charmtag.pairs import truthnorm as T


def _rep(tmp_path):
    tag = "fd_nd|trk0.85"
    R = {"pair_per_dstar_reco": {}, "reco_per_event": {}}
    for ch, v in (("D0D0", 0.020), ("DpD0", 0.035), ("DpDp", 0.010)):
        R["pair_per_dstar_reco"][f"{ch}|{tag}"] = [v, 0.1 * v]
        R["pair_per_dstar_reco"][f"{ch}|fd_nd|trk0.80"] = [1.05 * v, 0.1 * v]
    R["reco_per_event"].update({f"Dstar|{tag}": [1e-3, 1e-5], f"D0|{tag}": [7e-3, 1e-5], f"Dplus|{tag}": [6e-3, 1e-5]})
    rep = {"results": R, "nominal_tag": tag, "generator": {"n_events": 5e8},
           "variants_rel_to_nominal_pair_per_dstar": {"D0D0": {"fd_nd|trk0.80": 0.05, "fd_nd|trk0.85|norho": 0.3}}}
    p = tmp_path / "truth_norm.json"
    p.write_text(json.dumps(rep))
    return str(p)


def test_truth_alternative(tmp_path):
    p = _rep(tmp_path)
    cp = {"D0D0": 1.37, "DpD0": 1.20, "DpDp": 0.98}
    cs = {10: 0.93, 4: 0.79}
    a = T.truth_alternative(p, 3000.0, cp, cs, apply_cw=True)
    assert a["n_pair"]["D0D0"] == pytest.approx(3000 * 0.020 * 1.37)
    assert a["n_real"][10] == pytest.approx(3000 * 7.0 * 0.93)
    assert a["rel_stat"]["DpD0"] == pytest.approx(0.1)
    assert a["rel_syst_variants_max"]["D0D0"] == pytest.approx(0.05)   # norho pominięte
    b = T.truth_alternative(p, 3000.0, cp, cs, apply_cw=False)
    assert b["n_pair"]["DpDp"] == pytest.approx(3000 * 0.010)
    md = T.compare_md({"truth": a, "truth_nocw": b})
    assert "truth_nocw" in md and "D0D0" in md
