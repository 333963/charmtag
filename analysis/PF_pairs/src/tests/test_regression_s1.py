"""Test regresji liczb S1 (ε_B ratio próbek kontrolnych SB vs prawda MC, D0) — etap AA.

Źródło liczb: doc S (charm-S-bkggen.md) i niezależne odtworzenie S7 (charm-S7-przeglad.md §1).
Test czyta wynik `scripts/studies/s1_control.py` (JSON). Bieg na podpróbkach S_bkggen/sub trwa kilka–kilkanaście minut,
więc test jest pomijany, gdy JSON-u nie ma: uruchom
    CHARM_SUB=<...>/S_bkggen/sub CHARM_OUT=runs/s1_control python scripts/studies/s1_control.py
albo podaj gotowy plik w CHARM_S1_JSON.
"""
import json
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CANDIDATES = [os.environ.get("CHARM_S1_JSON", ""), os.path.join(HERE, "..", "runs", "s1_control", "s1_control.json")]

# (klucz, WP) → (wartość odniesienia, tolerancja). WP50: statystyka ±3–4 % na ½ MC (S7 §2.5) → tolerancja 0,05.
REF = {
    ("near_nocut", 80): (1.08, 0.02), ("near_nocut", 70): (1.16, 0.02), ("near_nocut", 50): (1.46, 0.05),
    ("far_nocut", 80): (0.96, 0.02), ("far_nocut", 70): (1.00, 0.02), ("far_nocut", 50): (1.19, 0.05),
    ("near_truthbkg_nocut", 80): (1.04, 0.02), ("near_truthbkg_nocut", 70): (1.07, 0.02), ("near_truthbkg_nocut", 50): (1.10, 0.05),
}


def _load():
    for p in CANDIDATES:
        if p and os.path.exists(p):
            return json.load(open(p)), p
    pytest.skip("brak s1_control.json — uruchom scripts/studies/s1_control.py (patrz docstring)")


def test_s1_effB_ratios_within_tolerance():
    d, p = _load()
    sb = d["sideband"]
    bad = []
    for (k, wp), (ref, tol) in REF.items():
        v = sb[k][f"effB_ratio_{wp}"]
        if abs(v - ref) > tol:
            bad.append(f"{k}@{wp}: {v:.3f} vs {ref:.2f}±{tol}")
    assert not bad, f"regresja S1 ({p}): " + "; ".join(bad)


def test_s1_reflection_template_has_errors_and_is_not_clipped():
    d, _ = _load()
    r = d["reflection"]
    assert "n_template_net_err" in r and "var_tmpl" in r, "s1_control bez poprawki AA (wariancja szablonu)"
    # S7: szablon z ½ MC jest statystycznie pusty (138 ± 184); obcięty sztucznie ~920 — liczby muszą się różnić
    assert r["n_template_net_clipped"] >= r["n_template_net"]
    assert r["n_template_net_err"] > 0
    # χ² z błędami Poissona ≥ χ² starą metodą (σ = h zaniżało σ → zaniżało χ²... poprawka podnosi χ²)
    assert r["chi2ndf"] >= 0.5 * r["chi2ndf_old_method"]
