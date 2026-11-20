"""Normalizacja par z prawdy Pythii 8 × mapa wydajności (truth_norm.json z Q_lxplus/truth/truth_combine.py) — etap T.

Addytywnie do build.py: `charmtag pairs --norm truth` (config: pairs.truth_norm_json, pairs.truth_apply_cw).
N_par(ch)  = N_D*(dane) × [par/D*]_reco,Pythia(nominal_tag) × c_w(ch)   (c_w opcjonalnie, jak w 'direct')
N_real(dt) = N_D*(dane) × [D/D*]_reco,Pythia(nominal_tag)     × c_w(dt)
Uwaga: poziom reco = preselekcja Donly × BR (bez ε(WP) i SF — te wchodzą przez Y kategorii, jak w 'direct').
"""
from __future__ import annotations

import json
import math

SP = {10: "D0", 4: "Dplus"}
CH = ("D0D0", "DpD0", "DpDp")


def load_truth(path):
    with open(path) as f:
        return json.load(f)


def _v(x):
    v = x[0] if isinstance(x, (list, tuple)) else x
    return float("nan") if v is None else float(v)


def _e(x):
    e = x[1] if isinstance(x, (list, tuple)) and len(x) > 1 else None
    return float("nan") if e is None else float(e)


def truth_ratios(rep, tag=None):
    """par/D* i D/D* (reco) dla danego wariantu (domyślnie nominal_tag) + względne błędy jackknife."""
    R = rep["results"]
    tag = tag or rep["nominal_tag"]
    ppd = {ch: R["pair_per_dstar_reco"][f"{ch}|{tag}"] for ch in CH}
    ds = R["reco_per_event"][f"Dstar|{tag}"]
    spd = {dt: _v(R["reco_per_event"][f"{sp}|{tag}"]) / _v(ds) for dt, sp in SP.items()}
    return ({ch: _v(v) for ch, v in ppd.items()},
            {ch: (_e(v) / _v(v) if _v(v) else float("nan")) for ch, v in ppd.items()},
            spd)


def truth_alternative(path, n_dstar, cp, cs, apply_cw=True, tag=None):
    rep = load_truth(path)
    ppd, rel, spd = truth_ratios(rep, tag)
    wp = {ch: (cp[ch] if apply_cw else 1.0) for ch in CH}
    ws = {dt: (cs[dt] if apply_cw else 1.0) for dt in SP}
    var = rep.get("variants_rel_to_nominal_pair_per_dstar", {})
    syst = {ch: max([abs(v) for k, v in var.get(ch, {}).items() if v is not None and not k.endswith("norho")
                     and math.isfinite(v)] or [0.0]) for ch in CH}
    return {"n_pair": {ch: n_dstar * ppd[ch] * wp[ch] for ch in CH},
            "n_real": {dt: n_dstar * spd[dt] * ws[dt] for dt in SP},
            "rel_stat": rel, "rel_syst_variants_max": syst,
            "pair_per_dstar": ppd, "single_per_dstar": {str(k): v for k, v in spd.items()},
            "apply_cw": bool(apply_cw), "tag": tag or rep["nominal_tag"],
            "n_gen_events": rep.get("generator", {}).get("n_events"), "source": path}


def compare_md(alt):
    """Tabela porównania wariantów normalizacji (oczekiwane pary na pełnej płaszczyźnie)."""
    L = ["| wariant | D0D0 | DpD0 | DpDp | D0 | D⁺ |", "|---|---|---|---|---|---|"]
    for k, a in alt.items():
        L.append(f"| {k} | " + " | ".join(f"{a['n_pair'][ch]:.1f}" for ch in CH) + " | "
                 + " | ".join(f"{a['n_real'][dt]:.0f}" for dt in (10, 4)) + " |")
    return "\n".join(L)
