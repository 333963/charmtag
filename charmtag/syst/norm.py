"""Systematyki normalizacji oczekiwanej liczby par (μ ∝ 1/N_par: Δμ/μ = 1/(1+δ) − 1; Z nie zależy od normalizacji).

Źródła (δ = względna zmiana N_par per kanał):
- mc_stat: statystyka MC par (pair_inputs meta: norm_alternatives[norm_mode].rel_stat) + MC D* (1/√n_D*);
- ndstar_fit: N_D* z fitu Δm w danych (N_err/N), w pełni skorelowane między kanałami;
- tune_MPI / tune_CR / tune_Monash: [par/D*]_reco z Pythii — warianty vs A3 (truth/combined_syst_<wariant>/truth_norm.json);
- eff_ntrk: wydajność(n_trk): warianty mapy ε (bez n_trk, n_trk = 0,80/0,90·n_ch) z truth_norm.json A3;
- eff_map_sample: mapa ε z całego MC (pooled) vs fd_nd;
- diffraction_x1p5: przekroje dyfrakcyjne ×1,5 (z truth_norm.json per proces; fallback 09-28: +1,0/+1,9/+1,0%);
- minbias_mix: ND vs ND+SD+DD+CD (09-28: +3/+4/+2%, tylko gdy brak truth).
"""
from __future__ import annotations

import json
import math
import os

CH = ("D0D0", "DpD0", "DpDp")
SP_OF = {"D0D0": ("D0", "D0"), "DpD0": ("Dplus", "D0"), "DpDp": ("Dplus", "Dplus")}
FALLBACK_0928 = {"diffraction_x1p5": {"D0D0": 0.010, "DpD0": 0.019, "DpDp": 0.010},
                 "minbias_mix": {"D0D0": 0.03, "DpD0": 0.04, "DpDp": 0.02}}
TUNE_GROUPS = {"tune_MPI": ("MPIup", "MPIdown"), "tune_CR": ("CRup", "CRdown", "CRqcd"), "tune_Monash": ("Monash",)}


def dmu_rel(delta):
    """N_par → N_par(1+δ): μ' = μ/(1+δ)."""
    return 1.0 / (1.0 + delta) - 1.0


def _v(x):
    return float(x[0]) if isinstance(x, (list, tuple)) else float(x)


def _ppd(rep, ch, tag=None):
    tag = tag or rep["nominal_tag"]
    v = rep["results"]["pair_per_dstar_reco"].get(f"{ch}|{tag}")
    return None if v is None or v[0] is None else _v(v)


def diffraction_scale(rep, ch, k=1.5):
    """Prawda generatora: [par/D*] przy σ_dyfr × k vs nominalnie (udziały procesów z σ Pythii)."""
    R, g = rep["results"], rep.get("generator", {})
    sig = g.get("sigma_mb", {})
    tot = sum(sig.values()) or None
    if not tot:
        return None

    def ratio(kk):
        num = den = 0.0
        for pc, s in sig.items():
            w = s / tot * (kk if pc != "nd" else 1.0)
            p = R["truth_pairs_per_event"].get(f"{ch}|{pc}")
            d = R["truth_per_event"].get(f"Dstar|{pc}")
            if p is None or d is None or p[0] is None or d[0] is None:
                continue
            num += w * _v(p)
            den += w * _v(d)
        return num / den if den else None
    r1, rk = ratio(1.0), ratio(k)
    return None if not r1 or rk is None else rk / r1 - 1.0


def collect(meta: dict, truth_dir=None):
    """Zwraca {źródło: {"delta": {ch: δ}, "dmu_rel": {ch}, "how": str, "pending": bool}}."""
    out = {}
    mode = meta.get("norm_mode", "direct")
    alt = (meta.get("norm_alternatives") or {}).get(mode, {})
    rel = alt.get("rel_stat")
    nds = (meta.get("mc_rates_unweighted") or {}).get("n_dstar")
    if mode == "frozen" or not rel:
        rel = {"D0D0": 0.35, "DpD0": 0.27, "DpDp": 0.37}        # fizyka_blok1 09-28 (21/41/15 par, n_eff 8/13/8)
        how = "MC 09-26 (n_eff 8/13/8)"
    else:
        how = f"MC par ({mode})"
    ds = (1 / math.sqrt(nds)) if (nds and mode != "frozen") else 0.0
    d = {ch: math.hypot(rel[ch], ds) for ch in CH}
    out["mc_stat"] = {"delta": d, "how": how + (f" ⊕ MC D* {ds:.1%}" if ds else ""), "sym": True}
    f = meta.get("dstar_fit_data") or {}
    if f.get("N"):
        r = f["N_err"] / f["N"]
        out["ndstar_fit"] = {"delta": {ch: r for ch in CH}, "how": f"fit Δm danych N_D*={f['N']:.0f}±{f['N_err']:.0f}", "sym": True}
    # --- truth (lxplus)
    td = os.path.expanduser(truth_dir) if truth_dir else None
    nomf = os.path.join(td, "combined_A3", "truth_norm.json") if td else None
    nom = json.load(open(nomf)) if nomf and os.path.exists(nomf) else None
    if nom:
        var = nom.get("variants_rel_to_nominal_pair_per_dstar", {})
        pick = lambda ch, keys: max([abs(v) for k, v in var.get(ch, {}).items() if v is not None and math.isfinite(v)
                                     and any(s in k for s in keys) and not k.endswith("norho")] or [0.0])
        out["eff_ntrk"] = {"delta": {ch: pick(ch, ("no_ntrk", "trk0.80", "trk0.90")) for ch in CH},
                           "how": "mapa ε: bez n_trk / n_trk=0,80·0,90·n_ch (truth A3)", "sym": True}
        out["eff_map_sample"] = {"delta": {ch: pick(ch, ("pooled_all_mc|trk0.85",)) for ch in CH}, "how": "mapa ε z całego MC", "sym": True}
        dd = {ch: diffraction_scale(nom, ch) for ch in CH}
        if all(v is not None for v in dd.values()):
            out["diffraction_x1p5"] = {"delta": dd, "how": "σ_SD,DD,CD ×1,5 (prawda Pythii A3)", "sym": True}
        for grp, names in TUNE_GROUPS.items():
            vals, found = {ch: 0.0 for ch in CH}, []
            for nm in names:
                fn = os.path.join(td, f"combined_syst_{nm}", "truth_norm.json")
                if not os.path.exists(fn):
                    continue
                rv = json.load(open(fn))
                found.append(nm)
                for ch in CH:
                    a, b = _ppd(rv, ch), _ppd(nom, ch)
                    if a and b:
                        vals[ch] = max(vals[ch], abs(a / b - 1.0))
            out[grp] = {"delta": vals if found else {}, "how": f"Pythia {'/'.join(names)} vs A3" + ("" if found else " — CZEKA (lxplus gen_syst)"),
                        "pending": not found, "sym": True}
    else:
        out["eff_ntrk"] = {"delta": {}, "how": "CZEKA na truth_norm.json A3 (lxplus combine)", "pending": True}
        for grp, names in TUNE_GROUPS.items():
            out[grp] = {"delta": {}, "how": f"Pythia {'/'.join(names)} vs A3 — CZEKA (lxplus gen_syst)", "pending": True}
        out["diffraction_x1p5"] = {"delta": FALLBACK_0928["diffraction_x1p5"], "how": "09-28 (mieszanka ×1,5 vs mieszanka)", "sym": True}
        out["minbias_mix"] = {"delta": FALLBACK_0928["minbias_mix"], "how": "09-28 (ND vs ND+SD+DD+CD)", "sym": True}
    for k, v in out.items():
        # Δμ/μ ≈ ±δ (1. rząd; dokładnie μ/(1±δ) — asymetria raportowana osobno)
        v["dmu_rel"] = {ch: abs(x) for ch, x in v["delta"].items()}
        v["dmu_rel_asym"] = {ch: [dmu_rel(abs(x)), dmu_rel(-abs(x)) if abs(x) < 1 else None] for ch, x in v["delta"].items()}
    return out
