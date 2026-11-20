"""Cechy logiczne z surowych kolumn (wspólne dla trening/apply/kontrola D*)."""
from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd

from .io import ColumnMap

log = logging.getLogger("charmtag")

SCALARS = ["run", "event", "weight", "decay_type", "m", "pt", "eta", "phi", "charge", "dm", "Lxy", "chi2",
           "impact", "impact_sig", "impact_z0", "cos_theta_star", "cos_theta_p", "ptcone", "truth", "m_kpi1", "m_kpi2"]
PER_D = ["d_pt", "d_eta", "d_phi", "d_p", "d_d0", "d_z0", "d_dedx", "d_npts"]
# HYP (2026-11-08): kolumny opcjonalne (`columns_opt` w YAML, v1.2) — dokładane do ramki TYLKO gdy istnieją w pliku
# (v1.1 → ramki bit w bit jak dotąd)
OPT_SCALARS = ["hypo_swap"]
OPT_PER_D = ["d_trk"]


def raw_columns(cm: ColumnMap, nd: int) -> list[str]:
    cols = [cm.get(k) for k in SCALARS]
    for k in PER_D:
        cols += [cm.get(f"{k}{i}") for i in range(nd)]
    cols += [cm.get(k) for k in OPT_SCALARS] + [cm.get(f"{k}{i}") for k in OPT_PER_D for i in range(nd)]  # HYP: None gdy brak
    return [c for c in cols if c]


def to_logical(df: pd.DataFrame, cm: ColumnMap, nd: int) -> pd.DataFrame:
    out = {}
    for k in SCALARS:
        c = cm.get(k)
        out[k] = df[c].to_numpy() if c and c in df else np.full(len(df), np.nan)
    for k in PER_D:
        for i in range(nd):
            c = cm.get(f"{k}{i}")
            out[f"{k}{i}"] = df[c].to_numpy(dtype=float) if c and c in df else np.full(len(df), np.nan)
    for k in OPT_SCALARS:  # HYP: tylko gdy kolumna jest w pliku
        c = cm.get(k)
        if c and c in df:
            out[k] = df[c].to_numpy()
    for k in OPT_PER_D:
        for i in range(nd):
            c = cm.get(f"{k}{i}")
            if c and c in df:
                out[f"{k}{i}"] = df[c].to_numpy()
    return pd.DataFrame(out)


def _stack(L, k, nd):
    return np.column_stack([L[f"{k}{i}"].to_numpy(dtype=float) for i in range(nd)])


def compute(L: pd.DataFrame, sp: dict, units: dict | None = None) -> pd.DataFrame:
    """Dodaje cechy pochodne do ramki logicznej L i ją zwraca. units: mnożniki do GeV (v1.1: pT córek w MeV)."""
    units = units or {}
    nd = int(sp["n_daughters"])
    upt = float(units.get("d_pt", 1.0))
    if upt != 1.0:
        for i in range(nd):
            L[f"d_pt{i}"] = L[f"d_pt{i}"] * upt
            L[f"d_p{i}"] = L[f"d_p{i}"] * upt
    ukp = float(units.get("m_kpi", 1.0))
    for c in ("m_kpi1", "m_kpi2"):
        L[c] = L[c].astype(float) * ukp
    f = lambda c: L[c].to_numpy(dtype=float)  # noqa: E731
    L["abs_eta"] = np.abs(f("eta"))
    L["log_pt"] = np.log1p(f("pt"))
    L["abs_impact"] = np.abs(f("impact"))
    L["abs_impact_sig"] = np.abs(f("impact_sig"))
    L["abs_impact_z0"] = np.abs(f("impact_z0"))
    L["abs_cos_theta_star"] = np.abs(f("cos_theta_star"))
    L["pointing"] = -np.log10(np.clip(1.0 - f("cos_theta_p"), 1e-7, 2.0))
    pt = _stack(L, "d_pt", nd)
    eta = _stack(L, "d_eta", nd)
    phi = _stack(L, "d_phi", nd)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        L["d_min_pt"] = np.nanmin(pt, axis=1)
        L["d_pt_ratio"] = np.nanmin(pt, axis=1) / np.nanmax(pt, axis=1)
        dphi = np.abs(((phi - f("phi")[:, None]) + np.pi) % (2 * np.pi) - np.pi)
        deta = eta - f("eta")[:, None]
        L["d_max_dR"] = np.nanmax(np.sqrt(dphi ** 2 + deta ** 2), axis=1)
        ad0 = np.abs(_stack(L, "d_d0", nd))
        L["d_min_abs_d0"] = np.nanmin(ad0, axis=1)
        L["d_max_abs_d0"] = np.nanmax(ad0, axis=1)
        L["d_max_abs_z0"] = np.nanmax(np.abs(_stack(L, "d_z0", nd)), axis=1)
        L["d_max_abs_eta"] = np.nanmax(np.abs(eta), axis=1)
        p = _stack(L, "d_p", nd)
        p = np.where(np.isfinite(p), p, pt * np.cosh(eta))
        dedx = _stack(L, "d_dedx", nd)
        dedx = np.where(dedx > 0, dedx, np.nan)  # ≤0 = brak pomiaru
        npts = _stack(L, "d_npts", nd)
        k = int(sp.get("kaon_index", 0))
        pis = list(sp.get("pion_indices", [1]))
        L["dedx_K"] = dedx[:, k]
        L["p_K"] = p[:, k]
        L["log_p_K"] = np.log(p[:, k])
        L["npts_K"] = npts[:, k]
        L["dedx_pi_max"] = np.nanmax(dedx[:, pis], axis=1)
        L["log_p_pi_min"] = np.log(np.nanmin(p[:, pis], axis=1))
    return L


def feature_list(cfg: dict, variant: str) -> list[str]:
    f = cfg["features"]
    feats = list(f["frozen7"])
    if variant == "H":
        feats += list(f["extra_H"])
    elif variant != "frozen7":
        raise ValueError(f"nieznany wariant {variant}")
    bad = set(feats) & set(f.get("forbidden", []))
    if bad or any("ptcone" in x for x in feats):
        raise ValueError(f"cechy zakazane (izolacja): {bad}")
    return feats


def usable(L: pd.DataFrame, feats: list[str]) -> list[str]:
    ok = []
    for c in feats:
        if c not in L or not np.isfinite(L[c].to_numpy(dtype=float)).any():
            log.warning("cecha %s niedostępna (same NaN) → pominięta", c)
        else:
            ok.append(c)
    return ok
