"""Kinematyka córek — JEDNA implementacja czterowektorów, masy niezmienniczej i masy po zamianie hipotez K↔π
(REVIEW P1-8: wcześniej 4 kopie w bkgsrc/refldata/bkggen.data/candrot).

Wejście: kolumny `d_pt{i}` [GeV], `d_eta{i}`, `d_phi{i}` po `features.compute` (pT córek w h5 jest w MeV —
konwersję robi `features`/`units.d_pt`, nie ten moduł).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .constants import M_K, M_PI


def four_vec(pt, eta, phi, m):
    """(E, px, py, pz) dla pT, η, φ [GeV, rad] i masy m [GeV]; kształt (n, 4)."""
    pt, eta, phi = (np.asarray(v, float) for v in (pt, eta, phi))
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    e = np.sqrt((pt * np.cosh(eta)) ** 2 + m * m)
    return np.stack([e, px, py, pz], -1)


def inv_mass(*p4):
    """Masa niezmiennicza sumy czterowektorów (n, 4) → (n,)."""
    s = sum(p4)
    return np.sqrt(np.maximum(s[..., 0] ** 2 - s[..., 1] ** 2 - s[..., 2] ** 2 - s[..., 3] ** 2, 0.0))


def _col(L: pd.DataFrame, q: str, i: int) -> np.ndarray:
    return L[f"{q}{i}"].to_numpy(dtype=float)


def daughter_mass(L: pd.DataFrame, kaon_index: int = 0, pion_index: int = 1, swap: bool = False) -> np.ndarray:
    """Masa pary (K, π) z kolumn córek. `swap=True` → córka K dostaje m_π i odwrotnie (refleksja K↔π dla D0)."""
    mk, mp = (M_PI, M_K) if swap else (M_K, M_PI)
    a = four_vec(_col(L, "d_pt", kaon_index), _col(L, "d_eta", kaon_index), _col(L, "d_phi", kaon_index), mk)
    b = four_vec(_col(L, "d_pt", pion_index), _col(L, "d_eta", pion_index), _col(L, "d_phi", pion_index), mp)
    return inv_mass(a, b)


def swapped_mass(L: pd.DataFrame, kaon_index: int = 0, pion_index: int = 1) -> np.ndarray:
    """D0: masa przy zamianie hipotez K↔π."""
    return daughter_mass(L, kaon_index, pion_index, swap=True)


def nominal_mass(L: pd.DataFrame, kaon_index: int = 0, pion_index: int = 1) -> np.ndarray:
    """D0: masa z nominalnymi hipotezami (kontrola zgodności z `Dmeson_mass`)."""
    return daughter_mass(L, kaon_index, pion_index, swap=False)


def swapped_mass_sp(L: pd.DataFrame, sp: dict) -> np.ndarray:
    """Jak `swapped_mass`, indeksy córek z wpisu `species` configu (`kaon_index`, `pion_indices[0]`)."""
    return swapped_mass(L, int(sp.get("kaon_index", 0)), int((sp.get("pion_indices") or [1])[0]))
