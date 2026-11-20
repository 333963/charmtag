"""Stałe fizyczne i jednostki — JEDNO źródło dla całego pakietu (REVIEW P1-7).

Konwencja: wewnątrz pakietu masy w **GeV**. Moduł `pairs` liczy histogramy w MeV (krawędzie binów z configu
`pairs.edges` są w MeV, jak w historycznym fit2d M7) — konwersja wyłącznie przez `to_mev`/`to_gev` na brzegu.
Wartości PDG 2024.
"""
from __future__ import annotations

M_PI = 0.13957039      # GeV, π±
M_K = 0.493677         # GeV, K±
M_D0 = 1.86484         # GeV
M_DPLUS = 1.86966      # GeV
M_DSTAR = 2.01026      # GeV
DM_DSTAR = 0.14543     # GeV, m(D*⁺) − m(D0)

# decayType (ntuple CharmAnalysis) → m0 [GeV]
M0_GEV = {10: M_D0, 4: M_DPLUS, 1: M_DSTAR}
# nazwa gatunku → m0 [GeV]
M0_SPECIES = {"D0": M_D0, "Dplus": M_DPLUS, "Dp": M_DPLUS, "Dstar": M_DSTAR}

MEV_PER_GEV = 1000.0


def to_mev(x):
    """GeV → MeV (skalar, ndarray lub Series)."""
    return x * MEV_PER_GEV


def to_gev(x):
    """MeV → GeV."""
    return x / MEV_PER_GEV
