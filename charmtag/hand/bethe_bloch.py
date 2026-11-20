"""Minimalna parametryzacja Bethe-Blocha dla dE/dx (jednostki: wielokrotność MIP).

Priorytet: użyj charmpy.bethe_bloch jeśli jest importowalne (API wykrywane
dynamicznie), inaczej lokalny fallback. Fallback = uproszczony BB (bez poprawki
gęstości) znormalizowany tak, że cząstka o βγ≈3.5 daje 1.0 (MIP).
Wystarcza do PULLI K-vs-π: liczy się kształt względem βγ, nie wartość absolutna.
"""
from __future__ import annotations
import numpy as np

from ..constants import M_K, M_PI  # noqa: F401  (GeV; re-eksport dla hand.features)
M_P = 0.93827
M_E = 0.000511
I_SI = 173e-9   # GeV (potencjał jonizacji Si)
_K_CONST = 0.307075e-3 * 14 / 28.0855  # GeV cm^2/g * Z/A dla Si

def _bb_raw(betagamma: np.ndarray) -> np.ndarray:
    bg = np.clip(np.asarray(betagamma, dtype=float), 1e-3, 1e4)
    beta2 = bg**2 / (1.0 + bg**2)
    # Tmax dla ciężkiej cząstki ~ 2 m_e β²γ² (pomijamy człony m_e/M)
    tmax = 2 * M_E * bg**2
    arg = 2 * M_E * bg**2 * tmax / I_SI**2
    raw = _K_CONST / beta2 * (0.5 * np.log(arg) - beta2)
    # efekt gęstości (Fermi plateau) – przybliżenie: powyżej minimum (βγ≈3.5) wzrost ograniczony do ~+6%
    bg_min = 3.5
    mn = _K_CONST / (bg_min**2 / (1 + bg_min**2)) * (0.5 * np.log(4 * M_E**2 * bg_min**4 / I_SI**2) - bg_min**2 / (1 + bg_min**2))
    plateau = mn * (1.0 + 0.03 * np.tanh(np.log(np.clip(bg / bg_min, 1e-6, None))))
    return np.where(bg > bg_min, plateau, raw)

_MIP_NORM = float(_bb_raw(np.array([3.5]))[0])

def dedx_expected(p_gev, mass_gev) -> np.ndarray:
    """Oczekiwane dE/dx w jednostkach MIP dla pędu p [GeV] i masy [GeV]."""
    p = np.clip(np.asarray(p_gev, dtype=float), 1e-3, None)
    return _bb_raw(p / mass_gev) / _MIP_NORM

def dedx_sigma(dedx_exp, rel=0.12) -> np.ndarray:
    """Rozdzielczość dE/dx: ~12% względna (ATLAS pixel, 3–4 klastry) – stała hipoteza."""
    return rel * np.asarray(dedx_exp, dtype=float)

def pull(dedx_meas, p_gev, mass_gev, scale=1.0, rel=0.12) -> np.ndarray:
    """(dE/dx_zmierzone/scale - BB(p,m)) / σ.  `scale` = wartość MIP w jednostkach danych."""
    exp = dedx_expected(p_gev, mass_gev)
    return (np.asarray(dedx_meas, dtype=float) / scale - exp) / dedx_sigma(exp, rel)

# --- opcjonalny most do charmpy -------------------------------------------
def _try_charmpy():
    try:
        import charmpy.bethe_bloch as cb  # type: ignore
    except Exception:
        return None
    for name in ("dedx_expected", "expected_dedx", "bethe_bloch", "dEdx"):
        f = getattr(cb, name, None)
        if callable(f):
            return f
    return None

CHARMPY_BB = _try_charmpy()

def dedx_expected_any(p_gev, mass_gev):
    """Użyj charmpy jeśli dostępne; zakładamy sygnaturę f(p, m). W razie błędu fallback."""
    if CHARMPY_BB is not None:
        try:
            return np.asarray(CHARMPY_BB(np.asarray(p_gev, float), mass_gev), float)
        except Exception:
            pass
    return dedx_expected(p_gev, mass_gev)
