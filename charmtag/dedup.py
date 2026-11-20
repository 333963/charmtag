"""Etap DD: wykrywanie i usuwanie dokładnych duplikatów kandydatów (addytywnie).

Duplikat = ten sam (run, event, decay_type) i te same ślady córek W TEJ SAMEJ KOLEJNOŚCI (kolejność córek =
hipoteza K/π). Identyfikator śladu = bity (d0, z0·sinθ_PV) córki — parametry toru względem PV, NIEZALEŻNE od refitu
wierzchołka kandydata (DD2/PF: pędy córek w h5 są po refit, więc klucz z (pt, eta, phi) nie rozpoznaje tego samego
śladu w dwóch kandydatach; wspólne ślady między kandydatami: 100 % par zgodnych z tolerancją 1 % ma identyczne bity
(d0, z0), kolizje klucza 5e-5). Gdy brak kolumn d_d0/d_z0 → klucz zapasowy (pt, eta, phi, m) z ostrzeżeniem.
Wyniki: v9 (trackId) — identyczna partycja co (trackId0, trackId1, pdgId0); h5 v1.1 — 0 duplikatów na 4,7 M wierszy
(DD2, klucz d0/z0; DD z kluczem pt/eta/phi też 0). Moduł jest strażnikiem na wypadek innej produkcji
(v9: 10,8 % D0 = kopia D0 z D*).

Konfiguracja: `samples.dedup: false` (domyślnie, bez zmian) | true (usuń, zostaw pierwszy) | "report" (tylko log).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("charmtag.dedup")

_TRACK_COLS = [f"d_{v}{i}" for i in range(3) for v in ("d0", "z0")]          # klucz właściwy (DD2): ślad = (d0, z0)
_DAUGHTER_COLS = [f"d_{v}{i}" for i in range(3) for v in ("pt", "eta", "phi")]  # klucz zapasowy (DD): po refit


def _bits(x) -> np.ndarray:
    a = np.ascontiguousarray(np.asarray(x, dtype=np.float32))
    a = np.where(np.isnan(a), np.float32(-1e30), a).astype(np.float32)  # NaN (brak 3. córki) → stała
    return a.view(np.int32).astype(np.int64)


def dup_key(L: pd.DataFrame) -> np.ndarray:
    """Macierz klucza (n, k) int64. Właściwy: (run, event, decay_type, d0/z0 córek w kolejności). Zapasowy (brak d0/z0):
    (run, event, decay_type, pt/eta/phi córek, m) — rozpoznaje tylko kopie bit w bit (nie po refit) → log.warning."""
    cols = []
    for c in ("run", "event", "decay_type"):
        if c in L.columns:
            cols.append(np.asarray(L[c]).astype(np.int64))
    if all(c in L.columns for c in _TRACK_COLS):
        for i in range(3):
            d0 = _bits(L[f"d_d0{i}"].to_numpy()); z0 = _bits(L[f"d_z0{i}"].to_numpy())
            pt = L[f"d_pt{i}"].to_numpy() if f"d_pt{i}" in L.columns else np.ones(len(L))
            tid = (d0 << 32) ^ (z0 & 0xFFFFFFFF)
            cols.append(np.where(np.asarray(pt) > 0, tid, -1))   # brak 3. córki (pt NaN/0) → -1
    else:
        missing = [c for c in _TRACK_COLS if c not in L.columns]
        log.warning("dedup: brak kolumn %s — klucz zapasowy (pt, eta, phi, m) po refit: wykrywa tylko kopie bit w bit", missing)
        for c in _DAUGHTER_COLS:
            if c in L.columns:
                cols.append(_bits(L[c].to_numpy()))
        if "m" in L.columns:
            cols.append(_bits(L["m"].to_numpy()))
    if not cols:
        raise ValueError("dedup: brak kolumn do klucza")
    return np.stack(cols, 1)


def dup_groups(L: pd.DataFrame):
    """(group_id, multiplicity) per wiersz."""
    k = dup_key(L)
    _, inv, cnt = np.unique(k, axis=0, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).ravel()
    return inv, cnt[inv]


def dedup_mask(L: pd.DataFrame) -> np.ndarray:
    """True = zostaw (pierwsze wystąpienie w kolejności ramki)."""
    inv, _ = dup_groups(L)
    keep = np.zeros(len(L), bool)
    first = np.unique(inv, return_index=True)[1]
    keep[first] = True
    return keep


def report(L: pd.DataFrame, label: str = "") -> dict:
    if len(L) == 0:
        return {"n": 0, "dup_rows": 0, "dup_frac": 0.0}
    inv, mult = dup_groups(L)
    isdup = mult > 1
    r = {"n": int(len(L)), "dup_rows": int(isdup.sum()), "dup_frac": float(isdup.mean()),
         "n_unique": int(len(np.unique(inv))), "max_mult": int(mult.max())}
    if "m" in L.columns and "decay_type" in L.columns and isdup.any():
        r["mass_range_dup"] = [float(L["m"][isdup].min()), float(L["m"][isdup].max())]
    log.info("dedup %s: %d wierszy, duplikatów %d (%.2f%%), krotność max %d", label, r["n"], r["dup_rows"],
             100 * r["dup_frac"], r["max_mult"])
    return r


def apply(L: pd.DataFrame, mode, label: str = "") -> pd.DataFrame:
    """mode: False/None → bez zmian; "report" → tylko log; True → usuń duplikaty."""
    if not mode or len(L) == 0:
        return L
    r = report(L, label)
    L.attrs["dedup_report"] = r
    if mode == "report":
        return L
    if r["dup_rows"] == 0:
        return L
    keep = dedup_mask(L)
    out = L[keep].reset_index(drop=True)
    out.attrs.update(L.attrs)
    log.info("dedup %s: usunięto %d wierszy", label, int((~keep).sum()))
    return out
