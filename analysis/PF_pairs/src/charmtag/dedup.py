"""Etap DD: wykrywanie i usuwanie dokładnych duplikatów kandydatów (addytywnie).

Duplikat = ten sam (run, event, decay_type), te same córki (pt, eta, phi bitowo, float32) i ta sama masa
(= ta sama hipoteza K/π). Klucz zwalidowany na h5 z trackId (H2, v9): identyczna partycja co
(trackId0, trackId1, pdgId0). W h5 v1.1 duplikatów NIE stwierdzono (0 / 5,4 M wierszy, DD) — moduł jest
strażnikiem na wypadek innej produkcji (v9: 10,8 % D0 = kopia D0 z D*).

Konfiguracja: `samples.dedup: false` (domyślnie, bez zmian) | true (usuń, zostaw pierwszy) | "report" (tylko log).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("charmtag.dedup")

_DAUGHTER_COLS = [f"d_{v}{i}" for i in range(3) for v in ("pt", "eta", "phi")]


def _bits(x) -> np.ndarray:
    a = np.ascontiguousarray(np.asarray(x, dtype=np.float32))
    a = np.where(np.isnan(a), np.float32(-1e30), a).astype(np.float32)  # NaN (brak 3. córki) → stała
    return a.view(np.int32).astype(np.int64)


def dup_key(L: pd.DataFrame) -> np.ndarray:
    """Macierz klucza (n, k) int64; kolumny brakujące w ramce są pomijane (klucz słabszy → log)."""
    cols = []
    for c in ("run", "event", "decay_type"):
        if c in L.columns:
            cols.append(np.asarray(L[c]).astype(np.int64))
    missing = [c for c in _DAUGHTER_COLS if c not in L.columns]
    for c in _DAUGHTER_COLS:
        if c in L.columns:
            cols.append(_bits(L[c].to_numpy()))
    if "m" in L.columns:
        cols.append(_bits(L["m"].to_numpy()))
    if missing:
        log.warning("dedup: brak kolumn %s — klucz słabszy", missing)
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
