"""Adapter h5 v1.1 (surowe kolumny Dmeson_*) → kanoniczne kolumny ręcznego dyskryminanta (etap W-val `adapt_v11.py`,
wpięte przez VAL). Bez σ(Lxy)/σ(d0) w h5 → pseudo-istotności: Lxy/0,20 mm, |d0|/σ_d0(pT). daughter0 = K."""
from __future__ import annotations
import numpy as np
import pandas as pd

LXY_SIGMA_MM = 0.20
DEDX_SCALE = 0.87          # MIP w jednostkach dE/dx danych v1.1 (W-val)
NDAU = {"D0": 2, "Dplus": 3, "Dp": 3}
RAW = ["Dmeson_mass", "Dmeson_pt", "Dmeson_eta", "Dmeson_Lxy", "Dmeson_cosThetaP", "Dmeson_chi2", "Dmeson_cosThetaStar",
       "Dmeson_impactSignificance"]
DAU = ["pt", "eta", "d0", "dEdx", "npts_dEdx"]


def raw_columns(nd: int) -> list[str]:
    return RAW + [f"Dmeson_daughter{i}_{v}" for i in range(nd) for v in DAU]


def sigma_d0_mm(pt_gev):
    return np.sqrt(0.012 ** 2 + (0.09 / np.clip(pt_gev, 0.3, None)) ** 2)


def to_canonical(g: dict, species: str) -> pd.DataFrame:
    """g: nazwa kolumny h5 → tablica (float). Zwraca ramkę z kolumnami COLMAP (m, pt, eta, lxy_sig, cos_point, chi2, cos_star,
    ipsig_min/max, d0_min/max, dedx_K, p_K, dedx_pi1, p_pi1, [dedx_pi2, p_pi2], d_ipsig)."""
    nd = NDAU[species]
    f = lambda k: np.asarray(g[k], np.float64)
    out = {"m": f("Dmeson_mass"), "pt": f("Dmeson_pt"), "eta": f("Dmeson_eta"), "chi2": f("Dmeson_chi2"),
           "cos_point": f("Dmeson_cosThetaP"), "cos_star": np.abs(f("Dmeson_cosThetaStar")), "d_ipsig": np.abs(f("Dmeson_impactSignificance"))}
    out["lxy_sig"] = f("Dmeson_Lxy") / LXY_SIGMA_MM
    d0s, ips, ps, des = [], [], [], []
    for i in range(nd):
        pt = f(f"Dmeson_daughter{i}_pt") / 1000.0; eta = f(f"Dmeson_daughter{i}_eta")
        p = pt * np.cosh(np.clip(eta, -3, 3)); d0 = np.abs(f(f"Dmeson_daughter{i}_d0"))
        de = f(f"Dmeson_daughter{i}_dEdx"); npts = f(f"Dmeson_daughter{i}_npts_dEdx")
        valid = (pt > 0) & (np.abs(eta) < 5)
        d0s.append(np.where(valid, d0, np.nan)); ips.append(np.where(valid, d0 / sigma_d0_mm(pt), np.nan))
        ps.append(np.where(valid, p, np.nan)); des.append(np.where((npts > 0) & (de > 0) & valid, de, np.nan))
    D0 = np.vstack(d0s); IP = np.vstack(ips)
    with np.errstate(all="ignore"):
        out["d0_min"] = np.nanmin(D0, 0); out["d0_max"] = np.nanmax(D0, 0)
        out["ipsig_min"] = np.nanmin(IP, 0); out["ipsig_max"] = np.nanmax(IP, 0)
    out["p_K"], out["dedx_K"] = ps[0], des[0]
    out["p_pi1"], out["dedx_pi1"] = ps[1], des[1]
    if nd == 3:
        out["p_pi2"], out["dedx_pi2"] = ps[2], des[2]
    return pd.DataFrame(out)


def hand_scores_for_rows(cfg, key: str, species: str, rows: np.ndarray, variant: str = "B", rev: int = 1, chunk: int = 1_000_000):
    """Surowe kolumny dla wierszy `rows` (globalne indeksy jak w scores) z pliku źródłowego próbki `key` → score ręczny.
    Zwraca DataFrame(row, score_hand)."""
    from .. import config as C
    from ..io import open_events, list_columns
    from .model import score as hand_score
    nd = NDAU[species]; cols = raw_columns(nd)
    path = C.sample_path(cfg, key); ds_name = cfg["io"]["dataset"]
    avail = set(list_columns(path, ds_name)); miss = [c for c in cols if c not in avail]
    if miss:
        raise KeyError(f"hand: brak kolumn w {key}: {miss}")
    rows = np.unique(np.asarray(rows, np.int64)); parts = []
    with open_events(C.sample_paths(cfg, key) if C.is_multi(cfg, key) else [path], ds_name) as ds:
        n = ds.shape[0]
        for start in range(0, n, chunk):
            stop = min(n, start + chunk)
            sel = rows[(rows >= start) & (rows < stop)]
            if not len(sel):
                continue
            a = ds.fields(cols)[start:stop]
            g = {c: np.asarray(a[c])[sel - start] for c in cols}
            df = to_canonical(g, species)
            s = hand_score(df, species="Dp" if species == "Dplus" else species, variant=variant, rev=rev, dedx_scale=DEDX_SCALE)
            parts.append(pd.DataFrame({"row": sel, "score_hand": s.iloc[:, 0].to_numpy()}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["row", "score_hand"])


def hand_scores_from_npz(globpat: str, species: str, rows: np.ndarray, variant: str = "B", rev: int = 1):
    """Fallback bez źródłowego h5: surowe kolumny z podpróbek npz (etapy S/S6, klucz `_row` = globalny wiersz h5)."""
    import glob
    from .model import score as hand_score
    nd = NDAU[species]; cols = raw_columns(nd); rows = np.unique(np.asarray(rows, np.int64)); parts = []
    for f in sorted(glob.glob(globpat)):
        z = np.load(f, mmap_mode="r"); r = np.asarray(z["_row"]).astype(np.int64)
        pos = np.clip(np.searchsorted(r, rows), 0, len(r) - 1); hit = r[pos] == rows
        if not hit.any():
            continue
        g = {c: np.asarray(z[c])[pos[hit]] for c in cols}
        s = hand_score(to_canonical(g, species), species="Dp" if species == "Dplus" else species, variant=variant, rev=rev, dedx_scale=DEDX_SCALE)
        parts.append(pd.DataFrame({"row": rows[hit], "score_hand": s.iloc[:, 0].to_numpy()}))
    return pd.concat(parts, ignore_index=True).drop_duplicates("row") if parts else pd.DataFrame(columns=["row", "score_hand"])
