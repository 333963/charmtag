"""Cechy wejściowe ręcznego dyskryminanta + deterministyczne transformacje (warstwa 0).

Kolumny źródłowe rozwiązywane przez aliasy (COLMAP) – dopasowanie do schematu h5
v1.1 / cech H z charmtag bez zmiany kodu: wystarczy dopisać alias albo podać
`colmap=` do score(). Kanoniczne nazwy wewnętrzne po lewej.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from . import bethe_bloch as bb

# kanoniczna nazwa -> lista aliasów (pierwszy znaleziony w df wygrywa)
COLMAP: dict[str, list[str]] = {
    "m":          ["m", "mass", "D_m", "cand_m", "Dmeson_m", "D0_m", "Dplus_m", "DMeson_m"],
    "pt":         ["pt", "D_pt", "cand_pt", "Dmeson_pt", "DMeson_pt"],
    "eta":        ["eta", "D_eta", "cand_eta", "Dmeson_eta", "DMeson_eta"],
    "lxy":        ["Lxy", "lxy", "D_Lxy", "DMeson_Lxy", "Lxy_mm"],
    "lxy_sig":    ["Lxy_sig", "lxy_sig", "LxySig", "Lxy_over_sigma", "DMeson_LxySig", "D_LxySig", "Lxy_significance"],
    "lxy_err":    ["Lxy_err", "LxyErr", "sigma_Lxy", "DMeson_LxyErr"],
    "cos_point":  ["cos_pointing", "cosTheta_pointing", "cos_theta_pointing", "cosPointing", "DMeson_cosPointing", "cos_point", "costheta_pointing"],
    "chi2":       ["chi2_vtx", "chi2", "chi2ndf", "chi2_ndof", "vtx_chi2", "DMeson_chi2", "chi2_ndf", "vertex_chi2"],
    "cos_star":   ["cos_theta_star", "cosThetaStar", "abs_cos_theta_star", "costhetastar", "cos_star", "DMeson_costhetastar"],
    "ipsig_min":  ["ipsig_min", "d0sig_min", "min_d0sig", "IPsig_min", "min_IPsig", "trk_d0sig_min"],
    "ipsig_max":  ["ipsig_max", "d0sig_max", "max_d0sig", "IPsig_max", "max_IPsig", "trk_d0sig_max"],
    "d0_min":     ["d0_min", "abs_d0_min", "min_abs_d0", "trk_d0_min"],
    "d0_max":     ["d0_max", "abs_d0_max", "max_abs_d0", "trk_d0_max"],
    "dedx_K":     ["dedx_K", "K_dedx", "dEdx_K", "trk_K_dedx", "dedx_k"],
    "p_K":        ["p_K", "K_p", "trk_K_p", "p_k"],
    "dedx_pi1":   ["dedx_pi1", "pi1_dedx", "dEdx_pi1", "dedx_pi", "pi_dedx", "dEdx_pi"],
    "p_pi1":      ["p_pi1", "pi1_p", "p_pi", "pi_p"],
    "dedx_pi2":   ["dedx_pi2", "pi2_dedx", "dEdx_pi2"],
    "p_pi2":      ["p_pi2", "pi2_p"],
    "d_ipsig":    ["d_ipsig", "D_ipsig", "Dmeson_impactSignificance", "impactSignificance", "DMeson_impactSignificance"],
    "trk_sig_d0": ["track_d0sig", "trk_d0sig", "d0sig"],   # 2D [N, ntrk] – alternatywa dla min/max
    "trk_d0":     ["track_d0", "trk_d0", "d0"],
}

# kolejność cech PRZYGOTOWANYCH (wejście sieci) – stała, używana też przez ONNX
PREPARED = [
    "lxy_sig_t",    # 0  asinh(Lxy/σ / 3)
    "point_t",      # 1  -log10(1 - cosθ_p + 1e-6), clip [0,6]
    "chi2_t",       # 2  log1p(χ²)
    "ipsig_min_t",  # 3  log1p(min IPsig córek)
    "ipsig_max_t",  # 4  log1p(max IPsig córek)
    "d0_min_t",     # 5  asinh(min|d0| / 0.05 mm)
    "cos_star",     # 6  |cosθ*|
    "kness_K",      # 7  (|z_π| - |z_K|) toru K, bramkowane pędem
    "kness_pi1",    # 8  (|z_K| - |z_π|) toru π1, bramkowane
    "kness_pi2",    # 9  j.w. π2 (0 dla D0)
    "logp_K",       # 10 ln p_K [GeV]
    "logp_pi_min",  # 11 ln min p_π
    "abs_eta",      # 12 |η| kandydata
    "d_ipsig_t",    # 13 log1p(|IPsig| kandydata D względem PV) – sygnał wskazuje na PV (mały), tło duży
]
MASS_SAFE = ["lxy_sig_t", "point_t", "chi2_t", "ipsig_min_t", "ipsig_max_t", "d0_min_t",
             "kness_K", "kness_pi1", "kness_pi2", "d_ipsig_t"]   # wariant B

def resolve(df: pd.DataFrame, colmap: dict | None = None) -> dict[str, str | None]:
    cm = dict(COLMAP)
    if colmap:
        for k, v in colmap.items():
            cm[k] = [v] if isinstance(v, str) else list(v) + cm.get(k, [])
    cols = set(df.columns)
    return {k: next((a for a in al if a in cols), None) for k, al in cm.items()}

def _get(df, res, key, default=np.nan):
    c = res.get(key)
    if c is None:
        return np.full(len(df), default, dtype=float)
    return df[c].to_numpy(dtype=float)

def dedx_gate(p_gev, p0=1.3, width=0.15):
    """Separacja K/π z dE/dx pikselowego zanika powyżej ~1.2–1.5 GeV → bramka sigmoidalna."""
    return 1.0 / (1.0 + np.exp((np.asarray(p_gev, float) - p0) / width))

def kaonness(dedx, p, want_mass, other_mass, scale=1.0):
    """Dodatnie gdy tor zgodny z hipotezą `want_mass` i niezgodny z `other_mass`."""
    zw = bb.pull(dedx, p, want_mass, scale)
    zo = bb.pull(dedx, p, other_mass, scale)
    k = np.clip(np.abs(zo) - np.abs(zw), -5, 5)
    k = np.where(np.isfinite(k), k, 0.0)
    return k * dedx_gate(p)

def prepare(df: pd.DataFrame, colmap: dict | None = None, dedx_scale: float = 1.0) -> pd.DataFrame:
    """df (kolumny H / h5 v1.1) -> macierz PREPARED (14 cech). Brakujące kolumny → 0 (neutralne)."""
    r = resolve(df, colmap)
    n = len(df)
    out = {}

    lxy_sig = _get(df, r, "lxy_sig")
    if np.all(np.isnan(lxy_sig)) and r["lxy"] and r["lxy_err"]:
        lxy_sig = _get(df, r, "lxy") / np.clip(_get(df, r, "lxy_err"), 1e-6, None)
    out["lxy_sig_t"] = np.arcsinh(np.nan_to_num(lxy_sig) / 3.0)

    cp = np.clip(np.nan_to_num(_get(df, r, "cos_point"), nan=0.0), -1, 1)
    out["point_t"] = np.clip(-np.log10(1.0 - cp + 1e-6), 0.0, 6.0)

    out["chi2_t"] = np.log1p(np.clip(np.nan_to_num(_get(df, r, "chi2"), nan=1.0), 0, 1e4))

    ipmin, ipmax = _get(df, r, "ipsig_min"), _get(df, r, "ipsig_max")
    if np.all(np.isnan(ipmin)) and r["trk_sig_d0"]:
        arr = np.abs(np.stack(df[r["trk_sig_d0"]].to_numpy()).astype(float))
        ipmin, ipmax = np.nanmin(arr, 1), np.nanmax(arr, 1)
    out["ipsig_min_t"] = np.log1p(np.abs(np.nan_to_num(ipmin)))
    out["ipsig_max_t"] = np.log1p(np.abs(np.nan_to_num(ipmax)))

    d0min = _get(df, r, "d0_min")
    if np.all(np.isnan(d0min)) and r["trk_d0"]:
        d0min = np.nanmin(np.abs(np.stack(df[r["trk_d0"]].to_numpy()).astype(float)), 1)
    out["d0_min_t"] = np.arcsinh(np.abs(np.nan_to_num(d0min)) / 0.05)

    out["cos_star"] = np.clip(np.abs(np.nan_to_num(_get(df, r, "cos_star"), nan=0.5)), 0, 1)

    def _kn(dk, pk, want, other):
        d, p = _get(df, r, dk), _get(df, r, pk)
        ok = np.isfinite(d) & np.isfinite(p) & (d > 0) & (p > 0)
        res = np.zeros(n)
        if ok.any():
            res[ok] = kaonness(d[ok], p[ok], want, other, dedx_scale)
        return res
    out["kness_K"] = _kn("dedx_K", "p_K", bb.M_K, bb.M_PI)
    out["kness_pi1"] = _kn("dedx_pi1", "p_pi1", bb.M_PI, bb.M_K)
    out["kness_pi2"] = _kn("dedx_pi2", "p_pi2", bb.M_PI, bb.M_K)

    pK = _get(df, r, "p_K"); p1 = _get(df, r, "p_pi1"); p2 = _get(df, r, "p_pi2")
    out["logp_K"] = np.where(np.isfinite(pK) & (pK > 0), np.log(np.clip(pK, 1e-3, None)), np.log(2.0))
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); ppi = np.nanmin(np.stack([p1, p2]), 0)
    out["logp_pi_min"] = np.where(np.isfinite(ppi) & (ppi > 0), np.log(np.clip(ppi, 1e-3, None)), np.log(1.5))

    out["abs_eta"] = np.abs(np.nan_to_num(_get(df, r, "eta"), nan=1.0))

    dip = _get(df, r, "d_ipsig")
    out["d_ipsig_t"] = np.log1p(np.abs(np.nan_to_num(dip, nan=2.5)))  # brak kolumny → wartość neutralna (≈próg)

    P = pd.DataFrame(out, index=df.index)[PREPARED].astype(np.float32)
    return P.replace([np.inf, -np.inf], 0.0).fillna(0.0)
