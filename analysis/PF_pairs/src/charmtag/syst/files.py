"""Jackknife po plikach danych i podziały danych — ŚLEPO (etap U).

Pliki = części h5 data_none (1 plik = 1 job produkcji; na lxplus 39). Noga → plik przez `cand` (= wiersz globalny w
sklejonej liście plików, jak charmtag.io) i offsety z liczby wierszy plików (h5py).
Ślepo liczone:
  1. nadmierna dyspersja: zliczenia par w ramionach / sidebandach per plik vs oczekiwanie ∝ zdarzeniom pliku
     (χ²/ndf → φ; oczekiwane σ_jackknife/σ_fit ≈ √φ, w 09-28 zmierzone 1,10–1,23);
  2. jackknife czułości: bez pliku k → nuisance'y z obszaru ślepego, Asimov → Z_k, σ_μ,k (wpływ pliku na czułość);
  3. podziały (połowy, naprzemiennie): κ z obu połówek (pull), zliczenia w ramionach (pull), Z/σ_μ połówek.
Odsłonięte Δμ (jackknife μ, pulle połówek) — tylko `unblind=True` i CHARM_UNBLIND=1 (decyzja o odsłonięciu).
"""
from __future__ import annotations

import copy
import logging
import os
import re

import numpy as np
import pandas as pd

from ..pairs import CHANNELS
from ..pairs.hists import build_pair_hists
from ..pairs.model import Model, region_mask
from . import core

log = logging.getLogger("charmtag")


def file_table(cfg, key="data_none"):
    """[(plik, job, wiersze, zdarzenia)] w kolejności sklejania (jak charmtag.config.sample_paths)."""
    import h5py
    from .. import config as C
    ds = cfg.get("io", {}).get("dataset", "events")
    rows = []
    for p in C.sample_paths(cfg, key):
        with h5py.File(p, "r") as f:
            n = int(f[ds].shape[0])
            ev = f.attrs.get("n_events_read")
        m = [x for x in re.findall(r"_(\d{6})(?=_)", os.path.basename(p)) if x.startswith("0")]
        rows.append({"file": os.path.basename(p), "job": m[-1] if m else os.path.basename(p), "rows": n,
                     "events": int(ev) if ev is not None else None})
    return pd.DataFrame(rows)


def assign_files(L: pd.DataFrame, ft: pd.DataFrame):
    off = np.r_[0, np.cumsum(ft.rows.to_numpy())]
    k = np.searchsorted(off, L.cand.to_numpy(), "right") - 1
    if (k < 0).any() or (k >= len(ft)).any():
        raise ValueError("cand poza zakresem wierszy plików — inna lista plików niż przy `charmtag pairs`?")
    return L.assign(file_idx=k)


def per_file_hists(L, I, nfile, region="blind"):
    H = {ch: [] for ch in CHANNELS}
    for k in range(nfile):
        h, _ = build_pair_hists(L[L.file_idx == k], I, region)
        for ch in CHANNELS:
            H[ch].append(h[ch])
    return {ch: np.array(v) for ch, v in H.items()}


def scaled_inputs(I, f):
    """Wejścia dla podpróbki zawierającej ułamek f zdarzeń (normalizacje ∝ f; gęstość fałszywych nóg bez zmian)."""
    I2 = copy.copy(I)
    I2.ndata = I.ndata * f
    I2.n_pair = {k: v * f for k, v in I.n_pair.items()}
    I2.n_real = {k: v * f for k, v in I.n_real.items()}
    I2.sbsb = dict(I.sbsb)
    return I2


def overdispersion(Hf, I, w):
    out = {}
    for ch in CHANNELS:
        for reg in ("arms", "sbsb"):
            mk = region_mask(I, ch, reg)
            n = Hf[ch][:, :, mk].sum(axis=(1, 2))
            e = w * n.sum()
            ok = e > 0
            chi2 = float(np.sum((n[ok] - e[ok]) ** 2 / e[ok]))
            ndf = int(ok.sum() - 1)
            out[f"{ch}_{reg}"] = {"chi2": chi2, "ndf": ndf, "phi": chi2 / max(ndf, 1), "n_total": int(n.sum()),
                                  "max_pull": float(np.max(np.abs((n[ok] - e[ok]) / np.sqrt(e[ok])))) if ok.any() else None}
    return out


def _sens(M, H, ch):
    if H[ch].sum() < 20:
        return None
    cf, x, src = core.base_point(M, H, ch)
    a = core.asimov_fit(cf, cf.lam(x), x)
    kap = {nm: [float(x[k]), None] for k, nm in enumerate(cf.names) if nm.startswith("kap")}
    return {"Z": a["Z"], "err": a["err"], "kappa": kap, "x": x, "cf": cf}


def _kappa_errors(cf, H, ch, x):
    m = cf.fit(H[ch], region_mask(cf.I, ch, "blind"), fix_mu=1.0, start=x)
    return {nm: [float(m.values[nm]), float(m.errors[nm])] for nm in cf.names if nm.startswith("kap")}


def run_files(cfg, I, L: pd.DataFrame, ft: pd.DataFrame, jackknife=True, splits=True, channels=None, unblind=False):
    nf = len(ft)
    rep = {"n_files": nf, "files": ft.to_dict("records")}
    if nf < 2:
        rep["skipped"] = "jeden plik danych — jackknife/podziały tylko na lxplus (części h5)"
        log.warning("syst files: %s", rep["skipped"])
        return rep
    L = assign_files(L, ft)
    ev = ft.events.to_numpy(float) if ft.events.notna().all() else ft.rows.to_numpy(float)
    w = ev / ev.sum()
    Hf = per_file_hists(L, I, nf, "blind")
    Ht = {ch: v.sum(0) for ch, v in Hf.items()}
    rep["overdispersion"] = od = overdispersion(Hf, I, w)
    for k, v in od.items():
        log.info("syst files: dyspersja %-12s χ²/ndf = %.1f/%d (φ=%.2f → σ_jk/σ ≈ %.2f)", k, v["chi2"], v["ndf"], v["phi"], np.sqrt(max(v["phi"], 1)))
    chs = channels or list(CHANNELS)
    M = Model(I)
    nom = {ch: _sens(M, Ht, ch) for ch in chs}
    rep["nominal"] = {ch: ({"Z": v["Z"], "err": v["err"]} if v else None) for ch, v in nom.items()}
    if jackknife:
        jk = {ch: [] for ch in chs}
        for k in range(nf):
            Mk = Model(scaled_inputs(I, 1 - w[k]))
            Hk = {ch: Ht[ch] - Hf[ch][k] for ch in CHANNELS}
            for ch in chs:
                r = _sens(Mk, Hk, ch)
                jk[ch].append(None if r is None else [r["Z"], r["err"]])
            log.info("syst jackknife: bez pliku %s → %s", ft.job[k], {ch: (None if jk[ch][-1] is None else round(jk[ch][-1][0], 2)) for ch in chs})
        rep["jackknife_sensitivity"] = {}
        for ch in chs:
            z = np.array([v[0] for v in jk[ch] if v is not None])
            if not len(z):
                continue
            n = len(z)
            rep["jackknife_sensitivity"][ch] = {"Z_drop": z.tolist(), "Z_mean": float(z.mean()),
                                                "Z_jk_sd": float(np.sqrt((n - 1) / n * np.sum((z - z.mean()) ** 2))),
                                                "max_drop_file": str(ft.job[int(np.argmin(z))]), "Z_min": float(z.min())}
    if splits:
        rep["splits"] = {}
        half = np.arange(nf) < nf // 2
        for nm, sel in (("halves", half), ("alternating", np.arange(nf) % 2 == 0)):
            r = {}
            fa = w[sel].sum()
            for ch in chs:
                HA = {c: Hf[c][sel].sum(0) for c in CHANNELS}
                HB = {c: Hf[c][~sel].sum(0) for c in CHANNELS}
                MA, MB = Model(scaled_inputs(I, fa)), Model(scaled_inputs(I, 1 - fa))
                sa, sb = _sens(MA, HA, ch), _sens(MB, HB, ch)
                if not sa or not sb:
                    continue
                ka, kb = _kappa_errors(sa["cf"], HA, ch, sa["x"]), _kappa_errors(sb["cf"], HB, ch, sb["x"])
                pulls = {k: (ka[k][0] - kb[k][0]) / np.hypot(ka[k][1], kb[k][1]) for k in ka if ka[k][1] and kb[k][1]}
                mk = region_mask(I, ch, "arms")
                na, nb = HA[ch][:, mk].sum(), HB[ch][:, mk].sum()
                ea = fa * (na + nb)
                r[ch] = {"Z": [sa["Z"], sb["Z"]], "sigma_mu": [sa["err"], sb["err"]], "kappa_pulls": pulls,
                         "kappa_chi2_ndf": [float(sum(p * p for p in pulls.values())), len(pulls)],
                         "arms_counts": [int(na), int(nb)], "arms_pull": float((na - ea) / np.sqrt(ea * (1 - fa)))}
                log.info("syst split %s %s: Z %.2f/%.2f, κ χ²=%.1f/%d, ramiona pull %.2f", nm, ch, sa["Z"], sb["Z"],
                         *r[ch]["kappa_chi2_ndf"], r[ch]["arms_pull"])
            rep["splits"][nm] = {"fraction_A": float(fa), "channels": r}
    if unblind:
        if os.environ.get("CHARM_UNBLIND") != "1":
            raise SystemExit("odsłonięcie jackknife wymaga CHARM_UNBLIND=1 (decyzja o odsłonięciu)")
        rep["unblinded"] = unblinded_jackknife(cfg, I, L, nf, w, Hf_blind=None, chs=chs, ft=ft)
    return rep


def unblinded_jackknife(cfg, I, L, nf, w, Hf_blind=None, chs=None, ft=None):
    """ODSŁONIĘTE: μ bez pliku k (pełna płaszczyzna). Wywoływane tylko po decyzji o odsłonięciu."""
    from ..pairs.fit2d import fit_channel
    Hf = per_file_hists(L, I, nf, "full")
    Ht = {ch: v.sum(0) for ch, v in Hf.items()}
    M = Model(I)
    out = {}
    for ch in chs:
        full = np.ones(Ht[ch].shape[1:], bool)
        r0, _, _ = fit_channel(M, ch, Ht[ch], full, scan=False)
        mus = []
        for k in range(nf):
            Mk = Model(scaled_inputs(I, 1 - w[k]))
            rk, _, _ = fit_channel(Mk, ch, Ht[ch] - Hf[ch][k], full, scan=False)
            mus.append(rk["mu_hat"])
        mus = np.array(mus)
        n = len(mus)
        sjk = float(np.sqrt((n - 1) / n * np.sum((mus - mus.mean()) ** 2)))
        out[ch] = {"mu": r0["mu_hat"], "sigma_fit": r0["mu_err_hesse"], "sigma_jk": sjk, "ratio": sjk / r0["mu_err_hesse"],
                   "mu_drop": mus.tolist(), "max_shift_file": str(ft.job[int(np.argmax(np.abs(mus - r0["mu_hat"])))])}
    return out
