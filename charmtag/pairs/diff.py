"""Pary różniczkowe (etap PD, 2026-11-05): Δφ / Δy / znak (like-sign) par D–D — rozróżnienie DPS od SPS.

Rozszerzenie ADDYTYWNE, domyślnie wyłączone (`pairs.diff: null`). Włączone (`charmtag pairs --unblind --diff dphi` albo
`pairs.diff: {var: dphi, edges: [0, 1.5708, 3.1416]}`) dokłada do nóg kinematykę (η, φ, ładunek, ładunek córki 0 = K),
a pary — PO wecie wspólnych śladów, tą samą pętlą co histogramy nominalne — dzieli na biny zmiennej różniczkowej:
  dphi  |φ1−φ2| zawinięte do [0, π]           (DPS ≈ płasko; SPS: pik przy π, z „gluon splitting" także przy 0)
  dy    |η1−η2|  (≈ |Δy| dla pT ≫ m)           (DPS: szerokie; SPS: skupione)
  sign  like / opposite (znak kwarku c obu nóg: D⁺ → ładunek, D0 → −ładunek K córki 0; v1.1 bez ładunków córek →
        D0 nieokreślone → tylko D⁺D⁺)         (DPS: like ≈ opposite; SPS: prawie wyłącznie opposite)
Per bin: histogramy (n_cat, 25, 25) per kanał (`pairs_full_diff.npz`, klucz `<kanał>__<bin>`), zliczenia SB×SB, oczekiwania
(n_pair × f_i, f_i = udział binu w parach MC minbias Pythii — `expect: mc` — albo 1/n_bin — `expect: flat`), fit2d per bin
tym samym modelem (μ_i, Z_i, 68 %), łączenie po kanałach per bin i testy kształtu: profil μ_i przeliczony na dowolną hipotezę
udziałów h_i (Pythia = f_MC, DPS = płasko / 1:1) → Δχ² = 2·min_μ Σ_i ΔNLL_i(μ·h_i/f_i) − Σ_i 2·ΔNLL_i(min), p-value χ²(n_bin−1).
Nic w nominale (selekcja, szablony, κ, normalizacja całkowita) się nie zmienia; μ_i sumują się do nominału z dokładnością fitów.
"""
from __future__ import annotations

import json
import logging
import os
import time

import numpy as np
import pandas as pd
from scipy.stats import chi2 as _chi2

from . import CHANNELS
from .model import Inputs, Model

log = logging.getLogger("charmtag")

PI = float(np.pi)
VARS = ("sign", "dphi", "dy", "sign_dphi")
DEFAULT_EDGES = {"dphi": [0.0, PI / 2, PI], "dy": [0.0, 1.0, 5.0], "sign_dphi": [0.0, PI / 2, PI]}
DEFAULT_CHANNELS = {"dphi": ["D0D0", "DpD0"], "dy": ["D0D0", "DpD0"], "sign": None, "sign_dphi": None}   # None → wszystkie
# INT (literatura DPS, 2026-11-05): przy ~220 parach jedyna rozstrzygająca obserwabla to LS/OS (SPS → tylko OS, DPS → 1:1;
# LHCb LS/OS ≈ 0,10–0,13) → `sign` domyślny; Δφ (2 biny) < 1σ łącznie, sensowny tylko wewnątrz OS (`sign_dphi`).
DEFAULTS = {"var": "sign", "edges": None, "expect": "mc", "channels": "default", "f_min": 0.05, "mc_min_pairs": 10,
            "charge_source": None,                              # {klucz_próbki: plik/glob/lista h5 v1.2} → ładunek córki 0 (K) po `row`
            "d0_flavour": "closest_mass"}                       # closest_mass (K z hipotezy bliższej m0; build.d0_flavour_closest_mass) | hypothesis
KIN_COLS = ("eta", "phi", "charge", "dq0")                     # kolumny nóg dokładane przez scan_source(kin=True)
D_CHARGE_CANDS = ["Dmeson_daughter{i}_charge"]                 # v1.2; brak w v1.1 → NaN (nie dotyka cfg.columns → cache nietknięty)
LHCB_LS_OS = (0.10, 0.13)                                      # LHCb: stosunek LS/OS par D–D (odniesienie w raporcie)
GRID = np.round(np.arange(-1.5, 10.01, 0.25), 2)               # szerszy niż fit2d.GRID: μ_i w binie o małym f_i może być > 5
FINE = np.linspace(GRID[0], GRID[-1], int((GRID[-1] - GRID[0]) / 0.0025) + 1)


# ----------------------------------------------------------------------------- specyfikacja
def spec_from_cfg(p: dict, var_override=None):
    """`pairs.diff` z configu (+ nadpisanie zmiennej z CLI) → pełna specyfikacja albo None (wyłączone)."""
    d = p.get("diff")
    if var_override is None and not d:
        return None
    if d is True or d is None:
        d = {}
    if isinstance(d, str):
        d = {"var": d}
    s = {**DEFAULTS, **d}
    if var_override and var_override != "cfg":
        s["var"] = var_override
        if not (d or {}).get("edges") or (d or {}).get("var") != var_override:
            s["edges"] = None
    if s["var"] not in VARS:
        raise SystemExit(f"pairs.diff.var={s['var']!r} — dozwolone: {VARS}")
    if s["var"] == "sign":
        s["edges"] = None
    elif not s["edges"]:
        s["edges"] = list(DEFAULT_EDGES[s["var"]])
    else:
        s["edges"] = [float(x) for x in s["edges"]]
        if len(s["edges"]) < 2 or np.any(np.diff(s["edges"]) <= 0):
            raise SystemExit(f"pairs.diff.edges muszą rosnąć: {s['edges']}")
    if s["expect"] not in ("mc", "flat"):
        raise SystemExit("pairs.diff.expect: mc | flat")
    if s["d0_flavour"] not in ("closest_mass", "hypothesis"):
        raise SystemExit("pairs.diff.d0_flavour: closest_mass | hypothesis")
    if s["channels"] == "default":
        s["channels"] = DEFAULT_CHANNELS[s["var"]]
    if s["channels"] is None or s["channels"] == "all":
        s["channels"] = list(CHANNELS)
    bad = [c for c in s["channels"] if c not in CHANNELS]
    if bad:
        raise SystemExit(f"pairs.diff.channels: nieznane kanały {bad}")
    return s


def _edge_labels(e):
    return [f"[{a:.3g},{b:.3g})" if k < len(e) - 2 else f"[{a:.3g},{b:.3g}]" for k, (a, b) in enumerate(zip(e[:-1], e[1:]))]


def labels(spec) -> list[str]:
    if spec["var"] == "sign":
        return ["like", "opposite"]
    if spec["var"] == "sign_dphi":                                # like (łącznie) + opposite × Δφ
        return ["like"] + [f"opp:{x}" for x in _edge_labels(spec["edges"])]
    return _edge_labels(spec["edges"])


def nbins(spec) -> int:
    return len(labels(spec))


def flat_fractions(spec) -> np.ndarray:
    """Hipoteza „DPS": Δφ płasko (udział ∝ szerokość binu), znak 1:1 (sign_dphi: like ½, opposite ½ płasko w Δφ);
    dla dy — równy udział (tylko punkt odniesienia)."""
    if spec["var"] == "sign":
        return np.array([0.5, 0.5])
    e = np.asarray(spec["edges"], float)
    w = np.diff(e) / np.diff(e).sum()
    if spec["var"] == "sign_dphi":
        return np.r_[0.5, 0.5 * w]
    return w if spec["var"] == "dphi" else np.full(len(w), 1.0 / len(w))


# ----------------------------------------------------------------------------- zmienne par
def wrap_dphi(phi1, phi2):
    d = np.abs(np.asarray(phi1, float) - np.asarray(phi2, float)) % (2 * PI)
    return np.where(d > PI, 2 * PI - d, d)


def charm_sign(decay_type, charge, dq0):
    """Znak kwarku c nogi: +1 = c (D⁺, D0 → K⁻π⁺), −1 = c̄; 0 = nieokreślony (D0 bez ładunku K w v1.1)."""
    dt = np.asarray(decay_type)
    q = np.nan_to_num(np.asarray(charge, float), nan=0.0)
    k = np.nan_to_num(np.asarray(dq0, float), nan=0.0)
    k = np.where(np.abs(k) > 1.5, 0.0, k)                     # sentinel −999 (v1.2) → nieokreślony
    s = np.where(dt == 4, np.sign(q), np.where(dt == 10, -np.sign(k), 0.0))
    return s.astype(int)


def pair_vars(P: pd.DataFrame) -> dict:
    """Zmienne pary z ramki par (kolumny nóg z sufiksami 1/2): dphi, dy, like_sign (−1 nieokreślony, 0 opposite, 1 like)."""
    need = [f"{c}{s}" for c in KIN_COLS for s in (1, 2)]
    miss = [c for c in need if c not in P]
    if miss:
        raise SystemExit(f"pairs diff: brak kolumn kinematyki nóg {miss} — nogi zbudowane bez `scan_source(kin=True)`")
    dphi = wrap_dphi(P.phi1.to_numpy(), P.phi2.to_numpy())
    dy = np.abs(P.eta1.to_numpy(float) - P.eta2.to_numpy(float))
    s1 = charm_sign(P.decay_type1.to_numpy(), P.charge1.to_numpy(), P.dq01.to_numpy())
    s2 = charm_sign(P.decay_type2.to_numpy(), P.charge2.to_numpy(), P.dq02.to_numpy())
    like = np.where((s1 == 0) | (s2 == 0), -1, (s1 == s2).astype(int))
    return {"dphi": dphi, "dy": dy, "like_sign": like}


def diff_bin(P: pd.DataFrame, spec) -> np.ndarray:
    """Indeks binu różniczkowego każdej pary (−1 = poza zakresem / nieokreślony). Ostatnia krawędź domknięta (π, Δy_max)."""
    if not len(P):
        return np.zeros(0, int)
    v = pair_vars(P)
    if spec["var"] == "sign":
        return np.where(v["like_sign"] < 0, -1, np.where(v["like_sign"] == 1, 0, 1))
    x = v["dphi" if spec["var"] == "sign_dphi" else spec["var"]]
    e = np.asarray(spec["edges"], float)
    b = np.digitize(x, e) - 1
    b = np.where(x == e[-1], len(e) - 2, b)
    b = np.where((b >= 0) & (b < len(e) - 1) & np.isfinite(x), b, -1)
    if spec["var"] == "sign_dphi":
        ls = v["like_sign"]
        return np.where(ls < 0, -1, np.where(ls == 1, 0, np.where(b >= 0, b + 1, -1)))
    return b


def hkey(ch, k) -> str:
    return f"{ch}__{k}"


# ----------------------------------------------------------------------------- oczekiwania per bin
def fractions_mc(n_mc: dict, spec, ch) -> tuple[np.ndarray, np.ndarray]:
    """Udział binów w parach MC minbias (Pythia) + błąd binomialny; brak MC → NaN."""
    n = np.asarray(n_mc.get(ch, []), float)
    if n.size != nbins(spec) or n.sum() <= 0:
        return np.full(nbins(spec), np.nan), np.full(nbins(spec), np.nan)
    f = n / n.sum()
    return f, np.sqrt(np.clip(f * (1 - f) / n.sum(), 0, None))


def fit_fractions(spec, n_mc: dict, ch) -> tuple[np.ndarray, str]:
    """Udziały użyte do normalizacji μ_i (SS_i = μ_i · n_pair · f_i): MC, a gdy za mało par MC / pusty bin → płasko (z notką)."""
    flat = flat_fractions(spec)
    if spec["expect"] == "flat":
        return flat, "flat"
    f, _ = fractions_mc(n_mc, spec, ch)
    ntot = float(np.nansum(np.asarray(n_mc.get(ch, []), float)))
    if not np.all(np.isfinite(f)) or ntot < spec["mc_min_pairs"] or np.any(f < spec["f_min"]):
        log.warning("pairs diff %s: udziały MC (%s, N=%.0f) nieużyteczne → normalizacja płaska", ch, np.round(f, 3), ntot)
        return flat, "flat (fallback: MC N=%.0f, f=%s)" % (ntot, np.round(f, 3).tolist())
    return f, "mc"


def diff_inputs_json(I: Inputs, spec, sbsb_bins: dict, n_mc: dict, counts: dict) -> dict:
    """Blok `pair_inputs_diff.json`: per bin n_pair, n_real, sbsb (per kanał), udziały MC/fit/flat. Czytelny dla fit2d --diff."""
    nb = nbins(spec)
    tot_sb = sum(sum(sbsb_bins[k].values()) for k in range(nb))
    bins = []
    for k in range(nb):
        # udział kombinatoryki (SB×SB) w binie → skala startowa SB/BS (κ wolne absorbuje resztę); dolna granica 1/(2·n_bin),
        # żeby pusty bin SB×SB (mało statystyki) nie zerował składowych SB/BS modelu
        g = max(sum(sbsb_bins[k].values()) / tot_sb, 0.5 / nb) if tot_sb > 0 else 1.0 / nb
        b = {"label": labels(spec)[k], "g_sbsb": g, "n_real": {str(dt): v * g for dt, v in I.n_real.items()},
             "n_pair": {}, "f_fit": {}, "f_mode": {}, "f_mc": {}, "f_mc_err": {}, "f_flat": float(flat_fractions(spec)[k]),
             "sbsb": {f"{ch}_{i}_{j}": int(v) for (ch, i, j), v in sbsb_bins[k].items()},
             "pairs_kept": int(counts.get("diff_kept", {}).get(str(k), 0))}
        for ch in CHANNELS:
            ff, mode = fit_fractions(spec, n_mc, ch)
            fm, fe = fractions_mc(n_mc, spec, ch)
            b["n_pair"][ch] = float(I.n_pair[ch] * ff[k])
            b["f_fit"][ch], b["f_mode"][ch] = float(ff[k]), mode
            b["f_mc"][ch], b["f_mc_err"][ch] = (float(fm[k]) if np.isfinite(fm[k]) else None), (float(fe[k]) if np.isfinite(fe[k]) else None)
        bins.append(b)
    return {"spec": spec, "labels": labels(spec), "nbins": nb, "n_mc_pairs": {ch: [int(x) for x in n_mc.get(ch, [])] for ch in CHANNELS},
            "n_mc_undefined": {ch: int(v) for ch, v in (counts.get("mc_undefined") or {}).items()},
            "pairs_undefined": int(counts.get("diff_undefined", 0)), "bins": bins}


def inputs_for_bin(I: Inputs, D: dict, k: int) -> Inputs:
    """Kopia wejść nominalnych z normalizacją i SB×SB binu k (szablony, Y, dataleg, nev — wspólne)."""
    b = D["bins"][k]
    sbsb = {}
    for key, v in b["sbsb"].items():
        ch, i, j = key.split("_")
        sbsb[(ch, int(i), int(j))] = int(v)
    return Inputs(edges=I.edges, ncat=I.ncat, sig=I.sig, sig_comp=I.sig_comp, Y=I.Y, dataleg=I.dataleg, nev_legs=I.nev_legs,
                  ndata=I.ndata, n_real={int(dt): float(v) for dt, v in b["n_real"].items()}, n_pair=dict(b["n_pair"]), sbsb=sbsb,
                  meta={**I.meta, "diff_bin": k, "diff_label": b["label"], "diff_spec": D["spec"]})


# ----------------------------------------------------------------------------- fit per bin + testy kształtu
def _interp_profile(scan: dict):
    mus = np.array(sorted(float(k) for k in scan))
    d = np.array([scan[f"{m:g}"] for m in mus])
    d = d - min(d.min(), 0.0)
    return mus, d


def profile_at(scan: dict, mu):
    """ΔNLL profilu w μ (liniowo między punktami siatki; poza siatką — ekstrapolacja ostatniego odcinka, ≥ 0)."""
    mus, d = _interp_profile(scan)
    mu = np.asarray(mu, float)
    out = np.interp(mu, mus, d)
    lo, hi = mu < mus[0], mu > mus[-1]
    if lo.any():
        out = np.where(lo, d[0] + (d[0] - d[1]) / (mus[0] - mus[1]) * (mu - mus[0]), out)
    if hi.any():
        out = np.where(hi, d[-1] + (d[-1] - d[-2]) / (mus[-1] - mus[-2]) * (mu - mus[-1]), out)
    return np.clip(out, 0.0, None)


def combine_bins(prof: dict, scale: dict | None = None):
    """Wspólne μ dla zbioru profili {klucz: scan}; `scale[klucz]` = h_i/f_i (hipoteza kształtu): ΔNLL_i(μ·scale_i).
    Zwraca μ̂, 68 %, q = 2·min ΣΔNLL (względem niezależnych minimów), Z względem μ=0."""
    tot = np.zeros_like(FINE)
    for key, scan in prof.items():
        s = 1.0 if scale is None else float(scale[key])
        tot += profile_at(scan, FINE * s)
    q = 2 * float(tot.min())
    tot0 = tot - tot.min()
    mu_c = float(FINE[np.argmin(tot0)])
    ok = FINE[tot0 <= 0.5]
    q0 = 2 * float(np.interp(0.0, FINE, tot0)) if mu_c > 0 else 0.0
    return {"mu": mu_c, "mu_68": [float(ok.min()), float(ok.max())], "q_shape": q, "Z": float(np.sqrt(max(q0, 0.0)))}


def shape_tests(prof_ch: dict, D: dict, ch, spec) -> dict:
    """Testy kształtu w kanale: hipotezy udziałów h (pythia = f_MC, dps = płasko) vs μ_i wolne; Δχ² ~ χ²(n_bin−1)."""
    nb = D["nbins"]
    out = {}
    hyps = {"dps_flat": flat_fractions(spec)}
    fm = np.array([D["bins"][k]["f_mc"][ch] if D["bins"][k]["f_mc"][ch] is not None else np.nan for k in range(nb)])
    if np.all(np.isfinite(fm)) and fm.min() > 0:
        hyps["pythia_mc"] = fm
    ff = np.array([D["bins"][k]["f_fit"][ch] for k in range(nb)])
    for name, h in hyps.items():
        sc = {k: float(h[k] / ff[k]) for k in prof_ch}
        r = combine_bins(prof_ch, sc)
        ndf = max(len(prof_ch) - 1, 1)
        out[name] = {"h": [float(x) for x in h], "mu_common": r["mu"], "mu_common_68": r["mu_68"], "chi2": r["q_shape"], "ndf": ndf,
                     "p": float(_chi2.sf(r["q_shape"], ndf))}
    return out


def run_diff(I: Inputs, D: dict, Hd: dict, out: str, ntoy: int = 0, seed: int = 2026, channels=None) -> dict:
    """Fit odsłonięty per bin (plan B, pełna płaszczyzna), łączenie po kanałach per bin, testy kształtu per kanał i łącznie."""
    from . import fit2d as F
    spec = D["spec"]
    nb = D["nbins"]
    channels = list(spec["channels"] if channels is None else channels)
    rep = {"mode": "unblind-diff", "spec": spec, "labels": D["labels"], "channels_fitted": channels, "channels": {}, "bins": {}, "tests": {}}
    models = [Model(inputs_for_bin(I, D, k)) for k in range(nb)]
    prof = {}                                                   # (ch, k) → profil
    for ch in channels:
        rc = {"bins": {}}
        for k in range(nb):
            key = hkey(ch, k)
            n = Hd.get(key)
            lab = D["labels"][k]
            if n is None or n.sum() < 20 or D["bins"][k]["n_pair"][ch] <= 0:
                rc["bins"][lab] = {"skipped": f"za mało par ({0 if n is None else n.sum():.0f}) albo S_oczek = 0"}
                log.warning("fit2d DIFF %s %s: %s", ch, lab, rc["bins"][lab]["skipped"])
                continue
            t0 = time.time()
            cf = models[k].channel(ch)
            m1 = cf.fit(n)
            m0 = cf.fit(n, fix_mu=0.0, start=np.array(m1.values))
            r = {"mu_hat": float(m1.values["mu"]), "mu_err_hesse": float(m1.errors["mu"]), "Z": F.z_from(m1, m0), "valid": bool(m1.valid),
                 "S_expected_full_plane": float(sum(v["SS"] for v in cf.exp.values())), "n_pairs": float(n.sum()),
                 "f_fit": D["bins"][k]["f_fit"][ch], "f_mode": D["bins"][k]["f_mode"][ch], "f_mc": D["bins"][k]["f_mc"][ch],
                 "f_mc_err": D["bins"][k]["f_mc_err"][ch], "f_flat": D["bins"][k]["f_flat"]}
            dev, ndf = cf.deviance(m1, n)
            r.update({"dev": dev, "ndf": ndf, "dev_ndf": dev / max(ndf, 1)})
            r["kappa"] = {nm: [float(m1.values[nm]), float(m1.errors[nm])] for nm in cf.names if nm.startswith("kap")}
            r["profile_dNLL"], r["mu_profile_68"] = F.profile(cf, n, np.ones(n.shape[1:], bool), m1, grid=GRID)
            if ntoy:
                r["toys_at_fit"] = F.toys(models[k], ch, cf.lam(np.array(m1.values)), ntoy, seed + k, mu_true=r["mu_hat"])
            prof[(ch, k)] = r["profile_dNLL"]
            rc["bins"][lab] = r
            log.info("fit2d DIFF %s %s: μ=%.2f [%.2f, %.2f] Z=%.2f S_oczek=%.1f f_%s=%.3f (%.0f s)", ch, lab, r["mu_hat"],
                     *r["mu_profile_68"], r["Z"], r["S_expected_full_plane"], r["f_mode"].split()[0], r["f_fit"], time.time() - t0)
        pc = {k: prof[(ch, k)] for k in range(nb) if (ch, k) in prof}
        if len(pc) >= 2:
            rc["tests"] = shape_tests(pc, D, ch, spec)
            rc["data_fractions"] = data_fractions(rc["bins"], D, ch)
            if nb == 2 and len(pc) == 2:
                rc["ratio_0_over_1"] = ratio_of_bins(rc["bins"][D["labels"][0]], rc["bins"][D["labels"][1]])
        rep["channels"][ch] = rc
    # per bin: łączenie po kanałach (jak fit2d.combine) → μ_bin, Z_bin
    for k in range(nb):
        pk = {ch: prof[(ch, k)] for ch in channels if (ch, k) in prof}
        if pk:
            rep["bins"][D["labels"][k]] = {**combine_bins(pk), "channels": list(pk)}
    # łącznie: (a) μ per bin wspólne w kanałach [nb parametrów] vs (b) jedno μ globalne z kształtem h [1 parametr]
    if len(prof) >= 2:
        tests = {}
        for name in ("pythia_mc", "dps_flat"):
            sc = {}
            for (ch, k) in prof:
                t = rep["channels"][ch].get("tests", {}).get(name)
                if t is None:
                    break
                sc[(ch, k)] = t["h"][k] / D["bins"][k]["f_fit"][ch]
            else:
                if sc:
                    r = combine_bins(prof, sc)
                    # model zagnieżdżony: μ_k per bin wspólne w kanałach, TA SAMA hipoteza kształtu h (SS = μ_k·N_par·h_k) → nb parametrów
                    ks = sorted({k for (_, k) in prof})
                    q_bins = sum(combine_bins({key: prof[key] for key in prof if key[1] == k}, sc)["q_shape"] for k in ks)
                    dchi2 = max(r["q_shape"] - q_bins, 0.0)
                    tests[name] = {"mu_global": r["mu"], "mu_global_68": r["mu_68"], "Z_global": r["Z"], "chi2_global": r["q_shape"],
                                   "chi2_bins_common": q_bins, "dchi2": dchi2, "ndf": max(len(ks) - 1, 1),
                                   "p_flat_in_bins": float(_chi2.sf(dchi2, max(len(ks) - 1, 1)))}
        rep["tests"] = tests
        tot = combine_bins(prof)                                  # bez kształtu: μ wspólne we wszystkich (ch, bin) — kontrola ≈ nominał
        rep["combined_all_bins"] = {"mu": tot["mu"], "mu_68": tot["mu_68"], "Z": tot["Z"], "note": "μ wspólne (kanały × biny) przy f_fit — kontrola vs nominał"}
    if spec["var"].startswith("sign"):
        rep["like_sign"] = ls_summary(rep, D)
    F._js(os.path.join(out, "fit2d_diff.json"), rep)
    return rep


def _z_from_p(p):
    """p-value testu χ² → równoważna istotność Z ≥ 0 (χ² o 1 st. sw.: Z = √χ²)."""
    return float(np.sqrt(_chi2.isf(min(max(p, 1e-300), 1.0), 1)))


def ls_summary(rep: dict, D: dict) -> dict:
    """INT: f_LS = N_LS/(N_LS+N_OS) w sygnale danych (z μ_i i oczekiwań Pythii per klasa), Pythia f_LS^MC ± binomialnie,
    LHCb LS/OS 0,10–0,13 → f_LS 0,09–0,12, Z odrzucenia „f_LS = Pythia" (z testu kształtu pythia_mc: Z = √χ² przy 1 st. sw.)."""
    labs = D["labels"]
    like = labs[0]
    out = {"lhcb_ls_os": list(LHCB_LS_OS), "lhcb_f_ls": [r / (1 + r) for r in LHCB_LS_OS], "channels": {}}
    for ch, rc in rep["channels"].items():
        df = rc.get("data_fractions", {})
        if like not in df or df[like][0] is None:
            continue
        b = rc["bins"]
        n_ls = b[like]["mu_hat"] * D["bins"][0]["n_pair"][ch]
        n_os = sum(b[lab]["mu_hat"] * D["bins"][k]["n_pair"][ch] for k, lab in enumerate(labs) if k > 0 and "mu_hat" in b.get(lab, {}))
        fmc = D["bins"][0]["f_mc"][ch]
        e = {"N_LS_sig": float(n_ls), "N_OS_sig": float(n_os), "f_LS": df[like], "f_LS_pythia": [fmc, D["bins"][0]["f_mc_err"][ch]],
             "ls_os": (float(n_ls / n_os) if n_os > 0 else None)}
        t = rc.get("tests", {}).get("pythia_mc")
        if t:
            e["Z_reject_pythia"] = _z_from_p(t["p"])
        t = rc.get("tests", {}).get("dps_flat")
        if t:
            e["Z_reject_dps"] = _z_from_p(t["p"])
        out["channels"][ch] = e
    for name in ("pythia_mc", "dps_flat"):
        t = rep.get("tests", {}).get(name)
        if t:
            out[f"Z_reject_{name}_combined"] = _z_from_p(t["p_flat_in_bins"])
    return out


def data_fractions(bins_rep: dict, D: dict, ch) -> dict:
    """Udziały binów w sygnale danych: f_i^dane = μ_i f_i / Σ μ_j f_j (błąd z Hesse, niezależne biny, propagacja liniowa)."""
    labs = [lab for lab in D["labels"] if "mu_hat" in bins_rep.get(lab, {})]
    s = np.array([bins_rep[lab]["mu_hat"] * bins_rep[lab]["f_fit"] for lab in labs])
    e = np.array([bins_rep[lab]["mu_err_hesse"] * bins_rep[lab]["f_fit"] for lab in labs])
    tot = s.sum()
    out = {}
    if tot <= 0:
        return {lab: [None, None] for lab in labs}
    for i, lab in enumerate(labs):
        f = s[i] / tot
        grad = np.array([-(s[i] / tot ** 2)] * len(s))
        grad[i] = (tot - s[i]) / tot ** 2
        out[lab] = [float(f), float(np.sqrt(np.sum((grad * e) ** 2)))]
    return out


def ratio_of_bins(r0: dict, r1: dict) -> dict:
    """μ_0/μ_1 (normalizacja f_fit) z błędem z Hesse; przy μ_1 ≤ 0 → None."""
    if r1["mu_hat"] <= 0 or r0["mu_hat"] <= 0:
        return {"ratio": None}
    ratio = r0["mu_hat"] / r1["mu_hat"]
    rel = np.sqrt((r0["mu_err_hesse"] / r0["mu_hat"]) ** 2 + (r1["mu_err_hesse"] / r1["mu_hat"]) ** 2)
    return {"ratio": float(ratio), "err": float(ratio * rel),
            "ratio_data_fractions": float((r0["mu_hat"] * r0["f_fit"]) / (r1["mu_hat"] * r1["f_fit"]))}


# ----------------------------------------------------------------------------- raport
def write_report(out: str, I: Inputs, D: dict, rep: dict, nominal: dict | None = None, path: str | None = None) -> str:
    """tfin_diff.md: tabela per kanał × bin, łączenie per bin, testy kształtu (Pythia / DPS), kontrola sumy vs nominał."""
    spec, labs = D["spec"], D["labels"]
    name = {"dphi": "Δφ", "dy": "Δy≈|Δη|", "sign": "znak (like/opposite)"}[spec["var"]]
    L = [f"# Pary D–D różniczkowo — {name} (ODSŁONIĘTE, etap PD)", "",
         f"- Zmienna: `{spec['var']}`, biny: {', '.join(labs)}; normalizacja μ_i: `{spec['expect']}` (SS_i = μ_i · N_par · f_i); "
         f"kanały fitowane: {', '.join(rep['channels_fitted'])}.",
         "- Pary MC minbias (Pythia) per bin: " + "; ".join(f"{ch} {D['n_mc_pairs'].get(ch)}" for ch in CHANNELS)
         + f"; nieokreślone w danych (poza binami): {D['pairs_undefined']}.",
         "- Nominał (łącznie, bez binów): " + (f"wspólne μ = {nominal['combined']['mu_common']:.2f} "
                                              f"[{nominal['combined']['mu_common_68'][0]:.2f}, {nominal['combined']['mu_common_68'][1]:.2f}], "
                                              f"Z = {nominal['combined']['Z']:.2f}σ" if nominal else "— (brak fit2d_unblind.json)") + ".", ""]
    L += ["| kanał | bin | N par | S oczek. | f_fit (tryb) | f_MC ± | μ_i | 68 % | Z_i | dev/ndf | f_i dane ± |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for ch, rc in rep["channels"].items():
        df = rc.get("data_fractions", {})
        for lab in labs:
            r = rc["bins"].get(lab, {})
            if "mu_hat" not in r:
                L.append(f"| {ch} | {lab} | — | — | — | — | {r.get('skipped', '—')} | | | | |")
                continue
            fmc = f"{r['f_mc']:.2f} ± {r['f_mc_err']:.2f}" if r.get("f_mc") is not None else "—"
            fd = df.get(lab, [None, None])
            fds = f"{fd[0]:.2f} ± {fd[1]:.2f}" if fd[0] is not None else "—"
            L.append(f"| {ch} | {lab} | {r['n_pairs']:.0f} | {r['S_expected_full_plane']:.1f} | {r['f_fit']:.2f} ({r['f_mode'].split()[0]}) | {fmc} | "
                     f"{r['mu_hat']:.2f} ± {r['mu_err_hesse']:.2f} | {r['mu_profile_68'][0]:.2f}–{r['mu_profile_68'][1]:.2f} | {r['Z']:.2f} | "
                     f"{r['dev_ndf']:.2f} | {fds} |")
        if "ratio_0_over_1" in rc and rc["ratio_0_over_1"].get("ratio") is not None:
            q = rc["ratio_0_over_1"]
            L.append(f"| {ch} | μ({labs[0]})/μ({labs[1]}) | | | | | {q['ratio']:.2f} ± {q['err']:.2f} | | | | |")
        for tn, t in rc.get("tests", {}).items():
            L.append(f"| {ch} | test `{tn}` h={np.round(t['h'], 2).tolist()} | | | | | μ={t['mu_common']:.2f} | "
                     f"{t['mu_common_68'][0]:.2f}–{t['mu_common_68'][1]:.2f} | | χ²/ndf={t['chi2']:.2f}/{t['ndf']} | p={t['p']:.2f} |")
    L.append("")
    if rep.get("bins"):
        L += ["Łączenie po kanałach per bin (suma profili, jak nominalne `combine`):", "", "| bin | kanały | μ | 68 % | Z |", "|---|---|---|---|---|"]
        for lab, b in rep["bins"].items():
            L.append(f"| {lab} | {', '.join(b['channels'])} | {b['mu']:.2f} | {b['mu_68'][0]:.2f}–{b['mu_68'][1]:.2f} | {b['Z']:.2f} |")
        L.append("")
    for tn, t in rep.get("tests", {}).items():
        L.append(f"- Test łączny `{tn}`: μ globalne = {t['mu_global']:.2f} [{t['mu_global_68'][0]:.2f}, {t['mu_global_68'][1]:.2f}] (Z = {t['Z_global']:.2f}σ); "
                 f"Δχ²(μ per bin wspólne w kanałach vs jedno μ z kształtem) = {t['dchi2']:.2f} / ndf {t['ndf']} → p = {t['p_flat_in_bins']:.2f}.")
    ls = rep.get("like_sign")
    if ls:
        L += ["", f"Like-sign (INT): f_LS = N_LS/(N_LS+N_OS) w sygnale; LHCb LS/OS {ls['lhcb_ls_os'][0]:.2f}–{ls['lhcb_ls_os'][1]:.2f} "
              f"→ f_LS {ls['lhcb_f_ls'][0]:.3f}–{ls['lhcb_f_ls'][1]:.3f}; DPS czyste → 0,5; SPS czyste → 0.", "",
              "| kanał | N_LS sig | N_OS sig | f_LS dane ± | f_LS Pythia ± | Z odrzucenia Pythii | Z odrzucenia DPS |", "|---|---|---|---|---|---|---|"]
        for ch, e in ls["channels"].items():
            fp = f"{e['f_LS_pythia'][0]:.3f} ± {e['f_LS_pythia'][1]:.3f}" if e["f_LS_pythia"][0] is not None else "—"
            L.append(f"| {ch} | {e['N_LS_sig']:.1f} | {e['N_OS_sig']:.1f} | {e['f_LS'][0]:.3f} ± {e['f_LS'][1]:.3f} | {fp} | "
                     f"{e.get('Z_reject_pythia', float('nan')):.2f} | {e.get('Z_reject_dps', float('nan')):.2f} |")
        for k in ("Z_reject_pythia_mc_combined", "Z_reject_dps_flat_combined"):
            if k in ls:
                L.append(f"- łącznie (kanały): {k} = {ls[k]:.2f}σ")
    c = rep.get("combined_all_bins")
    if c:
        L.append(f"- Kontrola: μ wspólne we wszystkich (kanał × bin) = {c['mu']:.2f} [{c['mu_68'][0]:.2f}, {c['mu_68'][1]:.2f}], Z = {c['Z']:.2f}σ "
                 f"(oczekiwane ≈ nominał dla kanałów fitowanych).")
    L += ["", "Interpretacja: `pythia_mc` = zgodność kształtu z Pythią (SPS+MPI), `dps_flat` = zgodność z czystym DPS (Δφ płasko / like:opposite 1:1). "
          "Małe p odrzuca hipotezę. f_MC z N par MC minbias — jego błąd nie jest wliczony w μ_i (patrz kolumna f_MC ±)."]
    path = path or os.path.join(out, "tfin_diff.md")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    return path


def load_diff(out: str) -> dict:
    p = os.path.join(out, "pair_inputs_diff.json")
    if not os.path.exists(p):
        raise SystemExit(f"brak {p}: najpierw `charmtag pairs --unblind --diff ...`")
    return json.load(open(p))
