#!/usr/bin/env python3
"""PREZ-2 / pz_ext (19.11.2026, projekt A) — ekstrakcja pojedynczych D (D0→Kπ, D+→Kππ) z MC h5 v1.1 (ta sama produkcja i config co nominalny
XGB H: runs/lxplus_v11/models/<gatunek>_H/config_used.yaml) do porównania klasyfikatorów drzewo vs sieć × wejścia low-level vs high-level.

Próbki (TYLKO MC, żadnych danych):
  --sample fd : sygnał = truth ∈ signal.truth_values (0/1; D0: swap closest_mass → truth+100 wyłączone), fd_nd (signal_mc)
  --sample mb : tło     = truth 8–99 (kombinatoryka), mb_nd (incl_mc)
Selekcja jak nominał: in_range (zakres mas, pT 3–30) ∧ weto D* (D+). Zapisywany cały zakres mas (okno sygnału = trening/ewaluacja,
pasma poza oknem = tylko test rzeźbienia masy). Podpróbka deterministyczna (keep_frac z liczby kandydatów → ~--cap wierszy).

Wyjście npz: Xh [n,18] (cechy H w kolejności meta.json nominału), T [n,nd,11] (ślady: low-level), G [n,1] (L_xy), m, pt, abs_eta, weight,
run, event, truth, fold (charmtag.samples.fold_of, 5 foldów = nominał), s_nom (score nominalnego XGB H, model foldu kandydata → OOF dla fd),
inwin (|m − m0| < signal_window), meta (json).
Ślad (low-level, 11): log pT, η, Δφ(ślad, D), Δη(ślad, D), log p, dE/dx (brak → NaN), brak dE/dx (0/1), N pkt dE/dx, d0, z0, rola (K = 1, π = 0).
Ładunek: kandydaci są budowani z poprawną kombinacją ładunków (K∓π±(π±)), więc ładunek śladu względem D = rola (K/π) — bez osobnej kolumny.
φ: absolutny azymut jest symetrią detektora → podajemy Δφ względem kierunku D.
Użycie: pz_ext.py --cfg config_used.yaml --h5-dir DIR --models DIR --species D0|Dplus --sample fd|mb --cap N --out OUT.npz [--max-rows N]"""
import argparse
import json
import os
import sys
import time

import numpy as np

T0 = time.time()
log = lambda *a: print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)  # noqa: E731
NTF = 11
TRK_NAMES = ["log_pt", "eta", "dphi_D", "deta_D", "log_p", "dedx", "dedx_missing", "npts", "d0", "z0", "is_kaon"]
BKG_TRUTH = (8, 99)


def track_tensor(L, sp):
    """[n, nd, 11] w kolejności: K, potem piony malejąco po pT."""
    nd = int(sp["n_daughters"]); k = int(sp.get("kaon_index", 0)); pis = list(sp.get("pion_indices", [1]))
    g = lambda c, i: L[f"{c}{i}"].to_numpy(dtype=float)  # noqa: E731
    n = len(L)
    phiD, etaD = L["phi"].to_numpy(dtype=float), L["eta"].to_numpy(dtype=float)
    pt_pi = np.stack([g("d_pt", j) for j in pis], 1)
    order = np.argsort(-np.nan_to_num(pt_pi, nan=-1), axis=1)                    # piony malejąco po pT
    slots = np.empty((n, nd), np.int64); slots[:, 0] = k
    slots[:, 1:] = np.asarray(pis)[order]
    T = np.full((n, nd, NTF), np.nan, np.float32)
    for s in range(nd):
        idx = slots[:, s]
        col = lambda c: np.choose(idx, [g(c, i) for i in range(nd)])  # noqa: E731
        pt, eta, phi = col("d_pt"), col("d_eta"), col("d_phi")
        p = col("d_p"); p = np.where(np.isfinite(p) & (p > 0), p, pt * np.cosh(eta))
        dedx = col("d_dedx"); miss = ~(dedx > 0)
        T[:, s, 0] = np.log(np.clip(pt, 1e-3, None)); T[:, s, 1] = eta
        T[:, s, 2] = (phi - phiD + np.pi) % (2 * np.pi) - np.pi; T[:, s, 3] = eta - etaD
        T[:, s, 4] = np.log(np.clip(p, 1e-3, None)); T[:, s, 5] = np.where(miss, np.nan, dedx); T[:, s, 6] = miss
        T[:, s, 7] = col("d_npts"); T[:, s, 8] = col("d_d0"); T[:, s, 9] = col("d_z0"); T[:, s, 10] = float(s == 0)
    return T


def nominal_score(Xh, fold, models_dir):
    import xgboost as xgb
    s = np.full(len(Xh), np.nan, np.float32)
    if not models_dir or not os.path.isdir(models_dir):
        return s
    for k in np.unique(fold):
        f = os.path.join(models_dir, f"fold{int(k)}.json")
        if not os.path.exists(f):
            continue
        b = xgb.Booster(); b.load_model(f)
        m = fold == k
        s[m] = b.predict(xgb.DMatrix(Xh[m])).astype(np.float32)
    return s


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfg", required=True); ap.add_argument("--h5-dir"); ap.add_argument("--models", default="")
    ap.add_argument("--species", required=True, choices=["D0", "Dplus"]); ap.add_argument("--sample", required=True, choices=["fd", "mb"])
    ap.add_argument("--cap", type=int, default=2_000_000); ap.add_argument("--out", required=True); ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--workdir", default=os.environ.get("TMPDIR", "/tmp"))
    a = ap.parse_args(argv)
    from charmtag import config as C, samples as S
    import yaml
    with open(a.cfg) as f:
        site = yaml.safe_load(f)
    ov = {"paths": {"out_dir": a.workdir}, "io": {"dedup": False}}
    if a.h5_dir:
        ov["paths"]["h5_dir"] = a.h5_dir
    if a.max_rows:
        ov["io"]["max_rows"] = a.max_rows
    cfg = C.load(None, overrides=S_merge(site, ov))
    sp = cfg["species"][a.species]
    key = "signal_mc" if a.sample == "fd" else "incl_mc"
    feats = list(cfg["features"]["frozen7"]) + list(cfg["features"]["extra_H"])
    meta_nom = os.path.join(a.models, "meta.json") if a.models else ""
    if meta_nom and os.path.exists(meta_nom):
        feats = json.load(open(meta_nom))["features"]
    n_rng, _ = S.count_species(cfg, key, a.species)
    keep = min(1.0, 1.2 * a.cap / max(0.9 * n_rng, 1)) if a.sample == "mb" else 1.0   # mb: ~wszystko kombinatoryka; fd: cały sygnał, potem --cap
    log(f"{key}/{a.species}: {n_rng} kandydatów w zakresie; keep_frac = {keep:.4f}")
    rng = np.random.default_rng(int(cfg["io"]["seed"]) + (7 if a.sample == "mb" else 3) + (0 if a.species == "D0" else 11))
    tv = cfg["signal"]["truth_values"]
    out = {k: [] for k in ("Xh", "T", "G", "m", "pt", "abs_eta", "weight", "run", "event", "truth")}
    nsel = 0
    for L, _ in S.iter_species(cfg, key, a.species):
        t = L["truth"].to_numpy(dtype=float)
        sel = S.in_range(L, sp).to_numpy() & S.dstar_veto_mask(L, sp)
        sel &= np.isin(t, tv) if a.sample == "fd" else ((t >= BKG_TRUTH[0]) & (t <= BKG_TRUTH[1]))
        if keep < 1:
            sel &= rng.random(len(L)) < keep
        if not sel.any():
            continue
        L = L[sel].reset_index(drop=True)
        out["Xh"].append(L[feats].to_numpy(dtype=np.float32)); out["T"].append(track_tensor(L, sp))
        out["G"].append(L[["Lxy"]].to_numpy(dtype=np.float32))
        for c in ("m", "pt", "abs_eta", "weight", "run", "event", "truth"):
            out[c].append(L[c].to_numpy())
        nsel += len(L)
        log(f"  +{len(L)} (razem {nsel})")
    O = {k: np.concatenate(v) if v else np.zeros(0) for k, v in out.items()}
    if len(O["m"]) > a.cap:                                                       # deterministyczne przycięcie do --cap
        idx = np.sort(np.random.default_rng(5).choice(len(O["m"]), a.cap, replace=False))
        O = {k: v[idx] for k, v in O.items()}
    O["run"] = np.nan_to_num(O["run"].astype(float)).astype(np.int64); O["event"] = O["event"].astype(np.int64)
    O["weight"] = np.nan_to_num(O["weight"].astype(np.float64), nan=1.0)
    O["fold"] = S.fold_of(O["run"], O["event"], int(cfg["training"]["n_folds"]))
    O["inwin"] = np.abs(O["m"] - sp["m0"]) < sp["signal_window"]
    O["s_nom"] = nominal_score(O["Xh"], O["fold"], a.models)
    meta = {"species": a.species, "sample": a.sample, "key": key, "n": int(len(O["m"])), "n_inwin": int(O["inwin"].sum()), "keep_frac": keep,
            "features_hi": feats, "track_features": TRK_NAMES, "m0": sp["m0"], "signal_window": sp["signal_window"], "mass_range": sp["mass_range"],
            "n_folds": int(cfg["training"]["n_folds"]), "reweight_bins": cfg["training"]["reweight_bins"], "nominal_models": a.models,
            "s_nom_available": bool(np.isfinite(O["s_nom"]).any())}
    np.savez_compressed(a.out, meta=json.dumps(meta), **{k: v for k, v in O.items()})
    log("OK", a.out, meta["n"], "w oknie", meta["n_inwin"])


def S_merge(site, ov):
    from charmtag.config import _merge
    return _merge(site, ov)


if __name__ == "__main__":
    sys.exit(main())
