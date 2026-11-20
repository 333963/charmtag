"""Budowa próbek: kandydaci danego gatunku z pliku, sygnał MC (truth), tło = sidebandy danych, weto D*, wagi."""
from __future__ import annotations

import hashlib
import logging
import os

import numpy as np
import pandas as pd

from . import config as C
from . import features as F
from .io import list_columns, read_rows, resolve_columns

log = logging.getLogger("charmtag")

from .constants import M_K, M_PI


def colmap(cfg, path, nd=3):
    avail = list_columns(path, cfg["io"]["dataset"])
    cm = resolve_columns(avail, cfg["columns"], nd)
    opt = cfg.get("columns_opt") or {}
    if opt:  # HYP: kolumny opcjonalne (v1.2: hypo_swap, trackId) — dopisywane tylko znalezione, bez wpisu do `missing`
        co = resolve_columns(avail, opt, nd)
        cm.mapping.update({k: v for k, v in co.mapping.items() if v and k not in cm.mapping})
    return cm


def _spec(cfg, key, species, decay_type=None, nd=None):
    sp = cfg["species"][species]
    dt = int(decay_type if decay_type is not None else sp["decay_type"])
    nd = int(nd or sp["n_daughters"])
    path = C.sample_path(cfg, key)
    cm = colmap(cfg, path, nd)
    for need in ("decay_type", "m", "pt", "eta", "event"):
        if not cm.has(need):
            raise KeyError(f"{path}: brak kolumny dla '{need}' (kandydaci: {cfg['columns'][need]}) — popraw columns w YAML")
    return sp, dt, nd, path, cm


def iter_species(cfg: dict, key: str, species: str, decay_type: int | None = None, nd: int | None = None):
    """Generator ramek logicznych (z cechami) porcjami — dla dużych plików (apply)."""
    sp, dt, nd, path, cm = _spec(cfg, key, species, decay_type, nd)
    dcol = cm["decay_type"]
    cols = F.raw_columns(cm, nd)
    chunk = int(cfg["io"]["chunk_rows"])
    from .io import n_rows
    src = C.source(cfg, key)  # etap P: 1 plik albo lista (row = wiersz globalny w sklejonej liście)
    n = n_rows(src, cfg["io"]["dataset"])
    if cfg["io"].get("max_rows"):
        n = min(n, int(cfg["io"]["max_rows"]))
    for start in range(0, n, chunk):
        raw = _rr_slice(src, cols, cfg["io"]["dataset"], start, min(n, start + chunk))
        hyp_blk = None
        if nd == 2 and cfg["signal"].get("swap_resolution") == "arbiter" and cm.has("event"):  # HYP: bloki zdarzeń (pary hipotez)
            from .hyp import event_blocks
            hyp_blk, _ = event_blocks(raw[cm["run"]].to_numpy() if cm.has("run") else np.zeros(len(raw), np.int64),
                                      raw[cm["event"]].to_numpy(), start)
        sel_dt = raw[dcol].to_numpy() == dt
        raw = raw[sel_dt]
        L = F.to_logical(raw.reset_index(drop=True), cm, nd)
        if not cm.has("run"):
            L["run"] = 0
        L["weight"] = L["weight"].fillna(1.0)
        L = F.compute(L, {**sp, "n_daughters": nd}, cfg.get("units"))
        L["row"] = raw.index.to_numpy() + start
        if hyp_blk is not None:
            L["evblk"] = hyp_blk[sel_dt]
        if nd == 2 and dt == sp["decay_type"]:
            mark_swaps(L, sp, cfg)
            if cfg["signal"].get("swap_resolution") == "arbiter":  # HYP: P(poprawna hipoteza) + wybór wiersza w parze
                from .hyp import apply_arbiter
                L = apply_arbiter(L, cfg)
        yield L, cm


def _rr_slice(path, cols, dataset, start, stop):
    import h5py
    cols = list(dict.fromkeys(cols))
    if isinstance(path, (list, tuple)):  # etap P: wiele plików
        from .io import open_events
        with open_events(path, dataset) as ds:
            arr = ds.fields(cols)[start:stop]
            return pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
    with h5py.File(path, "r") as f:
        if dataset in f and isinstance(f[dataset], h5py.Dataset) and f[dataset].dtype.names:
            arr = f[dataset].fields(cols)[start:stop]
            return pd.DataFrame({c: np.asarray(arr[c]) for c in cols})
        found = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and obj.ndim == 1:
                found.setdefault(name.split("/")[-1], name)

        f.visititems(visit)
        return pd.DataFrame({c: np.asarray(f[found[c]][start:stop]) for c in cols})


def _shrink(L):
    for c in L.columns:
        if L[c].dtype == np.float64 and c not in ("m", "dm"):
            L[c] = L[c].astype(np.float32)
    return L


def load_species(cfg: dict, key: str, species: str, cache: bool = True, decay_type: int | None = None,
                 nd: int | None = None, select=None, keep_frac: float = 1.0, tagx: str = "") -> pd.DataFrame:
    """Ramka logiczna + cechy; select(L)->maska, keep_frac = losowa podpróbka (deterministyczna)."""
    sp, dt, nd, path, cm = _spec(cfg, key, species, decay_type, nd)
    cdir = os.path.join(cfg["paths"]["out_dir"], "cache")
    src = C.source(cfg, key)
    sig = (path, os.path.getmtime(path)) if isinstance(src, str) else tuple((q, os.path.getmtime(q)) for q in src)
    _hx = ()
    if cfg["signal"].get("swap_resolution") == "arbiter":  # HYP: inny cache (hypo_swap/trackId/p_hyp); stare tagi bez zmian
        _hx = ("hyp", sorted((k, str(v)) for k, v in (cfg.get("columns_opt") or {}).items()),
               str((cfg["signal"].get("arbiter") or {}).get("model_dir")))
    tag = hashlib.md5(repr((*sig, dt, nd, cfg["io"].get("max_rows"), keep_frac, tagx,
                            sorted((k, str(v)) for k, v in cfg["columns"].items()), *_hx)).encode()).hexdigest()[:10]
    cfile = os.path.join(cdir, f"{key}_{dt}_{tag}.pkl")
    if cache and os.path.exists(cfile):
        return pd.read_pickle(cfile)
    rng = np.random.default_rng(cfg["io"]["seed"] + dt)
    parts = []
    for L, _ in iter_species(cfg, key, species, decay_type, nd):
        m = np.ones(len(L), bool) if select is None else np.asarray(select(L), bool)
        if keep_frac < 1:
            m &= rng.random(len(L)) < keep_frac
        parts.append(_shrink(L[m].reset_index(drop=True)))
    L = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    L["sample"] = key
    L.attrs["missing_columns"] = cm.missing
    _dd = cfg.get("io", {}).get("dedup", False)  # etap DD: dokładne duplikaty kandydatów (false | "report" | true; domyślnie false)
    if _dd:
        from .dedup import apply as _dedup_apply
        L = _dedup_apply(L, _dd, f"{key}/{species}")
    if cache:
        os.makedirs(cdir, exist_ok=True)
        L.to_pickle(cfile)
    log.info("%s: %d kandydatów decayType=%d (brak kolumn: %s)", key, len(L), dt, cm.missing)
    return L


def count_species(cfg, key, species, select_raw=None):
    """Szybkie liczenie kandydatów danego decayType (2 kolumny)."""
    sp, dt, nd, path, cm = _spec(cfg, key, species)
    df = read_rows(C.source(cfg, key), [cm["decay_type"], cm["m"]], cfg["io"]["dataset"], cfg["io"]["chunk_rows"],
                   row_filter=lambda d: d[cm["decay_type"]].to_numpy() == dt, max_rows=cfg["io"].get("max_rows"))
    m = df[cm["m"]].to_numpy()
    lo, hi = sp["mass_range"]
    return int(((m >= lo) & (m < hi)).sum()), int((np.abs(m - sp["m0"]) > sp["sideband_min"]).sum())


def _m2(pts, etas, phis, masses):
    px = py = pz = e = 0.0
    for pt, eta, phi, m in zip(pts, etas, phis, masses):
        px = px + pt * np.cos(phi)
        py = py + pt * np.sin(phi)
        pz = pz + pt * np.sinh(eta)
        e = e + np.sqrt((pt * np.cosh(eta)) ** 2 + m * m)
    return np.sqrt(np.maximum(e * e - px * px - py * py - pz * pz, 0))


def dstar_veto_mask(L: pd.DataFrame, sp: dict) -> np.ndarray:
    """True = przechodzi weto. D+→Kππ: odrzuć gdy min(m − m_Kπ_j) < dm_max (m_Kπ z kolumn mKpi1/2,
    a gdy ich brak — z kinematyki córek)."""
    v = sp.get("dstar_veto") or {}
    if not v.get("apply") or sp["n_daughters"] < 3:
        return np.ones(len(L), bool)
    m = L["m"].to_numpy(dtype=float)
    if "m_kpi1" in L and np.isfinite(L["m_kpi1"]).any():
        dm = np.fmin(m - L["m_kpi1"].to_numpy(dtype=float), m - L["m_kpi2"].to_numpy(dtype=float))
        return ~(dm < v["dm_max"])
    k = sp["kaon_index"]
    pis = sp["pion_indices"]
    g = lambda q, i: L[f"{q}{i}"].to_numpy(dtype=float)  # noqa: E731
    if not np.isfinite(g("d_pt", k)).any() or not np.isfinite(g("d_phi", k)).any():
        log.warning("weto D*: brak kinematyki córek → pominięte")
        return np.ones(len(L), bool)
    keep = np.ones(len(L), bool)
    for j in pis:
        m2 = _m2([g("d_pt", k), g("d_pt", j)], [g("d_eta", k), g("d_eta", j)],
                 [g("d_phi", k), g("d_phi", j)], [M_K, M_PI])
        keep &= ~((m - m2) < v["dm_max"])
    return keep


def mark_swaps(L: pd.DataFrame, sp: dict, cfg: dict):
    """D0: para torów z obiema hipotezami truth-matched → dalsza od m0 = swap (truth += 100).
    HYP: `swap_resolution: arbiter` → prawda z `Dmeson_hypo_swap` (v1.2, `columns_opt.hypo_swap`): truth += 100 dla hipotezy
    zamienionej (także w wierszach bez bliźniaka); bez kolumny → closest_mass jak dotąd (ostrzeżenie)."""
    mode = cfg["signal"].get("swap_resolution", "closest_mass")
    if mode == "arbiter":
        from .hyp import truth_from_hypo_swap
        if truth_from_hypo_swap(L, cfg):
            return
        if np.isin(L["truth"].to_numpy(dtype=float), cfg["signal"]["truth_values"]).any():
            log.warning("swap_resolution=arbiter: brak kolumny hypo_swap (v1.1?) → prawda swapu z closest_mass")
        mode = "closest_mass"
    if mode != "closest_mass":
        return
    t = L["truth"].to_numpy(dtype=float)
    tv = np.isin(t, cfg["signal"]["truth_values"])
    if not tv.any():
        return
    idx = np.flatnonzero(tv)
    a = np.round(L["d_pt0"].to_numpy(dtype=float)[idx], 4)
    b = np.round(L["d_pt1"].to_numpy(dtype=float)[idx], 4)
    key = pd.DataFrame({"run": L["run"].to_numpy()[idx], "ev": L["event"].to_numpy()[idx],
                        "lo": np.minimum(a, b), "hi": np.maximum(a, b),
                        "dm": np.abs(L["m"].to_numpy(dtype=float)[idx] - sp["m0"]), "i": idx})
    g = key.groupby(["run", "ev", "lo", "hi"])["dm"]
    rank = g.rank(method="first")
    n = g.transform("size")
    sw = key["i"].to_numpy()[((n > 1) & (rank > 1)).to_numpy()]
    L.loc[sw, "truth"] = L.loc[sw, "truth"] + 100


def in_range(L, sp):
    lo, hi = sp["mass_range"]
    return (L["m"] >= lo) & (L["m"] < hi) & (L["pt"] >= 3) & (L["pt"] < 30)


def signal_mask(L, cfg):
    return L["truth"].isin(cfg["signal"]["truth_values"]).to_numpy()


def sideband_mask(L, sp):
    return (np.abs(L["m"] - sp["m0"]) > sp["sideband_min"]).to_numpy()


def fold_of(run, event, n_folds: int) -> np.ndarray:
    """Deterministyczny fold z (run, event_number): splitmix64."""
    x = (np.asarray(run, dtype=np.uint64) * np.uint64(0x9E3779B97F4A7C15)) ^ np.asarray(event, dtype=np.int64).astype(np.uint64)
    with np.errstate(over="ignore"):
        x = x + np.uint64(0x9E3779B97F4A7C15)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        x = x ^ (x >> np.uint64(31))
    return (x % np.uint64(n_folds)).astype(np.int8)


def reweight(bkg: pd.DataFrame, sig: pd.DataFrame, bins: dict) -> np.ndarray:
    """Wagi tła tak, by rozkład (pT,|η|) zgadzał się z sygnałem; normalizacja sum_w(B)=sum_w(S)."""
    bp, be = np.asarray(bins["pt"]), np.asarray(bins["abs_eta"])
    hs, _, _ = np.histogram2d(sig["pt"], sig["abs_eta"], [bp, be], weights=sig["weight"])
    hb, _, _ = np.histogram2d(bkg["pt"], bkg["abs_eta"], [bp, be], weights=bkg["weight"])
    hs /= max(hs.sum(), 1e-12)
    hb /= max(hb.sum(), 1e-12)
    ratio = np.divide(hs, hb, out=np.zeros_like(hs), where=hb > 0)
    ip = np.clip(np.digitize(bkg["pt"], bp) - 1, 0, len(bp) - 2)
    ie = np.clip(np.digitize(bkg["abs_eta"], be) - 1, 0, len(be) - 2)
    w = bkg["weight"].to_numpy() * ratio[ip, ie]
    return w * (sig["weight"].sum() / max(w.sum(), 1e-12))


def training_set(cfg: dict, species: str):
    """(sig, bkg) — sygnał MC truth w zakresie, tło = sidebandy danych (none) po wetach."""
    sp = cfg["species"][species]
    rng = np.random.default_rng(cfg["io"]["seed"])
    S = load_species(cfg, "signal_mc", species, tagx="sig",
                     select=lambda L: in_range(L, sp).to_numpy() & signal_mask(L, cfg) & dstar_veto_mask(L, sp))
    # etap S5: training.bkg_source (domyślnie sideband = bez zmian). Źródła inne niż sideband są w pakiecie charmbkg (projekt B);
    # domyślna ścieżka go nie importuje, więc `charmtag train` i testy projektu A działają bez charmbkg (CI, 20.11.2026).
    if ((cfg.get("training") or {}).get("bkg_source") or "sideband") != "sideband":
        from .bkgsrc import source as _bkg_source, training_background
        _bkg_source(cfg)   # walidacja nazwy źródła (ValueError dla nieznanej) — jak dotąd
        B = training_background(cfg, species)
    else:
        _, n_sb = count_species(cfg, "data_none", species)
        cap = cfg["io"]["max_bkg_train"]
        frac = min(1.0, 1.3 * cap / max(n_sb, 1)) if cap else 1.0
        B = load_species(cfg, "data_none", species, keep_frac=frac, tagx="sb",
                         select=lambda L: in_range(L, sp).to_numpy() & sideband_mask(L, sp) & dstar_veto_mask(L, sp))
    for name, X, cap in (("sig", S, cfg["io"]["max_sig_train"]), ("bkg", B, cfg["io"]["max_bkg_train"])):
        if cap and len(X) > cap:
            idx = np.sort(rng.choice(len(X), int(cap), replace=False))
            if name == "sig":
                S = X.iloc[idx].reset_index(drop=True)
            else:
                B = X.iloc[idx].reset_index(drop=True)
    if len(S) == 0 or len(B) == 0:
        raise RuntimeError(f"pusta próbka treningowa: sygnał={len(S)} tło={len(B)}")
    B["w_train"] = reweight(B, S, cfg["training"]["reweight_bins"])
    S["w_train"] = S["weight"].to_numpy(dtype=float)
    S["y"], B["y"] = 1, 0
    return S, B
