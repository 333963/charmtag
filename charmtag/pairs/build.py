"""`charmtag pairs`: nogi i pary z wyników charmtag (scores/<sp>_<var>/<próbka>.h5) + śladów z h5 v1.1.

Reguły jak M7 v3 (unblinding_freeze_v20260927_01/code/m7_inputs_v3.py, unblind_lib.py):
- kategorie = przedziały wydajności sygnału MC: [0,30%), [30,45%), [45,65%) (próg score z MC truth, OOF),
  D0 dodatkowo z tagiem D* (|Δm − 145,43| < 1 MeV, ta sama hipoteza K/π co w D*): tag & ≥WP65 → 3, tag & WP90–65 → 4;
- D0: jedna hipoteza K/π na parę śladów, wybrana hashem klucza (ślepo na masę) spośród kandydatów ≥ WP90;
- D⁺: weto D* (veto_ok z charmtag), score ≥ WP65;
- pary ze wspólnym śladem usuwane (ślad = bity (d0, z0·sinθ_PV) córki — parametry toru niezależne od refitu kandydata; PF);
- ślepota: centrum (obie nogi w oknie |m−m0| < 75/70 MeV) wycinane przed histogramowaniem, chyba że --unblind.
Normalizacja (scenariusz R, jak 09-27): N_par = N_D*(dane, odporny fit Δm `dstar_refit.fit_single`, kształt z MC) × [par/D*]_MC minbias ND × c_w,
N_prawdziwych D = N_D* × [D/D*]_MC × c_w, gdzie c_w = (ważone/nieważone) z wag krotności 2D (n_trk, pT wiodący) z
mb_rates_2d.json (h5 v1.1 nie ma n_trk ani pT wiodącego śladu → wagi nie do odtworzenia; przenosimy ich efekt jako czynnik).
Lumi: dane = część 428770 w h5 (joby 25–47 z 47); N_D* liczone na tej samej próbce → lumi skraca się w μ.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os

import numpy as np
import pandas as pd

from .. import config as C
from .. import constants as _C
from .. import store
from ..io import list_columns, resolve_columns
from . import CHANNELS
from .model import Inputs, M0, SBW
from .frozen import sbsb_counts, shared_track
from .hists import build_pair_hists

log = logging.getLogger("charmtag")

DEFAULTS = {
    "variant": "H",
    "cats": [[0.0, 0.30], [0.30, 0.45], [0.45, 0.65]],
    "d0_arbitration_wp": 0.90,
    "tag_dm": [_C.DM_DSTAR, 0.0010],          # GeV: środek, pół-szerokość
    "tag_sf": 0.70,                       # SF wydajności tagu D* dane/MC (09-26: 0,70; nadpisać wynikiem etapu N)
    "sf_json": None,                      # skumulowane SF(WP) per gatunek (format sf_for_pairs_v2.json); None → 1
    "edges": {"10": [1750.0, 2000.0, 25], "4": [1700.0, 2200.0, 25]},
    "dstar_norm": {"d0_window": 0.040, "fit_range": [0.1396, 0.165]},
    "mb_rates_json": None,                # ścieżka do mb_rates_2d.json (zamrożone wejścia 09-27, niepubliczne) — ustaw w configu `pairs:`
    "freeze_dir": None,                   # katalog zamrożonych wejść (opcja --closure) — ustaw w configu `pairs:`
    "lumi_fraction": 23 / 47,             # joby 25–47 z 47 (data_428770_POLOWA.txt)
    "R_reference": 1.70,
    "chunk_rows": 2_000_000,
    "track_id": "d0z0",                   # PF 2026-11-04: ID śladu z bitów (d0, z0sinθ_PV), niezależne od refitu; "pt_eta" = stare (weto wspólnych śladów nieskuteczne)
    "d0_arbitration": "hash",             # HYP 2026-11-08: hash (ślepo, jak M7) | arbiter (wyższe p_hyp z charmtag.hyp; scores z apply+arbiter)
    "diff": None,                         # PD 2026-11-05: pary różniczkowe (sign | dphi | dy | sign_dphi), domyślnie WYŁĄCZONE — patrz pairs/diff.py
}


def pcfg(cfg):
    p = dict(DEFAULTS)
    p.update(cfg.get("pairs", {}) or {})
    return p


def pdir(cfg):
    return os.path.join(cfg["paths"]["out_dir"], "pairs")


def edges(p):
    return {int(k): np.linspace(v[0], v[1], int(v[2]) + 1) for k, v in p["edges"].items()}


# ----------------------------------------------------------------------------- ślady z h5
def _trk_id(pt, eta):
    """Stary identyfikator (bity pT, η córki). UWAGA (PF, 2026-11-04): pędy córek w h5 są PO refit wierzchołka kandydata,
    więc ten sam ślad w dwóch kandydatach ma inne bity → weto wspólnych śladów nie działało (0 % w D⁺, tylko swapy w D0)."""
    pt = np.asarray(pt, np.float32)
    eta = np.asarray(eta, np.float32)
    tid = (pt.view(np.uint32).astype(np.int64) << 31) | (eta.view(np.uint32).astype(np.int64) >> 1)   # 63 bity, ≥ 0
    return np.where(pt > 0, tid, -1)


def _trk_id_d0z0(pt, d0, z0):
    """Identyfikator śladu z parametrów toru względem PV (d0, z0·sinθ) — niezależnych od refitu kandydata (PF):
    na danych v1.1 100 % par z tolerancją 1 % w (pT,η,φ) ma identyczne bity (d0, z0), a 3,7 % par dzieli ślad przy różnicy pT > 1 %."""
    pt = np.asarray(pt, np.float32)
    d0 = np.asarray(d0, np.float32)
    z0 = np.asarray(z0, np.float32)
    # PF3-FIN (05.11): `<< 32` na int64 dawało UJEMNE ID dla d0 < 0 (bit znaku float32) → shared_track/arbitraż (t >= 0)
    # pomijały ≈ połowę torów (weto działało w połowie: D⁺D⁺ 68 % zamiast 85 %). Teraz 63 bity, zawsze ≥ 0: pełne 32 bity d0
    # (ze znakiem) << 31 | z0 bez najmłodszego bitu mantysy (kolizja wymaga bitowo identycznego d0 i z0 różnego o 1 ulp).
    tid = (d0.view(np.uint32).astype(np.int64) << 31) | (z0.view(np.uint32).astype(np.int64) >> 1)
    return np.where(pt > 0, tid, -1)


def track_ids(cfg, g):
    """Lista 3 tablic ID śladów córek wg `pairs.track_id`: d0z0 (domyślnie, poprawne) | pt_eta (stare, do odtworzenia starych liczb)."""
    mode = pcfg(cfg)["track_id"]
    if mode == "pt_eta":
        log.warning("pairs: track_id=pt_eta — ID śladów z (pT, η) po refit: weto wspólnych śladów NIESKUTECZNE (tylko do odtworzenia starych liczb)")
        return [_trk_id(g[f"d_pt{i}"], g[f"d_eta{i}"]) for i in range(3)]
    if mode != "d0z0":
        raise SystemExit(f"pairs: nieznany pairs.track_id={mode!r} (d0z0 | pt_eta)")
    missing = [c for i in range(3) for c in (f"d_d0{i}", f"d_z0{i}") if c not in g]
    if missing:
        raise SystemExit("pairs: brak kolumn parametrów torów córek %s w źródle h5 — weto wspólnych śladów wymaga d0/z0 (PF). "
                         "Sprawdź `columns.d_d0/d_z0` w configu; cichy powrót do (pT, η) jest wyłączony "
                         "(jawnie: pairs.track_id: pt_eta)." % missing)
    return [_trk_id_d0z0(g[f"d_pt{i}"], g[f"d_d0{i}"], g[f"d_z0{i}"]) for i in range(3)]


def scan_source(cfg, key, rows: np.ndarray, chunk=2_000_000, dstar=True, count_events=True, kin=False, charge_src=None):
    """Jeden przebieg po pliku źródłowym: ślady dla `rows`, kandydaci D* (Δm < 0.165), liczba zdarzeń z kandydatami.
    PD: `kin=True` dokłada do nóg η, φ, ładunek kandydata i `dq0` = ładunek córki 0 (K; kolumna v1.2 `Dmeson_daughter0_charge`,
    brak → NaN); `charge_src` = lista plików h5 v1.2 z tymi samymi wierszami (60 kolumn bit w bit) → `dq0` czytane stamtąd
    (join po `row`, kontrola run/event/decayType/masa w wybranych wierszach — niezgodność = błąd)."""
    path = C.sample_path(cfg, key)
    ds_name = cfg["io"]["dataset"]
    cm = resolve_columns(list_columns(path, ds_name), cfg["columns"], 3)
    need = (["run", "event", "decay_type", "m", "dm", "truth", "pt"] + [f"d_pt{i}" for i in range(3)] + [f"d_eta{i}" for i in range(3)]
            + [f"d_d0{i}" for i in range(3)] + [f"d_z0{i}" for i in range(3)])  # PF: d0/z0 → ID śladu niezależne od refitu
    if kin:
        need += ["eta", "phi", "charge"]
    cols = {k: cm.get(k) for k in need if cm.get(k)}
    if kin:
        from .diff import D_CHARGE_CANDS
        dq = resolve_columns(list_columns(path, ds_name), {"d_charge": cfg["columns"].get("d_charge", D_CHARGE_CANDS)}, 1)
        if dq.get("d_charge0"):
            cols["dq0"] = dq["d_charge0"]
        miss = [k for k in ("eta", "phi", "charge") if k not in cols]
        if miss:
            raise SystemExit(f"pairs diff: brak kolumn kinematyki {miss} w {path} (cfg.columns)")
    rows = np.unique(np.asarray(rows, np.int64))
    trk, dst, evs = [], [], []
    from ..io import open_events  # etap P: 1 plik albo lista (row = wiersz globalny)
    cs = _ChargeSource(charge_src, ds_name, cfg) if (kin and charge_src) else None
    with open_events(C.sample_paths(cfg, key) if C.is_multi(cfg, key) else [path], ds_name) as ds:
        n = ds.shape[0]
        if cfg["io"].get("max_rows"):
            n = min(n, int(cfg["io"]["max_rows"]))
        if cs is not None and cs.ds.shape[0] != ds.shape[0]:
            raise SystemExit(f"pairs diff charge_source[{key}]: {cs.ds.shape[0]} wierszy vs {ds.shape[0]} w źródle — to nie te same pliki/kolejność")
        for start in range(0, n, chunk):
            stop = min(n, start + chunk)
            a = ds.fields(list(cols.values()))[start:stop]
            g = {k: np.asarray(a[v]) for k, v in cols.items()}
            t = track_ids(cfg, g)
            run = g["run"].astype(np.int64) if "run" in g else np.zeros(stop - start, np.int64)
            sel = rows[(rows >= start) & (rows < stop)] - start
            if len(sel):
                T1 = {"row": sel + start, "trk0": t[0][sel], "trk1": t[1][sel], "trk2": t[2][sel]}
                if kin:
                    T1.update({"eta": g["eta"][sel].astype(float), "phi": g["phi"][sel].astype(float), "charge": g["charge"][sel].astype(float),
                               "dq0": g["dq0"][sel].astype(float) if "dq0" in g else np.full(len(sel), np.nan)})
                    if cs is not None:
                        T1["dq0"] = cs.dq0(start, stop, sel, g, key)
                trk.append(pd.DataFrame(T1))
            if dstar:
                dsm = (g["decay_type"] == cfg["dstar_control"]["dstar_decay_type"]) & (g["dm"] > 0) & (g["dm"] < 0.165)
                if dsm.any():
                    dst.append(pd.DataFrame({"run": run[dsm], "event": g["event"][dsm].astype(np.int64), "m": g["m"][dsm].astype(float),
                                             "dm": g["dm"][dsm].astype(float), "pt": g["pt"][dsm].astype(float),
                                             "truth": g["truth"][dsm].astype(float) if "truth" in g else np.nan,
                                             "trk0": t[0][dsm], "trk1": t[1][dsm], "trk2": t[2][dsm]}))
            if count_events:
                evs.append(np.unique((run << 40) ^ g["event"].astype(np.int64)))
    T = pd.concat(trk, ignore_index=True) if trk else pd.DataFrame(columns=["row", "trk0", "trk1", "trk2"] + (["eta", "phi", "charge", "dq0"] if kin else []))
    D = pd.concat(dst, ignore_index=True) if dst else pd.DataFrame(columns=["run", "event", "m", "dm", "pt", "truth", "trk0", "trk1", "trk2"])
    nev = int(len(np.unique(np.concatenate(evs)))) if evs else 0
    if cs is not None:
        cs.close()
        log.info("pairs diff: ładunek córki 0 z charge_source[%s] (%d plików, %d wierszy sprawdzonych, 0 niezgodnych)", key, len(cs.paths), cs.checked)
    return T, D, nev


class _ChargeSource:
    """PD: h5 v1.2 jako źródło `Dmeson_daughter0_charge` dla nóg z v1.1 (te same wiersze). Kontrola zgodności w wybranych wierszach."""
    CHECK = ("run", "event", "decay_type", "m")

    def __init__(self, paths, ds_name, cfg):
        from ..io import open_events
        from .diff import D_CHARGE_CANDS
        self.paths = list(paths)
        self.ds = open_events(self.paths, ds_name)
        avail = list_columns(self.paths, ds_name)
        self.cm = resolve_columns(avail, {k: cfg["columns"][k] for k in self.CHECK}, 3)
        dq = resolve_columns(avail, {"d_charge": cfg["columns"].get("d_charge", D_CHARGE_CANDS)}, 1)
        if not dq.get("d_charge0"):
            raise SystemExit(f"pairs diff charge_source: brak kolumny ładunku córki 0 ({D_CHARGE_CANDS}) w {self.paths[0]} — to nie v1.2?")
        self.col = dq["d_charge0"]
        self.checked = 0

    def dq0(self, start, stop, sel, g, key):
        cols = [self.cm[k] for k in self.CHECK if self.cm.get(k)] + [self.col]
        a = self.ds.fields(cols)[start:stop]
        for k in self.CHECK:
            if self.cm.get(k) and k in g:
                x, y = np.asarray(a[self.cm[k]])[sel], g[k][sel]
                bad = ~(x == y) if not np.issubdtype(x.dtype, np.floating) else ~(np.isclose(x.astype(float), y.astype(float), rtol=0, atol=1e-9) | (np.isnan(x.astype(float)) & np.isnan(y.astype(float))))
                if bad.any():
                    raise SystemExit(f"pairs diff charge_source[{key}]: {int(bad.sum())} z {len(sel)} wierszy ma inne `{k}` niż źródło v1.1 "
                                     f"(wiersze {start}–{stop}) — join po `row` nieważny (inne pliki/kolejność?)")
        self.checked += len(sel)
        v = np.asarray(a[self.col])[sel].astype(float)
        return np.where(np.abs(v) > 1.5, np.nan, v)          # sentinel −999 → NaN

    def close(self):
        self.ds.close()


def charge_files(cfg, spec, key):
    """`pairs.diff.charge_source`: {klucz: plik|glob|lista} albo katalog z plikami o nazwach jak paths.samples → lista plików v1.2 lub None."""
    src = (spec or {}).get("charge_source")
    if not src:
        return None
    import glob as _glob
    if isinstance(src, str):
        d = os.path.expanduser(src)
        if not os.path.isdir(d):
            raise SystemExit(f"pairs.diff.charge_source={src!r}: katalog nie istnieje (albo podaj słownik {{klucz_próbki: plik/glob}})")
        h5 = os.path.abspath(cfg["paths"]["h5_dir"])
        out = []
        for q in C.sample_paths(cfg, key):          # ta sama ścieżka względem h5_dir (podkatalogi próbek!), zapas: sama nazwa pliku
            rel = os.path.relpath(os.path.abspath(q), h5)
            cand = [os.path.join(d, rel)] if not rel.startswith("..") else []
            cand.append(os.path.join(d, os.path.basename(q)))
            hit = next((c for c in cand if os.path.exists(c)), None)
            if hit is None:
                log.warning("pairs diff charge_source[%s]: brak %s w %s → ładunek córki 0 nieokreślony (v1.1)", key, rel, d)
                return None
            out.append(hit)
        return out
    it = src.get(key)
    if not it:
        return None
    items = list(it) if isinstance(it, (list, tuple)) else [it]
    out = []
    for q in items:
        q = os.path.expanduser(str(q))
        q = q if os.path.isabs(q) else os.path.join(cfg["paths"]["h5_dir"], q)
        out += [x for x in sorted(_glob.glob(q)) if not x.endswith(".part")] if any(c in q for c in "*?[") else [q]
    if not out or not all(os.path.exists(q) for q in out):
        raise SystemExit(f"pairs.diff.charge_source[{key}]: brak plików {out or items}")
    return out


# ----------------------------------------------------------------------------- nogi
def thresholds(sig, effs):
    s = np.sort(sig["score"].to_numpy())
    return {f"{e:g}": float(np.quantile(s, 1 - e)) for e in effs}


def cat_of(sc, thr, cats):
    lab = np.full(len(sc), -1)
    for i, (e0, e1) in enumerate(cats):
        hi = thr[f"{e0:g}"] if e0 > 0 else 9.0
        lo = thr[f"{e1:g}"]
        lab[(sc >= lo) & (sc < hi)] = i
    return lab


def arbitrate(df):
    """Jedna hipoteza na (zdarzenie, nieuporządkowana para śladów), wybór hashem klucza (ślepo na masę) — jak m7_templates."""
    lo = np.minimum(df.trk0, df.trk1)
    hi = np.maximum(df.trk0, df.trk1)
    key = df.source_file.astype(str) + ":" + df.entry.astype(str) + ":" + lo.astype(str) + ":" + hi.astype(str)
    df = df.assign(_key=key.values).sort_values(["_key", "cand"])
    df["_n"] = df.groupby("_key").cand.transform("size")
    df["_r"] = df.groupby("_key").cumcount()
    pick = df._key.map(lambda k: int(hashlib.sha256(k.encode()).hexdigest()[:6], 16)) % df._n
    return df[df._r == pick].drop(columns=["_key", "_n", "_r"])


def arbitrate_arbiter(df):
    """HYP: jedna hipoteza na (zdarzenie, nieuporządkowana para śladów) — wiersz z wyższym P(poprawna hipoteza) `p_hyp`
    (arbiter charmtag.hyp: dE/dx/p bez masy → bez rzeźbienia). Remis → niższy `cand`."""
    if "p_hyp" not in df:
        raise SystemExit("pairs.d0_arbitration=arbiter: brak kolumny p_hyp w scores D0 — `charmtag apply` z signal.swap_resolution=arbiter "
                         "i signal.arbiter.model_dir")
    lo = np.minimum(df.trk0, df.trk1)
    hi = np.maximum(df.trk0, df.trk1)
    key = df.source_file.astype(str) + ":" + df.entry.astype(str) + ":" + lo.astype(str) + ":" + hi.astype(str)
    df = df.assign(_key=key.values, _p=-df.p_hyp.fillna(-1.0).astype(float)).sort_values(["_key", "_p", "cand"])
    return df[df.groupby("_key").cumcount() == 0].drop(columns=["_key", "_p"])


def d0_arbitrate(L, p):
    mode = p.get("d0_arbitration", "hash")
    if mode == "hash":
        return arbitrate(L)
    if mode == "arbiter":
        return arbitrate_arbiter(L)
    raise SystemExit(f"pairs: nieznany pairs.d0_arbitration={mode!r} (hash | arbiter)")


def d0_flavour_closest_mass(L: pd.DataFrame) -> pd.DataFrame:
    """PD: ładunek córki „K" (`dq0`) do znaku zapachu D0 brany z hipotezy K/π bliższej m0 spośród rzędów TEJ SAMEJ pary torów
    (obie hipotezy = osobne wiersze h5). Swap (zła hipoteza) odwraca znak K → bez tego LS/OS w D0 rozmyte o ~30 % na nogę
    (udział swapów w szablonach 27–39 %). Zmienia TYLKO `dq0` (klasyfikacja LS/OS); arbitraż (hash, ślepy na masę), masa, kategoria
    i wszystko w nominale bez zmian. To samo w MC minbias (`legs_for_truth`) → udział f_MC per bin z tą samą resztkową dylucją."""
    if not len(L):
        return L
    lo, hi = np.minimum(L.trk0, L.trk1), np.maximum(L.trk0, L.trk1)
    key = (L.source_file.astype(str) + ":" + L.entry.astype(str) + ":" + lo.astype(str) + ":" + hi.astype(str)).to_numpy()
    dm = pd.Series(np.abs(L.m.to_numpy(float) - M0[10]), index=L.index)
    best = dm.groupby(key).idxmin()                              # indeks wiersza najbliższego m0 w grupie
    dq_best = pd.Series(L.dq0.loc[best.to_numpy()].to_numpy(), index=best.index)
    return L.assign(dq0=dq_best.reindex(key).to_numpy())


def tag_keys(D, p):
    c, w = p["tag_dm"]
    t = D[np.abs(D.dm - c) < w]
    return set(zip(t.run.astype(np.int64), t.event.astype(np.int64), t.trk0.astype(np.int64), t.trk1.astype(np.int64)))


def make_legs(sc: pd.DataFrame, T: pd.DataFrame, dt: int, thr: dict, p: dict, tags=None, E=None):
    """sc: wyniki charmtag (row, run, event, m [GeV], score, veto_ok, truth); zwraca nogi w kategoriach (m w MeV)."""
    L = sc.merge(T, on="row", how="inner")
    L = L.assign(source_file=L["run"].astype(np.int64).astype(str), entry=L["event"].astype(np.int64), cand=L["row"].astype(np.int64),
                 decay_type=dt, m=L["m"].astype(float) * 1000.0)
    if E is not None:
        L = L[(L.m >= E[dt][0]) & (L.m < E[dt][-1])]
    cats = p["cats"]
    if dt == 10:
        t90, t65 = thr[f"{p['d0_arbitration_wp']:g}"], thr["0.65"]
        L = L[L.score >= t90]
        if "dq0" in L and (p.get("diff") or {}).get("d0_flavour", "closest_mass") == "closest_mass":
            L = d0_flavour_closest_mass(L)                      # PD: znak D0 z hipotezy bliższej m0 (tylko klasyfikacja LS/OS)
        L = d0_arbitrate(L, p)                                  # HYP: hash (domyślnie) | arbiter
        tg = np.zeros(len(L), bool)
        if tags:
            tg = np.fromiter(((int(r), int(e), int(a), int(b)) in tags for r, e, a, b in
                              zip(L.run, L.event, L.trk0, L.trk1)), bool, len(L))
        s = L.score.to_numpy()
        c = cat_of(s, thr, cats)
        L = L.assign(tag=tg, cat=np.where(tg & (s >= t65), 3, np.where(tg & (s >= t90), 4, np.where(tg, -1, c))),
                     slice_=c)
    else:
        L = L[L.veto_ok.astype(bool) & (L.score >= thr["0.65"])]
        L = L.assign(tag=False, cat=cat_of(L.score.to_numpy(), thr, cats), slice_=-1)
        L["slice_"] = L["cat"]
    return L[L.cat >= 0].reset_index(drop=True)


# ----------------------------------------------------------------------------- MC: szablony i wydajności
def signal_templates(legs, n_correct_all, dt, E, p, ncat):
    """Szablon = poprawne + swap (po arbitrażu); Y = N(kat)/N(poprawnych przed cięciem). D0: SF tagu jak m7_inputs_v3."""
    sig, comp, Y = {}, {}, {}
    tv = legs.truth.to_numpy()
    corr, swap = np.isin(tv, [0, 1]), np.isin(tv, [100, 101])
    h = lambda m: np.histogram(legs.m[m], E[dt])[0] / max(n_correct_all, 1)
    cat, sl = legs.cat.to_numpy(), legs.slice_.to_numpy()
    for c in range(ncat):
        if dt == 10 and c < 3:
            back = (cat == 3) & (sl == c)
            hc = h((cat == c) & corr) + (1 - p["tag_sf"]) * h(back & corr)
            hs = h((cat == c) & swap) + (1 - p["tag_sf"]) * h(back & swap)
        elif dt == 10:
            hc, hs = p["tag_sf"] * h((cat == c) & corr), p["tag_sf"] * h((cat == c) & swap)
        else:
            hc, hs = h((cat == c) & corr), h((cat == c) & swap)
        t = hc + hs
        sig[(dt, c)] = t / max(t.sum(), 1e-300)
        comp[(dt, c)] = (hc / max(t.sum(), 1e-300), hs / max(t.sum(), 1e-300))
        Y[(dt, c)] = float(t.sum())
    return sig, comp, Y


def sf_cat(p, dt, c):
    if not p.get("sf_json"):
        return 1.0
    from .frozen import _sf_cum
    SFJ = json.load(open(os.path.expanduser(p["sf_json"])))
    if dt == 10 and c == 3:
        return _sf_cum(SFJ, 10, 0.65)
    if dt == 10 and c == 4:
        return (_sf_cum(SFJ, 10, 0.9) * 0.9 - _sf_cum(SFJ, 10, 0.65) * 0.65) / 0.25
    e0, e1 = p["cats"][c]
    return (_sf_cum(SFJ, dt, e1) * e1 - (_sf_cum(SFJ, dt, e0) * e0 if e0 > 0 else 0.0)) / (e1 - e0)


# ----------------------------------------------------------------------------- normalizacja
def n_dstar_data(D, p, D_mc=None):
    """N_D* w danych do normalizacji par. Nominał: `dstar_refit.fit_single` (2 Gaussy z kształtem z MC truth D*, skala σ wolna ≥0,8,
    tło u^a(1+c1t+c2t²) [NDS; pełny model z e^{-bu} jako kolejna próba], zakres [m_π+0,9; 160] MeV) — odporny na artefakt
    limitu σ≤2,5 MeV starego `massfit.fit_dm` (V).
    Zwraca dict z N/N_err nominału + `old` (fit_dm, do porównania) + `method`. Fallback na fit_dm tylko gdy refit niepoprawny."""
    from ..massfit import fit_dm
    w = p["dstar_norm"]["d0_window"]
    sel = np.abs(D.m - D.dm - M0[10] / 1000) < w
    dm_mev = D.dm[sel].to_numpy(float) * 1000.0
    old = fit_dm(D.dm[sel].to_numpy(), *p["dstar_norm"]["fit_range"])
    try:
        from .. import dstar_refit as R
        shape = None
        if D_mc is not None and len(D_mc):
            tm = D_mc[(np.abs(D_mc.m - D_mc.dm - M0[10] / 1000) < w) & D_mc.truth.isin([0, 1])]
            if len(tm) >= 200:
                shape = R.shape_from_mc(tm.dm.to_numpy(float) * 1000.0)
        if shape is None:
            shape = {"mu": 145.43, "sc": 0.6, "k": 3.0, "f": 0.7}
            log.warning("pairs: brak MC D* do kształtu Δm — kształt domyślny (mu 145,43; σ 0,6 MeV; k 3; f 0,7)")
        edges = np.arange(R.LO, R.HI + 1e-9, 0.1)
        h, _ = np.histogram(dm_mev, edges)
        f = R.fit_single(h, edges, shape)
        f = {**f, "method": "dstar_refit.fit_single", "shape_mc": shape, "old": old}
        if not f["valid"] or f["N"] <= 0:
            raise RuntimeError(f"refit niepoprawny (valid={f['valid']}, N={f['N']:.1f})")
        if old["valid"] and old["N"] > 0 and abs(f["N"] / old["N"] - 1) > 0.25:
            log.warning("pairs: N_D* refit/fit_dm = %.3f — duża rozbieżność (artefakt σ≤2,5 MeV w fit_dm?)", f["N"] / old["N"])
        return f
    except Exception as e:  # noqa: BLE001
        log.warning("pairs: dstar_refit nie zadziałał (%r) — N_D* z massfit.fit_dm (σ≤2,5 MeV: znany artefakt V, μ obciążone)", e)
        return {**old, "method": "massfit.fit_dm (fallback)", "old": old, "fallback_reason": repr(e)}


def check_event_key_unique(L: pd.DataFrame, what: str = "MC"):
    """Twardy test klucza zdarzenia (source_file, entry) = (run, event): wiersze (cand = wiersz h5) jednego zdarzenia muszą tworzyć
    rozłączne przedziały. W MC fd (i v2.0) event_number powtarza się między plikami AOD (errata v2.0) → pary z różnych zdarzeń
    (audyt 16.11: ×70–140 fałszywych par). mb_nd v1.1/v1.2 ma unikalne numery. Przerywa, gdy przedziały kluczy się nakładają."""
    if not len(L):
        return
    g = L.groupby(["source_file", "entry"], sort=False)["cand"].agg(["min", "max"]).sort_values("min")
    lo, hi = g["min"].to_numpy(), g["max"].to_numpy()
    bad = int((np.maximum.accumulate(hi)[:-1] >= lo[1:]).sum())
    if bad:
        raise SystemExit(f"{what}: klucz (run, event) nie jest unikalny — {bad} zdarzeń nakłada się wierszami z innymi "
                         f"(fd MC / v2.x? użyj bloku zdarzenia lub event_uid). Pary byłyby fałszywe.")


def mc_rates(legs_c: dict, D, p, diff=None):
    """Nieważone tempa w MC minbias ND (h5): poprawne D, pary poprawnych (bez wspólnych śladów), D* truth (to samo okno).
    PD: `diff` → dodatkowo `n_pairs_diff[ch]` = pary MC per bin różniczkowy (ta sama zmienna co w danych) i `n_pairs_diff_undefined`."""
    w = p["dstar_norm"]["d0_window"]
    dsel = (np.abs(D.m - D.dm - M0[10] / 1000) < w) & D.truth.isin([0, 1]) & (D.dm < p["dstar_norm"]["fit_range"][1])
    n_ds = int(dsel.sum())
    L = pd.concat([legs_c[10], legs_c[4]], ignore_index=True)
    check_event_key_unique(L, "mc_rates")
    P = L.merge(L, on=["source_file", "entry"], suffixes=("1", "2"))
    P = P[(P.decay_type1 < P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.cand1 < P.cand2))]
    P = P[~shared_track(P)]
    out = {"n_dstar": n_ds, "n_single": {dt: int(len(legs_c[dt])) for dt in (10, 4)}, "n_pairs": {}}
    for ch, (a, b) in CHANNELS.items():
        out["n_pairs"][ch] = int(((P.decay_type1 == min(a, b)) & (P.decay_type2 == max(a, b))).sum())
    if diff is not None:
        from .diff import diff_bin, nbins
        db = diff_bin(P, diff)
        out["n_pairs_diff"], out["n_pairs_diff_undefined"] = {}, {}
        for ch, (a, b) in CHANNELS.items():
            s = ((P.decay_type1 == min(a, b)) & (P.decay_type2 == max(a, b))).to_numpy()
            out["n_pairs_diff"][ch] = [int((s & (db == k)).sum()) for k in range(nbins(diff))]
            out["n_pairs_diff_undefined"][ch] = int((s & (db < 0)).sum())
    return out


def mc_pair_rel_stat(n_pairs: dict, mb: dict, p: dict) -> dict:
    """Względna stat. MC par w trybie direct. n_pair ∝ Σw (wagi krotności 2D) → 1/√n_eff, n_eff = (Σw)²/Σw² (VER_budget 16.11:
    1/√n zaniżało wiersz mc_stat 0,188 → 0,308 wspólnie). Źródło n_eff: `pairs.mc_pairs_neff` {kanał: n_eff} w configu albo
    `pairs_<ch>.n_eff` w mb_rates_json; bez nich stare 1/√n z ostrzeżeniem."""
    ov = p.get("mc_pairs_neff") or {}
    out = {}
    for ch in CHANNELS:
        ne = ov.get(ch, (mb.get(f"pairs_{ch}") or {}).get("n_eff"))
        if ne is None:
            log.warning("pairs: brak n_eff dla %s (pairs.mc_pairs_neff / mb_rates_json) → mc_stat = 1/√n (zaniżone przy MC ważonym)", ch)
            ne = n_pairs[ch]
        out[ch] = 1 / np.sqrt(max(float(ne), 1.0))
    return out


def weight_factors(mb):
    """c_w = (X/D*)_ważone / (X/D*)_nieważone z mb_rates_2d.json (wagi krotności 2D, 09-26)."""
    ds_w, ds_u = mb["dstar_true_per_event"]["w"], mb["dstar_true_per_event"]["unw"]
    cs = {10: (mb["D0_correct_per_event"]["w"] / ds_w) / (mb["D0_correct_per_event"]["unw"] / ds_u),
          4: (mb["Dplus_correct_per_event"]["w"] / ds_w) / (mb["Dplus_correct_per_event"]["unw"] / ds_u)}
    cp = {ch: (mb[f"pairs_{ch}"]["per_event_w"] / ds_w) / (mb[f"pairs_{ch}"]["per_event_unw"] / ds_u) for ch in CHANNELS}
    ratios_w = {"pair_per_dstar": {ch: mb[f"pairs_{ch}"]["per_event_w"] / ds_w for ch in CHANNELS},
                "single_per_dstar": {10: mb["D0_correct_per_event"]["w"] / ds_w, 4: mb["Dplus_correct_per_event"]["w"] / ds_w}}
    return cs, cp, ratios_w


# ----------------------------------------------------------------------------- komenda
def cmd_pairs(cfg, unblind=False, norm_mode="direct", diff=None):
    """PD: `diff` = zmienna różniczkowa z CLI ('cfg' → z `pairs.diff`) albo None (wyłączone — wyjścia jak dotąd)."""
    p = pcfg(cfg)
    var = p["variant"]
    E = edges(p)
    out = pdir(cfg)
    os.makedirs(out, exist_ok=True)
    from .diff import spec_from_cfg
    spec = spec_from_cfg(p, diff)
    if spec is not None and not unblind:
        raise SystemExit("pairs --diff: pary różniczkowe wymagają jawnego --unblind (pomiar jest rozszerzeniem odsłoniętego wyniku)")
    kin = spec is not None
    if kin:
        p = {**p, "diff": spec}                                   # znormalizowana specyfikacja (make_legs: d0_flavour)
    sdir = lambda sp: os.path.join(cfg["paths"].get("scores_dir") or os.path.join(cfg["paths"]["out_dir"], "scores"), f"{sp}_{var}")  # S5 nom_mix (z Q_lxplus/s5_bkg/src, brakowało w L): honoruj paths.scores_dir
    cols = ["row", "run", "event", "m", "pt", "score", "veto_ok", "truth"]
    if p.get("d0_arbitration", "hash") == "arbiter":
        cols = cols + ["p_hyp"]                                   # HYP (store.load pomija brakujące → błąd w arbitrate_arbiter)
    ncat = {10: 5, 4: 3}
    meta = {"variant": var, "config": {k: v for k, v in p.items()}, "unblind": bool(unblind)}
    if spec is not None:
        meta["diff"] = spec

    # --- progi kategorii z MC sygnału (OOF), poprawne kandydaty w zakresie mas
    thr, scs = {}, {}
    for sp, dt in (("D0", 10), ("Dplus", 4)):
        s = store.load(os.path.join(sdir(sp), "signal_mc.h5"), cols)
        s = s[np.isin(s.truth, [0, 1]) & (s.m * 1000 >= E[dt][0]) & (s.m * 1000 < E[dt][-1])]
        if dt == 4:
            s = s[s.veto_ok.astype(bool)]
        thr[dt] = thresholds(s, [0.3, 0.45, 0.65, 0.9])
        scs[dt] = s
    meta["thresholds"] = {str(k): v for k, v in thr.items()}
    log.info("pairs: progi score (ε_S 30/45/65/90%%): D0 %s, D+ %s", thr[10], thr[4])

    def legs_for(key, keep_truth=False):
        S = {}
        for sp, dt in (("D0", 10), ("Dplus", 4)):
            x = store.load(os.path.join(sdir(sp), f"{key}.h5"), cols)
            lo = thr[dt]["0.9"] if dt == 10 else thr[dt]["0.65"]
            S[dt] = x[(x.score >= lo) & (x.m * 1000 >= E[dt][0]) & (x.m * 1000 < E[dt][-1])]
        rows = np.concatenate([S[10].row.to_numpy(), S[4].row.to_numpy()])
        T, D, nev = scan_source(cfg, key, rows, chunk=p["chunk_rows"], kin=kin, charge_src=charge_files(cfg, spec, key) if kin else None)
        tags = tag_keys(D, p)
        legs = {dt: make_legs(S[dt], T, dt, thr[dt], p, tags if dt == 10 else None, E) for dt in (10, 4)}
        return legs, D, nev

    # --- MC sygnału: szablony, Y
    legs_mc, D_mc, _ = legs_for("signal_mc")
    sig, comp, Y = {}, {}, {}
    for dt in (10, 4):
        a, b, c = signal_templates(legs_mc[dt][np.isin(legs_mc[dt].truth, [0, 1, 100, 101])], len(scs[dt]), dt, E, p, ncat[dt])
        sig.update(a)
        comp.update(b)
        Y.update({k: v * sf_cat(p, *k) for k, v in c.items()})
    meta["Y_raw"] = {f"{k[0]}_{k[1]}": v for k, v in Y.items()}
    meta["swapfrac"] = {f"{k[0]}_{k[1]}": float(v[1].sum()) for k, v in comp.items()}
    log.info("pairs: Y (z SF) %s; udział swapów %s", meta["Y_raw"], meta["swapfrac"])

    # --- dane
    legs_d, D_d, nev_d = legs_for("data_none")
    L = pd.concat([legs_d[10], legs_d[4]], ignore_index=True)
    sb = L[np.abs(L.m - L.decay_type.map(M0)) > L.decay_type.map(SBW)]
    dataleg = {(dt, c): np.histogram(sb.m[(sb.decay_type == dt) & (sb.cat == c)], E[dt])[0].astype(float) for dt in (10, 4) for c in range(ncat[dt])}
    sbsb = sbsb_counts(sb, ncat)
    fD = n_dstar_data(D_d, p, D_mc=D_mc)
    meta.update({"n_events_with_candidates": nev_d, "legs_per_cat": {f"{dt}_{c}": int(((L.decay_type == dt) & (L.cat == c)).sum())
                                                                    for dt in (10, 4) for c in range(ncat[dt])}, "dstar_fit_data": fD})
    log.info("pairs: dane %d zdarzeń z kandydatami; N_D* (%s) = %.0f ± %.0f (fit_dm: %.0f)", nev_d, fD.get("method"), fD["N"], fD["N_err"],
             fD.get("old", {}).get("N", float("nan")))

    # --- normalizacja z MC minbias ND
    if not p.get("mb_rates_json"):
        raise SystemExit("pairs: brak `pairs.mb_rates_json` w configu (plik mb_rates_2d.json z zamrożonych wejść; patrz docs/double_charm.md)")
    mb = json.load(open(os.path.expanduser(p["mb_rates_json"])))
    cs, cp, rw = weight_factors(mb)
    n_pair, n_real = {}, {}
    alt = {}
    # (a) 'frozen': ważone stosunki X/D* z 09-26 (selekcja bazowa ntupli CharmAnalysis; zakłada tę samą wydajność bazową)
    alt["frozen"] = {"n_pair": {ch: fD["N"] * rw["pair_per_dstar"][ch] for ch in CHANNELS},
                     "n_real": {dt: fD["N"] * rw["single_per_dstar"][dt] for dt in (10, 4)}}
    # (b) 'direct': nieważone stosunki X/D* w MC minbias ND z h5 v1.1 (ta sama selekcja co dane) × c_w
    legs_mb, D_mb, _ = legs_for_truth(cfg, p, "incl_mc", E, sdir, cols, kin=spec)
    r = mc_rates(legs_mb, D_mb, p, diff=spec)
    meta["mc_rates_unweighted"] = r
    ok_direct = r["n_dstar"] >= 50 and all(v >= 5 for v in r["n_pairs"].values())
    if r["n_dstar"] > 0:
        alt["direct"] = {"n_pair": {ch: fD["N"] * r["n_pairs"][ch] / r["n_dstar"] * cp[ch] for ch in CHANNELS},
                         "n_real": {dt: fD["N"] * r["n_single"][dt] / r["n_dstar"] * cs[dt] for dt in (10, 4)},
                         "rel_stat": mc_pair_rel_stat(r["n_pairs"], mb, p)}
    # (c) 'truth' (etap T): prawda Pythii 8 × mapa ε (truth_norm.json z lxplus), addytywnie
    tn = p.get("truth_norm_json")
    if tn and os.path.exists(os.path.expanduser(tn)):
        from .truthnorm import truth_alternative
        alt["truth"] = truth_alternative(os.path.expanduser(tn), fD["N"], cp, cs, apply_cw=p.get("truth_apply_cw", True))
        alt["truth_nocw"] = truth_alternative(os.path.expanduser(tn), fD["N"], cp, cs, apply_cw=False)
    if norm_mode == "truth" and "truth" not in alt:
        raise SystemExit("pairs --norm truth: brak pliku pairs.truth_norm_json (%s)" % tn)
    if norm_mode == "direct" and not ok_direct:
        log.warning("pairs: za mało par/D* w MC minbias h5 (%s, D*=%d) → normalizacja 'frozen'", r["n_pairs"], r["n_dstar"])
        norm_mode = "frozen"
    n_pair, n_real = alt[norm_mode]["n_pair"], alt[norm_mode]["n_real"]
    meta["norm_alternatives"] = alt
    meta.update({"norm_mode": norm_mode, "c_w_pair": cp, "c_w_single": {str(k): v for k, v in cs.items()},
                 "lumi_fraction": p["lumi_fraction"], "R_reference": p["R_reference"]})
    I = Inputs(edges={k: v for k, v in E.items()}, ncat=ncat, sig=sig, sig_comp=comp, Y=Y, dataleg=dataleg, nev_legs=float(nev_d),
               ndata=float(nev_d), n_real=n_real, n_pair=n_pair, sbsb=sbsb, meta=meta)
    I.to_files(os.path.join(out, "pair_inputs"))
    region = "full" if unblind else "blind"
    if spec is None:
        H, counts = build_pair_hists(L, I, region)
    else:
        H, counts, Hd = build_pair_hists(L, I, region, diff=spec)            # H identyczne jak bez diff (test); Hd per bin
    np.savez(os.path.join(out, f"pairs_{region}.npz"), **H)
    if spec is not None:
        from .diff import diff_inputs_json
        np.savez(os.path.join(out, f"pairs_{region}_diff.npz"), **Hd)
        counts["mc_undefined"] = r.get("n_pairs_diff_undefined", {})
        Dj = diff_inputs_json(I, spec, sbsb_counts(sb, ncat, diff=spec), r.get("n_pairs_diff", {}), counts)
        with open(os.path.join(out, "pair_inputs_diff.json"), "w") as f:
            json.dump(Dj, f, indent=1, default=float)
        log.info("pairs diff[%s]: biny %s; pary danych per bin %s (nieokreślone %d); pary MC per bin %s (nieokreślone %s)", spec["var"],
                 Dj["labels"], counts["diff_kept"], counts["diff_undefined"], Dj["n_mc_pairs"], Dj["n_mc_undefined"])
    with open(os.path.join(out, "pairs_counts.json"), "w") as f:  # PF: liczniki weta wspólnych śladów (pairs_before_veto, pairs_shared_track)
        json.dump({**counts, "track_id": p["track_id"]}, f, indent=1)
    if p.get("mix_control"):  # etap S5: kontrola tła par mixingiem nóg (opcja, domyślnie wyłączona)
        from .mixctl import run as _mix_run
        _mix_run(L, I, out, p)
    legs_out = legs_for_export(L, unblind)
    if _has_pa():
        legs_out.to_parquet(os.path.join(out, "legs_data.parquet" if unblind else "legs_data_sb.parquet"), index=False)
    log.info("pairs: region=%s, par w histogramach %d; oczekiwane pary (pełna płaszczyzna przed kat.) %s", region, counts["pairs_kept"],
             {k: round(v, 1) for k, v in n_pair.items()})
    return I, H, counts


def legs_for_export(L: pd.DataFrame, unblind: bool) -> pd.DataFrame:
    """Nogi danych do zapisu. ŚLEPO: tylko nogi z sidebandów (|m−m0| > 75/70 MeV) — z pełnej listy nóg centrum SS dałoby się
    odtworzyć jednym joinem po zdarzeniu (recenzja kodu 2.0.1, P0-1). Pełna lista dopiero z --unblind."""
    cols = ["source_file", "entry", "cand", "decay_type", "m", "cat", "trk0", "trk1", "trk2"]
    cols += [c for c in ("eta", "phi", "charge", "dq0") if c in L]          # PD: kinematyka nóg (tylko gdy pairs --diff)
    if unblind:
        return L[cols]
    sb = np.abs(L.m - L.decay_type.map(M0)) > L.decay_type.map(SBW)
    return L.loc[sb, cols]


def legs_for_truth(cfg, p, key, E, sdir, cols, kin=None):
    """MC minbias: poprawne D (truth 0/1, bez cięcia score) ze śladami — do temp par/D*. PD: `kin` = spec diff → także η/φ/ładunki."""
    S = {}
    cm_d0 = bool(kin) and (kin or {}).get("d0_flavour", "closest_mass") == "closest_mass"
    for sp, dt in (("D0", 10), ("Dplus", 4)):
        x = store.load(os.path.join(sdir(sp), f"{key}.h5"), cols)
        tv = [0, 1, 100, 101] if (cm_d0 and dt == 10) else [0, 1]   # PD: swapy potrzebne do wyboru hipotezy bliższej m0 (potem odrzucane)
        x = x[np.isin(x.truth, tv) & (x.m * 1000 >= E[dt][0]) & (x.m * 1000 < E[dt][-1])]
        if dt == 4:
            x = x[x.veto_ok.astype(bool)]
        S[dt] = x
    T, D, nev = scan_source(cfg, key, np.concatenate([S[10].row.to_numpy(), S[4].row.to_numpy()]), chunk=p["chunk_rows"], kin=bool(kin),
                            charge_src=charge_files(cfg, kin, key) if kin else None)
    legs = {}
    for dt in (10, 4):
        z = S[dt].merge(T, on="row")
        z = z.assign(source_file=z.run.astype(np.int64).astype(str), entry=z.event.astype(np.int64), cand=z.row.astype(np.int64),
                     decay_type=dt, m=z.m * 1000.0)
        if cm_d0 and dt == 10:
            z = d0_flavour_closest_mass(z)
            z = z[np.isin(z.truth, [0, 1])]
        legs[dt] = z
    return legs, D, nev


def _has_pa():
    try:
        import pyarrow  # noqa: F401
        return True
    except Exception:
        return False
