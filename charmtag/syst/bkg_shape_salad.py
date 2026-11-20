"""B4 — systematyka ``bkg_shape_salad``: kształt tła pod pikiem (ε_B) dla pojedynczych D i par D–D (U).

Źródła:
* korekta centralna ``rel_salad_vs_SB`` z JSON-u S8 (``charmtag s8salad`` → ``<out_dir>/syst/bkg_shape_salad_<sp>.json``,
  klucze per ``wp|bin``) — ile ε_B pod pikiem różni się od proxy z asymetrycznych SB (D⁺: ~+25 %);
* niepewność względna ε_B (B1): D⁺ ±5/10/25 % @WP80/70/50; D0 ±5/10/15 % (z ``mc_truth_min: 2``; ±20 % @50 bez).

Wyjście w formacie jak Y (``syst_bkg_model_<sp>.json``) / U (wiersze per źródło): ``per["wp|bin"]`` i blok ``pairs``
z niepewnością przewidywanego B pod pikiem dla kanałów D0D0 / DplusD0 / DplusDplus (nogi niezależne → kwadratura).
Użycie w U: wiersz „bkg_shape_salad” = ±rel_unc_pair × B_pred(Asimov) → Δμ przez pochodną μ(B) z fit2d (nie dotyczy
kształtu pol2 w ficie — ten jest nieobciążony, Y/X); dla Zproxy/ε_B z SB-asym zastosować ``rel_salad_vs_SB`` jako korektę.

CLI (bez haka w ``cli.py``): ``python -m charmtag.syst.bkg_shape_salad --species Dplus [--s8-json F] [--out DIR]``.
Hak opcjonalny (1 linia w ``_SUBCOMMANDS``): ``"bkgsyst-salad": "charmtag.syst.bkg_shape_salad:main"``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from collections.abc import Iterable

SOURCE = "bkg_shape_salad"
VERSION = "B4-2026-11-05"
WPS = (80, 70, 50)
BINS = ("3-5", "5-8", "8-30")

# B1 (D⁺ §4; D0 §6f). Klucz "D0_truth8" = D0 z MC truth ≥ 8 (bez poprawki mc_truth_min: 2).
DEFAULT_REL_UNC = {
    "Dplus": {80: 0.05, 70: 0.10, 50: 0.25},
    "D0": {80: 0.05, 70: 0.10, 50: 0.15},
    "D0_truth8": {80: 0.05, 70: 0.10, 50: 0.20},
}
PAIR_CHANNELS = {"D0D0": ("D0", "D0"), "DplusD0": ("Dplus", "D0"), "DplusDplus": ("Dplus", "Dplus")}


def _norm_species(sp: str) -> str:
    s = sp.replace("+", "plus").replace("⁺", "plus")
    return "Dplus" if s.lower().startswith("dp") else "D0"


def default_unc(species: str, wp: int, d0_truth_min2: bool = True) -> float:
    sp = _norm_species(species)
    key = sp if (sp == "Dplus" or d0_truth_min2) else "D0_truth8"
    table = DEFAULT_REL_UNC[key]
    if wp in table:
        return table[wp]
    # interpolacja liniowa w WP między znanymi punktami, poza zakresem — najbliższy
    ks = sorted(table)
    if wp >= ks[-1]:
        return table[ks[-1]]
    if wp <= ks[0]:
        return table[ks[0]]
    for a, b in zip(ks[:-1], ks[1:]):  # noqa: RUF007  (py3.9 na lxplus)
        if a <= wp <= b:
            f = (wp - a) / (b - a)
            return table[a] + f * (table[b] - table[a])
    return table[ks[0]]


def _is_bin(k) -> bool:
    return "-" in str(k) and not str(k).lstrip("wpWP").isdigit()


def _wp(k) -> str:
    return str(k).lower().replace("wp", "").strip()


def load_s8(path: str | None) -> dict:
    """Czyta JSON S8 i zwraca {"wp|bin": eps_B(SALAD, SR) / eps_B(SB danych)} (stosunek, 1.0 = brak korekty).

    Prawdziwy wynik ``charmtag s8salad`` (S8/B4, pełna lumi 05.11): ``effB_salad_SR[bin][wp]``, ``effB_data_SB[bin][wp]``
    oraz ``rel_salad_vs_SB[bin][wp] = effB_salad_SR/effB_data_SB − 1`` (bin zewnętrzny, WP wewnętrzny, forma „−1”).
    Stosunek liczony wprost z ``effB_*`` gdy są; w przeciwnym razie z ``rel_salad_vs_SB`` w dowolnym układzie:
    płaski ``wp|bin``, zagnieżdżony ``[wp][bin]`` lub ``[bin][wp]``, lista rekordów {wp, bin, rel_salad_vs_SB};
    wartość |x| < 0.9 traktowana jako forma „−1” (→ 1+x), inaczej jako stosunek. Klucz "all" pomijany."""
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        js = json.load(f)
    out: dict[str, float] = {}

    def put(wp, b, x, ratio: bool):
        if "all" in (str(b), str(wp)):
            return
        x = float(x)
        out[f"{_wp(wp)}|{b}"] = x if ratio else 1.0 + x

    sr, sb = js.get("effB_salad_SR") if isinstance(js, dict) else None, js.get("effB_data_SB") if isinstance(js, dict) else None
    if isinstance(sr, dict) and isinstance(sb, dict):
        for k1, v1 in sr.items():
            if not isinstance(v1, dict):
                continue
            for k2, x in v1.items():
                try:
                    y = sb[k1][k2]
                except (KeyError, TypeError):
                    continue
                if y:
                    b, wp = (k1, k2) if _is_bin(k1) else (k2, k1)
                    put(wp, b, float(x) / float(y), True)
        if out:
            return out
    rel = js.get("rel_salad_vs_SB", js) if isinstance(js, dict) else js
    def as_ratio(x) -> bool:
        return abs(float(x)) >= 0.9

    if isinstance(rel, dict):
        for k, v in rel.items():
            if isinstance(v, dict) and all(isinstance(x, (int, float)) for x in v.values()):
                for k2, x in v.items():  # rel[wp][bin] lub rel[bin][wp]
                    b, wp = (k, k2) if _is_bin(k) else (k2, k)
                    put(wp, b, x, as_ratio(x))
            elif isinstance(v, dict) and "rel_salad_vs_SB" in v:
                out[str(k)] = float(v["rel_salad_vs_SB"]) if as_ratio(v["rel_salad_vs_SB"]) else 1.0 + float(v["rel_salad_vs_SB"])
            elif isinstance(v, (int, float)) and "|" in str(k):
                wp, b = str(k).split("|", 1)
                put(wp, b, v, as_ratio(v))
    elif isinstance(rel, list):
        for rec in rel:
            try:
                put(int(rec["wp"]), rec["bin"], rec["rel_salad_vs_SB"], as_ratio(rec["rel_salad_vs_SB"]))
            except (KeyError, TypeError, ValueError):
                continue
    return out


def load_s8_closure(path: str | None) -> dict:
    """{"wp|bin": eps_B(SALAD, MC SB ważone) / eps_B(SB danych)} — zamknięcie SALAD w pasmach, z ``effB_salad_SB``/``effB_data_SB``
    prawdziwego JSON-u S8 (1.0 = idealne). Pusty dict, gdy brak pól."""
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        js = json.load(f)
    a, b = (js.get("effB_salad_SB"), js.get("effB_data_SB")) if isinstance(js, dict) else (None, None)
    out: dict[str, float] = {}
    if not (isinstance(a, dict) and isinstance(b, dict)):
        return out
    for k1, v1 in a.items():
        if not isinstance(v1, dict):
            continue
        for k2, x in v1.items():
            try:
                y = b[k1][k2]
            except (KeyError, TypeError):
                continue
            bn, wp = (k1, k2) if _is_bin(k1) else (k2, k1)
            if y and "all" not in (str(bn), str(wp)):
                out[f"{_wp(wp)}|{bn}"] = float(x) / float(y)
    return out


def build(
    species: str,
    s8_json: str | None = None,
    wps: Iterable[int] = WPS,
    bins: Iterable[str] = BINS,
    d0_truth_min2: bool = True,
) -> dict:
    sp = _norm_species(species)
    s8 = load_s8(s8_json)
    clo = load_s8_closure(s8_json)
    per = {}
    for wp in wps:
        for b in bins:
            key = f"{wp}|{b}"
            rel = s8.get(key)
            u0 = default_unc(sp, int(wp), d0_truth_min2)
            c = clo.get(key)
            # 05.11 (S8 pełna lumi): niepewność = max(B1 default, |zamknięcie SALAD w SB − 1|) — sk_mlp nie domyka
            # D⁺ 5–8 (1,16) i D0 8–30 @50 (0,64); B1 default liczone na hgb/pol2.
            u = max(u0, abs(c - 1.0)) if c is not None else u0
            per[key] = {
                "rel_salad_vs_SB": rel,
                "rel_unc": u,
                "closure_SB": c,
                "method": "S8 json" if rel is not None else "brak w S8 → korekta 1.0",
                "unc_source": "B1 default" if u == u0 else "|closure_SB-1| (S8) > B1 default",
            }
    return {
        "source": SOURCE,
        "version": VERSION,
        "species": sp,
        "d0_truth_min2": bool(d0_truth_min2) if sp == "D0" else None,
        "s8_json": s8_json if s8 else None,
        "per": per,
        "legend": {
            "rel_salad_vs_SB": "eps_B(SALAD, pod pikiem) / eps_B(SB asymetryczne) — korekta centralna proxy; None = brak pomiaru (1.0)",
            "rel_unc": "niepewnosc wzgledna eps_B pod pikiem: max(B1 closure adwersarialny + fit pol2 danych, |closure_SB-1|); symetryczna",
            "closure_SB": "eps_B(SALAD, MC SB wazone)/eps_B(SB danych) z JSON S8 — zamkniecie w pasmach (None = brak)",
        },
    }


def pair_unc(syst_a: dict, syst_b: dict, wp: int, bin_a: str = "8-30", bin_b: str = "8-30") -> float:
    """Niepewność względna przewidywanego B pod pikiem pary (nogi niezależne): sqrt(u_a² + u_b²)."""
    ua = syst_a["per"][f"{wp}|{bin_a}"]["rel_unc"]
    ub = syst_b["per"][f"{wp}|{bin_b}"]["rel_unc"]
    return math.sqrt(ua * ua + ub * ub)


def pairs_block(syst_d0: dict, syst_dp: dict, wps: Iterable[int] = WPS, bins: Iterable[str] = BINS) -> dict:
    by = {"D0": syst_d0, "Dplus": syst_dp}
    out = {}
    for ch, (a, b) in PAIR_CHANNELS.items():
        out[ch] = {}
        for wp in wps:
            for bn in bins:
                out[ch][f"{wp}|{bn}"] = {
                    "rel_unc_B": pair_unc(by[a], by[b], int(wp), bn, bn),
                    "rel_corr_B": (by[a]["per"][f"{wp}|{bn}"]["rel_salad_vs_SB"] or 1.0)
                    * (by[b]["per"][f"{wp}|{bn}"]["rel_salad_vs_SB"] or 1.0),
                }
    return {
        "source": SOURCE,
        "version": VERSION,
        "channels": out,
        "legend": {
            "rel_unc_B": "niepewnosc wzgledna B_pred pod pikiem pary = kwadratura nog (ten sam bin pT obu nog)",
            "rel_corr_B": "iloczyn korekt rel_salad_vs_SB nog, gdy B_pred liczone z eps_B(SB asym); 1.0 gdy brak S8",
        },
        "note": "Dotyczy przewidywanego B (Asimov/Zproxy/mixing-norm), nie ksztaltu tla w fit2d (tilt+rho — U).",
    }


def write(out_dir: str, species: str, s8_json: str | None = None, **kw) -> str:
    os.makedirs(os.path.join(out_dir, "syst"), exist_ok=True)
    js = build(species, s8_json, **kw)
    path = os.path.join(out_dir, "syst", f"{SOURCE}_syst_{js['species']}.json")
    with open(path, "w") as f:
        json.dump(js, f, indent=1)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", required=True, choices=["D0", "Dplus"])
    ap.add_argument("--s8-json", default=None, help="S8 bkg_shape_salad_<sp>.json (rel_salad_vs_SB)")
    ap.add_argument("--out", default=".", help="out_dir (zapis do <out>/syst/)")
    ap.add_argument("--d0-truth8", action="store_true", help="D0 bez mc_truth_min 2 → ±20 %% @50")
    a = ap.parse_args(argv)
    p = write(a.out, a.species, a.s8_json, d0_truth_min2=not a.d0_truth8)
    print(p)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
