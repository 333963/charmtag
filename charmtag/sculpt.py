"""Metryki rzeźbienia masy do charmtag validate (etap X-val, addytywnie; nic w pipeline nie zmienia).

Dlaczego: charmtag `rwin(..., "window")` liczy okno vs CAŁE SB w zakresie 1.70–2.20 GeV. Dla D⁺ SB są
asymetryczne (lewy 100 MeV, prawy 260 MeV), a ε_B(m) klasyfikatora maleje monotonicznie z m
(×2.5 na 0.5 GeV) → R_win ≈ 1.24 nawet gdy pod pikiem NIE ma żadnej struktury. Poniżej metryki,
które rozdzielają gładkie nachylenie (nieszkodliwe dla fitu exp/pol2) od lokalnego rzeźbienia.
"""
from __future__ import annotations
import numpy as np


def rwin_sym(m, passed, m0, win, sb_min, sb_max, w=None, lo=None, hi=None):
    """R = [N_okno/N_SB]_po / [...]_przed z SYMETRYCZNYM SB (sb_min<|m−m0|<sb_max po obu stronach;
    sb_max przycinane do min(m0−lo, hi−m0), żeby granica zakresu nie łamała symetrii).
    Nachylenie liniowe ε(m) znosi się w pierwszym rzędzie → R≈1 ⇔ brak lokalnej struktury."""
    m = np.asarray(m, float); passed = np.asarray(passed, bool)
    if lo is not None: sb_max = min(sb_max, m0 - lo)
    if hi is not None: sb_max = min(sb_max, hi - m0)
    w = np.ones(len(m)) if w is None else np.asarray(w, float)
    dm = np.abs(m - m0); A = dm < win; B = (dm > sb_min) & (dm < sb_max)
    nA, nB, kA, kB = w[A].sum(), w[B].sum(), w[A & passed].sum(), w[B & passed].sum()
    if min(nA, nB, kA, kB) <= 0:
        return float("nan"), float("nan")
    r = (kA / kB) / (nA / nB)
    return float(r), float(r * np.sqrt((1 - kA / nA) / kA + (1 - kB / nB) / kB))


def slope_eff(m, passed, m0, sb_min, lo, hi, nb=25, min_n=50):
    """Nachylenie względnej wydajności cięcia ε(m)/ε̄ [1/GeV] z fitu liniowego w binach m poza |m−m0|<sb_min."""
    m = np.asarray(m, float); passed = np.asarray(passed, bool)
    e = np.linspace(lo, hi, nb + 1); c = 0.5 * (e[1:] + e[:-1])
    n, _ = np.histogram(m, e); k, _ = np.histogram(m[passed], e)
    ok = (n > min_n) & (np.abs(c - m0) > sb_min)
    if ok.sum() < 3:
        return float("nan")
    eff = k[ok] / n[ok]
    if not np.isfinite(eff.mean()) or eff.mean() <= 0:   # nic nie przeszło cięcia → nachylenie nieokreślone
        return float("nan")
    return float(np.polyfit(c[ok] - m0, eff / eff.mean(), 1)[0])


def pol2_closure(m, passed, m0, win, sb_min, lo, hi, nb=50):
    """Test istotny dla fitu masy: histogram m PO cięciu (próbka bez sygnału: MC tła / kontrola),
    fit exp(pol2) tylko w SB (|m−m0|>sb_min), ekstrapolacja do okna |m−m0|<win → (N_pred/N_true, σ_stat).
    1 ⇔ tło pod pikiem przewidywalne z sidebandów modelem exp(pol2) (jak w charmtag fit)."""
    m = np.asarray(m, float); passed = np.asarray(passed, bool)
    e = np.linspace(lo, hi, nb + 1); c = 0.5 * (e[1:] + e[:-1]); k, _ = np.histogram(m[passed], e)
    sb = (np.abs(c - m0) > sb_min) & (k > 0); wn = np.abs(c - m0) < win
    if sb.sum() < 6 or k[wn].sum() == 0:
        return float("nan"), float("nan")
    x = c - m0
    p = np.polyfit(x[sb], np.log(k[sb]), 2, w=np.sqrt(k[sb]))
    pred, true = np.exp(np.polyval(p, x[wn])).sum(), k[wn].sum()
    return float(pred / true), float(np.sqrt(true) / true)


def sculpt_report(m, passed, sp, lo=None, hi=None, sb_max=0.20):
    """Komplet dla słownika gatunku charmtag (m0, signal_window, sideband_min, mass_range)."""
    m0, win, sbm = sp["m0"], sp["signal_window"], sp["sideband_min"]
    lo = sp["mass_range"][0] if lo is None else lo; hi = sp["mass_range"][1] if hi is None else hi
    r, re = rwin_sym(m, passed, m0, win, sbm, sb_max, lo=lo, hi=hi)
    pc, pe = pol2_closure(m, passed, m0, win, sbm, lo, hi)
    return {"R_win_sym": r, "R_win_sym_err": re, "slope_perGeV": slope_eff(m, passed, m0, sbm, lo, hi),
            "pol2_closure": pc, "pol2_closure_err": pe}


# ---------------------------------------------------------------- etap VAL: haki dla charmtag validate
def eps_regions(m, passed, m0, win, sb_min, sb_max, lo, hi, inner=None, w=None):
    """ε cięcia w regionach: SR (|m−m0|<win), SB asymetryczne (|m−m0|>sb_min w [lo,hi] — stary proxy ε_B(SB)),
    SB symetryczne (sb_min<|m−m0|<min(sb_max, m0−lo, hi−m0)), lewy/prawy SB, pasmo wewnętrzne (inner, opcjonalnie)."""
    m = np.asarray(m, float); passed = np.asarray(passed, bool)
    w = np.ones(len(m)) if w is None else np.asarray(w, float)
    sbm = min(sb_max, m0 - lo, hi - m0)
    dm = np.abs(m - m0)
    reg = {"SR": dm < win, "SB_asym": (dm > sb_min) & (m >= lo) & (m < hi), "SB_sym": (dm > sb_min) & (dm < sbm),
           "SB_left": (dm > sb_min) & (m < m0) & (m >= lo), "SB_right": (dm > sb_min) & (m > m0) & (m < hi)}
    if inner is not None:
        reg["inner_sym"] = (dm > inner[0]) & (dm < inner[1])
    out = {}
    for k, r in reg.items():
        n = w[r].sum(); kk = w[r & passed].sum()
        out[k] = {"eps": float(kk / n) if n > 0 else float("nan"), "n": int(r.sum()), "n_pass": int((r & passed).sum())}
    return out


def validate_block(m, passed_by_wp: dict, sp, lo=None, hi=None, sb_max=0.20, has_signal=False, w=None):
    """Słownik per WP do validation/<species>.json. has_signal=True (dane): R_win_sym liczone na paśmie wewnętrznym
    (sp['inner_band'], symetrycznie), pol2_closure pomijane (okno zawiera sygnał)."""
    m0, win, sbm = sp["m0"], sp["signal_window"], sp["sideband_min"]
    lo = sp["mass_range"][0] if lo is None else lo; hi = sp["mass_range"][1] if hi is None else hi
    out = {}
    for wp, passed in passed_by_wp.items():
        e = eps_regions(m, passed, m0, win, sbm, sb_max, lo, hi, inner=sp.get("inner_band"), w=w)
        row = {"eps": e, "slope_perGeV": slope_eff(m, passed, m0, sbm, lo, hi)}
        if has_signal and sp.get("inner_band"):
            a, b = sp["inner_band"]
            # pasmo wewnętrzne (a<|dm|<b) vs SB symetryczne, bez centrum (sygnał w |dm|<a wyłączony) — liczone wprost
            m_ = np.asarray(m, float); dm = np.abs(m_ - m0); sbM = min(sb_max, m0 - lo, hi - m0)
            A = (dm > a) & (dm < b); B = (dm > sbm) & (dm < sbM); p = np.asarray(passed, bool)
            ww = np.ones(len(m_)) if w is None else np.asarray(w, float)
            nA, nB, kA, kB = ww[A].sum(), ww[B].sum(), ww[A & p].sum(), ww[B & p].sum()
            if min(nA, nB, kA, kB) > 0:
                r = (kA / kB) / (nA / nB); re = r * np.sqrt((1 - kA / nA) / kA + (1 - kB / nB) / kB)
            else:
                r, re = float("nan"), float("nan")
            row["R_inner_sym"] = float(r); row["R_inner_sym_err"] = float(re)
        else:
            r, re = rwin_sym(m, passed, m0, win, sbm, sb_max, w=w, lo=lo, hi=hi)
            pc, pe = pol2_closure(m, passed, m0, win, sbm, lo, hi)
            row.update({"R_win_sym": float(r), "R_win_sym_err": float(re), "pol2_closure": pc, "pol2_closure_err": pe})
            if e["SR"]["eps"] == e["SR"]["eps"]:  # ε_B pod pikiem (prawda w MC tła) vs proxy
                row["eps_B_SR_true"] = e["SR"]["eps"]
                row["proxy_ratio_SB_asym"] = e["SB_asym"]["eps"] / e["SR"]["eps"] if e["SR"]["eps"] else float("nan")
                row["proxy_ratio_SB_sym"] = e["SB_sym"]["eps"] / e["SR"]["eps"] if e["SR"]["eps"] else float("nan")
        out[str(int(round(wp * 100)))] = row
    return out


LEGEND = ("R_win_window (stare) = okno/CAŁE SB; przy nachyleniu ε_B(m) (slope≠0, własność tła) i asymetrycznych SB "
          "daje >1 BEZ rzeźbienia. Rzeźbienie ⇔ R_win_sym≠1 (SB symetryczne) lub pol2_closure≠1 (MC tła). "
          "ε_B pod pikiem: proxy_ratio_SB_sym≈1, proxy_ratio_SB_asym<1 (asymetryczne SB zaniżają ε_B ~20 % dla D⁺).")
