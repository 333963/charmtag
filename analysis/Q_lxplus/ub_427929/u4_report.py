#!/usr/bin/env python3
"""UB427 raport (18.11.2026, projekt A): μ par D–D w runie 427929 (klasa MB) vs nominał 428770.

Wejścia: wyjście `charmtag pairs --unblind` + `fit2d --unblind` dla 427929 (pair_inputs.*, pairs_full.npz, fit2d_unblind.json),
to samo dla 428770 (odsłonięcie 05.11), cw.json (u4_cw.py), statystyki filtra. μ wspólne i przedział 68 % tą samą funkcją co
jackknife odsłonięcia (`ub_jk/ubjk_run.replica_fit`: minimum sumy profili kanałów). Porównanie: Δ = μ427 − μ428 w σ
(stat. 427929 z profilu ⊕ stat. 428770 = +0,48/−0,46 po jackknife zdarzeń; systematyki MC wspólne → nie wchodzą do Δ).
  u4_report.py --pairs PAIRS427 --nominal PAIRS428 [--cw cw.json] [--fstats DIR] --out DIR"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ub_jk"))
STAT428 = (0.46, 0.48)          # μ = 1,50 +0,48/−0,46 (stat, jackknife zdarzeń ×1,18; UNB_WYNIK §1b)


N_NDS_428770 = (4436.390108391369, 670.8561890721969)   # NDS rerun 12.11 (nominał μ = 1,50 liczony na tej bazie)


def rescale(R, n_new, n_new_err):
    """μ ∝ 1/N_D* dokładnie (NDS ndstar.md §2): przeskalowanie μ, błędów i przedziału z N_D* biegu par na N_D* metody NDS."""
    k = R["N_dstar"] / n_new
    f = R["fit"]
    for ch, v in f.items():
        if isinstance(v, dict) and "mu_hat" in v:
            v["mu_hat"] *= k; v["mu_err_hesse"] *= k
    c = f["combined"]
    c["mu"] *= k
    if c.get("mu_68"):
        c["mu_68"] = [x * k if x is not None else None for x in c["mu_68"]]
        c["sigma_68"] = (c["sigma_68"] * k) if c.get("sigma_68") else None
    R.update({"N_dstar_pairs_run": R["N_dstar"], "N_dstar": n_new, "N_dstar_err": n_new_err, "rescale_k": k})
    return R


def fit_dir(d, U):
    _, CHANNELS, Model, load_inputs, _, _ = U._charmtag()
    I = load_inputs(os.path.join(d, "pair_inputs"))
    H = dict(np.load(os.path.join(d, "pairs_full.npz")))
    r = U.replica_fit(Model(I), H, list(CHANNELS), combined=True, interval=True)
    js = json.load(open(os.path.join(d, "pair_inputs.json")))
    meta = js.get("meta", {})
    fD = meta.get("dstar_fit_data", {})
    return {"fit": r, "N_dstar": fD.get("N"), "N_dstar_err": fD.get("N_err"), "nev": meta.get("n_events_with_candidates"),
            "n_pair_exp": js.get("n_pair"), "c_w_pair": meta.get("c_w_pair"), "pairs_hist": {ch: int(H[ch].sum()) for ch in CHANNELS},
            "fit2d": json.load(open(os.path.join(d, "fit2d_unblind.json"))) if os.path.exists(os.path.join(d, "fit2d_unblind.json")) else None}, CHANNELS


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True); ap.add_argument("--nominal", required=True); ap.add_argument("--cw"); ap.add_argument("--fstats")
    ap.add_argument("--ndstar", help="ndstar_427929.json (u4_ndstar.py) → μ 427929 przeskalowane na N_D* metody NDS")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    import ubjk_run as U
    R7, CH = fit_dir(a.pairs, U)
    R8, _ = fit_dir(a.nominal, U)
    R8 = rescale(R8, *N_NDS_428770)                                  # 05.11: 5584 (fallback fit_dm) → 4436 (NDS) ⇒ μ = 1,50
    if a.ndstar:
        nd = json.load(open(a.ndstar))
        if not nd["fit"].get("valid"):
            raise SystemExit(f"STOP: fit N_D* 427929 (NDS) niepoprawny: {nd['fit']}")
        R7 = rescale(R7, nd["fit"]["N"], nd["fit"]["N_err"]); R7["ndstar_nds"] = nd
    c7, c8 = R7["fit"]["combined"], R8["fit"]["combined"]
    lo7, hi7 = c7.get("mu_68") or [None, None]
    d = c7["mu"] - c8["mu"]
    s7 = (c7["mu"] - lo7) if d > 0 else (hi7 - c7["mu"])               # strona 427929 w kierunku 428770
    s8 = STAT428[1] if d > 0 else STAT428[0]
    r7 = R7["N_dstar_err"] / R7["N_dstar"]; r8 = N_NDS_428770[1] / N_NDS_428770[0]
    tot = math.sqrt(s7 ** 2 + s8 ** 2 + (c7["mu"] * r7) ** 2 + (c8["mu"] * r8) ** 2) if (s7 and s8) else None
    rep = {"r427929": R7, "r428770": R8, "delta": d, "sigma_427929_side": s7, "sigma_428770_side": s8,
           "delta_sigma": d / math.hypot(s7, s8) if (s7 and s8) else None,
           "delta_sigma_with_ndstar_stat": d / tot if tot else None, "ndstar_rel_stat": {"427929": r7, "428770": r8}}
    if a.cw and os.path.exists(a.cw):
        cw = json.load(open(a.cw)); rep["cw"] = {k: cw.get(k) for k in ("base_428770", "R", "R_vs_alldata", "c_w_nominal", "c_w_427929", "repro_428770_vs_nominal")}
        rep["cw_neff"] = {k: v.get("n_eff_pairs") for k, v in cw.get("samples", {}).items()}
        rep["cw_class"] = {k: v.get("info") for k, v in cw.get("samples", {}).items()}
    if a.fstats:
        S = []
        for f in sorted(glob.glob(os.path.join(a.fstats, "*.json"))):
            S += json.load(open(f))
        tot = {k: int(sum(s.get(k, 0) for s in S)) for k in ("rows_in", "rows_out", "D0_in", "D0_out", "Dstar_in", "Dstar_out", "Dplus_in", "Dplus_out",
                                                            "evtable_in", "evtable_class", "evtable_grl")}
        tot["files"] = len(S)
        rep["filter"] = tot
    json.dump(U._clean(rep), open(os.path.join(a.out, "ub427.json"), "w"), indent=1)
    open(os.path.join(a.out, "UB427.md"), "w").write(to_md(rep, CH))
    print(open(os.path.join(a.out, "UB427.md")).read())
    return 0


def _row(nm, R, CH):
    f = R["fit"]; c = f["combined"]; lo, hi = c.get("mu_68") or [float("nan")] * 2
    return (f"| {nm} | {R['N_dstar']:.0f} ± {R['N_dstar_err']:.0f} | {R['nev']} | " + " / ".join(str(R["pairs_hist"][ch]) for ch in CH) + " | "
            + " | ".join(f"{f[ch]['mu_hat']:.2f} ± {f[ch]['mu_err_hesse']:.2f}" for ch in CH)
            + f" | **{c['mu']:.2f}** [{lo:.2f}; {hi:.2f}] | {c.get('Z', float('nan')):.2f} |")


def to_md(rep, CH):
    L = ["# UB427 — μ par D–D w czystym runie minimum-bias 427929 (klasa MB: GRL ∧ HLT_mb_mbts_L1MBTS_1) vs 428770", "",
         "Ta sama selekcja co nominał (h5 v2.0 z cięciem L_xy(D0, D*) ≥ 0,25 mm i D0 1,75–2,00 GeV ≡ v1.2/v1.1), te same modele H, progi, szablony,",
         "Y×SF i model fit2d; N_D* metodą NDS (fit_single z poprawką NLL, kształt MC V; μ ∝ 1/N_D* przeskalowane z N_D* biegu par, jak nominał 5584 → 4436); c_w = c_w(nominał) × R_X, R_X = c_X(427929)/c_X(428770) z wag w2 do rozkładu (n_trk, lead_pt) każdego runu.",
         "Runy NIE łączone. μ wspólne = minimum sumy profili kanałów (jak `ub_jk`), 68 % z profilu (bez korekty jackknife zdarzeń).", "",
         "| run | N_D* | zdarzenia z kand. | pary w hist. D0D0 / DpD0 / DpDp | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | Z |", "|---|---|---|---|---|---|---|---|---|",
         _row("427929 (MB)", rep["r427929"], CH), _row("428770 (nominał, odtw.)", rep["r428770"], CH), "",
         f"Δμ = μ427929 − μ428770 = {rep['delta']:+.2f}; σ = {rep['sigma_427929_side']:.2f} (427929, profil) ⊕ {rep['sigma_428770_side']:.2f} (428770, jackknife) "
         f"→ **{rep['delta_sigma']:+.2f}σ**; z błędem stat. N_D* obu runów ({100 * rep['ndstar_rel_stat']['427929']:.0f} % / {100 * rep['ndstar_rel_stat']['428770']:.0f} %): "
         f"**{rep['delta_sigma_with_ndstar_stat']:+.2f}σ**." if rep.get("delta_sigma") is not None else "Δμ: przedział 68 % 427929 nieokreślony.", ""]
    cw = rep.get("cw")
    if cw:
        L += ["## c_w (wagi zdarzeń) dla 427929", "",
              "| X | c_w nominał | R_X = c(427929)/c(" + str(cw["base_428770"]) + ") | R_X vs wszystkie dane 428770 | c_w 427929 | odtw. 428770 / nominał |", "|---|---|---|---|---|---|"]
        for x in ("D0", "Dplus", "D0D0", "DpD0", "DpDp"):
            L.append(f"| {x} | {cw['c_w_nominal'][x]:.3f} | {cw['R'][x]:.3f} | {cw['R_vs_alldata'][x]:.3f} | {cw['c_w_427929'][x]:.3f} | {cw['repro_428770_vs_nominal'][x]:.3f} |")
        L += ["", f"n_eff par MC przy wagach: {json.dumps({k: {c: round(v, 1) for c, v in (d or {}).items()} for k, d in rep['cw_neff'].items()}, ensure_ascii=False)}",
              f"Zdarzenia w mapach: {json.dumps(rep['cw_class'], ensure_ascii=False)}", ""]
    ft = rep.get("filter")
    if ft:
        L += ["## Filtr h5 v2.0 → selekcja nominału", "",
              f"Pliki {ft['files']}; wiersze {ft['rows_in']} → {ft['rows_out']}; D0 {ft['D0_in']} → {ft['D0_out']} (L_xy ≥ 0,25 i okno), D* {ft['Dstar_in']} → {ft['Dstar_out']}, "
              f"D⁺ {ft['Dplus_in']} → {ft['Dplus_out']}; zdarzenia (evtable) {ft['evtable_in']}, z GRL {ft['evtable_grl']}, w klasie MB {ft['evtable_class']} "
              f"({100 * ft['evtable_class'] / max(ft['evtable_in'], 1):.1f} %).", ""]
    L += ["Uwagi: σ 427929 bez współczynnika korelacji zdarzeń (w 428770 jackknife dał ×1,18); systematyki normalizacji MC (mc_stat, model wag, tune) są w obu runach",
          "wspólne i nie wchodzą do Δ; c_w 427929 zależy od ekstrapolacji w2 do rozkładu krotności 427929 (n_eff par w MC w tabeli).", ""]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
