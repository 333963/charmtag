#!/usr/bin/env python3
"""Budżet μ dla 427929 jako pomiaru głównego (19.11.2026, projekt A, studium W2). Lokalnie, tylko z istniejących wyników.
Wiersz po wierszu budżetu 428770 (pairs/syst/syst_blind.json + VER_budget/RAPORT.md + TUNE_pythia/budget_with_tune.md) → status dla 427929:
  „bez zmian” (własność MC/metody), „skaluje” (przeliczone z wielkości 427929), „427929” (zmierzone w 427929), „lxplus” (czeka, wartość tymczasowa).
Konwencja jak charmtag.syst.report: wiersze względne (Δμ/μ), „wspólne” = Σ w·δ, w ∝ 1/σ²_stat — tu σ z fitu 427929 (Hesse), nie Asimov 428770.
  python3 budzet_427929.py  (z katalogu r427929) → BUDZET_427929.md, budzet_427929.json"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))       # charm_work
CH = ["D0D0", "DpD0", "DpDp"]


def J(*p):
    return json.load(open(os.path.join(*p)))


def main():
    ub = J(HERE, "ub427.json"); jk = J(HERE, "jk_ev200", "ubjk.json") if os.path.exists(os.path.join(HERE, "jk_ev200", "ubjk.json")) else None
    cw = J(HERE, "cw", "cw.json"); nds = J(HERE, "ndstar_427929.json")
    sb = J(CW, "PF2_lxplus", "runs", "lxplus", "pairs", "syst", "syst_blind.json")
    r7 = ub["r427929"]; fit = r7["fit"]; comb = fit["combined"]
    mu = comb["mu"]; lo68, hi68 = comb["mu_68"]
    # --- statystyka: profil × czynnik jackknife po zdarzeniach (własny dla 427929)
    f_jk, src_jk = 1.18, "założone z 428770 (jk_ev200)"
    if jk:
        try:
            f_jk = float(jk["per_mode"]["fixed"]["combined"]["ratio_68"]); src_jk = "jackknife 427929 (200 grup zdarzeń, tryb fixed)"
        except (KeyError, TypeError):
            pass
    sig_up, sig_dn = (hi68 - mu) * f_jk, (mu - lo68) * f_jk
    s_hesse = {c: fit[c]["mu_err_hesse"] for c in CH}
    w = np.array([1 / s_hesse[c] ** 2 for c in CH]); w = w / w.sum()
    # --- wiersze 428770 (syst_blind) po nazwie źródła
    rows428 = {r["desc"]: r for r in sb["rows"]}

    def ch_vals(src):
        d = rows428[src]["dmu"]
        return [float(d.get(c) or 0.0) for c in CH]
    # --- c_w 427929: metoda stosunku (użyta) vs bezpośrednia v2.0 (wagi do 427929 MB)
    cw_used = {c: cw["c_w_427929"][c] for c in CH}; cw_dir = cw["samples"]["427929_mb"]["c_w"]
    neff7 = cw["samples"]["427929_mb"]["n_eff_pairs"]
    # --- SF przeniesienie (HYBSF v3, nomcat; D⁺ z DpDp: SF pary 427929 / 428770; D⁰: tag D* 0,98 ± 0,08 w obu runach → 0)
    sfr_dpdp = 1.038 / 1.159; sf_leg = math.sqrt(sfr_dpdp)
    sf_tr = [0.0, abs(1 / sf_leg - 1), abs(1 / sfr_dpdp - 1)]
    # --- N_D*: stat 427929 ⊕ syst NDS 428770 przeniesiony względnie (do czasu ub427nds)
    nd_stat = nds["fit"]["N_err"] / nds["fit"]["N"]; nd_sys_nds, nd_shape = 605.0 / 4436.39, 336.0 / 4436.39
    nds_lx = os.path.join(HERE, "nds", "ndstar_variants_427929.json")
    nd_src = "syst NDS 428770 przeniesiony względnie (13,6 % ⊕ kształt MC 7,6 %) — CZEKA ub427nds"
    if os.path.exists(nds_lx):
        s = J(nds_lx)["summary"]["rel"]; nd_sys_nds, nd_shape = s["nds_syst"], s["shape_mc"]; nd_src = "ub427nds (warianty NDS na 427929)"
    nd_syst = math.hypot(nd_sys_nds, nd_shape); nd_tot = math.hypot(nd_stat, nd_syst)
    R = []   # (grupa, źródło, [δ per kanał], status, uwaga)

    def add(g, s, v, st, note=""):
        R.append({"group": g, "source": s, "ch": [float(x) for x in v], "status": st, "note": note})
    shapes = ["szerokość piku D0 ±6%", "szerokość piku D⁺ ±3%", "skala masy D0 ±1 MeV", "skala masy D⁺ ±1 MeV", "tło fałszywych nóg exp(poly2)→exp(poly3)",
              "udział swapów K↔π D0 ±20%", "refleksja Ds→KKπ w D⁺ (5%, płaska)", "refleksja Λc→pKπ w D⁺ (1%, płaska)"]
    for s in shapes:
        add("Kształty", s, ch_vals(s), "bez zmian", "zachowawczo; S/B par w 427929 ≈ 4× lepsze (tło przy tych samych D*) → wiersze tła raczej mniejsze")
    add("SF", "SF D0 (score) ±1σ", ch_vals("SF D0 (score) ±1σ"), "bez zmian", "te same SF co nominał (z 428770)")
    add("SF", "SF D⁺ (score) ±1σ", ch_vals("SF D⁺ (score) ±1σ"), "bez zmian", "")
    add("SF", "SF tagu D* ±14%", ch_vals("SF tagu D* ±14%"), "bez zmian", "")
    add("SF", "przeniesienie SF 428770 → 427929", sf_tr, "427929",
        f"HYBSF v3 nomcat: SF pary DpDp 1,038/1,159 = {sfr_dpdp:.3f} (0,6σ od 1) → noga D⁺ {sf_leg:.3f}; D⁰ z tagu D* 0,98 ± 0,08 w obu runach → 0")
    add("Normalizacja", "statystyka MC par (n_eff, baza nominału)", [0.35, 0.27, 0.37], "bez zmian",
        "c_w(427929) = c_w(nominał)·R_X ⇒ dziedziczy szum bazy nominału (VER_budget A: 1/√n_eff)")
    add("Normalizacja", "model wag (w2↔w1)", [abs(cw_dir[c] - 1) for c in CH], "skaluje",
        f"w 427929 wagi ≈ 1 (⟨n_trk⟩ {cw['samples']['427929_mb']['mean_ntrk']:.1f} vs 39,9 w 428770); cały efekt wag |c_w,427929 − 1| jako syst (428770: 0,052/0,189/0,307)")
    add("Normalizacja", "κ_cls (korelacja ε nóg)", [0.061, 0.029, 0.058], "bez zmian", "własność klasyfikatora w MC fd_nd; jednostronnie w dół")
    # 19.11 (decyzja): tune asymetrycznie — dół: MPI ⊕ CR↓ (0,040/0,108 wiersza CR) ⊕ Monash; góra: MPI ⊕ fragmentacja ALICE (1/k − 1, frag/FRAG.md;
    #   zastępuje górną stronę CRqcd −11 % → −29 %, bez podwójnego liczenia) ⊕ Monash↑ (0)
    _mpi, _cr, _mon = [0.040, 0.041, 0.037], [0.112, 0.109, 0.096], [0.020, 0.019, 0.017]
    TUNE_DN = [math.sqrt(_mpi[i] ** 2 + (_cr[i] * 0.040 / 0.108) ** 2 + _mon[i] ** 2) for i in range(3)]
    TUNE_UP = [math.hypot(_mpi[i], 1 / k - 1) for i, k in enumerate((0.647, 0.725, 0.812))]
    add("Normalizacja", "tune Pythii (MPI ⊕ CR ⊕ Monash) — strona dolna μ", TUNE_DN, "19.11 asym.",
        "dół: MPI ⊕ CR↓ ⊕ Monash; górna strona (fragmentacja ALICE, k = 0,65/0,73/0,81) liczona osobno → syst +góra/−dół")
    add("Normalizacja", "wydajność(n_trk)", ch_vals("wydajność(n_trk)"), "bez zmian", "")
    add("Normalizacja", "mapa ε: próbka MC", ch_vals("mapa ε: próbka MC"), "bez zmian", "")
    add("Normalizacja", "dyfrakcja ×1,5", ch_vals("dyfrakcja ×1,5"), "bez zmian", "w MBTS domieszka dyfrakcji większa, ale wariant ×1,5 dał 0,000 (D przy pT > 3 GeV z ND)")
    add("Normalizacja", "trigger MBTS (427929)", [0.004] * 3, "427929", "PS = 1 w 100 % LB (MULT G1a), plateau sptrk 0,996–0,998 (G1b); D* i pary ten sam trigger → znosi się")
    add("Normalizacja", "pile-up (nogi z dwóch zderzeń)", [0.001] * 3, "427929", "npu > 0 w 0,65 % zdarzeń; ~0,0065·r_D²/r_par ≈ 0,1 %")
    add("κ, kategorie", "κ: plan B ↔ jedno κ na kanał", ch_vals("κ: plan B ↔ jedno κ na kanał"), "bez zmian", "")
    add("κ, kategorie", "wydajności kategorii score ±5% (niezależnie)", ch_vals("wydajności kategorii score ±5% (niezależnie)"), "bez zmian", "")
    for r in R:
        r["common"] = float(np.dot(w, r["ch"]))
    groups = {}
    for r in R:
        groups.setdefault(r["group"], []).append(r)
    gsum = {g: {"ch": [math.sqrt(sum(r["ch"][i] ** 2 for r in rs)) for i in range(3)], "common": math.sqrt(sum(r["common"] ** 2 for r in rs))} for g, rs in groups.items()}
    tot = {"ch": [math.sqrt(sum(v["ch"][i] ** 2 for v in gsum.values())) for i in range(3)], "common": math.sqrt(sum(v["common"] ** 2 for v in gsum.values()))}
    # strona górna: ten sam budżet z wierszem tune zastąpionym TUNE_UP
    _tr = [r for r in R if r["source"].startswith("tune Pythii")][0]
    _up_common = float(np.dot(w, TUNE_UP))
    _norm_up = math.sqrt(gsum["Normalizacja"]["common"] ** 2 - _tr["common"] ** 2 + _up_common ** 2)
    tot_up = math.sqrt(tot["common"] ** 2 - gsum["Normalizacja"]["common"] ** 2 + _norm_up ** 2)
    tot_dn = tot["common"]
    syst_abs = tot["common"] * mu; nd_abs = nd_tot * mu
    # --- Z i porównania
    Z = comb["Z"]; Z_jk = Z / f_jk
    s_up_tot_at1 = math.sqrt(sig_up ** 2 + tot_up ** 2 + nd_tot ** 2)        # test μ = 1: skale liczone przy μ = 1 (zachowawczo)
    pull_pythia = (1 - mu) / s_up_tot_at1
    m8, s8 = 1.495, 0.47                                                              # 428770: μ, stat (jk) — UNB_WYNIK §1b
    sh = {"tune": 0.117, "SF": 0.132, "mc_stat": 0.308, "kcls": 0.042, "shapes": 0.071, "kappa": 0.094}       # wspólne syst (względne) — znoszą się w różnicy
    run_spec8 = math.hypot(0.163, 0.151)                                              # 428770: model wag ⊕ stat N_D*
    run_spec7 = math.sqrt(dict((r["source"], r["common"]) for r in R)["model wag (w2↔w1)"] ** 2 + nd_stat ** 2
                          + dict((r["source"], r["common"]) for r in R)["przeniesienie SF 428770 → 427929"] ** 2)
    sd = math.sqrt(s8 ** 2 + ((hi68 - lo68) / 2 * f_jk) ** 2 + (m8 * run_spec8) ** 2 + (mu * run_spec7) ** 2)
    diff = mu - m8
    out = {"mu": mu, "stat": [sig_up, sig_dn], "f_jk": f_jk, "f_jk_src": src_jk, "w": w.tolist(), "rows": R, "groups": gsum, "total": tot,
           "syst_abs": syst_abs, "syst_rel_up": tot_up, "syst_rel_dn": tot_dn, "tune_up": TUNE_UP, "tune_dn": TUNE_DN, "ndstar": {"stat": nd_stat, "syst_nds": nd_sys_nds, "shape_mc": nd_shape, "total": nd_tot, "abs": nd_abs, "src": nd_src,
                                            "N": nds["fit"]["N"], "N_err": nds["fit"]["N_err"], "scale": nds["fit"]["scale"]},
           "Z": Z, "Z_jk": Z_jk, "pull_vs_pythia": pull_pythia, "diff_vs_428770": diff, "diff_sigma": diff / sd, "sigma_diff": sd,
           "cw_ratio_vs_direct": {c: cw_dir[c] / cw_used[c] for c in CH}, "neff_427929": neff7}
    json.dump(out, open(os.path.join(HERE, "budzet_427929.json"), "w"), indent=1, default=float)
    f = lambda x: f"{x:.3f}"  # noqa: E731
    md = ["# Budżet μ — run 427929 (czysty minimum bias) jako pomiar główny", "",
          "Wariant na wypadek potwierdzenia obciążenia 428770 triggerem J12 (UBTRIG). Wejścia: UB427 (`UB427.md`, 18.11), jackknife 427929 (`jk_ev200/ubjk.md`), "
          "budżet 428770 (`pairs/syst/syst_blind.json`, `VER_budget/RAPORT.md`, `TUNE_pythia/budget_with_tune.md`). Skrypt `budzet_427929.py`.", "",
          f"**μ(427929) = {mu:.2f} +{sig_up:.2f}/−{sig_dn:.2f} (stat) +{tot_up * mu:.2f}/−{tot_dn * mu:.2f} (syst) ± {nd_abs:.2f} (N_D*)**; Z = {Z:.2f}σ (profil), ≈ {Z_jk:.2f}σ z czynnikiem jackknife {f_jk:.2f}.", "",
          f"- stat: profil 68 % [{lo68:.3f}; {hi68:.3f}] × {f_jk:.2f} ({src_jk}).",
          f"- syst (bez N_D*): wspólne względne +{tot_up:.3f} / −{tot_dn:.3f} × μ (asymetria: tune ↑ = fragmentacja ALICE, decyzja z 19.11).",
          f"- N_D* = {nds['fit']['N']:.0f} ± {nds['fit']['N_err']:.0f} (stat {nd_stat:.1%}) ⊕ syst {nd_syst:.1%} [{nd_src}] → {nd_tot:.1%}.",
          f"- Zgodność z Pythią (μ = 1): {pull_pythia:+.2f}σ w stronę deficytu (skale syst przy μ = 1, zachowawczo).",
          f"- Różnica z 428770 (μ = 1,50): {diff:+.2f} ± {sd:.2f} → {diff / sd:+.2f}σ (wspólne syst MC — tune, SF, mc_stat, κ, kształty — znoszą się; "
          "zostają stat obu runów, stat N_D* obu runów, model wag 428770 i wiersze specyficzne 427929).", "",
          "## Wiersze (Δμ/μ)", "", f"Wagi „wspólne” (1/σ²_Hesse 427929): {', '.join(f'{c} {x:.3f}' for c, x in zip(CH, w))}.", "",
          "| grupa | źródło | D0D0 | D⁺D0 | D⁺D⁺ | wspólne | status | uwaga |", "|---|---|---|---|---|---|---|---|"]
    for r in R:
        md.append(f"| {r['group']} | {r['source']} | " + " | ".join(f(x) for x in r["ch"]) + f" | {f(r['common'])} | {r['status']} | {r['note']} |")
    for g, v in gsum.items():
        md.append(f"| **suma** | {g} | " + " | ".join(f(x) for x in v["ch"]) + f" | **{f(v['common'])}** | | |")
    md.append(f"| **suma** | wszystko (bez N_D*) | " + " | ".join(f(x) for x in tot["ch"]) + f" | **{f(tot['common'])}** | | |")
    md += ["", "## Uwagi do decyzji", "",
           f"1. **N_D* w 427929 jest najsłabszym ogniwem.** Fit Δm daje skalę szerokości piku {nds['fit']['scale']:.2f} × MC (428770: 1,02, ten sam detektor); "
           f"per bin pT: 3–5 GeV {nds['per_bin']['3-5']['N']:.0f} ± {nds['per_bin']['3-5']['N_err']:.0f} (skala {nds['per_bin']['3-5']['scale']:.2f}), "
           f"5–8 {nds['per_bin']['5-8']['N']:.0f} ± {nds['per_bin']['5-8']['N_err']:.0f} ({nds['per_bin']['5-8']['scale']:.2f}), "
           f"8–30 {nds['per_bin']['8-30']['N']:.0f} ± {nds['per_bin']['8-30']['N_err']:.0f} ({nds['per_bin']['8-30']['scale']:.2f}). "
           "Szeroki pik przy dużym tle (3–5 GeV) sugeruje tło wciągnięte w sygnał → N_D* zawyżone → μ zaniżone. Etap `ub427nds` (paczka `Q_lxplus/ub427_nds`) liczy "
           "skan wariantów NDS i fity per bin ze stałą skalą σ; do tego czasu syst N_D* przeniesiony z 428770 jest prawdopodobnie zaniżony.",
           "2. **c_w: metoda stosunku vs bezpośrednio.** Bezpośrednie c_w(427929) z MC v2.0 (wagi ≈ 1) / użyte c_w = "
           + ", ".join(f"{c} {out['cw_ratio_vs_direct'][c]:.3f}" for c in CH)
           + " — dokładnie to samo, co nieodtworzenie c_w nominału 428770 na v2.0 (`repro_428770_vs_nominal`); to szum bazy normalizacji wspólny dla obu runów, "
           "pokryty wierszem n_eff — nie dodaję drugi raz. Wariant bezpośredni dałby μ ≈ 0,45 (D⁺D0 ×0,72). Normalizacja z prawdy (`gen_norm`, PASS 19.11) ma zastąpić bazę w obu runach.",
           "3. Wiersze kształtów/κ przeniesione zachowawczo (Asimov 428770); w 427929 S/B par ≈ 4× lepsze, więc realnie mniejsze.",
           f"4. MC pary przy wagach 427929: n_eff {', '.join(f'{c} {neff7[c]:.1f}' for c in CH)} (428770: 16,0 / 14,6 / 7,5 na v2.0) — po przejściu na normalizację z prawdy "
           "statystyka MC w 427929 będzie wyraźnie lepsza niż w 428770.", ""]
    open(os.path.join(HERE, "BUDZET_427929.md"), "w").write("\n".join(md))
    print("\n".join(md[:14]))


if __name__ == "__main__":
    main()
