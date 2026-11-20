#!/usr/bin/env python3
"""TRIG_A/syst_trig.py (19.11.2026, projekt A) — wiersz systematyki „trigger / selekcja zdarzeń” wg PRErejestracji TRIG_A/SYST_TRIG.md.

Czyta ubtrig.json (UBTRIG, μ per klasa triggera w 428770) i ub427.json (UB427, μ w 427929 vs 428770), opcjonalnie ubjk.json z UB427JK
(czynnik f_jk = σ_jk/σ_profil dla 427929). Liczy d₁ (między runami), d₂ (klasy w runie), d, Z_trig, d_emu (TRIGCAL), status i wiersz.
Bez ważnego UBTRIG → status „niekompletne/BRAK”, nic do budżetu.
  syst_trig.py --ubtrig ubtrig.json --ub427 ub427.json [--ub427jk ubjk.json] [--out DIR]
  syst_trig.py --selftest [--ub427 ub427.json --ubtrig-bad ubtrig_nan.json]"""
import argparse
import json
import math
import os
import sys

MU = 1.50                      # zapis bieżący (N_D* = 4436)
MU_428_REF = 1.50              # odtworzenie nominału w UB427 (tol 0,01)
F_JK_DEFAULT = 1.18            # σ_jk/σ_profil w 428770 (jk_ev200) — założone dla 427929 do czasu UB427JK
F_EV_CLASS = 1.18              # czynnik zdarzeń dla klas UBTRIG (podzbiory 428770)
KAPPA_EMU, KAPPA_LO, KAPPA_HI = 1.24, 1.146, 1.388   # TRIGCAL, kalibracje χ²/ndf < 2 (em 27,5/4; all 40/9), k = 0..2
SIG_WAG = 0.24                 # „model wag” w μ (VER_budget, zapis 16.11)
CLASSES = ["J12only", "MBTS_j20", "MBTS_noj20", "other"]
CLEAN = ["MBTS_noj20", "other"]


def _fin(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def wmean(vals, errs):
    w = [1 / e ** 2 for e in errs]
    m = sum(wi * v for wi, v in zip(w, vals)) / sum(w)
    return m, 1 / math.sqrt(sum(w))


def d_emu(r, k):
    x = r * (k - 1)
    return x / (1 + x)


def check_ubtrig(U):
    """→ (ok, powody)."""
    why = []
    df = U.get("dstar_fraction", {})
    rs = [df.get(c, {}).get("r") for c in CLASSES]
    if not all(_fin(r) for r in rs) or not all(df.get(c, {}).get("determined") for c in CLASSES):
        why.append("r_c nieokreślone/NaN: " + str({c: df.get(c, {}).get("r") for c in CLASSES}))
    elif not 0.95 <= sum(rs) <= 1.05:
        why.append(f"Σ r_c = {sum(rs):.3f} ∉ [0,95; 1,05]")
    if not U.get("nominal_ok"):
        why.append("nominał nieodtworzony")
    if not all((U.get("sum_equals_pairs_full") or {}).values()) or not U.get("sum_equals_pairs_full"):
        why.append("Σ klas ≠ pairs_full")
    if U.get("legs_unmatched", 1) != 0:
        why.append(f"nogi bez klasy: {U.get('legs_unmatched')}")
    F = U.get("fits", {})
    ok_c = [c for c in CLASSES if _fin((F.get(c, {}).get("combined") or {}).get("mu")) and _fin((F.get(c, {}).get("combined") or {}).get("sigma_68"))]
    if "J12only" not in ok_c or len(ok_c) < 2:
        why.append(f"brak fitu wspólnego μ dla J12only i ≥1 innej klasy (dopasowane: {ok_c})")
    if not _fin((F.get("nominal", {}).get("combined") or {}).get("mu")):
        why.append("brak μ nominalnego")
    return (not why), why


def check_ub427(V):
    why = []
    m8 = V.get("r428770", {}).get("fit", {}).get("combined", {}).get("mu")
    if not _fin(m8) or abs(m8 - MU_428_REF) > 0.01:
        why.append(f"μ₄₂₈₇₇₀ = {m8} ≠ {MU_428_REF}")
    for k in ("sigma_427929_side", "sigma_428770_side"):
        if not _fin(V.get(k)):
            why.append(f"brak {k}")
    if not _fin(V.get("r427929", {}).get("fit", {}).get("combined", {}).get("mu")):
        why.append("brak μ₄₂₇₉₂₉")
    return (not why), why


def compute_d1(V, f_jk):
    m7 = V["r427929"]["fit"]["combined"]["mu"]; m8 = V["r428770"]["fit"]["combined"]["mu"]
    e7, e8 = V.get("ndstar_rel_stat", {}).get("427929", 0.0), V.get("ndstar_rel_stat", {}).get("428770", 0.0)
    s = math.sqrt((f_jk * V["sigma_427929_side"]) ** 2 + V["sigma_428770_side"] ** 2 + (m7 * e7) ** 2 + (m8 * e8) ** 2)
    return {"mu_427929": m7, "mu_428770": m8, "f_jk": f_jk, "d": (m8 - m7) / m8, "sigma": s / m8, "Z_delta": (m8 - m7) / s}


def class_mu(U, c):
    cm = U["fits"][c]["combined"]; r_rel = U["fits"][c].get("r_rel_err") or 0.0
    return cm["mu"], math.hypot(F_EV_CLASS * cm["sigma_68"], cm["mu"] * r_rel)


def compute_d2(U, rest_classes):
    F = U["fits"]
    have = [c for c in rest_classes if _fin((F.get(c, {}).get("combined") or {}).get("mu"))]
    if not have:
        return None
    mj, sj = class_mu(U, "J12only")
    mr, sr = wmean(*zip(*[class_mu(U, c) for c in have]))
    mn = F["nominal"]["combined"]["mu"]; r = U["dstar_fraction"]["J12only"]["r"]
    return {"mu_J12": mj, "err_J12": sj, "mu_rest": mr, "err_rest": sr, "rest": have, "mu_nom": mn, "r_J12": r,
            "d": r * (mj - mr) / mn, "sigma": r * math.hypot(sj, sr) / mn}


def decide(d, s, de, flag_incons):
    Z = d / s
    zemu = (d - de) / s
    row = {}
    if abs(Z) < 2 and abs(zemu) <= 2:
        st = "brak istotnego efektu"; row = {"trigger": -MU * de, "jednostronny": True, "wag": SIG_WAG, "zastepuje_wag": False}
    elif abs(Z) < 2:
        st = "emulacja odrzucona"
        row = {"trigger": (-MU * max(abs(d), s)) if d >= 0 else MU * s, "jednostronny": d >= 0, "wag": SIG_WAG, "zastepuje_wag": False}
    elif Z < 3 and Z >= 2:
        st = "wskazówka"; row = {"selekcja_zdarzen": -max(MU * d, math.hypot(SIG_WAG, MU * de)), "jednostronny": True, "zastepuje_wag": True}
    elif Z >= 3:
        st = "efekt (zmiana nominału tylko osobną decyzją)"
        row = {"selekcja_zdarzen": -max(MU * d, math.hypot(SIG_WAG, MU * de)), "jednostronny": True, "zastepuje_wag": True,
               "wariant_mu_corr": [MU * (1 - d), MU * s]}
    else:
        st = "nadwyżka w czystych zdarzeniach"; row = {"trigger": MU * abs(d), "jednostronny": False, "wag": SIG_WAG, "zastepuje_wag": False}
    if flag_incons:
        st += "; d₁ i d₂ niezgodne"
    return {"Z_trig": Z, "Z_vs_emu": zemu, "status": st, "row": row}


def run(U, V, f_jk=None, f_jk_how=None):
    rep = {"prereg": "TRIG_A/SYST_TRIG.md (19.11.2026)", "mu": MU}
    okv, whyv = check_ub427(V) if V is not None else (False, ["brak ub427.json"])
    oku, whyu = check_ubtrig(U) if U is not None else (False, ["brak ubtrig.json"])
    rep["ub427_ok"], rep["ub427_why"], rep["ubtrig_ok"], rep["ubtrig_why"] = okv, whyv, oku, whyu
    fj = f_jk if f_jk is not None else F_JK_DEFAULT
    rep["f_jk"], rep["f_jk_how"] = fj, f_jk_how or "założone (428770 jk_ev200)"
    if okv:
        rep["d1"] = compute_d1(V, fj)
    r = U["dstar_fraction"]["J12only"]["r"] if oku else 0.769
    rep["r_J12"] = r
    rep["d_emu"] = {"kappa": KAPPA_EMU, "d": d_emu(r, KAPPA_EMU), "d_lo": d_emu(r, KAPPA_LO), "d_hi": d_emu(r, KAPPA_HI)}
    if not oku:
        rep["status"] = "BRAK — UBTRIG nieważny; decyzja odłożona (niekompletne)" + ("" if okv else "; UB427 też nieważny")
        rep["row"] = None
        return rep
    rep["d2"] = compute_d2(U, [c for c in CLASSES if c != "J12only"])
    rep["d2_clean"] = compute_d2(U, CLEAN)
    if not okv:
        rep["status"] = "BRAK — UB427 nieważny; decyzja odłożona (niekompletne)"; rep["row"] = None
        return rep
    d1, d2 = rep["d1"], rep["d2"]
    d, s = wmean([d1["d"], d2["d"]], [d1["sigma"], d2["sigma"]])
    incons = abs(d1["d"] - d2["d"]) / math.hypot(d1["sigma"], d2["sigma"]) > 2
    if incons:
        s = math.hypot(s, abs(d1["d"] - d2["d"]) / 2)
    rep["d"], rep["sigma_d"], rep["d1_d2_inconsistent"] = d, s, incons
    rep.update(decide(d, s, rep["d_emu"]["d"], incons))
    return rep


def to_md(rep):
    f = lambda x: "–" if x is None else f"{x:+.3f}"  # noqa: E731
    L = ["# SYST_TRIG — wiersz „trigger / selekcja zdarzeń” (reguły: TRIG_A/SYST_TRIG.md, prerejestracja 19.11)", "",
         f"Status: **{rep['status']}**", ""]
    if not rep["ubtrig_ok"]:
        L.append("UBTRIG nieważny: " + "; ".join(rep["ubtrig_why"]))
    if not rep["ub427_ok"]:
        L.append("UB427 nieważny: " + "; ".join(rep["ub427_why"]))
    L += ["", "| pomiar | μ (wejście) | d = (μ_nom − μ_czyste)/μ_nom | σ_d | Z |", "|---|---|---|---|---|"]
    if "d1" in rep:
        x = rep["d1"]
        L.append(f"| d₁ UB427 (427929 vs 428770; f_jk = {x['f_jk']:.2f}, {rep['f_jk_how']}) | {x['mu_427929']:.2f} vs {x['mu_428770']:.2f} | {f(x['d'])} | {x['sigma']:.3f} | {x['d'] / x['sigma']:+.2f} |")
    for k, nm in (("d2", "d₂ UBTRIG (J12only vs reszta)"), ("d2_clean", "d₂ kontrola (J12only vs MBTS_noj20 ⊕ other)")):
        x = rep.get(k)
        if x:
            L.append(f"| {nm}; r_J12 = {x['r_J12']:.3f} | {x['mu_J12']:.2f} ± {x['err_J12']:.2f} vs {x['mu_rest']:.2f} ± {x['err_rest']:.2f} (μ_nom {x['mu_nom']:.2f}) | {f(x['d'])} | {x['sigma']:.3f} | {x['d'] / x['sigma']:+.2f} |")
    if "d" in rep:
        L.append(f"| **d (średnia)** | | **{f(rep['d'])}** | {rep['sigma_d']:.3f} | **{rep['Z_trig']:+.2f}** |")
    e = rep["d_emu"]
    L += ["", f"Emulacja J12 (TRIGCAL): κ = {e['kappa']} → d_emu = {e['d']:.3f} (zakres {e['d_lo']:.3f}–{e['d_hi']:.3f}) przy r_J12 = {rep['r_J12']:.3f}"
          + (f"; zgodność danych z emulacją: {rep['Z_vs_emu']:+.2f}σ" if "Z_vs_emu" in rep else "")]
    if rep.get("row"):
        L += ["", f"Wiersz (μ = {rep['mu']}): `{json.dumps(rep['row'], ensure_ascii=False)}`"]
    return "\n".join(L) + "\n"


def selftest(ub427_path=None, bad_path=None):
    def fake_u(mj, sj, mr, sr, r=0.769):
        fits = {"nominal": {"combined": {"mu": 1.189, "sigma_68": 0.317}}}
        for c, (m, s) in zip(CLASSES, ((mj, sj), (mr, sr), (mr, 2 * sr), (mr, 1.5 * sr))):
            fits[c] = {"combined": {"mu": m, "sigma_68": s}, "r_rel_err": 0.03}
        return {"dstar_fraction": {c: {"r": x, "determined": True} for c, x in zip(CLASSES, (r, 0.135, 0.010, 0.084))},
                "nominal_ok": True, "sum_equals_pairs_full": {"D0D0": True, "DpD0": True, "DpDp": True}, "legs_unmatched": 0, "fits": fits}
    V0 = {"r427929": {"fit": {"combined": {"mu": 1.40}}}, "r428770": {"fit": {"combined": {"mu": 1.50}}}, "sigma_427929_side": 0.30,
          "sigma_428770_side": 0.46, "ndstar_rel_stat": {"427929": 0.2, "428770": 0.15}}
    res, ok = [], True
    cases = [("brak efektu, zgodne z emulacją", fake_u(1.25, 0.35, 1.05, 0.5), V0, "brak istotnego efektu"),
             ("silny efekt", fake_u(1.60, 0.05, 0.60, 0.05), {**V0, "r427929": {"fit": {"combined": {"mu": 0.55}}}, "sigma_427929_side": 0.05,
                                                             "sigma_428770_side": 0.05, "ndstar_rel_stat": {"427929": 0.01, "428770": 0.01}}, "efekt"),
             ("odwrotny znak", fake_u(0.80, 0.05, 1.90, 0.06), {**V0, "r427929": {"fit": {"combined": {"mu": 2.6}}}, "sigma_427929_side": 0.1,
                                                              "sigma_428770_side": 0.05, "ndstar_rel_stat": {"427929": 0.01, "428770": 0.01}}, "nadwyżka"),
             ("NaN r_c", {**fake_u(1.2, 0.3, 1.1, 0.4), "dstar_fraction": {c: {"r": None, "determined": False} for c in CLASSES}}, V0, "BRAK")]
    if ub427_path:
        V = json.load(open(ub427_path))
        cases.append(("prawdziwy UB427 + wymyślony UBTRIG (efekt jak w emulacji)", fake_u(1.30, 0.40, 1.00, 0.55), V, None))
        if bad_path:
            cases.append(("prawdziwy UB427 + błędny UBTRIG 18.11 (NaN)", json.load(open(bad_path)), V, "BRAK"))
    for nm, U, V, want in cases:
        r = run(U, V)
        good = want is None or r["status"].startswith(want)
        ok &= good
        res.append((nm, r["status"], r.get("d"), r.get("sigma_d"), r.get("Z_trig"), r.get("row"), good))
        print(f"[{'OK' if good else 'BŁĄD'}] {nm}: {r['status']}; d = {r.get('d')}, σ = {r.get('sigma_d')}, Z = {r.get('Z_trig')}; wiersz {r.get('row')}")
        if "d1" in r:
            print(f"       d₁ = {r['d1']['d']:+.3f} ± {r['d1']['sigma']:.3f}" + (f", d₂ = {r['d2']['d']:+.3f} ± {r['d2']['sigma']:.3f}" if r.get("d2") else ""))
    print("SELFTEST", "OK" if ok else "BŁĄD")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ubtrig"); ap.add_argument("--ub427"); ap.add_argument("--ub427jk"); ap.add_argument("--out")
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--ubtrig-bad")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest(a.ub427, a.ubtrig_bad)
    U = json.load(open(a.ubtrig)) if a.ubtrig and os.path.exists(a.ubtrig) else None
    V = json.load(open(a.ub427)) if a.ub427 and os.path.exists(a.ub427) else None
    fj, how = None, None
    if a.ub427jk and os.path.exists(a.ub427jk):
        J = json.load(open(a.ub427jk))
        try:
            pm = J["per_mode"][J.get("headline_mode", "fixed")]
            c = pm["combined"]
            fj = float(c["ratio_68"]); how = f"zmierzone (UB427JK, ±{100 * c.get('rel_precision_sigma_jk', 0):.0f} %)"
        except Exception as ex:  # noqa: BLE001
            print("UWAGA: nie odczytano f_jk z", a.ub427jk, ex, file=sys.stderr)
    rep = run(U, V, fj, how)
    md = to_md(rep)
    print(md)
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        json.dump(rep, open(os.path.join(a.out, "SYST_TRIG.json"), "w"), indent=1, default=float)
        open(os.path.join(a.out, "SYST_TRIG_wynik.md"), "w").write(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
