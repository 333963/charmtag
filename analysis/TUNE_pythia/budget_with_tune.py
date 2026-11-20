"""Budżet systematyk μ par D–D z tune Pythii (MPI, CR z CRqcd, Monash) w grupie „Normalizacja” (decyzja z 12.11.2026).

Nie zmienia nominalnego fitu, Z ani μ centralnego — tylko przelicza budżet systematyk.
Reguły odtworzone 1:1 z charmtag.syst (kopia Q_lxplus/syst/charmtag_syst, identyczna logicznie z pf_pairs/src i AA_merge):
- norm.py collect() l. 97–110: δ_grupy = max po wariantach |[par/D*]_war / [par/D*]_A3 − 1|;
  l. 119: dmu_rel = |δ| (1. rząd; dmu_rel() = 1/(1+δ) − 1 tylko do asymetrii, l. 120).
- report.py build() l. 40–41, 76–77: „wspólne” = Σ_ch w_ch·Δμ_ch / Σ w_ch, w = 1/σ²_stat (brak w kanale = 0);
  l. 86–93: sumy grup = kwadratura wierszy (bez pending i bez kappa_partner*), suma „wspólne” = kwadratura wierszy „wspólne”;
  l. 112–113: „wszystko” = kwadratura sum grup.
Uruchom z charm_work:  ../.venv_h5/bin/python TUNE_pythia/budget_with_tune.py
"""
import json
import math
import os
import sys

sys.path.insert(0, "Q_lxplus/syst")
from charmtag_syst import norm as N  # noqa: E402

SYST = "PF2_lxplus/runs/lxplus/pairs/syst/syst_blind.json"
TUNE = "TUNE_pythia/out/tune_syst.json"
TRUTH = "../runs/lxplus/truth"
OUTJ = "TUNE_pythia/out/budget_with_tune.json"
OUTM = "TUNE_pythia/budget_with_tune.md"
CH = N.CH
TUNE_ROWS = ("tune_MPI", "tune_CR", "tune_Monash")
# --- baza normalizacji (decyzja z 12.11.2026) ---
# Źródło: charm_work/NDS_ndstar/rerun_2026-11-12/rerun.md (TL;DR) — powtórka skanu NDS po poprawce NLL w dstar_refit.fit_single:
#   N_D* = 4436, δN/N = 21,9 % (MINOS sym.), μ_wsp. = 1,495 (≈1,50) +0,41/−0,39 (stat, profil) ± 0,33 (N_D*).
#   Od 13.11: stat × 1,18 (jackknife po 200 grupach zdarzeń, unblind/jk_ev200/ubjk.md) → +0,48/−0,46.
MU_COMMON = 1.495         # μ_wsp. centralne (rerun 12.11)
NDSTAR_REL = 0.219        # δN_D*/N_D* (rerun 12.11, MINOS sym. 680/4436), podawane osobno w wyniku
NDSTAR_N = 4436           # N_D* (rerun 12.11) — tylko do opisu
# Poprzednia baza (UNB_WYNIK §1b, aktualizacja 09.11): N_D* = 4443, 20,6 %, μ = 1,49 ± 0,31 (N_D*) — zachowana jako „previous_base”.
PREV_MU_COMMON = 1.49
PREV_NDSTAR_REL = 0.206
PREV_NDSTAR_N = 4443
TOL = 0.002

S = json.load(open(SYST))
T = json.load(open(TUNE))
sig = S["sigma_stat"]
w = {ch: 1.0 / sig[ch] ** 2 for ch in CH}
W = sum(w.values())


def common(vals):                                   # report.py l. 76–77
    if not any(vals.get(ch) is not None for ch in CH):
        return None
    return sum(w[ch] * (vals.get(ch) or 0.0) for ch in CH) / W


def counted(r):                                     # report.py l. 86
    return not r["pending"] and not r["source"].startswith("kappa_partner")


def sums(rows):                                     # report.py l. 86–93, 110–113
    groups = []
    for r in rows:
        if r["group"] not in groups:
            groups.append(r["group"])
    tot = {g: {ch: 0.0 for ch in CH} for g in groups}
    totc = {g: 0.0 for g in groups}
    for r in rows:
        if not counted(r):
            continue
        for ch in CH:
            if r["dmu"].get(ch) is not None:
                tot[r["group"]][ch] += r["dmu"][ch] ** 2
        if r["dmu_common"] is not None:
            totc[r["group"]] += r["dmu_common"] ** 2
    out = {g: {**{ch: math.sqrt(tot[g][ch]) for ch in CH}, "common": math.sqrt(totc[g])} for g in groups}
    out["wszystko"] = {**{ch: math.sqrt(sum(tot[g][ch] for g in groups)) for ch in CH},
                       "common": math.sqrt(sum(totc.values()))}
    return out


# ---- 1. odtworzenie starych wartości „wspólne” i sum
chk = {"rows_common_maxdiff": 0.0, "sums": {}}
for r in S["rows"]:
    if r["dmu_common"] is not None:
        chk["rows_common_maxdiff"] = max(chk["rows_common_maxdiff"], abs(common(r["dmu"]) - r["dmu_common"]))
old = sums(S["rows"])
ref = {g: {**v, "common": None} for g, v in S["totals"].items()}
ref["wszystko"] = dict(S["total_all"])
ref["wszystko"]["common"] = S["total_common"]
ref_common_md = {"Kształty": 0.071, "SF": 0.132, "Normalizacja": 0.216, "κ, kategorie": 0.094}   # syst_blind.md (json nie zapisuje sum „wspólne” grup)
maxd = 0.0
for g in old:
    d = {}
    for ch in CH:
        d[ch] = old[g][ch] - ref[g][ch]
    rc = ref[g]["common"] if ref[g]["common"] is not None else ref_common_md[g]
    d["common"] = old[g]["common"] - rc
    chk["sums"][g] = d
    maxd = max(maxd, max(abs(x) for x in d.values()))
chk["sums_maxdiff"] = maxd
if maxd > TOL or chk["rows_common_maxdiff"] > 1e-9:
    print(json.dumps(chk, indent=1, ensure_ascii=False))
    sys.exit("STOP: nie odtworzono starych sum co do 0,002 — reguła sumowania inna niż w report.py")

# ---- 2. wiersze tune: ścieżka kodu norm.collect na truth_norm.json (lokalna kopia) + kontrola z tune_syst.json
coll = N.collect({}, TRUTH) if os.path.exists(os.path.join(TRUTH, "combined_A3", "truth_norm.json")) else None
tune = {}
for g in TUNE_ROWS:
    tj = T["groups"][g]
    if coll and not coll[g].get("pending"):
        dm = dict(coll[g]["dmu_rel"])
        src = "norm.collect(truth_norm.json)"
    else:
        dm = {ch: abs(tj[ch]) for ch in CH}
        src = "tune_syst.json (groups)"
    diff = max(abs(dm[ch] - abs(tj[ch])) for ch in CH)
    # asymetria (informacyjnie): Δμ/μ = 1/(1+δ) − 1 per wariant, skrajne w górę i w dół
    vs = [v for v in T["variants"].values() if v["group"] == g]
    up = {ch: max([N.dmu_rel(v["delta"][ch]) for v in vs] + [0.0]) for ch in CH}
    dn = {ch: min([N.dmu_rel(v["delta"][ch]) for v in vs] + [0.0]) for ch in CH}
    tune[g] = {"dmu": dm, "dmu_common": common(dm), "source_of_numbers": src, "maxdiff_vs_tune_syst_json": diff,
               "asym_up": {**up, "common": common(up)}, "asym_down": {**dn, "common": common(dn)},
               "variants": [k for k, v in T["variants"].items() if v["group"] == g]}

# ---- 3. nowe wiersze i sumy
rows_new = []
for r in S["rows"]:
    r = dict(r)
    if r["source"] in tune:
        r["dmu"] = tune[r["source"]]["dmu"]
        r["dmu_common"] = tune[r["source"]]["dmu_common"]
        r["pending"] = False
    rows_new.append(r)
new = sums(rows_new)
tune_only = {**{ch: math.sqrt(sum(tune[g]["dmu"][ch] ** 2 for g in TUNE_ROWS)) for ch in CH},
             "common": math.sqrt(sum(tune[g]["dmu_common"] ** 2 for g in TUNE_ROWS))}

# ---- 4. bez wiersza N_D* (ndstar_fit) — kwadratowo odjęty
nd = next(r for r in S["rows"] if r["source"] == "ndstar_fit")
ndv = {**nd["dmu"], "common": nd["dmu_common"]}
K = list(CH) + ["common"]
qsub = lambda a, b: math.sqrt(max(a * a - b * b, 0.0))  # noqa: E731
res = {}
for lab, s in (("przed", old), ("po", new)):
    res[lab] = {"Normalizacja": {k: s["Normalizacja"][k] for k in K},
                "wszystko": {k: s["wszystko"][k] for k in K},
                "Normalizacja_bez_NDS": {k: qsub(s["Normalizacja"][k], ndv[k]) for k in K},
                "wszystko_bez_NDS": {k: qsub(s["wszystko"][k], ndv[k]) for k in K}}

# ---- 5. przeskalowanie do μ_wsp. (MU_COMMON; wariant historyczny PREV_MU_COMMON)
def scaled(mu, nds_rel, n_dstar):
    s = {"mu_common": mu, "ndstar_rel": nds_rel, "ndstar_N": n_dstar,
         "tune_abs": tune_only["common"] * mu,
         "syst_bez_NDS_z_tune_abs": res["po"]["wszystko_bez_NDS"]["common"] * mu,
         "syst_bez_NDS_bez_tune_abs": res["przed"]["wszystko_bez_NDS"]["common"] * mu,
         "syst_z_NDS_z_tune_abs": res["po"]["wszystko"]["common"] * mu,
         "NDS_abs": nds_rel * mu,
         "tune_CR_up_abs": tune["tune_CR"]["asym_up"]["common"] * mu,
         "tune_CR_down_abs": tune["tune_CR"]["asym_down"]["common"] * mu}
    # wariant asymetryczny (informacyjnie): syst bez N_D* i bez tune ⊕ tune ↑ / ↓ (per grupa skrajne 1/(1+δ)−1)
    for _d in ("up", "down"):
        _t = math.sqrt(sum(tune[g][f"asym_{_d}"]["common"] ** 2 for g in TUNE_ROWS))
        s[f"tune_asym_{_d}_rel"] = _t
        s[f"syst_bez_NDS_asym_{_d}_abs"] = math.hypot(qsub(old["wszystko"]["common"], ndv["common"]), _t) * mu
    return s


sc = scaled(MU_COMMON, NDSTAR_REL, NDSTAR_N)
sc_prev = scaled(PREV_MU_COMMON, PREV_NDSTAR_REL, PREV_NDSTAR_N)
sc_prev["note"] = "zapis historyczny z 12.11 przed zmianą bazy N_D* (UNB_WYNIK §1b, aktualizacja 09.11)"

out = {"rules": {"dmu_norm": "Δμ/μ = |δ| (norm.py l. 119); δ_grupy = max |a/b−1| po wariantach (l. 97–110)",
                 "common": "Σ w·Δμ / Σ w, w = 1/σ²_stat, brak = 0 (report.py l. 40–41, 76–77)",
                 "sums": "kwadratura, bez pending i kappa_partner* (report.py l. 86–93, 112–113)"},
       "check_old": chk, "sigma_stat": sig, "tune_rows": tune, "tune_quadrature": tune_only,
       "ndstar_fit_row": ndv, "budget": res,
       "base_source": "NDS_ndstar/rerun_2026-11-12/rerun.md (TL;DR), decyzja z 12.11.2026",
       "scaled_common": sc, "previous_base": sc_prev}
os.makedirs(os.path.dirname(OUTJ), exist_ok=True)
json.dump(out, open(OUTJ, "w"), indent=1, ensure_ascii=False)

# ---- 6. markdown
f3 = lambda v: f"{v:.3f}"  # noqa: E731
f2 = lambda v: f"{v:.2f}".replace(".", ",")  # noqa: E731
f3c = lambda v: f"{v:.3f}".replace(".", ",")  # noqa: E731
f1p = lambda v: f"{100 * v:.1f}".replace(".", ",")  # noqa: E731
fmu = lambda v: f"{round(v + 1e-9, 2):.2f}".replace(".", ",")  # noqa: E731  (1.495 → 1,50; float 1.495 < 1.495 dziesiętnie)
CHN = {"D0D0": "D⁰D⁰", "DpD0": "D⁺D⁰", "DpDp": "D⁺D⁺", "common": "wspólne"}
L = ["# Budżet systematyk μ par D–D z tune Pythii (MPI, CR z CRqcd, Monash)", "",
     f"Skrypt `TUNE_pythia/budget_with_tune.py` (z `charm_work/`). Wejścia: `{SYST}`, `{TUNE}`"
     + (f", `truth_norm.json` z `{TRUTH}` (ścieżka `norm.collect`)." if coll else "."),
     "Wartości względne (Δμ/μ, tj. Δμ przy μ=1). Nominalny fit, Z i μ centralne bez zmian.", "",
     "## Reguły (z kodu `charmtag.syst`)",
     "- Normalizacja: `norm.py` l. 97–110 — δ grupy tune = max po wariantach |[par/D*]_war/[par/D*]_A3 − 1|; "
     "l. 119 — **Δμ/μ = |δ|** (1. rząd). `dmu_rel(δ) = 1/(1+δ) − 1` (l. 26–28) służy tylko do pola asymetrii (l. 120), nie do tabeli.",
     "- „wspólne”: `report.py` l. 40–41, 76–77 — Σ w·Δμ / Σ w, w = 1/σ²_stat (Asimov ślepy: "
     + ", ".join(f"{CHN[ch]} {sig[ch]:.3f}" for ch in CH) + "), brak w kanale = 0.",
     "- Sumy: l. 86–93 — kwadratura wierszy (pomijane: pending oraz `kappa_partner*`); suma „wspólne” = kwadratura wierszy „wspólne”; "
     "l. 112–113 — „wszystko” = kwadratura grup.", "",
     "## Kontrola: odtworzenie starych sum",
     f"Max |Δ| „wspólne” w wierszach: {chk['rows_common_maxdiff']:.1e}; max |Δ| sum grup/„wszystko”: {chk['sums_maxdiff']:.1e} (próg {TOL}).", "",
     "## Wiersze tune (Normalizacja)", "",
     "| grupa | warianty | Δμ D⁰D⁰ | Δμ D⁺D⁰ | Δμ D⁺D⁺ | Δμ wspólne | asym. wspólne (↑ / ↓) |", "|---|---|---|---|---|---|---|"]
for g in TUNE_ROWS:
    t = tune[g]
    L.append(f"| {g} | {'/'.join(t['variants'])} | " + " | ".join(f3(t["dmu"][ch]) for ch in CH) +
             f" | {f3(t['dmu_common'])} | +{f3(t['asym_up']['common'])} / {f3(t['asym_down']['common'])} |")
L.append("| **tune razem (kwadr.)** | | " + " | ".join(f3(tune_only[ch]) for ch in CH) + f" | {f3(tune_only['common'])} | |")
L += ["", "Asymetria: Δμ/μ = 1/(1+δ) − 1. CRqcd obniża oczekiwaną liczbę par (δ ≈ −10…−11 %), więc przesuwa μ **tylko w górę** "
      f"(CR: +{f3(tune['tune_CR']['asym_up']['common'])} w górę vs {f3(tune['tune_CR']['asym_down']['common'])} w dół, wspólne). "
      "Tabela (zgodnie z kodem) używa symetrycznego |δ|.", "",
      "## Budżet przed / tune / po (względne)", "",
      "| suma | " + " | ".join(CHN[k] + " przed | tune | po" for k in K) + " |", "|---|" + "---|" * (3 * len(K))]
for lab, key in (("Normalizacja", "Normalizacja"), ("wszystko", "wszystko"),
                 ("Normalizacja bez N_D*", "Normalizacja_bez_NDS"), ("wszystko bez N_D*", "wszystko_bez_NDS")):
    L.append(f"| {lab} | " + " | ".join(f"{f3(res['przed'][key][k])} | {f3(tune_only[k])} | **{f3(res['po'][key][k])}**" for k in K) + " |")
L += ["", f"Wiersz N_D* z tabeli ślepej (`ndstar_fit`, {nd['note']}): {f3(ndv['common'])} we wszystkich kanałach — odjęty kwadratowo w wierszach „bez N_D*”, "
      f"bo wynik końcowy podaje N_D* osobno (powtórka NDS 12.11: N_D* = {NDSTAR_N}, {f1p(NDSTAR_REL)} %).", "",
      f"## Przeskalowanie do μ_wsp. = {f3c(MU_COMMON)} (≈ {fmu(MU_COMMON)})", "",
      f"Baza N_D* = {NDSTAR_N} ({f1p(NDSTAR_REL)} %) z `NDS_ndstar/rerun_2026-11-12/rerun.md` (decyzja z 12.11.2026).", "",
      f"| wielkość | względna | × {f3c(MU_COMMON)} |", "|---|---|---|",
      f"| (a) syst. tune (MPI ⊕ CR ⊕ Monash) | {f3(tune_only['common'])} | ± {f2(sc['tune_abs'])} |",
      f"| (b) syst. łączna bez N_D*, z tune | {f3(res['po']['wszystko_bez_NDS']['common'])} | ± {f2(sc['syst_bez_NDS_z_tune_abs'])} |",
      f"| (c) syst. łączna bez N_D*, bez tune | {f3(res['przed']['wszystko_bez_NDS']['common'])} | ± {f2(sc['syst_bez_NDS_bez_tune_abs'])} |",
      f"| info: syst. łączna z N_D* (0,105), z tune | {f3(res['po']['wszystko']['common'])} | ± {f2(sc['syst_z_NDS_z_tune_abs'])} |",
      f"| info: N_D* (NDS, rerun 12.11) | {f3(NDSTAR_REL)} | ± {f2(sc['NDS_abs'])} |",
      f"| info: tune CR asym. (↑ / ↓) | +{f3(tune['tune_CR']['asym_up']['common'])} / {f3(tune['tune_CR']['asym_down']['common'])} | "
      f"+{f2(sc['tune_CR_up_abs'])} / {f2(sc['tune_CR_down_abs'])} |", "",
      "Per kanał tylko wartości względne (tabela wyżej).", "",
      f"Poprzednia baza (historycznie, zapis z 12.11 przed zmianą): μ = {f2(PREV_MU_COMMON)}, N_D* = {PREV_NDSTAR_N} ({f1p(PREV_NDSTAR_REL)} %) → "
      f"tune ± {f2(sc_prev['tune_abs'])}, syst bez N_D* ± {f2(sc_prev['syst_bez_NDS_z_tune_abs'])} "
      f"(asym. +{f2(sc_prev['syst_bez_NDS_asym_up_abs'])}/−{f2(sc_prev['syst_bez_NDS_asym_down_abs'])}), N_D* ± {f2(sc_prev['NDS_abs'])}.", "",
      "## Zapis wyniku (baza N_D* od 12.11)", "",
      f"> **μ_wsp. = {fmu(MU_COMMON)} +0,48/−0,46 (stat; profil +0,41/−0,39 × 1,18 z jackknife po zdarzeniach, 13.11) ± {f2(sc['syst_bez_NDS_z_tune_abs'])} (syst, bez N_D*) ± {f2(sc['NDS_abs'])} (N_D*)**", "",
      f"Składnik tune w ± {f2(sc['syst_bez_NDS_z_tune_abs'])} jest symetryzowany (|δ|); w rzeczywistości CRqcd przesuwa μ tylko w górę "
      f"(CR: +{f2(sc['tune_CR_up_abs'])} / {f2(sc['tune_CR_down_abs'])}). Z = 4,06σ (4,11σ) bez zmian — normalizacja nie zmienia Z.",
      f"Wariant asymetryczny (informacyjnie, tune jako skrajne 1/(1+δ)−1 per grupa: +{f3(sc['tune_asym_up_rel'])} / −{f3(sc['tune_asym_down_rel'])}): "
      f"syst. bez N_D* = +{f2(sc['syst_bez_NDS_asym_up_abs'])} / −{f2(sc['syst_bez_NDS_asym_down_abs'])}.", ""]
open(OUTM, "w").write("\n".join(L) + "\n")
print("\n".join(L))
