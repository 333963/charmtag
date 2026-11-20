#!/usr/bin/env python3
"""TRUTHNORM — raport (19.11.2026). Czyta res/tn_*.json → TRUTHNORM.md, truthnorm.json."""
import glob
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = {os.path.basename(p)[3:-5]: json.load(open(p)) for p in glob.glob(os.path.join(HERE, "res", "tn_*.json"))}
CH = ("D0D0", "DpD0", "DpDp")
CLS = ["J12only", "MBTS_j20", "MBTS_noj20", "other"]
ND_05 = 5583.649730277749; ND_NOW = 4436.0; SC = ND_05 / ND_NOW       # skala zapisu bieżącego 428770 (N_D* NDS 12.11)
JK428, JK427 = 1.18, 1.13                                            # jackknife zdarzeń (jk_ev200 428770; ub427jk)
SYST_REL_428 = 0.42 / 1.50                                           # VER_budget/RAPORT_gen_norm.md (prawda mb_nd), względnie
ND_REL_428 = 0.33 / 1.50
_B = json.load(open(os.path.join(HERE, "..", "r427929", "budzet_427929.json")))   # 19.11: syst asym. (tune ↑ = fragmentacja ALICE), N_D* z ub427nds
SYST_UP_427, SYST_DN_427, ND_REL_427 = _B["syst_rel_up"], _B["syst_rel_dn"], _B["ndstar"]["total"]
SYST_REL_427 = SYST_DN_427


def comb(x):
    c = x["fit"]["combined"]; lo, hi = c.get("mu_68") or [float("nan")] * 2
    return c["mu"], lo, hi, c["Z"]


def row(name, x, scale=1.0):
    mu, lo, hi, z = comb(x)
    ch = " | ".join(f"{x['fit'][c]['mu_hat'] * scale:.2f} ± {x['fit'][c]['mu_err_hesse'] * scale:.2f}" for c in CH)
    return f"| {name} | {ch} | **{mu * scale:.2f}** [{lo * scale:.2f}; {hi * scale:.2f}] | {z:.2f} |", (mu * scale, lo * scale, hi * scale, z)


out, md = {}, ["# TRUTHNORM — μ par D–D z normalizacją z prawdy mb_nd (decyzja z 19.11; GEN_mcstat/PROCEDURA.md §6)", "",
               f"n_pair(ch) × f_ch, f = {R['nom428']['f']} (`GEN_mcstat/out/truth_mbnd.json`); tło bez zmian; pełny ponowny fit (`truthnorm_run.py`, "
               "ub_jk.replica_fit, charmtag PF_pairs/src). Kontrola f ≡ 1 = odtworzenie nominału.", "",
               "## (a) 428770 i (c) 427929", "",
               "| bieg | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %, profil] | Z |", "|---|---|---|---|---|---|"]
for key, nm in (("nom428_f1", "428770 kontrola f ≡ 1 (próbka 05.11)"), ("nom428", "428770 prawda mb_nd (próbka 05.11, N_D* 5584)"),
                ("nom427_f1", "427929 kontrola f ≡ 1"), ("nom427", "427929 prawda mb_nd")):
    if key in R:
        s, v = row(nm, R[key]); md.append(s); out[key] = v
if "nom428" in R:
    s, v = row("428770 prawda mb_nd — skala zapisu (N_D* 4436)", R["nom428"], SC); md.append(s); out["nom428_scaled"] = v
md.append("")
ok428 = abs(out.get("nom428_f1", [0])[0] - 1.1889) < 2e-3 if "nom428_f1" in out else None
ok427 = abs(out.get("nom427_f1", [0])[0] - 0.567) < 5e-3 if "nom427_f1" in out else None
md.append(f"Kontrole: 428770 f ≡ 1 → {out.get('nom428_f1', [float('nan')])[0]:.4f} (ref. 1,1889: {'OK' if ok428 else 'FAIL'}); "
          f"427929 f ≡ 1 → {out.get('nom427_f1', [float('nan')])[0]:.4f} (ref. jackknife 0,567, nogi przebudowane; UB427 0,53: {'OK' if ok427 else 'FAIL'}).")
md += ["", "## Zapis wyników (stat × jackknife zdarzeń; syst bez N_D* i N_D* względnie z budżetów)", ""]
if "nom428_scaled" in out:
    mu, lo, hi, z = out["nom428_scaled"]
    md.append(f"- **428770:** μ = {mu:.2f} +{(hi - mu) * JK428:.2f}/−{(mu - lo) * JK428:.2f} (stat, ×{JK428}) ± {mu * SYST_REL_428:.2f} (syst, prawda mb_nd) "
              f"± {mu * ND_REL_428:.2f} (N_D*); Z = {z:.2f} (profil; ≈ {z / JK428:.2f} z korelacjami zdarzeń). Trigger: osobno (TRIG_A/SYST_TRIG.md).")
    out["write_428770"] = {"mu": mu, "stat_up": (hi - mu) * JK428, "stat_dn": (mu - lo) * JK428, "syst": mu * SYST_REL_428, "nds": mu * ND_REL_428, "Z": z}
if "nom427" in out:
    mu, lo, hi, z = out["nom427"]
    tot = math.hypot((hi - mu) * JK427, mu * math.hypot(SYST_UP_427, ND_REL_427))      # strona w kierunku μ = 1 (góra)
    md.append(f"- **427929 (wynik główny, decyzja z 19.11):** μ = {mu:.2f} +{(hi - mu) * JK427:.2f}/−{(mu - lo) * JK427:.2f} (stat, ×{JK427}) "
              f"+{mu * SYST_UP_427:.2f}/−{mu * SYST_DN_427:.2f} (syst) ± {mu * ND_REL_427:.2f} (N_D*); "
              f"Z = {z:.2f} (≈ {z / JK427:.2f}); vs Pythia (μ = 1): {(mu - 1) / tot:+.2f}σ (σ_tot górne ≈ {tot:.2f}, syst względne z BUDZET_427929 przy μ jak w biegu).")
    out["write_427929"] = {"mu": mu, "stat_up": (hi - mu) * JK427, "stat_dn": (mu - lo) * JK427, "syst_up": mu * SYST_UP_427, "syst_dn": mu * SYST_DN_427,
                           "syst": mu * SYST_DN_427, "nds": mu * ND_REL_427, "Z": z, "vs_pythia_sigma": (mu - 1) / tot}
# (b) klasy
md += ["", "## (b) 428770 per klasa triggera (prawda mb_nd; r_c, nogi i tło jak UBTRIG)", "",
       "| klasa | r_c | pary D0D0/DpD0/DpDp (= UBTRIG?) | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | Z |", "|---|---|---|---|---|---|---|---|"]
vals, errs, used = [], [], []
for i, c in enumerate(CLS):
    x = R.get(f"cls{i}")
    if not x:
        md.append(f"| {c} | – | – | – | – | – | brak | – |"); continue
    s, v = row(c, x)
    same = x["pairs_class"] == x["pairs_class_ubtrig"]
    cells = s.split(" | ")
    md.append(f"| {c} | {x['r']:.3f} ± {x['r_err']:.3f} | {'/'.join(str(x['pairs_class'][h]) for h in CH)} ({'tak' if same else 'NIE'}) | " + " | ".join(cells[1:]))
    mu, lo, hi, z = v; sig = 0.5 * (hi - lo)
    vals.append(mu); errs.append(math.hypot(sig, mu * x["r_err"] / x["r"])); used.append(c); out[f"cls_{c}"] = v
if "J12only" in used and len(used) > 1:
    j = used.index("J12only"); rest = [k for k in range(len(used)) if k != j]
    w = [1 / errs[k] ** 2 for k in rest]; mr = sum(wk * vals[k] for wk, k in zip(w, rest)) / sum(w); er = 1 / math.sqrt(sum(w))
    d = (vals[j] - mr) / math.hypot(errs[j], er)
    out["J12_vs_rest"] = {"mu_J12": vals[j], "err_J12": errs[j], "mu_rest": mr, "err_rest": er, "diff_sigma": d}
    md += ["", f"J12only vs reszta (próbka 05.11): {vals[j]:.2f} ± {errs[j]:.2f} vs {mr:.2f} ± {er:.2f} → {d:+.2f}σ (bez czynnika jackknife; ≈ {d / JK428:+.2f}σ z nim). "
           f"Na skali zapisu (×{SC:.3f}): J12only {vals[j] * SC:.2f}, reszta {mr * SC:.2f}."]
md += ["", "Uwagi: (1) 427929 — f_ch z 428770 (metoda stosunku R_X; przybliżenie: R_X liczony na kandydatach reko, przyjęty równy dla prawdy). "
       "(2) Syst względne przeniesione z budżetów (RAPORT_gen_norm: 0,42/1,50; BUDZET_427929: 0,387 bez N_D*, N_D* 0,284 — ten ostatni czeka na ub427nds). "
       "(3) Nominał w UNB_WYNIK.md NIE zmieniony.", ""]
json.dump(out, open(os.path.join(HERE, "truthnorm.json"), "w"), indent=1)
open(os.path.join(HERE, "TRUTHNORM.md"), "w").write("\n".join(md))
print("\n".join(md))
