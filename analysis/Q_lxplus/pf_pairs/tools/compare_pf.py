#!/usr/bin/env python3
"""PF2: porównanie powtórki par (poprawne weto wspólnych śladów) z biegiem nominalnym (stare weto), ½ lumi (1,65σ) i wrześniem (2,63σ).
Wejścia: <pf>/pairs/{pairs_counts.json, pair_inputs.json, fit2d_blind.json, fit2d_toys.md|json?}, <pf>/pairs/syst/*.json, <pf>/mix/pairs/mix_control.json,
<nominal>/pairs/{pair_inputs.json, fit2d_blind.json}. ŚLEPO: tylko Z_Asimov, σ_μ, dev/ndf, zliczenia SB×SB/ramion — bez μ z centrum danych.
Użycie: compare_pf.py --pf /eos/.../pf_pairs --nominal /eos/.../runs/lxplus_v11 --out compare_pf  → compare_pf.md/.json"""
import argparse, glob, json, os
REF = {"half_lumi_mac": {"Z": 1.65, "sigma_mu": {"D0D0": 1.05, "DpD0": 0.95, "DpDp": 3.87}, "dev_ndf": {"D0D0": 0.68, "DpD0": 1.44, "DpDp": 2.24}},
       "september_v10": {"Z": 2.63, "sigma_mu": {"D0D0": 0.67, "DpD0": 0.50, "DpDp": 1.23}},
       "full_lumi_tilt_rho_V": {"Z": 2.27, "dev_ndf": {"D0D0": 0.79, "DpD0": 1.27, "DpDp": 1.28}}}
CH = ("D0D0", "DpD0", "DpDp")
def jload(p):
    try: return json.load(open(p))
    except Exception: return None
def g(d, *ks, default=None):
    for k in ks:
        if not isinstance(d, dict) or k not in d: return default
        d = d[k]
    return d
def fit_summary(fb):
    out = {}
    if not fb: return out
    for ch in CH:
        c = g(fb, "channels", ch) or {}
        out[ch] = {"Z_asimov_arms": g(c, "asimov_from_arms", "Z"), "sigma_mu": g(c, "asimov_from_arms", "sigma_mu"),
                   "dev_ndf_blind": c.get("blind_region_dev_ndf"), "arms_dev_per_bin": c.get("arms_dev_per_bin"),
                   "rho": c.get("rho") or g(c, "bkg", "rho"), "tilt": c.get("tilt") or g(c, "bkg", "tilt"), "skipped": c.get("skipped")}
    out["Z_quadrature"] = fb.get("Z_asimov_quadrature") or g(fb, "asimov_nominal", "Z_quadrature")
    return out
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pf", required=True); ap.add_argument("--nominal", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); P = os.path.join(a.pf, "pairs"); N = os.path.join(a.nominal, "pairs")
    res = {"pf": a.pf, "nominal": a.nominal, "reference": REF}
    res["counts_pf"] = jload(os.path.join(P, "pairs_counts.json"))
    res["counts_mix"] = jload(os.path.join(a.pf, "mix", "pairs", "pairs_counts.json"))
    for tag, d in (("pf", P), ("nominal", N)):
        pi = jload(os.path.join(d, "pair_inputs.json")) or {}
        res[f"sbsb_{tag}"] = pi.get("sbsb"); res[f"n_pair_{tag}"] = pi.get("n_pair"); res[f"norm_mode_{tag}"] = g(pi, "meta", "norm_mode")
        res[f"fit2d_{tag}"] = fit_summary(jload(os.path.join(d, "fit2d_blind.json")))
    toys = jload(os.path.join(P, "toys", "fit2d_toys.json")) or jload(os.path.join(P, "fit2d_toys.json"))
    res["toys_pf"] = toys
    res["syst_pf"] = {os.path.basename(f): jload(f) for f in glob.glob(os.path.join(P, "syst", "*.json"))}
    res["mix_pf"] = jload(os.path.join(a.pf, "mix", "pairs", "mix_control.json"))
    # --- markdown
    L = ["# PF2 — powtórka par z poprawnym wetem wspólnych śladów (ŚLEPO)", ""]
    c = res["counts_pf"] or {}
    if c:
        b, s = c.get("pairs_before_veto"), c.get("pairs_shared_track")
        L += [f"- track_id: `{c.get('track_id')}`; pary przed wetem {b}, usunięte (wspólny ślad) {s}"
              + (f" = {100*s/b:.1f} %" if b else "") + f"; w histogramach {c.get('pairs_kept')} (region {c.get('region')})", ""]
    L += ["| kanał | SB×SB nominal | SB×SB PF | zmiana | Z_A nominal | Z_A PF | σ_μ nominal | σ_μ PF | dev/ndf nominal | dev/ndf PF |", "|---|---|---|---|---|---|---|---|---|---|"]
    for ch in CH:
        sn, sp = (res["sbsb_nominal"] or {}).get(ch), (res["sbsb_pf"] or {}).get(ch)
        sn = sum(sn.values()) if isinstance(sn, dict) else sn; sp = sum(sp.values()) if isinstance(sp, dict) else sp
        fn, fp = res["fit2d_nominal"].get(ch, {}), res["fit2d_pf"].get(ch, {})
        f = lambda x: "—" if x is None else (f"{x:.2f}" if isinstance(x, float) else str(x))
        chg = f"{100*(sp/sn-1):+.0f} %" if (sn and sp) else "—"
        L.append(f"| {ch} | {f(sn)} | {f(sp)} | {chg} | {f(fn.get('Z_asimov_arms'))} | {f(fp.get('Z_asimov_arms'))} | {f(fn.get('sigma_mu'))} | {f(fp.get('sigma_mu'))} | {f(fn.get('dev_ndf_blind'))} | {f(fp.get('dev_ndf_blind'))} |")
    L += ["", f"Z łączne (kwadratura): nominal {res['fit2d_nominal'].get('Z_quadrature')} → PF {res['fit2d_pf'].get('Z_quadrature')}; "
          f"odniesienia: ½ lumi 1,65σ, wrzesień 2,63σ, pełna lumi tilt+rho (V) 2,27σ — wszystkie ze starym wetem.", ""]
    rho = {ch: res["fit2d_pf"].get(ch, {}).get("rho") for ch in CH}
    L += [f"- ρ/tilt (PF): {rho} — oczekiwane ρ≈0 po wecie (MC: 0,24→0,01 w D⁺D⁺); jeśli ρ nadal ≠0 → model tilt+rho zostaje.",
          f"- toys: {'brak' if not toys else 'patrz toys_pf w JSON (pokrycie, pull, σ_μ mediana)'}",
          f"- syst (U) po PF: {len(res['syst_pf'])} plików; mixing (S5): {'brak' if not res['mix_pf'] else 'mix_control.json'}", "",
          "Szablon do docu T: kolumny powyżej + wiersz „wrzesień/½ lumi (stare weto)” dla kontekstu. Odsłonięcie: decyzja po tej tabeli."]
    open(a.out + ".md", "w").write("\n".join(L) + "\n")
    json.dump(res, open(a.out + ".json", "w"), indent=1, default=str)
    print("\n".join(L))
if __name__ == "__main__":
    main()
