"""Tabela systematyk (markdown + json) z części: core_<kanał>.json, norm.json, files.json, wp.json (etap U)."""
from __future__ import annotations

import json
import math
import os


from ..pairs import CHANNELS
from . import SOURCES_DOC

GROUPS = [
    ("Kształty", ["width_D0", "width_Dp", "mscale_D0", "mscale_Dp", "bkg_poly3", "swap_D0", "refl_Ds", "refl_Lc"]),
    ("SF", ["sf_D0", "sf_Dp", "sf_tag"]),
    ("Normalizacja", ["mc_stat", "ndstar_fit", "tune_MPI", "tune_CR", "tune_Monash", "eff_ntrk", "eff_map_sample",
                      "diffraction_x1p5", "minbias_mix"]),
    ("κ, kategorie", ["kappa_channel", "kappa_partner03", "kappa_partner05", "cat_Y"]),
]
NORM_DOC = {"mc_stat": "statystyka MC par (+D*)", "ndstar_fit": "N_D* z fitu Δm (dane)", "tune_MPI": "tune: MPI pT0Ref ±",
            "tune_CR": "tune: CR range / tryb", "tune_Monash": "tune: Monash vs A3", "eff_ntrk": "wydajność(n_trk)",
            "eff_map_sample": "mapa ε: próbka MC", "diffraction_x1p5": "dyfrakcja ×1,5", "minbias_mix": "skład minbias ND vs mieszanka"}
SCALE_ONLY = {"sf_D0", "sf_Dp", "sf_tag"}       # tylko skala sygnału → nie zmienia obserwowanego Z (test μ=0)


def _ld(p):
    return json.load(open(p)) if os.path.exists(p) else None


def load_parts(d):
    core = {ch: _ld(os.path.join(d, f"core_{ch}.json")) for ch in CHANNELS}
    return core, _ld(os.path.join(d, "norm.json")), _ld(os.path.join(d, "files.json")), _ld(os.path.join(d, "wp.json"))


def build(d, title="Systematyki par D–D (ŚLEPO)"):
    core, norm, files, wp = load_parts(d)
    chs = [ch for ch in CHANNELS if core.get(ch)]
    sig = {ch: core[ch]["nominal"]["err"] for ch in chs}
    Z0 = {ch: core[ch]["nominal"]["Z"] for ch in chs}
    w = {ch: 1 / sig[ch] ** 2 for ch in chs}
    W = sum(w.values()) or 1.0
    rows, tot = [], {g: {ch: 0.0 for ch in chs} for g, _ in GROUPS}
    totc = {g: 0.0 for g, _ in GROUPS}
    shapeZ = {ch: 0.0 for ch in chs}            # kwadratura Δμ źródeł zmieniających tło/kształt (wpływ na Z)
    for g, names in GROUPS:
        for nm in names:
            dm, note, pend = {}, "", False
            if g == "Normalizacja":
                if not norm or nm not in norm:
                    continue
                v = norm[nm]
                pend = bool(v.get("pending"))
                dm = {ch: v["dmu_rel"].get(ch) for ch in chs}
                note = v.get("how", "")
                desc = NORM_DOC.get(nm, nm)
            else:
                desc = SOURCES_DOC.get(nm, nm)
                for ch in chs:
                    s = core[ch]["sources"].get(nm)
                    if not s or "dmu_sym" not in s:
                        continue
                    if "up" in s and "down" in s:
                        dm[ch] = s["dmu_sym"]
                    else:
                        dm[ch] = s["dmu_sym"]
                    if nm == "bkg_poly3" and "blind_dev" in s:
                        bd = s["blind_dev"]
                        note += f"{ch}: Δdev {bd['ddev']:.0f}/{bd['dndf']}; "
                    if nm.startswith("kappa") and "down" in s:
                        note += f"{ch}: B→alt {s['up']['dmu']:+.2f}, alt→B {s['down']['dmu']:+.2f}; "
            if not dm and not pend:
                continue
            vals = {ch: dm.get(ch) for ch in chs}
            comb = None
            if any(vals[ch] is not None for ch in chs):                 # brak w kanale = źródło nie dotyczy (0)
                comb = sum(w[ch] * (vals[ch] or 0.0) for ch in chs) / W   # w pełni skorelowane w obrębie źródła
            # ΔZ łączne (kwadratura kanałów), z wariantów Asimov
            dz = None
            if g != "Normalizacja" and nm not in SCALE_ONLY:
                zz = []
                for ch in chs:
                    s = core[ch]["sources"].get(nm)
                    zz.append(Z0[ch] + (s["dZ_min"] if s and "dZ_min" in s else 0.0))
                dz = math.sqrt(sum(z * z for z in zz)) - math.sqrt(sum(Z0[ch] ** 2 for ch in chs))
            if not pend and not nm.startswith("kappa_partner"):
                for ch in chs:
                    if vals[ch] is not None:
                        tot[g][ch] += vals[ch] ** 2
                        if g != "Normalizacja" and nm not in SCALE_ONLY:
                            shapeZ[ch] += vals[ch] ** 2
                if comb is not None:
                    totc[g] += comb ** 2
            rows.append((g, nm, desc, vals, comb, dz, note.strip(), pend))
    # κ partnera: wariant interpretacyjny (T12), nie dodajemy do sumy; plan B ↔ kanał wchodzi
    L = [f"# {title}", ""]
    L.append("Δμ przy μ=1 (bezwzględne = względne), Asimov z nuisance'ów obszaru ślepego; „wspólne” = średnia ważona 1/σ²_stat "
             "(w pełni skorelowane w obrębie źródła). ΔZ = zmiana Z_Asimov w kwadraturze kanałów (najgorszy kierunek). "
             "Normalizacja: Δμ/μ = ±δ(N_par), nie zmienia obserwowanego Z.")
    L.append("")
    L.append("| grupa | źródło | Δμ D0D0 | Δμ D⁺D0 | Δμ D⁺D⁺ | Δμ wspólne | ΔZ | uwagi |")
    L.append("|---|---|---|---|---|---|---|---|")
    f = lambda v: "—" if v is None else f"{v:.3f}"
    for g, nm, desc, vals, comb, dz, note, pend in rows:
        L.append(f"| {g} | {desc}{' (czeka)' if pend else ''} | " + " | ".join(f(vals.get(ch)) for ch in CHANNELS) +
                 f" | {f(comb)} | {('skala' if (g == 'Normalizacja' or nm in SCALE_ONLY) else '—') if dz is None else f'{dz:+.2f}'} | {note} |")
    L.append("| **σ_stat** | Asimov | " + " | ".join(f(sig.get(ch)) for ch in CHANNELS) + f" | {f(1 / math.sqrt(W))} | | Z_A: " +
             ", ".join(f"{ch} {Z0[ch]:.2f}" for ch in chs) + f"; kwadr. {math.sqrt(sum(z * z for z in Z0.values())):.2f} |")
    for g, _ in GROUPS:
        L.append(f"| **suma** | {g} | " + " | ".join(f(math.sqrt(tot[g][ch]) if ch in tot[g] else None) for ch in CHANNELS) +
                 f" | {f(math.sqrt(totc[g]))} | | |")
    allc = {ch: math.sqrt(sum(tot[g][ch] for g, _ in GROUPS)) for ch in chs}
    L.append("| **suma** | wszystko | " + " | ".join(f(allc.get(ch)) for ch in CHANNELS) + f" | {f(math.sqrt(sum(totc.values())))} | | |")
    # wpływ na Z: źródła tła/kształtu poszerzają σ pod H0 → Z_ch/√(1+(σ_sh/σ_stat)²)
    Zs = {ch: Z0[ch] / math.sqrt(1 + shapeZ[ch] / sig[ch] ** 2) for ch in chs}
    zq0, zq1 = math.sqrt(sum(z * z for z in Z0.values())), math.sqrt(sum(z * z for z in Zs.values()))
    L.append("")
    L.append(f"**Wpływ na Z (Asimov, μ=1):** bez systematyk {zq0:.2f}σ → z systematykami tła/kształtu/κ (kwadratura σ) {zq1:.2f}σ "
             f"(" + ", ".join(f"{ch} {Z0[ch]:.2f}→{Zs[ch]:.2f}" for ch in chs) + "). SF i normalizacja zmieniają tylko skalę μ.")
    prof = {ch: core[ch].get("profiled_shape_syst") for ch in chs}
    if all(prof.values()):
        L.append("- Profil (PairFit syst=True: szerokość, skala, tag, poly3) na Asimovie: " +
                 ", ".join(f"{ch} σ_μ {sig[ch]:.2f}→{prof[ch]['err']:.2f}, Z {Z0[ch]:.2f}→{prof[ch]['Z']:.2f}" for ch in chs))
    if files and "skipped" not in files:
        L.append("")
        L.append("**Pliki danych (ślepo):**")
        for k, v in files.get("overdispersion", {}).items():
            L.append(f"- dyspersja {k}: χ²/ndf {v['chi2']:.1f}/{v['ndf']} → σ_jk/σ ≈ {math.sqrt(max(v['phi'], 1)):.2f}")
        for ch, v in files.get("jackknife_sensitivity", {}).items():
            L.append(f"- jackknife czułości {ch}: Z bez pliku {v['Z_min']:.2f}–{max(v['Z_drop']):.2f} (sd_jk {v['Z_jk_sd']:.2f}), "
                     f"najbardziej wpływowy {v['max_drop_file']}")
        for nm, s in files.get("splits", {}).items():
            for ch, v in s["channels"].items():
                L.append(f"- podział {nm} {ch}: Z {v['Z'][0]:.2f}/{v['Z'][1]:.2f}, κ χ² {v['kappa_chi2_ndf'][0]:.1f}/{v['kappa_chi2_ndf'][1]}, "
                         f"ramiona pull {v['arms_pull']:+.2f}")
    elif files:
        L.append(f"\n**Pliki danych:** {files['skipped']}")
    if wp:
        L.append("")
        L.append("**Warianty WP/kategorii (ślepo; Δμ_dane po odsłonięciu):**")
        for nm, v in wp["variants"].items():
            if "missing" in v:
                L.append(f"- {nm}: brak ({v['missing']})")
                continue
            L.append(f"- {nm}: Z kwadr. {v['Z_quadrature']:.2f} (nominał {wp['nominal']['Z_quadrature']:.2f}); oczek. rozrzut stat. Δμ " +
                     ", ".join(f"{ch} {v[ch]['exp_stat_spread_dmu']:.2f}" for ch in CHANNELS if ch in v and "exp_stat_spread_dmu" in v[ch]))
    bss = _ld(os.path.join(d, "bkg_shape_salad_pairs.json"))   # B4: informacyjnie (B płynie w fit2d → nie wchodzi do Δμ)
    if bss and bss.get("channels"):
        L.append("")
        L.append("**ε_B pod pikiem (B4, `bkg_shape_salad`, informacyjnie):** niepewność względna przewidywanego B_pred pary "
                 "(Asimov/Zproxy/mixing-norm; w fit2d tło płynie, więc Δμ nie jest liczone) — " +
                 "; ".join(f"{ch} @80/70/50 (8–30): " + "/".join(f"{bss['channels'][ch].get(f'{wp}|8-30', {}).get('rel_unc_B', float('nan')):.2f}"
                                                              for wp in (80, 70, 50)) for ch in bss["channels"]))
    md = "\n".join(L) + "\n"
    js = {"rows": [{"group": g, "source": nm, "desc": desc, "dmu": vals, "dmu_common": comb, "dZ": dz, "note": note, "pending": pend}
                   for g, nm, desc, vals, comb, dz, note, pend in rows],
          "sigma_stat": sig, "Z_asimov": Z0, "Z_with_shape_syst": Zs, "totals": {g: {ch: math.sqrt(v) for ch, v in t.items()} for g, t in tot.items()},
          "total_all": allc, "total_common": math.sqrt(sum(totc.values()))}
    with open(os.path.join(d, "syst_blind.md"), "w") as fh:
        fh.write(md)
    with open(os.path.join(d, "syst_blind.json"), "w") as fh:
        json.dump(js, fh, indent=1, default=float)
    return md, js
