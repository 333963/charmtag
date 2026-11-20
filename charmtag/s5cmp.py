"""`charmtag s5compare`: porównanie wariantów tła (etap S5) — „ile zyskuje sygnał” względem nominalnego H.

  charmtag s5compare -c cfg --run nominal=<katalog_biegu> --run farveto=<katalog_wariantu> ... --nominal nominal
                     --out <prefiks_wyjścia> [--species D0 Dplus] [--wp 0.70]
Katalog wariantu może mieć s5_variant.json: {"models_dir","scores_dir","validation_dir","pairs_dir","fit_tag","label"}
(wariant „tylko fit”, np. szablon refleksji z danych na nominalnych score). Brak pliku → wszystko w katalogu.
Metryki per gatunek × bin pT (dane = data_none):
  ε_B OOF (własne tło treningowe — NIE porównywalne między źródłami), ε_B(dane SB) (wspólna definicja: |m−m0|>sideband_min),
  ε_B(MC SR) (incl_mc truth≥8 w oknie sygnału) i Zproxy_MC = ε_S/√ε_B(MC SR), R_win(dane, pasmo wewn.), R_lewy/prawy SB,
  kontrola D* (ε dane/MC @WP, D0), fit masy: N ± σ, Z = N/σ_N (główna miara), S/√(S+B) w ±2σ, χ²/ndf;
  zysk = Z/Z_nom i (Z/Z_nom)² = ekwiwalent lumi. Pary: Z_Asimov (fit2d ślepo, nuisance z ramion) per kanał i kwadratura.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys

import numpy as np

log = logging.getLogger("charmtag")


def _j(p):
    try:
        with open(p) as f:
            return json.load(f)
    except Exception:
        return None


def variant_paths(d):
    v = _j(os.path.join(d, "s5_variant.json")) or {}
    return {"dir": d, "models": v.get("models_dir") or os.path.join(d, "models"),
            "scores": v.get("scores_dir") or os.path.join(d, "scores"),
            "validation": v.get("validation_dir") or os.path.join(d, "validation"),
            "fits": os.path.join(d, "fits"), "pairs": v.get("pairs_dir") or os.path.join(d, "pairs"),
            "fit_tag": v.get("fit_tag", ""), "label": v.get("label", ""), "bkg_source": v.get("bkg_source", "")}


def _scan(path, cols, chunk=5_000_000):
    import h5py
    with h5py.File(path, "r") as f:
        n = f[cols[0]].shape[0]
        for s in range(0, n, chunk):
            yield {c: f[c][s:s + chunk] for c in cols if c in f}


def eff_tables(cfg, species, V, wps):
    """ε_B(dane SB) i ε_B(MC SR) per bin pT przy progach wariantu; z plików scores (porcjami)."""
    from .model import pass_wp
    sp = cfg["species"][species]
    meta = _j(os.path.join(V["models"], f"{species}_H", "meta.json"))
    if not meta:
        return {}
    thr = meta["thresholds"]
    bins = cfg["training"]["pt_bins"]
    out = {}
    for key, kind in (("data_none", "dataSB"), ("incl_mc", "mcSR")):
        p = os.path.join(V["scores"], f"{species}_H", f"{key}.h5")
        if not os.path.exists(p):
            continue
        acc = {}
        for c in _scan(p, ["m", "pt", "score", "veto_ok", "truth"]):
            dm = np.abs(c["m"] - sp["m0"])
            ok = c["veto_ok"].astype(bool)
            if kind == "dataSB":
                ok &= dm > sp["sideband_min"]
            else:
                t = c["truth"]
                ok &= (dm < sp["signal_window"]) & (t >= 8) & (t < 100)
            for lo, hi in zip(bins[:-1], bins[1:]):
                b = ok & (c["pt"] >= lo) & (c["pt"] < hi)
                a = acc.setdefault(f"{lo}-{hi}", {"n": 0, **{w: 0 for w in wps}})
                a["n"] += int(b.sum())
                for w in wps:
                    a[w] += int((b & pass_wp(c["score"], c["pt"], thr, w)).sum())
        for k, a in acc.items():
            row = out.setdefault(k, {})
            for w in wps:
                e = a[w] / a["n"] if a["n"] else float("nan")
                row[f"effB_{kind}_{int(w*100)}"] = e
                row[f"effB_{kind}_{int(w*100)}_err"] = float(np.sqrt(max(e * (1 - e), 0) / a["n"])) if a["n"] else float("nan")
                if kind == "mcSR":
                    row[f"Zproxy_MC_{int(w*100)}"] = w / np.sqrt(e) if e > 0 else float("nan")
    return out


def species_table(cfg, species, V, wps, fit_wps, ds="data_none", with_eff=True):
    bins = cfg["training"]["pt_bins"]
    keys = [f"{lo}-{hi}" for lo, hi in zip(bins[:-1], bins[1:])]
    T = {k: {} for k in keys}
    mo = _j(os.path.join(V["models"], f"{species}_H", "metrics_oof.json")) or {}
    for k in keys:
        for w in wps:
            T[k][f"effB_oof_{int(w*100)}"] = (mo.get(k) or {}).get(f"effB_{int(w*100)}")
    if with_eff:
        for k, r in eff_tables(cfg, species, V, wps).items():
            T.setdefault(k, {}).update(r)
    va = ((_j(os.path.join(V["validation"], f"{species}.json")) or {}).get("variants") or {}).get("H") or {}
    for k in keys:
        rw = (va.get("R_win_data_inner_bins") or {}).get(k) or {}
        for w in wps:
            x = rw.get(f"{int(w*100)}") or {}
            T[k][f"Rwin_{int(w*100)}"] = x.get("R")
            T[k][f"Rwin_{int(w*100)}_err"] = x.get("err")
            d, m = (va.get("dstar_data") or {}).get(k) or {}, (va.get("dstar_mc") or {}).get(k) or {}
            if f"eff_{int(w*100)}" in d:
                T[k][f"dstar_data_{int(w*100)}"] = d[f"eff_{int(w*100)}"]
                T[k][f"dstar_data_{int(w*100)}_err"] = d.get(f"err_{int(w*100)}")
                T[k][f"dstar_mc_{int(w*100)}"] = m.get(f"eff_{int(w*100)}")
                mc = m.get(f"eff_{int(w*100)}")
                T[k][f"dstar_ratio_{int(w*100)}"] = d[f"eff_{int(w*100)}"] / mc if mc else None
    glob_ = {"R_sb_lowhigh": va.get("R_sb_data_lowhigh"), "R_win_mcbkg": va.get("R_win_mcbkg_window")}
    fj = _j(os.path.join(V["fits"], f"{species}_H{V['fit_tag']}.json")) or {}
    for k in keys:
        e = (fj.get("bins") or {}).get(k) or {}
        for w in fit_wps:
            f = e.get(f"{ds}_wp{int(w*100)}") or {}
            if "N_sig" in f and f.get("N_sig_err"):
                T[k][f"N_{int(w*100)}"] = f["N_sig"]
                T[k][f"Nerr_{int(w*100)}"] = f["N_sig_err"]
                T[k][f"Z_{int(w*100)}"] = f["N_sig"] / f["N_sig_err"]
                T[k][f"ssb_{int(w*100)}"] = f.get("signif")
                T[k][f"chi2ndf_{int(w*100)}"] = f["chi2"] / max(f["ndf"], 1)
                T[k][f"valid_{int(w*100)}"] = f.get("valid")
                if "refl_ratio_used" in f:
                    T[k][f"refl_r_{int(w*100)}"] = f["refl_ratio_used"]
                    T[k][f"refl_r_mc_{int(w*100)}"] = f.get("refl_ratio_mc")
    return T, glob_


def pairs_table(V):
    r = _j(os.path.join(V["pairs"], "fit2d_blind.json"))
    if not r:
        return None
    out = {"Z_quadrature": r.get("Z_asimov_quadrature"), "channels": {}}
    for ch, v in (r.get("channels") or {}).items():
        if "skipped" in v:
            out["channels"][ch] = {"skipped": v["skipped"]}
            continue
        a = v.get("asimov_from_arms") or {}
        out["channels"][ch] = {"Z": a.get("Z"), "sigma_mu": a.get("sigma_mu"), "dev_ndf": v.get("blind_region_dev_ndf")}
    mix = _j(os.path.join(V["pairs"], "mix_control.json"))
    if mix:
        out["mix_control"] = {ch: {k: c.get(k) for k in ("sbsb_chi2_ndf", "corr_m1m2_same_sbsb", "corr_m1m2_mix_sbsb",
                                                         "arms_ratio_same_over_mix", "centre_mix_pred_uncorrelated")}
                              for ch, c in mix["channels"].items()}
    return out


def _f(x, n=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return f"{x:.{n}g}" if abs(x) < 1e4 else f"{x:.0f}"


def compare(cfg, runs: dict, nominal: str, species=("D0", "Dplus"), wp=0.70, wps=(0.80, 0.70, 0.50), with_eff=True):
    fit_wps = sorted(set(list(wps) + [wp]), reverse=True)
    res = {"nominal": nominal, "wp": wp, "variants": {}}
    for name, d in runs.items():
        V = variant_paths(d)
        res["variants"][name] = {"paths": V, "species": {}, "pairs": pairs_table(V)}
        for sp in species:
            T, G = species_table(cfg, sp, V, list(wps), fit_wps, with_eff=with_eff)
            res["variants"][name]["species"][sp] = {"bins": T, "global": G}
    nom = res["variants"].get(nominal)
    w = int(wp * 100)
    L = ["# S5 — warianty tła: ile zyskuje sygnał (vs nominalny H)", "",
         f"WP główny {w}% (ε_S MC); dane = data_none; Z = N/σ_N z fitu masy; zysk lumi = (Z/Z_nom)². "
         "ε_B(OOF) liczone na własnym tle treningowym (nieporównywalne między źródłami); ε_B(dane SB) i ε_B(MC SR) — wspólne definicje.", ""]
    for sp in species:
        L += [f"## {sp}", "", f"| wariant | pT | ε_B(dane SB)@{w} | ε_B(MC SR)@{w} | Zproxy_MC | R_win@{w} | D* dane/MC@{w} | N ± σ | Z | Z/Z_nom | lumi× | χ²/ndf |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for name, R in res["variants"].items():
            for k, r in R["species"][sp]["bins"].items():
                n0 = ((nom or {}).get("species", {}).get(sp, {}).get("bins", {}).get(k) or {}).get(f"Z_{w}") if nom else None
                z = r.get(f"Z_{w}")
                g = z / n0 if (z is not None and n0) else None
                r[f"gain_Z_{w}"] = g
                r[f"gain_lumi_{w}"] = g * g if g is not None else None
                ds = r.get(f"dstar_ratio_{w}")
                L.append(f"| {name} | {k} | {_f(r.get(f'effB_dataSB_{w}'))} | {_f(r.get(f'effB_mcSR_{w}'))} | {_f(r.get(f'Zproxy_MC_{w}'))} | "
                         f"{_f(r.get(f'Rwin_{w}'))} | {_f(ds)} | {_f(r.get(f'N_{w}'), 5)} ± {_f(r.get(f'Nerr_{w}'), 3)} | {_f(z)} | {_f(g)} | "
                         f"{_f(r[f'gain_lumi_{w}'])} | {_f(r.get(f'chi2ndf_{w}'))} |")
        # najlepszy WP per bin (informacyjnie — wybór na danych, nie do raportu fizycznego)
        L += ["", f"Z przy innych WP ({', '.join(str(int(x*100)) for x in fit_wps)}):", ""]
        for name, R in res["variants"].items():
            L.append(f"- {name}: " + "; ".join(f"{k}: " + "/".join(_f(r.get(f'Z_{int(x*100)}')) for x in fit_wps)
                                               for k, r in R["species"][sp]["bins"].items()))
        L.append("")
    L += ["## Pary D–D (fit2d ślepo, Z_Asimov z nuisance'ów z ramion)", "", "| wariant | D0D0 | D⁺D0 | D⁺D⁺ | kwadratura | Z/Z_nom |", "|---|---|---|---|---|---|"]
    zn = ((nom or {}).get("pairs") or {}).get("Z_quadrature") if nom else None
    for name, R in res["variants"].items():
        P = R["pairs"]
        if not P:
            L.append(f"| {name} | — | — | — | — | — |")
            continue
        c = P["channels"]
        zq = P.get("Z_quadrature")
        L.append(f"| {name} | {_f((c.get('D0D0') or {}).get('Z'))} | {_f((c.get('DpD0') or {}).get('Z'))} | {_f((c.get('DpDp') or {}).get('Z'))} | "
                 f"{_f(zq)} | {_f(zq / zn if (zq and zn) else None)} |")
    for name, R in res["variants"].items():
        mc = (R["pairs"] or {}).get("mix_control")
        if mc:
            L += ["", f"Kontrola mixingiem nóg ({name}): " + "; ".join(
                f"{ch}: SB×SB χ²/ndf {_f(v.get('sbsb_chi2_ndf'))}, corr same {_f(v.get('corr_m1m2_same_sbsb'), 2)} / mix {_f(v.get('corr_m1m2_mix_sbsb'), 2)}, "
                f"ramiona same/mix {_f(v.get('arms_ratio_same_over_mix'))}" for ch, v in mc.items())]
    return res, "\n".join(L) + "\n"


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "s5compare":
        argv = argv[1:]
    ap = argparse.ArgumentParser(prog="charmtag s5compare", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--run", action="append", required=True, help="nazwa=katalog")
    ap.add_argument("--nominal", default=None)
    ap.add_argument("--out", required=True, help="prefiks wyjścia (…/compare → compare.md, compare.json)")
    ap.add_argument("--species", nargs="*", default=["D0", "Dplus"])
    ap.add_argument("--wp", type=float, default=0.70)
    ap.add_argument("--no-eff", action="store_true", help="bez czytania plików scores (szybko)")
    a = ap.parse_args(argv)
    from . import config as C
    cfg = C.load(a.config)  # s5compare: ścieżki biegów z --run, paths.* niewymagane
    runs = dict(r.split("=", 1) for r in a.run)
    res, md = compare(cfg, runs, a.nominal or next(iter(runs)), tuple(a.species), a.wp, with_eff=not a.no_eff)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out + ".json", "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    with open(a.out + ".md", "w") as f:
        f.write(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
