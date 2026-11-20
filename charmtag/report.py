"""Raport: PNG + JSON + krótki markdown."""
from __future__ import annotations

import json
import os

import numpy as np

from . import store


def _plt():
    """Leniwy import matplotlib (REVIEW P2-18: start CLI bez ~1 s na import)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(_RC)
    return plt

C1, C2, C3, INK, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#333333", "#dddddd"
_RC = {"axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.edgecolor": "#888888",
       "axes.titlesize": 10, "font.size": 9, "lines.linewidth": 2, "figure.dpi": 110}


def _load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def plot_scores(oof, pt_bins, path, title):
    plt = _plt()
    fig, axs = plt.subplots(1, len(pt_bins) - 1, figsize=(4 * (len(pt_bins) - 1), 3.2))
    for ax, lo, hi in zip(np.atleast_1d(axs), pt_bins[:-1], pt_bins[1:]):
        b = (oof["pt"] >= lo) & (oof["pt"] < hi)
        for y, c, lab in ((1, C1, "sygnał MC (truth)"), (0, C2, "tło: sidebandy danych (ważone)")):
            s = oof[b & (oof["y"] == y)]
            ax.hist(s["score"], 50, range=(0, 1), weights=s["w_train"], density=True, histtype="step", color=c, label=lab, lw=1.6)
        ax.set_yscale("log")
        ax.set_title(f"{title}, pT {lo}–{hi} GeV")
        ax.set_xlabel("score (OOF)")
    np.atleast_1d(axs)[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_roc(oofs: dict, pt_bins, path, title):
    plt = _plt()
    fig, axs = plt.subplots(1, len(pt_bins) - 1, figsize=(4 * (len(pt_bins) - 1), 3.2))
    cols = [C1, C2, C3]
    for ax, lo, hi in zip(np.atleast_1d(axs), pt_bins[:-1], pt_bins[1:]):
        for (name, oof), c in zip(oofs.items(), cols):
            b = oof[(oof["pt"] >= lo) & (oof["pt"] < hi)]
            o = np.argsort(-b["score"].to_numpy())
            y, w = b["y"].to_numpy()[o], b["w_train"].to_numpy()[o]
            es = np.cumsum(w * (y == 1)) / max((w * (y == 1)).sum(), 1e-12)
            eb = np.cumsum(w * (y == 0)) / max((w * (y == 0)).sum(), 1e-12)
            ax.plot(es, np.maximum(eb, 1e-6), color=c, label=name)
        for wp in (0.5, 0.7, 0.8):
            ax.axvline(wp, color="#999999", lw=0.8, ls=":")
        ax.set_yscale("log")
        ax.set_xlabel("ε sygnału")
        ax.set_ylabel("ε tła")
        ax.set_title(f"{title}, pT {lo}–{hi} GeV")
    np.atleast_1d(axs)[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_fits(fitres, path, wp_tag, datasets):
    plt = _plt()
    bins = list(fitres["bins"])
    fig, axs = plt.subplots(len(datasets), len(bins), figsize=(4 * len(bins), 3.2 * len(datasets)), squeeze=False)
    for i, ds in enumerate(datasets):
        for j, key in enumerate(bins):
            ax = axs[i][j]
            f = fitres["bins"][key].get(f"{ds}_{wp_tag}", {})
            if "_hist" not in f:
                ax.text(0.5, 0.5, f.get("error", "brak"), ha="center", transform=ax.transAxes)
                continue
            xe = np.array(f["_edges"])
            xc = 0.5 * (xe[1:] + xe[:-1])
            h = np.array(f["_hist"])
            ax.errorbar(xc, h, np.sqrt(np.maximum(h, 1)), fmt="o", ms=2.5, color=INK, lw=0.8, label="dane")
            ax.plot(xc, f["_pred"], color=C1, label="fit")
            ax.plot(xc, f["_comp_bkg"], color=C2, ls="--", lw=1.4, label="tło")
            if f.get("_comp_refl"):
                ax.plot(xc, np.array(f["_comp_bkg"]) + np.array(f["_comp_refl"]), color=C3, ls=":", lw=1.4,
                        label="tło+refleksja")
            ax.set_title(f"{ds} {wp_tag} pT {key}: N={f['N_sig']:.0f}±{f['N_sig_err']:.0f}", fontsize=8)
            ax.set_xlabel("m [GeV]")
        axs[i][0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def make(cfg, species):
    out = cfg["paths"]["out_dir"]
    rd = os.path.join(out, "report")
    os.makedirs(rd, exist_ok=True)
    pt_bins = cfg["training"]["pt_bins"]
    oofs = {}
    for v in ("H", "frozen7"):
        p = os.path.join(out, "models", f"{species}_{v}", "oof.h5")
        if os.path.exists(p):
            oofs[v] = store.load(p)
    pngs = []
    if oofs:
        for v, oof in oofs.items():
            p = os.path.join(rd, f"{species}_{v}_scores.png")
            plot_scores(oof, pt_bins, p, f"{species} {v}")
            pngs.append(p)
        p = os.path.join(rd, f"{species}_roc.png")
        plot_roc(oofs, pt_bins, p, species)
        pngs.append(p)
    fitj = {}
    for v in ("H", "frozen7"):
        fr = _load(os.path.join(out, "fits", f"{species}_{v}.json"))
        if fr:
            fitj[v] = fr
            for tag in ("nocut", f"wp{int(cfg['fit']['wp']*100)}"):
                p = os.path.join(rd, f"{species}_{v}_fit_{tag}.png")
                plot_fits(fr, p, tag, cfg["fit"]["datasets"])
                pngs.append(p)
    val = _load(os.path.join(out, "validation", f"{species}.json"))
    md = markdown(cfg, species, val, fitj)
    with open(os.path.join(rd, f"{species}.md"), "w") as f:
        f.write(md)
    summ = {"species": species, "validation": val,
            "fits": {v: {k: {kk: {x: y for x, y in vv.items() if not x.startswith("_")} if isinstance(vv, dict) else vv
                             for kk, vv in e.items()} for k, e in fr["bins"].items()} for v, fr in fitj.items()},
            "pngs": pngs}
    with open(os.path.join(rd, f"{species}_summary.json"), "w") as f:
        json.dump(summ, f, indent=1, default=str)
    return md, pngs


def _f(x, n=3):
    try:
        return f"{x:.{n}g}"
    except Exception:  # noqa: BLE001
        return "—"


def sculpt_section(r, wps):
    """VAL: tabela rzeźbienia (X) i ε_B pod pikiem vs proxy (S8) + ręczna kontrola (W)."""
    mc, da = r.get("sculpt_mcbkg"), r.get("sculpt_data")
    if not mc and not da:
        return []
    L = ["", "### Rzeźbienie masy i ε_B pod pikiem (VAL)", ""]
    if mc:
        L += ["MC tła (bez sygnału), globalnie:", "",
              "| WP | R_win_window (stare) | R_win_sym | pol2_closure | slope ε_B(m) [1/GeV] | ε_B(SR) prawda | proxy SB_asym/prawda | proxy SB_sym/prawda |",
              "|---|---|---|---|---|---|---|---|"]
        for k, x in mc["all"].items():
            old = (r.get("R_win_mcbkg_window") or {}).get(k, {})
            L.append(f"| {k} | {_f(old.get('R'))}±{_f(old.get('err'),2)} | {_f(x.get('R_win_sym'))}±{_f(x.get('R_win_sym_err'),2)} | "
                     f"{_f(x.get('pol2_closure'))}±{_f(x.get('pol2_closure_err'),2)} | {_f(x.get('slope_perGeV'),2)} | "
                     f"{_f(x.get('eps_B_SR_true'),4)} | {_f(x.get('proxy_ratio_SB_asym'))} | {_f(x.get('proxy_ratio_SB_sym'))} |")
        if mc.get("bins"):
            L += ["", "per bin pT (R_win_sym / pol2_closure):", ""]
            for b, blk in mc["bins"].items():
                L.append(f"- {b}: " + ", ".join(f"@{k} {_f(x.get('R_win_sym'))} / {_f(x.get('pol2_closure'))}" for k, x in blk.items()))
    if da:
        L += ["", "Dane (pasmo wewnętrzne vs SB symetryczne; sygnał wyłączony):", "",
              "| WP | R_inner_sym | slope ε(m) SB [1/GeV] | ε(SB_asym) | ε(SB_sym) | ε(SB lewy) | ε(SB prawy) |", "|---|---|---|---|---|---|---|"]
        for k, x in da["all"].items():
            e = x["eps"]
            L.append(f"| {k} | {_f(x.get('R_inner_sym'))}±{_f(x.get('R_inner_sym_err'),2)} | {_f(x.get('slope_perGeV'),2)} | "
                     f"{_f(e['SB_asym']['eps'],4)} | {_f(e['SB_sym']['eps'],4)} | {_f(e['SB_left']['eps'],4)} | {_f(e['SB_right']['eps'],4)} |")
    hc = r.get("hand_control_data")
    if hc and "error" not in hc:
        L += ["", f"Ręczny dyskryminant (W, wariant {hc['variant']}, bez cech skorelowanych z m; próg dobrany do tej samej ε w SB_sym), n={hc['n']}:", "",
              "| WP | R_inner_sym XGB | R_inner_sym ręczny | slope XGB | slope ręczny |", "|---|---|---|---|---|"]
        for k in hc["hand"]:
            a, b = hc["xgb_same_rows"][k], hc["hand"][k]
            L.append(f"| {k} | {_f(a.get('R_inner_sym'))}±{_f(a.get('R_inner_sym_err'),2)} | {_f(b.get('R_inner_sym'))}±{_f(b.get('R_inner_sym_err'),2)} | "
                     f"{_f(a.get('slope_perGeV'),2)} | {_f(b.get('slope_perGeV'),2)} |")
    elif hc:
        L += ["", f"Ręczny dyskryminant: {hc['error']}"]
    eb = r.get("eps_B_SR")                                       # B4: nominał ε_B pod pikiem (pol2/SALAD) vs proxy SB
    if eb and "error" not in eb:
        L += ["", f"ε_B(SR) — nominał B4 (`{eb.get('method')}`; proxy SB symetryczne, asymetryczne = legacy):", "",
              "| WP | ε_B(SR) nominał | metoda | stat. wzgl. | SB_sym/nom | SB_asym(legacy)/nom | χ²/ndf fitu SB (all/pass) |",
              "|---|---|---|---|---|---|---|"]
        for k, x in eb["wps"].items():
            p2 = x.get("pol2", {})
            L.append(f"| {k} | {_f(x.get('nominal'),4)} | {x.get('method_used')} | {_f(x.get('stat_rel'),2)} | {_f(x.get('rel_sym_vs_nom'))} | "
                     f"{_f(x.get('rel_asym_vs_nom'))} | {_f(p2.get('chi2_ndf_all'),2)}/{_f(p2.get('chi2_ndf_pass'),2)} |")
    elif eb:
        L += ["", f"ε_B(SR) nominał B4: {eb['error']}"]
    sal = r.get("eps_B_SR_salad")
    if sal:
        L += ["", "ε_B(SR) SALAD: " + (sal.get("skipped") or sal.get("error") or str(sal))]
    L += ["", "Legenda: " + (r.get("sculpt_legend") or ""), ""]
    return L


def markdown(cfg, species, val, fitj):
    wps = cfg["training"]["working_points"]
    L = [f"# charmtag — {species}", ""]
    if val:
        for v, r in val["variants"].items():
            L += [f"## Wariant {v} ({len(r['features'])} cech)", "", "Cechy: " + ", ".join(r["features"]), ""]
            L += ["| pT [GeV] | AUC | " + " | ".join(f"ε_B@{int(w*100)}" for w in wps) + " | " +
                  " | ".join(f"R_win(dane)@{int(w*100)}" for w in wps) + " |",
                  "|---|---|" + "---|" * (2 * len(wps))]
            for key, m in r["metrics_oof"].items():
                rw = (r.get("R_win_data_inner_bins") or {}).get(key, {})
                L.append(f"| {key} | {m['auc']:.4f} | " + " | ".join(_f(m[f'effB_{int(w*100)}']) for w in wps) + " | " +
                         " | ".join(_f(rw.get(f"{int(w*100)}", {}).get("R")) for w in wps) + " |")
            L.append("")
            L.append("R_sb(dane) lewy/prawy sideband: " + ", ".join(f"@{k}: {_f(x['R'])}±{_f(x['err'],2)}" for k, x in (r.get("R_sb_data_lowhigh") or {}).items()))
            L.append("")
            L.append("R_win(dane) = [N(pasmo wewn.)/N(sidebandy)] po cięciu / przed (per bin pT; >1 może znaczyć ogon sygnału/refleksje w paśmie). Globalnie: " + ", ".join(f"@{k}: {_f(x['R'])}" for k, x in (r.get("R_win_data_inner") or {}).items()) + ". R_win(MC tło, okno sygnału): " +
                     ", ".join(f"@{k}: {_f(x['R'])}±{_f(x['err'],2)}" for k, x in (r.get("R_win_mcbkg_window") or {}).items()))
            L += sculpt_section(r, wps)   # VAL
            if r.get("dstar_data") and "error" not in r["dstar_data"]:
                L += ["", "Kontrola D0 z D* (ε przy WP, dane / MC):", "",
                      "| pT | " + " | ".join(f"ε@{int(w*100)} dane | MC" for w in wps) + " |",
                      "|---|" + "---|" * len(wps)]
                for key in [k for k in r["dstar_data"] if "-" in k]:
                    a, b = r["dstar_data"][key], (r.get("dstar_mc") or {}).get(key, {})
                    L.append(f"| {key} | " + " | ".join(
                        f"{_f(a.get(f'eff_{int(w*100)}'))}±{_f(a.get(f'err_{int(w*100)}'),2)} / {_f(b.get(f'eff_{int(w*100)}'))}"
                        for w in wps) + " |")
            if r.get("importance"):
                L += ["", "Ważność (gain): " + ", ".join(f"{k} {v:.3f}" for k, v in list(r["importance"].items())[:8])]
            L.append("")
        if val.get("H_vs_frozen7"):
            L += ["## H vs frozen7 (wzrost odrzucania 1/ε_B)", ""]
            for key, g in val["H_vs_frozen7"].items():
                L.append(f"- pT {key}: " + ", ".join(f"WP{k.split('_')[-1]}: ×{_f(x)}" for k, x in g.items()))
            L.append("")
    for v, fr in fitj.items():
        tag = f"wp{int(fr['wp']*100)}"
        L += [f"## Fit masy — {species} {v} ({tag})", "", "| pT | próbka | N_D ± σ | σ_eff [MeV] (σ1 MC) | χ²/ndf | ε_fit dane |",
              "|---|---|---|---|---|---|"]
        for key, e in fr["bins"].items():
            for ds in cfg["fit"]["datasets"]:
                f = e.get(f"{ds}_{tag}", {})
                if "N_sig" not in f:
                    L.append(f"| {key} | {ds} | {f.get('error', '—')} | | | |")
                    continue
                ef = e.get(f"{ds}_eff_fit", {})
                L.append(f"| {key} | {ds} | {f['N_sig']:.1f} ± {f['N_sig_err']:.1f} | {1e3*f['sigma_eff']:.1f} "
                         f"({1e3*f['sigma_mc']:.1f}) | {f['chi2']:.0f}/{f['ndf']} | {_f(ef.get('eff'))}±{_f(ef.get('err'),2)} |")
        L.append("")
    return "\n".join(L)
