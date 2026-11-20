#!/usr/bin/env python3
"""PREZ-2 / pz_report (19.11.2026, projekt A) — porównanie 4 klasyfikatorów pojedynczych D (XGB/sieć × high/low-level) + nominalny XGB H (produkcja).

Na OOF w oknie sygnału (wagi treningowe: tło w (pT, |η|) jak sygnał, klasy zrównoważone):
  AUC (± rozrzut foldów/√k), ε_B przy ε_S = 90 % (WP90 nominału), ε_S przy ε_B nominału w jego WP90 (stały poziom tła), zysk odrzucenia ε_B(nom)/ε_B,
  rzeźbienie: ε_B w oknie / ε_B poza oknem (tło mb, przy progu ε_S = 90 %) — 1 = brak zależności od masy.
Wykresy (PNG + PDF): roc_<sp>, scores_<sp>, importance_<sp>, sculpt_<sp>; tabela TABLE.md + table.json.
Użycie: pz_report.py --in DIR_D0 --in DIR_Dplus --out DIR    (DIR_* zawiera scores_xgb.npz, scores_net.npz, imp_*.json, info_*.json)"""
import argparse
import json
import os
import sys

import numpy as np

VARS = [("xgb", "hi", "BDT (XGBoost), high-level", "C0", "-"), ("net", "hi", "NN (MLP), high-level", "C1", "-"),
        ("xgb", "lo", "BDT (XGBoost), low-level", "C0", "--"), ("net", "lo", "NN (Deep Sets), low-level", "C1", "--")]
NOM = ("nom", "nom", "BDT H (nominal, trained vs data SB)", "k", ":")
SPN = {"D0": r"$D^0\to K^-\pi^+$", "Dplus": r"$D^+\to K^-\pi^+\pi^+$"}
WP_EFF = 0.90


def atlas(ax, sp, extra=""):
    ax.text(0.04, 0.95, "ATLAS", transform=ax.transAxes, fontweight="bold", fontstyle="italic", fontsize=12, va="top")
    ax.text(0.20, 0.95, "", transform=ax.transAxes, fontsize=12, va="top")
    ax.text(0.04, 0.88, r"Pythia8 A3, $\sqrt{s}$ = 13.6 TeV, " + SPN.get(sp, sp) + extra, transform=ax.transAxes, fontsize=9, va="top")


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=150 if ext == "png" else None, bbox_inches="tight")


def roc(s, y, w):
    ok = np.isfinite(s); s, y, w = s[ok], y[ok], w[ok]
    o = np.argsort(-s, kind="stable"); s, y, w = s[o], y[o], w[o]
    tp = np.cumsum(w * (y == 1)); fp = np.cumsum(w * (y == 0))
    return np.r_[0, fp / fp[-1]], np.r_[0, tp / tp[-1]], np.r_[np.inf, s]


def auc(s, y, w):
    f, t, _ = roc(s, y, w)
    return float(np.trapezoid(t, f) if hasattr(np, "trapezoid") else np.trapz(t, f))


def thr_at_eff(s, w, eff):
    """Próg, przy którym ważona wydajność (s ≥ thr) = eff."""
    o = np.argsort(-s); c = np.cumsum(w[o]) / w.sum()
    return float(s[o][min(np.searchsorted(c, eff), len(o) - 1)])


def effw(s, w, t):
    return float(w[s >= t].sum() / max(w.sum(), 1e-12))


def metrics(s, y, w, inwin, fold, thr_nom_b=None):
    """thr_nom_b: ε_B nominału w WP90 (stały poziom tła) → ε_S modelu przy tym ε_B."""
    m = inwin & np.isfinite(s)
    S, B = m & (y == 1), m & (y == 0)
    a = auc(s[m], y[m], w[m])
    fa = [auc(s[m & (fold == k)], y[m & (fold == k)], w[m & (fold == k)]) for k in np.unique(fold)]
    t90 = thr_at_eff(s[S], w[S], WP_EFF); eb90 = effw(s[B], w[B], t90)
    out = {"auc": a, "auc_err": float(np.std(fa) / np.sqrt(len(fa))), "thr_eS90": t90, "eB_at_eS90": eb90}
    if thr_nom_b is not None:
        tb = thr_at_eff(s[B], w[B], thr_nom_b); out["eS_at_eB_nom"] = effw(s[S], w[S], tb)
    outb = (~inwin) & (y == 0) & np.isfinite(s)                                     # rzeźbienie: tło poza oknem (wagi = 1, kształt masy)
    inb = inwin & (y == 0) & np.isfinite(s)
    if outb.sum() > 100:
        e_in = float((s[inb] >= t90).mean()); e_out = float((s[outb] >= t90).mean())
        out["sculpt_ratio"] = e_in / e_out if e_out > 0 else np.nan
    return out


def load_sp(d):
    X = {f: np.load(os.path.join(d, f"scores_{f}.npz")) for f in ("xgb", "net") if os.path.exists(os.path.join(d, f"scores_{f}.npz"))}
    base = next(iter(X.values()))
    for f, z in X.items():
        assert np.array_equal(z["y"], base["y"]) and np.array_equal(z["m"], base["m"]), "niezgodne wiersze między rodzinami"
    imp = {f: json.load(open(os.path.join(d, f"imp_{f}.json"))) for f in X if os.path.exists(os.path.join(d, f"imp_{f}.json"))}
    info = {f: json.load(open(os.path.join(d, f"info_{f}.json"))) for f in X if os.path.exists(os.path.join(d, f"info_{f}.json"))}
    S = {(f, i): X[f][f"s_{i}"] for f in X for i in ("hi", "lo") if f"s_{i}" in X[f].files}
    if np.isfinite(base["s_nom"]).any():
        S[("nom", "nom")] = base["s_nom"]
    return base, S, imp, info


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="ins", action="append", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(a.out, exist_ok=True)
    table, md = {}, ["# PREZ-2 — single-D classifier: BDT vs NN × high-level vs low-level inputs (MC, OOF, signal window)", "",
                     "Signal: forced-decay MC, truth-matched D; background: minimum-bias MC combinatorics (truth 8–99); |m − m₀| < window; "
                     "background reweighted to signal (pT, |η|), classes balanced. WP90 = ε_S 90 %. „ε_S @ ε_B(nom)” = signal efficiency at the background "
                     "efficiency of the nominal BDT H at its WP90. Sculpting = ε_B(in window)/ε_B(outside window) at WP90 (1 = mass-independent).", ""]
    for d in a.ins:
        base, S, imp, info = load_sp(d)
        sp = (info.get("xgb") or info.get("net"))["species"]
        y, w, inwin, fold, m = base["y"], base["w"], base["inwin"], base["fold"], base["m"]
        nomb = None
        if ("nom", "nom") in S:
            sn = S[("nom", "nom")]; mm = inwin & np.isfinite(sn)
            t = thr_at_eff(sn[mm & (y == 1)], w[mm & (y == 1)], WP_EFF); nomb = effw(sn[mm & (y == 0)], w[mm & (y == 0)], t)
        rows = {}
        for v in VARS + [NOM]:
            key = (v[0], v[1])
            if key in S:
                rows[v[2]] = metrics(S[key], y, w, inwin, fold, nomb)
        nparam = {}
        for f in info:
            for i, fl in info[f]["folds"].items():
                k0 = fl[0]
                nparam[(f, i)] = k0.get("n_params") or k0.get("n_trees")
        table[sp] = {"rows": rows, "n_sig_win": int(((y == 1) & inwin).sum()), "n_bkg_win": int(((y == 0) & inwin).sum()),
                     "eB_nom_WP90": nomb, "size": {f"{k[0]}_{k[1]}": v for k, v in nparam.items()}}
        md += [f"## {sp}  (signal {table[sp]['n_sig_win']}, background {table[sp]['n_bkg_win']} candidates in window)", "",
               "| model | size | AUC | ε_B @ ε_S 90 % | rejection gain vs nominal | ε_S @ ε_B(nom) | sculpting |", "|---|---|---|---|---|---|---|"]
        for v in VARS + [NOM]:
            r = rows.get(v[2])
            if not r:
                continue
            sz = nparam.get((v[0], v[1]), "–")
            sz = f"{sz} {'params' if v[0] == 'net' else 'trees'}" if sz != "–" else "–"
            gain = f"{nomb / r['eB_at_eS90']:.2f}" if nomb and r["eB_at_eS90"] > 0 else "–"
            md.append(f"| {v[2]} | {sz} | {r['auc']:.4f} ± {r['auc_err']:.4f} | {r['eB_at_eS90']:.4f} | {gain} | "
                      f"{r.get('eS_at_eB_nom', float('nan')):.3f} | {r.get('sculpt_ratio', float('nan')):.2f} |")
        md.append("")
        # --- ROC
        fig, ax = plt.subplots(figsize=(6, 5))
        for v in VARS + [NOM]:
            if (v[0], v[1]) not in S:
                continue
            mm = inwin & np.isfinite(S[(v[0], v[1])])
            f_, t_, _ = roc(S[(v[0], v[1])][mm], y[mm], w[mm])
            ax.plot(t_, 1 / np.clip(f_, 1e-6, None), color=v[3], ls=v[4], label=f"{v[2]} (AUC {rows[v[2]]['auc']:.3f})")
        ax.set_yscale("log"); ax.set_xlabel(r"Signal efficiency $\varepsilon_S$"); ax.set_ylabel(r"Background rejection $1/\varepsilon_B$")
        ax.set_xlim(0.3, 1.0); ax.axvline(WP_EFF, color="grey", lw=0.8, ls=":"); ax.legend(fontsize=7, loc="lower left"); atlas(ax, sp)
        ax.set_ylim(1, ax.get_ylim()[1] * 5); save(fig, a.out, f"roc_{sp}"); plt.close(fig)
        # --- rozkłady score
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        for ax, v in zip(axs.flat, VARS):
            if (v[0], v[1]) not in S:
                ax.set_visible(False); continue
            s = S[(v[0], v[1])]; e = np.linspace(0, 1, 51)
            for cls, lab, col in ((1, "signal", "C0"), (0, "background", "C3")):
                mm = inwin & (y == cls) & np.isfinite(s)
                ax.hist(s[mm], e, weights=w[mm] / w[mm].sum(), histtype="step", color=col, label=lab, lw=1.5)
            ax.axvline(rows[v[2]]["thr_eS90"], color="grey", ls=":", lw=1, label="WP90")
            ax.set_yscale("log"); ax.set_xlabel("Classifier score"); ax.set_ylabel("Normalised"); ax.set_title(v[2], fontsize=10)
            ax.legend(fontsize=8, loc="upper center"); atlas(ax, sp)
            ax.set_ylim(top=ax.get_ylim()[1] * 30)
        fig.tight_layout(); save(fig, a.out, f"scores_{sp}"); plt.close(fig)
        # --- ważność (high-level)
        names = (info.get("xgb") or info.get("net"))["meta_sig"]["features_hi"]
        fig, ax = plt.subplots(figsize=(7, 6))
        yy = np.arange(len(names)); hb = 0.38
        for off, f, lab, col in ((-hb / 2, "xgb", "BDT (XGBoost)", "C0"), (hb / 2, "net", "NN (MLP)", "C1")):
            pi = (imp.get(f, {}).get("hi") or {}).get("perm_dauc")
            if pi:
                err = imp[f]["hi"].get("perm_dauc_err", {})
                ax.barh(yy + off, [pi[n] for n in names], hb, xerr=[err.get(n, 0) for n in names], color=col, label=lab)
        ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=8); ax.invert_yaxis(); ax.set_xlabel(r"Permutation importance: $\Delta$AUC")
        ax.legend(fontsize=8, loc="lower right"); atlas(ax, sp, ", high-level inputs")
        ax.set_ylim(len(names) + 2.5, -1); save(fig, a.out, f"importance_{sp}"); plt.close(fig)
        sh = (imp.get("xgb", {}).get("hi") or {}).get("shap_mean_abs")
        if sh:
            o = sorted(names, key=lambda n: -sh[n])
            fig, ax = plt.subplots(figsize=(6, 6)); ax.barh(np.arange(len(o)), [sh[n] for n in o], color="C0")
            ax.set_yticks(np.arange(len(o))); ax.set_yticklabels(o, fontsize=8); ax.invert_yaxis(); ax.set_xlabel("mean |SHAP| (log-odds)")
            atlas(ax, sp, ", BDT high-level"); ax.set_ylim(len(o) + 2.5, -1); save(fig, a.out, f"shap_{sp}"); plt.close(fig)
        # --- rzeźbienie: masa tła przed/po cięciu WP90
        lo_, hi_ = np.nanmin(m), np.nanmax(m)
        e = np.linspace(lo_, hi_, 51)
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        for ax, v in zip(axs.flat, VARS):
            if (v[0], v[1]) not in S:
                ax.set_visible(False); continue
            s = S[(v[0], v[1])]; bb = (y == 0) & np.isfinite(s); t = rows[v[2]]["thr_eS90"]
            h0, _ = np.histogram(m[bb], e); h1, _ = np.histogram(m[bb & (s >= t)], e)
            ax.step(e[:-1], h0 / max(h0.sum(), 1), where="post", color="grey", label="background, all")
            ax.step(e[:-1], h1 / max(h1.sum(), 1), where="post", color="C3", label="background, pass WP90 (GOOD)")
            ms_ = (info.get("xgb") or info.get("net"))["meta_sig"]
            ax.axvspan(ms_["m0"] - ms_["signal_window"], ms_["m0"] + ms_["signal_window"], color="C0", alpha=0.08, label="signal window (training)")
            ax.set_xlabel("m [GeV]"); ax.set_ylabel("Normalised"); ax.set_title(f"{v[2]} — sculpting {rows[v[2]].get('sculpt_ratio', float('nan')):.2f}", fontsize=9)
            ax.legend(fontsize=7, loc="lower center"); atlas(ax, sp); ax.set_ylim(0, ax.get_ylim()[1] * 1.6)
        fig.tight_layout(); save(fig, a.out, f"sculpt_{sp}"); plt.close(fig)
    json.dump(table, open(os.path.join(a.out, "table.json"), "w"), indent=1, default=float)
    open(os.path.join(a.out, "TABLE.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    sys.exit(main())
