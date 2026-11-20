#!/usr/bin/env python3
"""PREZ 19.11 — slajd 6 „Train on pairs of D mesons” + szkice topologii (projekt A).

(a) topologie par D–D z generatora Pythia 8 A3 (AFP_ddpred/gen/A3/gen_*.tsv.gz; prawda, bez detektora):
    pary D0/D⁺ w kanale analizy (mode = 1), pT 3–30 GeV, |η| < 2,5; OS (c + c̄) / LS (cc, c̄c̄ ≈ DPS) / pary niezależne
    (nogi z różnych zdarzeń = granica czystego DPS). Generator nie zapisuje procesu partonowego → podział OS/LS jako proxy.
(b) PAIRNET: architektura (Deep Sets), rozkład score i ROC na parach MC (pn_scores.npz, bieg v1 17.11, OOF), σ_μ/σ_nom z SBTRAIN v3.
(c) schematy: flavour creation (back-to-back), gluon splitting (kolinearne), DPS.
Użycie: make_p6.py --gen 'AFP_ddpred/gen/A3/gen_*.tsv.gz' --scores pn_scores.npz --sbtrain PAIRNET/runs/v3/sbtrain/sbtrain.json --out DIR"""
import argparse
import glob
import json
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle  # noqa: E402

plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "legend.fontsize": 10, "figure.dpi": 110})
CHN = {0: "D$^0$D$^0$", 1: "D$^+$D$^0$", 2: "D$^+$D$^+$"}
COL = {"os": "#1f5fa8", "ls": "#d1495b", "ind": "#555555", "all": "#000000"}


def atlas(ax, sim=True, extra=None, x=0.04, y=0.94):
    ax.text(x, y, "", transform=ax.transAxes, fontsize=13, va="top")
    if extra:
        ax.text(x, y - 0.07, extra, transform=ax.transAxes, fontsize=10.5, va="top")


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ (a) generator
def read_one(f):
    d = pd.read_csv(f, sep="\t", comment="#", header=None, names=["t", "ev", "pdg", "pt", "eta", "phi", "fb", "mode"], dtype={"t": str})
    d = d[d.t == "D"]
    pdg = d.pdg.astype(int).to_numpy(); apdg = np.abs(pdg)
    q = np.isin(apdg, [411, 421]) & (d["mode"].astype(int).to_numpy() == 1) & (d.pt.to_numpy() >= 3) & (d.pt.to_numpy() < 30) & (np.abs(d.eta.to_numpy()) < 2.5)
    return np.stack([d.ev.to_numpy()[q].astype(float), pdg[q].astype(float), d.pt.to_numpy()[q], d.eta.to_numpy()[q], d.phi.to_numpy()[q]], 1)


def pairs_same_event(L):
    """L: (fi, ev, pdg, pt, eta, phi) → tablica par (pt1, pt2, dphi, deta, ls, chan) z tego samego zdarzenia."""
    key = L[:, 0] * 1e9 + L[:, 1]
    o = np.argsort(key, kind="stable"); L = L[o]; key = key[o]
    _, st, cnt = np.unique(key, return_index=True, return_counts=True)
    rows = []
    for s, c in zip(st[cnt >= 2], cnt[cnt >= 2]):
        for i in range(s, s + c):
            for j in range(i + 1, s + c):
                rows.append(pair_row(L[i], L[j]))
    return np.array(rows)


def pair_row(a, b):
    pt1, pt2 = (a[3], b[3]) if a[3] >= b[3] else (b[3], a[3])
    dphi = abs((a[5] - b[5] + np.pi) % (2 * np.pi) - np.pi)
    ls = float(np.sign(a[2]) == np.sign(b[2]))
    ch = {(421, 421): 0, (411, 421): 1, (421, 411): 1, (411, 411): 2}[(int(abs(a[2])), int(abs(b[2])))]
    return (pt1, pt2, dphi, abs(a[4] - b[4]), ls, ch)


def pairs_indep(L, n, rng):
    i = rng.integers(0, len(L), n); j = rng.integers(0, len(L), n)
    k = (L[i, 0] != L[j, 0]) | (L[i, 1] != L[j, 1])
    return np.array([pair_row(L[a], L[b]) for a, b in zip(i[k], j[k])])


def plot_gen(P, I, out, nfiles):
    os_, ls_ = P[:, 4] == 0, P[:, 4] == 1
    info = {"n_pairs": int(len(P)), "n_os": int(os_.sum()), "n_ls": int(ls_.sum()), "f_LS": float(ls_.mean())}
    lab = f"Pythia 8 A3, $\\sqrt{{s}}$ = 13.6 TeV, generator level\nD$^0$/D$^+$ pairs, $p_T$ 3–30 GeV, $|\\eta|$ < 2.5"
    # Δφ
    eb = np.linspace(0, np.pi, 25)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for m, nm, c in ((os_, f"OS (c$\\bar{{c}}$), {os_.sum()} pairs", COL["os"]), (ls_, f"LS (cc, $\\bar{{c}}\\bar{{c}}$; DPS-like), {ls_.sum()} pairs", COL["ls"])):
        h, _ = np.histogram(P[m, 2], eb, density=True); e = np.sqrt(np.histogram(P[m, 2], eb)[0]) / max(m.sum(), 1) / np.diff(eb)
        ax.errorbar(0.5 * (eb[1:] + eb[:-1]), h, yerr=e, fmt="o", ms=4, color=c, label=nm)
    h, _ = np.histogram(I[:, 2], eb, density=True)
    ax.stairs(h, eb, color=COL["ind"], ls="--", lw=1.6, label="independent D's (pure DPS limit)")
    ax.set_xlabel(r"$|\Delta\phi|$ (D, D)"); ax.set_ylabel("normalised"); ax.set_xlim(0, np.pi); ax.set_ylim(0, ax.get_ylim()[1] * 1.45)
    ho, _ = np.histogram(P[os_, 2], eb, density=True)
    ax.annotate("gluon splitting\n(collinear)", xy=(0.07, ho[0]), xytext=(0.35, 0.12), fontsize=9.5, color="#333", arrowprops=dict(arrowstyle="->", color="#777"))
    ax.annotate("flavour creation\n(back-to-back)", xy=(3.07, ho[-1]), xytext=(2.0, 0.12), fontsize=9.5, color="#333", arrowprops=dict(arrowstyle="->", color="#777"))
    atlas(ax, extra=lab); ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.80), frameon=False)
    save(fig, out, "p6_gen_dphi")
    near, away = P[:, 2] < np.pi / 3, P[:, 2] > 2 * np.pi / 3
    info["os_frac_dphi_lt_pi3"] = float((near & os_).sum() / os_.sum()); info["os_frac_dphi_gt_2pi3"] = float((away & os_).sum() / os_.sum())
    info["ls_frac_dphi_lt_pi3"] = float((near & ls_).sum() / ls_.sum()); info["ls_frac_dphi_gt_2pi3"] = float((away & ls_).sum() / ls_.sum())
    # |Δη|
    eb = np.linspace(0, 5, 26)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for m, nm, c in ((os_, "OS (c$\\bar{c}$)", COL["os"]), (ls_, "LS (DPS-like)", COL["ls"])):
        h, _ = np.histogram(P[m, 3], eb, density=True); e = np.sqrt(np.histogram(P[m, 3], eb)[0]) / max(m.sum(), 1) / np.diff(eb)
        ax.errorbar(0.5 * (eb[1:] + eb[:-1]), h, yerr=e, fmt="o", ms=4, color=c, label=nm)
    h, _ = np.histogram(I[:, 3], eb, density=True); ax.stairs(h, eb, color=COL["ind"], ls="--", lw=1.6, label="independent D's (pure DPS limit)")
    ax.set_xlabel(r"$|\Delta\eta|$ (D, D)"); ax.set_ylabel("normalised"); ax.set_xlim(0, 5); ax.set_ylim(0, ax.get_ylim()[1] * 1.5)
    atlas(ax, extra=lab); ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.78), frameon=False)
    save(fig, out, "p6_gen_deta")
    for k, m in (("all", np.ones(len(P), bool)), ("os", os_), ("ls", ls_)):
        info[f"{k}_frac_deta_lt1"] = float((P[m, 3] < 1).mean())
    info["indep_frac_deta_lt1"] = float((I[:, 3] < 1).mean())
    # 2D pT1–pT2
    eb = np.linspace(3, 15, 25)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)
    for ax, m, nm in ((axs[0], os_, "OS (c$\\bar{c}$)"), (axs[1], ls_, "LS (DPS-like)")):
        h, xe, ye = np.histogram2d(P[m, 0], P[m, 1], [eb, eb])
        im = ax.pcolormesh(xe, ye, (h / h.sum()).T, norm=LogNorm(vmin=1e-4, vmax=0.1), cmap="viridis")
        ax.plot([3, 15], [3, 15], color="w", lw=0.8, ls=":")
        ax.set_xlabel(r"leading D $p_T$ [GeV]"); ax.set_title(nm, fontsize=12)
        atlas(ax, extra=f"Pythia 8 A3, generator level\n{m.sum()} pairs")
    axs[0].set_ylabel(r"subleading D $p_T$ [GeV]")
    fig.colorbar(im, ax=axs, label="fraction of pairs", pad=0.01)
    save(fig, out, "p6_gen_pt1pt2")
    info["os_mean_ptasym"] = float(np.mean((P[os_, 0] - P[os_, 1]) / (P[os_, 0] + P[os_, 1])))
    # 2D Δφ–Δη (OS)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    h, xe, ye = np.histogram2d(P[os_, 2], P[os_, 3], [np.linspace(0, np.pi, 19), np.linspace(0, 4, 17)])
    im = ax.pcolormesh(xe, ye, (h / h.sum()).T, cmap="magma_r")
    ax.set_xlabel(r"$|\Delta\phi|$"); ax.set_ylabel(r"$|\Delta\eta|$"); fig.colorbar(im, ax=ax, label="fraction of OS pairs")
    atlas(ax, extra="Pythia 8 A3, generator level, OS D pairs"); ax.set_ylim(0, 4)
    save(fig, out, "p6_gen_dphi_deta_os")
    info["nfiles"] = nfiles
    return info


# ------------------------------------------------------------------------------------------------ (b) PAIRNET
def plot_pairnet(scores, sbtrain, out):
    z = np.load(scores)
    y, smp, ch, win, nom = z["y"], z["smp"], z["chan"], z["inwin"], z["nom"]
    models = [("s_base_sum", "single-D BDT: $\\Sigma$ logit (nominal-like)", "#888888", "-"), ("s_xgb_leg", "pair BDT (2×18 leg features)", "#e08a00", "--"),
              ("s_pnet_leg", "pair NN (Deep Sets, legs)", "#1f5fa8", "-"), ("s_pnet_kin", "pair NN + pair kinematics", "#2a9d8f", ":")]
    info = {}
    # score sig vs bkg (pair NN), wszystkie kanały
    sig = (y == 1) & win; bkg = (y == 0) & (smp == 1) & win
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    s = z["s_pnet_leg"]; eb = np.linspace(0, 1, 41)
    ax.hist(s[sig], eb, density=True, histtype="stepfilled", alpha=0.35, color="#1f5fa8", label=f"signal: both legs true D ({sig.sum()})")
    ax.hist(s[bkg], eb, density=True, histtype="step", lw=1.8, color="#d1495b", label=f"background: ≥1 combinatorial leg ({bkg.sum()})")
    ax.set_yscale("log"); ax.set_xlabel("pair NN score (v1, background from MC)"); ax.set_ylabel("normalised"); ax.set_ylim(1e-2, ax.get_ylim()[1] * 30)
    atlas(ax, extra="minimum-bias + forced-decay MC, out-of-fold\nall channels, both legs |m − m$_0$| < 55 MeV"); ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.76), frameon=False)
    save(fig, out, "p6_pnet_score")
    # ROC per kanał
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
    for c, ax in zip((0, 1, 2), axs):
        m = (ch == c) & win
        sg, bg = m & (y == 1), m & (y == 0) & (smp == 1)
        info[f"chan{c}"] = {"n_sig": int(sg.sum()), "n_bkg": int(bg.sum())}
        for k, nm, col, ls in models:
            sc = z[k]
            thr = np.quantile(sc[sg], np.linspace(0.0, 0.99, 200))
            es = np.array([(sc[sg] >= t).mean() for t in thr]); eb_ = np.array([(sc[bg] >= t).mean() for t in thr])
            ax.plot(es, np.maximum(eb_, 1e-5), color=col, ls=ls, lw=1.8, label=nm)
            from sklearn.metrics import roc_auc_score  # noqa: E402
            yy = np.r_[np.ones(sg.sum()), np.zeros(bg.sum())]; ss = np.r_[sc[sg], sc[bg]]
            t50 = np.quantile(sc[sg], 0.5)
            info[f"chan{c}"][k] = {"auc": float(roc_auc_score(yy, ss)), "effB_at_effS50": float((sc[bg] >= t50).mean())}
        en_s, en_b = nom[sg].mean(), nom[bg].mean()
        ax.plot([en_s], [en_b], "k*", ms=13, label="nominal categories (both legs ≥ WP65)")
        info[f"chan{c}"]["nominal"] = {"effS": float(en_s), "effB": float(en_b)}
        ax.set_yscale("log"); ax.set_xlim(0.2, 1.0); ax.set_ylim(5e-4, 1); ax.set_xlabel(r"signal efficiency $\varepsilon_S$")
        atlas(ax, extra=f"{CHN[c]} pairs, MC (OOF)")
    axs[0].set_ylabel(r"background efficiency $\varepsilon_B$")
    hh, ll = axs[0].get_legend_handles_labels(); fig.legend(hh, ll, loc="lower center", ncol=5, fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.06))
    save(fig, out, "p6_pnet_roc")
    # σ_μ / σ_nom (SBTRAIN v3, sekcja C: tło z danych SB×SB, sygnał MC × SF)
    S = json.load(open(sbtrain))
    C = S["C"]["SF"]
    names = [("nomcat", "nominal\ncategories"), ("base_sum", "single-D BDT\nΣ logit"), ("pnet1", "pair NN\n(bkg: MC)"), ("xgbd", "pair BDT\n(bkg: data SB)"), ("pnetd", "pair NN\n(bkg: data SB)")]
    vals, errs = [], []
    for k, _ in names:
        vals.append(float(C["ratio"]["data"][k]["all"])); errs.append(float(C["ratio_err"][k]["all"]))
    info["sigma_ratio"] = dict(zip([k for k, _ in names], vals))
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    xs = np.arange(len(names)); cols = ["#888888", "#aaaaaa", "#7fa7d6", "#e08a00", "#1f5fa8"]
    ax.bar(xs, vals, yerr=errs, color=cols, width=0.62, capsize=4)
    for x_, v in zip(xs, vals):
        ax.text(x_, v + 0.035, f"{v:.2f}", ha="center", fontsize=11)
        ax.set_xticks(xs, [n for _, n in names], fontsize=10); ax.set_ylim(0, 1.5); ax.set_ylabel(r"expected $\sigma_\mu / \sigma_\mu^{\rm nominal}$")
    atlas(ax, sim=False, extra="428770: background from data sidebands (SB×SB),\nsignal MC × data SF; K = 6 categories")
    save(fig, out, "p6_pnet_sigma")
    return info


def plot_arch(out):
    fig, ax = plt.subplots(figsize=(11, 4.2)); ax.set_xlim(0, 11); ax.set_ylim(0, 4.2); ax.axis("off")

    def box(x, y, w, h, t, c):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc=c, ec="k", lw=1)); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=10)

    def arr(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14, lw=1.2, color="k"))
    box(0.1, 2.6, 2.0, 1.1, "D leg A\n18 features + species\n+ single-D BDT logit", "#e8eef7")
    box(0.1, 0.5, 2.0, 1.1, "D leg B\n18 features + species\n+ single-D BDT logit", "#e8eef7")
    box(2.75, 2.6, 2.4, 1.1, "shared leg encoder φ\nMLP 20→128→128→64", "#cfe0f3")
    box(2.75, 0.5, 2.4, 1.1, "shared leg encoder φ\n(same weights)", "#cfe0f3")
    box(5.7, 1.55, 1.9, 1.1, "symmetric pooling\nφ(A) + φ(B)\n|φ(A) − φ(B)|", "#fde9c9")
    box(8.2, 1.55, 2.6, 1.1, "head ρ: MLP\n128→128→128→1\n→ P(both legs true D)", "#d8efe9")
    arr(2.15, 3.15, 2.75, 3.15); arr(2.15, 1.05, 2.75, 1.05); arr(5.2, 3.15, 5.7, 2.4); arr(5.2, 1.05, 5.7, 1.8); arr(7.6, 2.1, 8.2, 2.1)
    ax.text(5.5, 3.95, "pair NN (Deep Sets): invariant under leg exchange; 5-fold out-of-fold; no pair kinematics in nominal variant", ha="center", fontsize=10.5)
    ax.text(8.2, 0.75, "optional: + pair kinematics\n(|Δφ|, |Δy|, ΔR, log p$_T$, p$_T$ asym.)", fontsize=9, color="#555")
    save(fig, out, "p6_pnet_arch")


# ------------------------------------------------------------------------------------------------ (c) schematy
def plot_schematics(out):
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))

    def arrow(ax, x0, y0, x1, y1, c, lab=None, lw=2.2, ls="-", wavy=False):
        if wavy:
            t = np.linspace(0, 1, 200); dx, dy = x1 - x0, y1 - y0; L = np.hypot(dx, dy); nx, ny = -dy / L, dx / L
            w = 0.05 * np.sin(2 * np.pi * 7 * t)
            ax.plot(x0 + dx * t + nx * w, y0 + dy * t + ny * w, color=c, lw=1.8)
        else:
            ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16, lw=lw, color=c, ls=ls))
        if lab:
            ax.text(x1 + 0.06 * np.sign(x1 - x0 + 1e-9), y1 + 0.06 * np.sign(y1 - y0 + 1e-9), lab, fontsize=13, color=c, ha="center", va="center")
    for ax in axs:
        ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.15, 1.15); ax.set_aspect("equal"); ax.axis("off")
        ax.add_patch(Circle((0, 0), 0.07, color="k"))
    # flavour creation
    ax = axs[0]; ax.set_title("flavour creation (gg → c$\\bar{c}$)\nback-to-back: |Δφ| ≈ π", fontsize=12)
    arrow(ax, 0, 0, 1.0, 0.25, "#1f5fa8", "c"); arrow(ax, 0, 0, -1.0, -0.25, "#d1495b", "$\\bar{c}$")
    ax.text(0, -0.85, "transverse plane", ha="center", fontsize=9.5, color="#666")
    # gluon splitting
    ax = axs[1]; ax.set_title("gluon splitting (g → c$\\bar{c}$)\ncollinear: |Δφ|, ΔR small", fontsize=12)
    arrow(ax, 0, 0, -1.0, -0.2, "#444", "g")
    arrow(ax, 0, 0, 0.45, 0.1, "#2a9d8f", wavy=True)
    arrow(ax, 0.45, 0.1, 1.05, 0.38, "#1f5fa8", "c"); arrow(ax, 0.45, 0.1, 1.07, 0.0, "#d1495b", "$\\bar{c}$")
    ax.text(0.2, 0.28, "g*", fontsize=12, color="#2a9d8f")
    # DPS
    ax = axs[2]; ax.set_title("double parton scattering (DPS)\ntwo independent c$\\bar{c}$ pairs → also cc (LS)", fontsize=12)
    arrow(ax, 0, 0, 0.95, 0.45, "#1f5fa8", "c"); arrow(ax, 0, 0, -0.95, -0.45, "#d1495b", "$\\bar{c}$")
    arrow(ax, 0, 0, -0.35, 1.0, "#1f5fa8", "c", ls="--"); arrow(ax, 0, 0, 0.35, -1.0, "#d1495b", "$\\bar{c}$", ls="--")
    ax.text(0.75, -0.75, "scattering 2", fontsize=9.5, color="#666"); ax.text(0.55, 0.62, "scattering 1", fontsize=9.5, color="#666")
    save(fig, out, "p6_schematic_topologies")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", required=True); ap.add_argument("--scores", required=True); ap.add_argument("--sbtrain", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--max-files", type=int, default=0); ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    files = sorted(glob.glob(a.gen))[: a.max_files or None]
    with Pool(a.threads) as P:
        parts = P.map(read_one, files)
    L = np.concatenate([np.c_[np.full(len(p), i, float), p] for i, p in enumerate(parts)])
    Pp = pairs_same_event(L)
    I = pairs_indep(L, 200000, np.random.default_rng(1))
    info = {"gen": plot_gen(Pp, I, a.out, len(files)), "pairnet": plot_pairnet(a.scores, a.sbtrain, a.out)}
    plot_arch(a.out); plot_schematics(a.out)
    json.dump(info, open(os.path.join(a.out, "p6_numbers.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
