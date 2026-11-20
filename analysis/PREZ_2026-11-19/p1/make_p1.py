#!/usr/bin/env python3
"""PREZ-1 (19.11.2026, projekt A): wykresy do prezentacji wg kartki — pkt 1 (definicje sygnał/tło), 2 (M_inv sygnał/tło, cutflow),
4 (score XGB H, GOOD/BAD, 2D), 5 (dane 428770 GOOD/BAD). Źródła:
  - scores nominalne v11_full (L_charmtag/runs/v11_full/scores/{D0_H,Dplus_H}/{signal_mc=fd_nd, incl_mc=mb_nd, data_none=428770}; pełna statystyka,
    po zakresie masy i pT 3–30; score OOF z 5 foldów),
  - surowe h5 v1.2 (FEAT_charmnet/data/{mb_nd 24, fd_nd 8} pliki) — początek cutflow i D* (decayType 1) w MC,
  - tabele B_trig bt_{data_428770,fd_nd}.npz — D⁰←D* z Δm i score nogi D⁰ (D* w danych i fd_nd).
Uruchomienie: .venv_h5/bin/python make_p1.py  (wyjście: ten katalog, PNG 150 dpi + PDF, numbers.json, cutflow.md)."""
import glob
import json
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
CW = os.path.abspath(os.path.join(OUT, "..", ".."))
SC = os.path.join(CW, "L_charmtag/runs/v11_full/scores")
RAW = os.path.join(CW, "FEAT_charmnet/data")
BT = os.path.join(CW, "B_trig/runs/lxplus")
PTB = [3.0, 5.0, 8.0, 30.0]
WPS = (0.80, 0.70, 0.50)
WP_GOOD = 0.70                                                   # nominalny WP fitu pojedynczych D (default.yaml fit.wp)
SPEC = {"D0": dict(dt=10, m0=1.86484, rng=(1.75, 2.00), win=0.055, lab=r"$D^0\to K^-\pi^+$", nb=50),
        "Dplus": dict(dt=4, m0=1.86966, rng=(1.70, 2.20), win=0.045, lab=r"$D^+\to K^-\pi^+\pi^+$", nb=50)}
# grupy kategorii prawdy (Dmeson_truth_category ↔ truth_channel CharmAnalysis; 100/101 = swap K↔π z charmtag.samples.mark_swaps)
GROUPS = [("Signal (prompt c + from b)", [0, 1], "#1f77b4"),
          (r"Reflection $K\leftrightarrow\pi$", [100, 101], "#9467bd"),
          ("Other charm / other decay mode", [2, 3, 6, 7], "#2ca02c"),
          ("Partially reco (≥1 daughter from c/b)", [8, 9, 10, 11], "#ff7f0e"),
          ("Combinatorial", [4, 5, 12, 13, 14, 15], "#bbbbbb")]
SMP = {"incl_mc": "MC minimum bias (mb_nd, Pythia A3)", "signal_mc": "MC forced decay (fd_nd)", "data_none": "Data 2022, run 428770"}
NUM = {}

plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "legend.fontsize": 9, "figure.dpi": 150})


def atlas(ax, sub=""):
    ax.text(0.03, 0.96, "", transform=ax.transAxes, va="top", fontsize=12)
    ax.text(0.03, 0.89, r"$\sqrt{s}$ = 13.6 TeV" + ("\n" + sub if sub else ""), transform=ax.transAxes, va="top", fontsize=9, linespacing=1.4)


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=150)
    plt.close(fig)


def load(sp, smp):
    h = h5py.File(os.path.join(SC, f"{sp}_H", f"{smp}.h5"), "r")
    d = {k: h[k][:] for k in ("m", "pt", "score", "truth", "veto_ok")}
    d["truth"] = np.nan_to_num(d["truth"], nan=-1).astype(int)
    return d


def wp_thresholds(sp):
    """Próg score dający ε_S = WP na sygnale MC (fd_nd, truth 0/1, zakres masy) per bin pT (definicja WP z treningu)."""
    s = load(sp, "signal_mc")
    q = np.isin(s["truth"], [0, 1]) & s["veto_ok"]
    th = {}
    for wp in WPS:
        th[wp] = [float(np.quantile(s["score"][q & (s["pt"] >= lo) & (s["pt"] < hi)], 1 - wp)) for lo, hi in zip(PTB[:-1], PTB[1:])]
    return th


def thr_of(pt, th):
    i = np.clip(np.digitize(pt, PTB) - 1, 0, len(PTB) - 2)
    return np.asarray(th)[i]


def group_masks(t):
    return [(n, np.isin(t, codes), c) for n, codes, c in GROUPS]


def stack_mass(ax, m, t, sp, w=None, title=""):
    S = SPEC[sp]; E = np.linspace(*S["rng"], S["nb"] + 1)
    gm = [(n, q, c) for n, q, c in group_masks(t) if q.any()]
    ax.hist([m[q] for _, q, _ in gm][::-1], E, stacked=True, color=[c for *_, c in gm][::-1], label=[n for n, *_ in gm][::-1], histtype="stepfilled", alpha=0.9)
    sig = np.isin(t, [0, 1]); fs = sig.sum() / max(len(t), 1)
    if 0 < fs < 0.05:                                                 # sygnał niewidoczny w stosie → dodatkowo nałożony ×k
        k = 10 ** int(np.floor(np.log10(0.25 / fs)))
    if 0 < fs < 0.05 and k > 1:
        hs, _ = np.histogram(m[sig], E); ax.stairs(hs * k, E, color="#1f77b4", lw=2, ls="--", label=f"Signal ×{k}")
    ax.set_xlabel(r"$m_{\rm inv}$ [GeV]"); ax.set_ylabel(f"Candidates / {1e3 * (E[1] - E[0]):.0f} MeV"); ax.set_xlim(*S["rng"])
    ax.set_ylim(0, ax.get_ylim()[1] * 1.75)
    h, l = ax.get_legend_handles_labels(); ax.legend(h[::-1], l[::-1], loc="upper right", frameon=False, fontsize=8.5)
    atlas(ax, f"{S['lab']}\n{title}")


def comp_counts(t, m, sp):
    S = SPEC[sp]; win = np.abs(m - S["m0"]) < S["win"]
    return {n: {"all": int(q.sum()), "window": int((q & win).sum())} for n, q, _ in group_masks(t)}


def main():
    TH = {sp: wp_thresholds(sp) for sp in SPEC}
    NUM["wp_thresholds"] = {sp: {str(k): v for k, v in TH[sp].items()} for sp in SPEC}
    cut_rows = {}
    for sp, S in SPEC.items():
        th = TH[sp][WP_GOOD]
        D = {smp: load(sp, smp) for smp in ("incl_mc", "signal_mc", "data_none")}
        for smp in ("incl_mc", "signal_mc", "data_none"):
            v = D[smp]["veto_ok"]
            for k in list(D[smp]):
                D[smp][k] = D[smp][k][v] if k != "veto_ok" else D[smp][k]
            D[smp]["good"] = D[smp]["score"] >= thr_of(D[smp]["pt"], th)
        # (2) masa: sygnał/tło stack, mb_nd i fd_nd
        for smp in ("incl_mc", "signal_mc"):
            d = D[smp]
            fig, ax = plt.subplots(figsize=(7, 5)); stack_mass(ax, d["m"], d["truth"], sp, title=SMP[smp] + "\n" + r"$3<p_T<30$ GeV, no MVA cut")
            save(fig, f"p2_mass_{sp}_{smp}")
            NUM[f"comp_{sp}_{smp}"] = comp_counts(d["truth"], d["m"], sp)
        # (4) score: sygnał (fd_nd truth 0/1) vs tło (mb_nd ≠ sygnał/refleksja), z progami WP
        sig = D["signal_mc"]; bkg = D["incl_mc"]
        qs = np.isin(sig["truth"], [0, 1]); qb = ~np.isin(bkg["truth"], [0, 1, 100, 101])
        fig, ax = plt.subplots(figsize=(7, 5)); E = np.linspace(0, 1, 51)
        ax.hist(sig["score"][qs], E, density=True, histtype="step", lw=2, color="#1f77b4", label="Signal (MC fd_nd, truth-matched)")
        ax.hist(bkg["score"][qb], E, density=True, histtype="step", lw=2, color="#d62728", label="Background (MC mb_nd, not signal)")
        for i, (lo, hi) in enumerate(zip(PTB[:-1], PTB[1:])):
            ax.axvline(th[i], ls=("--", "-.", ":")[i], color="k", lw=1.2, label=f"WP70 cut, {lo:.0f}<pT<{hi:.0f} GeV: {th[i]:.3f}")
        ax.set_yscale("log"); ax.set_xlabel("XGBoost H score"); ax.set_ylabel("Normalised"); ax.set_ylim(1e-2, 3e2)
        ax.text(0.86, 0.12, "BAD ← | → GOOD", transform=ax.transAxes, ha="center", fontsize=10)
        ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.62), frameon=False, fontsize=8.5); atlas(ax, S["lab"]); save(fig, f"p4_score_{sp}")
        NUM[f"score_{sp}"] = {"eff_S_good": float(sig["good"][qs].mean()), "eff_B_good": float(bkg["good"][qb].mean()),
                              "rejection": float(1 / max(bkg["good"][qb].mean(), 1e-12))}
        # (4) masa GOOD / BAD w MC (mb_nd i fd_nd)
        for smp in ("incl_mc", "signal_mc"):
            d = D[smp]
            fig, axs = plt.subplots(1, 2, figsize=(13, 5))
            for ax, (lab, q) in zip(axs, (("GOOD (score ≥ WP70)", d["good"]), ("BAD (score < WP70)", ~d["good"]))):
                stack_mass(ax, d["m"][q], d["truth"][q], sp, title=f"{SMP[smp]}\n{lab}")
            save(fig, f"p4_mass_goodbad_{sp}_{smp}")
            NUM[f"comp_good_{sp}_{smp}"] = comp_counts(d["truth"][d["good"]], d["m"][d["good"]], sp)
            NUM[f"comp_bad_{sp}_{smp}"] = comp_counts(d["truth"][~d["good"]], d["m"][~d["good"]], sp)
        # (4) 2D: score vs masa, sygnał i tło (mb_nd)
        d = D["incl_mc"]
        fig, axs = plt.subplots(1, 2, figsize=(13, 5))
        for ax, (lab, q) in zip(axs, (("Signal (truth-matched)", np.isin(d["truth"], [0, 1])), ("Background", ~np.isin(d["truth"], [0, 1, 100, 101])))):
            hh = ax.hist2d(d["m"][q], d["score"][q], bins=[np.linspace(*S["rng"], 51), np.linspace(0, 1, 41)], cmap="viridis", norm=matplotlib.colors.LogNorm())
            fig.colorbar(hh[3], ax=ax, label="Candidates")
            ax.axhline(np.median(th), color="w", ls="--", lw=1); ax.set_xlabel(r"$m_{\rm inv}$ [GeV]"); ax.set_ylabel("XGBoost H score")
            ax.set_title(f"{SMP['incl_mc']}: {lab}", fontsize=10); ax.text(0.03, 0.96, "", transform=ax.transAxes, va="top", color="w")
        save(fig, f"p4_2d_score_mass_{sp}_mc")
        # (5) dane: GOOD vs BAD (1D) i 2D
        d = D["data_none"]; E = np.linspace(*S["rng"], S["nb"] + 1)
        fig, axs = plt.subplots(1, 2, figsize=(13, 5))
        for ax, (lab, q, c) in zip(axs, (("GOOD (score ≥ WP70)", d["good"], "#1f77b4"), ("BAD (score < WP70)", ~d["good"], "#d62728"))):
            ax.hist(d["m"][q], E, histtype="step", lw=1.8, color=c, label=f"{lab}: {q.sum():,} cand.")
            ax.axvline(S["m0"], color="k", ls=":", lw=1); ax.set_xlabel(r"$m_{\rm inv}$ [GeV]"); ax.set_ylabel(f"Candidates / {1e3 * (E[1] - E[0]):.0f} MeV")
            ax.set_ylim(0, ax.get_ylim()[1] * 1.35); ax.legend(loc="upper right", frameon=False); atlas(ax, f"{S['lab']}\n{SMP['data_none']}")
        save(fig, f"p5_data_mass_goodbad_{sp}")
        fig, ax = plt.subplots(figsize=(7, 5))
        hh = ax.hist2d(d["m"], d["score"], bins=[np.linspace(*S["rng"], 51), np.linspace(0, 1, 41)], cmap="viridis", norm=matplotlib.colors.LogNorm())
        fig.colorbar(hh[3], ax=ax, label="Candidates"); ax.axhline(np.median(th), color="w", ls="--", lw=1)
        ax.set_xlabel(r"$m_{\rm inv}$ [GeV]"); ax.set_ylabel("XGBoost H score"); ax.set_title(f"{SMP['data_none']}: {S['lab']}", fontsize=10)
        ax.text(0.03, 0.96, "", transform=ax.transAxes, va="top", color="w"); save(fig, f"p5_data_2d_score_mass_{sp}")
        win = np.abs(d["m"] - S["m0"]) < S["win"]; sb = np.abs(d["m"] - S["m0"]) > 2 * S["win"]
        wfrac = S["win"] * 2 / (S["rng"][1] - S["rng"][0] - 4 * S["win"])
        NUM[f"data_{sp}"] = {"n_all": int(len(d["m"])), "n_good": int(d["good"].sum()), "n_bad": int((~d["good"]).sum()),
                             "good_window_over_sb_scaled": float((d["good"] & win).sum() / max((d["good"] & sb).sum() * wfrac, 1)),
                             "bad_window_over_sb_scaled": float((~d["good"] & win).sum() / max((~d["good"] & sb).sum() * wfrac, 1))}
        # cutflow (pełna statystyka, scores): wejście klasyfikatora → okno → WP
        rows = []
        for smp in ("incl_mc", "signal_mc"):
            z = load(sp, smp); s = np.isin(z["truth"], [0, 1]); b = ~np.isin(z["truth"], [0, 1, 100, 101]); w = np.abs(z["m"] - S["m0"]) < S["win"]
            steps = [("mass range + 3<pT<30 GeV (classifier input)", np.ones(len(s), bool))]
            if sp == "Dplus":
                steps.append((r"D* veto (min(m − m_Kπ) ≥ 160 MeV)", z["veto_ok"].astype(bool)))
            acc = steps[-1][1].copy(); steps.append((f"signal window |m − m0| < {1e3 * S['win']:.0f} MeV", acc & w))
            for wp in WPS:
                steps.append((f"score ≥ WP{int(wp * 100)}", acc & w & (z["score"] >= thr_of(z["pt"], TH[sp][wp]))))
            rows.append((smp, [(n, int((q & s).sum()), int((q & b).sum())) for n, q in steps]))
        cut_rows[sp] = rows
    # D*: Δm (bt: D⁰←D*), fd_nd z prawdą D* i D⁰, dane GOOD/BAD wg score nogi D⁰ (WP90 nogi jak w parach)
    WP90 = 0.38729336857795715
    dmE = np.linspace(139.5, 170, 62)
    for nm in ("fd_nd", "data_428770"):
        z = np.load(os.path.join(BT, f"bt_{nm}.npz"))
        q = (np.abs(z["m"] - SPEC["D0"]["m0"]) < 0.055) & (z["pt"] >= 3) & (z["pt"] < 30)
        dm = z["dm"][q] * 1e3; sc = z["score"][q]
        fig, ax = plt.subplots(figsize=(7, 5))
        if nm == "fd_nd":
            sigq = np.isin(z["truth"][q], [0, 1]) & np.isin(z["truth_ds"][q], [0, 1])
            ax.hist([dm[~sigq], dm[sigq]], dmE, stacked=True, color=["#bbbbbb", "#1f77b4"], label=["Background (not true D*)", "Signal (true D* and D⁰)"], histtype="stepfilled")
            NUM["dstar_fd_nd"] = {"n": int(len(dm)), "n_sig": int(sigq.sum())}
            sub = "MC forced decay (fd_nd)"
        else:
            g = sc >= WP90
            ax.hist(dm[g], dmE, histtype="step", lw=1.8, color="#1f77b4", label=f"D⁰ leg GOOD (score ≥ WP90): {g.sum():,}")
            ax.hist(dm[~g], dmE, histtype="step", lw=1.8, color="#d62728", label=f"D⁰ leg BAD: {(~g).sum():,}")
            NUM["dstar_data"] = {"n": int(len(dm)), "n_good": int(g.sum())}
            sub = SMP["data_none"]
        ax.set_xlabel(r"$\Delta m = m(K\pi\pi_s) - m(K\pi)$ [MeV]"); ax.set_ylabel("Candidates / 0.5 MeV"); ax.set_xlim(139.5, 170)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.6); ax.legend(loc="upper right", frameon=False)
        atlas(ax, r"$D^{*+}\to D^0(K\pi)\pi^+_s$" + "\n" + sub + "\n" + r"$|m_{K\pi}-m_{D^0}|<55$ MeV, $3<p_T<30$ GeV"); save(fig, f"p2_dstar_dm_{nm}")
    # cutflow wstępny z surowych h5 (lokalne pliki v1.2) + D* w mb_nd
    pre = {}
    for smp, pat in (("mb_nd", "mb_nd/*.h5"), ("fd_nd", "fd_nd/*.h5")):
        fs = sorted(glob.glob(os.path.join(RAW, pat)))
        acc = {sp: np.zeros((4, 2), np.int64) for sp in ("D0", "Dplus", "Dstar")}
        dstar_dm, dstar_t = [], []
        for f in fs:
            ev = h5py.File(f, "r")["events"]
            dt = ev["Dmeson_decayType"][:]; m = ev["Dmeson_mass"][:]; pt = ev["Dmeson_pt"][:]; t = np.nan_to_num(ev["Dmeson_truth_category"][:], nan=-1)
            for sp, code in (("D0", 10), ("Dplus", 4), ("Dstar", 1)):
                q0 = dt == code; s = np.isin(t, [0, 1])
                if sp == "Dstar":
                    q1 = q0; q2 = q1 & (pt >= 3) & (pt < 30)
                    dm_ = ev["Dmeson_delta_mass"][:]
                    dstar_dm.append(dm_[q2]); dstar_t.append(t[q2])
                else:
                    lo, hi = SPEC[sp]["rng"]; q1 = q0 & (m >= lo) & (m < hi); q2 = q1 & (pt >= 3) & (pt < 30)
                for i, q in enumerate((q0, q1, q2)):
                    acc[sp][i] += [int((q & s).sum()), int((q & ~s).sum())]
        pre[smp] = {"files": len(fs), **{sp: acc[sp][:3].tolist() for sp in acc}}
        if smp == "mb_nd":
            dm = np.concatenate(dstar_dm); t = np.concatenate(dstar_t)
            dm = dm * (1e3 if np.nanmax(dm) < 5 else 1.0)
            fig, ax = plt.subplots(figsize=(7, 5)); s = np.isin(t, [0, 1])
            ax.hist([dm[~s], dm[s]], dmE, stacked=True, color=["#bbbbbb", "#1f77b4"], label=["Background", "Signal (truth-matched D*)"], histtype="stepfilled")
            ax.set_xlabel(r"$\Delta m$ [MeV]"); ax.set_ylabel("Candidates / 0.5 MeV"); ax.set_xlim(139.5, 170); ax.set_ylim(0, ax.get_ylim()[1] * 1.6)
            ax.legend(loc="upper right", frameon=False); atlas(ax, r"$D^{*+}\to D^0\pi^+_s$" + "\n" + f"MC mb_nd ({len(fs)} files)" + "\n" + r"$3<p_T<30$ GeV")
            save(fig, "p2_dstar_dm_mb_nd")
            NUM["dstar_mb_nd_local"] = {"files": len(fs), "n": int(len(dm)), "n_sig": int(s.sum())}
    NUM["cutflow_raw_local"] = pre
    # cutflow.md
    md = ["# Cutflow (S = truth-matched c/b, B = reszta poza refleksjami K↔π)", "",
          "## Surowe ntuple (finder CharmAnalysis = preselekcja), lokalne pliki h5 v1.2 (mb_nd 24, fd_nd 8)", "",
          "| próbka | gatunek | krok | S | B |", "|---|---|---|---|---|"]
    names = ["candidates in ntuple (finder preselection)", "mass range", "3 < pT < 30 GeV"]
    for smp, P in pre.items():
        for sp in ("D0", "Dplus", "Dstar"):
            for i, (s, b) in enumerate(P[sp]):
                n = names[i] if not (sp == "Dstar" and i == 1) else "(no m cut; Δm range in fit)"
                md.append(f"| {smp} | {sp} | {n} | {s:,} | {b:,} |")
    md += ["", "## Pełna statystyka (scores v11_full; wejście klasyfikatora = zakres masy i pT już nałożone)", ""]
    for sp, rows in cut_rows.items():
        for smp, steps in rows:
            md += [f"### {sp}, {SMP[smp]}", "", "| krok | S | B | S/B |", "|---|---|---|---|"]
            md += [f"| {n} | {s:,} | {b:,} | {s / max(b, 1):.3g} |" for n, s, b in steps]
            md.append("")
    open(os.path.join(OUT, "cutflow.md"), "w").write("\n".join(md) + "\n")
    NUM["cutflow_full"] = {sp: {smp: steps for smp, steps in rows} for sp, rows in cut_rows.items()}
    json.dump(NUM, open(os.path.join(OUT, "numbers.json"), "w"), indent=1)
    print("OK", OUT)


if __name__ == "__main__":
    main()
