#!/usr/bin/env python3
"""UNB-DRY: raport po odsłonięciu par D–D. Wejście = katalog wyników unblind po fetch (EOS charm_q/pf_pairs/unblind → Mac).
Użycie:
  fill_unblind.py --ub <UB_DIR> --pf <PF_DIR> [--ysyst Y.json] [--salad S8.json] [--template unblind_report_template.md] --out <prefix> [--no-plots]
UB_DIR/pairs/: fit2d_unblind.json (wymagany), pairs_full.npz + pair_inputs.{json,npz} (do wykresów rzutów; bez nich — tylko tabele).
PF_DIR: compare_pf.json, pairs/syst/syst_blind.json (budżet U), pairs/fit2d_toys.json.
Nigdy nie zgaduje: brak pliku/klucza → „—”."""
import argparse, json, os, sys, time
import numpy as np
CH = ("D0D0", "DpD0", "DpDp"); NAME = {"D0D0": "D⁰D⁰", "DpD0": "D⁺D⁰", "DpDp": "D⁺D⁺"}
SEPT = {"D0D0": {"mu": 1.37, "Z": 1.9}, "DpD0": {"mu": 1.62, "Z": 2.9}, "DpDp": {"mu": 1.57, "Z": 1.0}, "combined_Z": 3.6, "mu_common": 1.50}
def jl(p):
    try: return json.load(open(p))
    except Exception: return None
def g(d, *ks):
    for k in ks:
        if not isinstance(d, dict) or k not in d: return None
        d = d[k]
    return d
def f(x, nd=2):
    try: return "—" if x is None else f"{float(x):.{nd}f}"
    except Exception: return str(x)
def pval(z):
    from math import erfc, sqrt
    return None if z is None else 0.5 * erfc(float(z) / sqrt(2))

def plots(ubp, out, src):
    """Rzuty m1/m2 (suma po kategoriach) dane vs model (tło μ=0 + sygnał) z ponownego fitu per kanał."""
    sys.path.insert(0, src)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from charmtag.pairs.model import Model, load_inputs
    I = load_inputs(os.path.join(ubp, "pair_inputs")); M = Model(I); H = dict(np.load(os.path.join(ubp, "pairs_full.npz")))
    files = {}
    for ch in CH:
        n = H[ch]; cf = M.channel(ch); full = np.ones(n.shape[1:], bool)
        m = cf.fit(n, full); x = np.array(m.values); lam = cf.lam(x)
        x0 = x.copy(); x0[cf.names.index("mu")] = 0.0; bkg = cf.lam(x0)
        a, b = cf.a, cf.b; e1, e2 = np.asarray(I.edges[str(a)] if str(a) in I.edges else I.edges[a]), np.asarray(I.edges[str(b)] if str(b) in I.edges else I.edges[b])
        e1 = np.linspace(*e1[:2], int(e1[2]) + 1) if len(e1) == 3 else e1; e2 = np.linspace(*e2[:2], int(e2[2]) + 1) if len(e2) == 3 else e2
        fig, axs = plt.subplots(2, 2, figsize=(10, 6.5), gridspec_kw={"height_ratios": [2, 1]}, sharex="col"); ax = axs[0]
        for k, (axis, e, lab) in enumerate(((2, e1, "m₁"), (1, e2, "m₂"))):
            d = n.sum(axis=(0, axis)); L = lam.sum(axis=(0, axis)); B = bkg.sum(axis=(0, axis)); c = 0.5 * (e[1:] + e[:-1])
            ax[k].errorbar(c, d, yerr=np.sqrt(np.clip(d, 1, None)), fmt="ko", ms=3, label="dane")
            ax[k].step(c, L, where="mid", color="C0", label="fit (S+B)"); ax[k].step(c, B, where="mid", color="C1", ls="--", label="tło (μ=0)")
            ax[k].set_ylabel("pary / bin"); ax[k].set_title(f"{NAME[ch]} — rzut {lab}")
            r = axs[1][k]; r.errorbar(c, d - B, yerr=np.sqrt(np.clip(d, 1, None)), fmt="ko", ms=3); r.step(c, L - B, where="mid", color="C3", label="sygnał z fitu")
            r.axhline(0, color="gray", lw=0.8); r.set_xlabel(f"{lab} [MeV]"); r.set_ylabel("dane − tło")
            if k == 0: r.legend(fontsize=8)
        ax[0].legend(fontsize=8); fig.tight_layout(); p = f"{out}_{ch}.png"; fig.savefig(p, dpi=130); plt.close(fig); files[ch] = os.path.basename(p)
    return files

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ub", required=True); ap.add_argument("--pf", required=True); ap.add_argument("--ysyst"); ap.add_argument("--salad")
    ap.add_argument("--template", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "unblind_report_template.md"))
    ap.add_argument("--out", required=True); ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--src", default=os.path.expanduser("~/mnt/charm/charm_work/Q_lxplus/pf_pairs/src") if os.path.isdir(os.path.expanduser("~/mnt/charm")) else os.path.expanduser("~/Documents/charm/charm_work/Q_lxplus/pf_pairs/src"))
    ap.add_argument("--pseudo", action="store_true", help="znacznik: dane pseudo (próba generalna) — wpisany w nagłówek")
    a = ap.parse_args()
    ubp = os.path.join(a.ub, "pairs") if os.path.isdir(os.path.join(a.ub, "pairs")) else a.ub
    U = jl(os.path.join(ubp, "fit2d_unblind.json"))
    if not U: sys.exit(f"brak {ubp}/fit2d_unblind.json")
    C = jl(os.path.join(a.pf, "compare_pf.json")) or {}; S = jl(os.path.join(a.pf, "pairs", "syst", "syst_blind.json")) or {}
    PI = jl(os.path.join(ubp, "pair_inputs.json")) or {}
    out = {"generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()), "ub": a.ub, "pf": a.pf,
           "header_note": "**PRÓBA GENERALNA — pseudo-dane z modelu, NIE wynik.**" if a.pseudo else "Wynik odsłonięty (decyzja: flaga UNBLIND_OK, bramki T-PREP).",
           "norm_mode": g(PI, "meta", "norm_mode") or C.get("norm_mode_pf") or "—",
           "ndstar": f(g(PI, "meta", "dstar_fit_data", "N"), 0), "ndstar_method": g(PI, "meta", "dstar_fit_data", "method") or "—",
           "ndstar_old": f(g(PI, "meta", "old", "dstar_fit_data", "N") or g(PI, "meta", "dstar_fit_data_old", "N"), 0)}
    rows, prow, srows, nrows = [], [], [], []
    stat = S.get("sigma_stat") or {}; tot = S.get("total_all") or {}
    for ch in CH:
        r = g(U, "channels", ch) or {}; p68 = r.get("mu_profile_68") or [None, None]
        zexp = g(C, "fit2d_pf", ch, "Z_asimov_arms"); zo = r.get("Z")
        sy = tot.get(ch)
        N_obs = PI and None
        rows.append(f"| {NAME[ch]} | {f(r.get('mu_hat'))} | [{f(p68[0])}, {f(p68[1])}] | ± {f(r.get('mu_err_hesse'))} | ± {f(sy)} | {f(zo)} | {f(zexp)} | {f(pval(zo),4)} | {f(r.get('dev_ndf'))} | {f(r.get('S_expected_full_plane'),1)} |")
        t = g(r, "toys_at_fit") or {}
        prow.append(f"| {NAME[ch]} | {f(t.get('bias'),3)} | {f(t.get('mu_corrected'))} | {f(t.get('pull_std') or t.get('pull_sigma'))} | {f(t.get('coverage68'))} | {f(g(r,'mu_no_tag',0))} ± {f(g(r,'mu_no_tag',1))} | {f(g(r,'crosscheck_kchannel','mu_hat'))} ({f(g(r,'crosscheck_kchannel','Z'))}) |")
        sp = g(U, "syst_profiled", ch) or {}
        srows.append(f"| {NAME[ch]} | {f(stat.get(ch))} | " + " | ".join(f(g(S, "totals", grp, ch)) for grp in ("Kształty", "SF", "Normalizacja", "κ, kategorie")) + f" | {f(sy)} | {f(sp.get('mu_hat'))} | {f(sp.get('Z'))} |")
        out[f"mu_{ch}"], out[f"Z_{ch}"] = f(r.get("mu_hat")), f(zo)
        nrows.append(f"| {NAME[ch]} | {f(g(C,'n_pair_pf',ch),1)} | {f(SEPT[ch]['mu'])} / {f(SEPT[ch]['Z'],1)} | {f(r.get('mu_hat'))} / {f(zo)} |")
    c = U.get("combined") or {}; c68 = c.get("mu_common_68") or [None, None]; cs = g(U, "syst_profiled", "combined") or {}
    out.update(main_rows="\n".join(rows), toy_rows="\n".join(prow), syst_rows="\n".join(srows), cmp_rows="\n".join(nrows),
               mu_common=f(c.get("mu_common")), mu_common_68=f"[{f(c68[0])}, {f(c68[1])}]", Z_comb=f(c.get("Z")), p_comb=f(pval(c.get("Z")), 5),
               Z_comb_syst=f(cs.get("Z")), mu_comb_syst=f(cs.get("mu_common")), Z_exp_comb=f(g(C, "fit2d_pf", "Z_quadrature")),
               Z_sept=f(SEPT["combined_Z"], 1), mu_sept=f(SEPT["mu_common"]), syst_common=f(S.get("total_common")))
    out["date"] = time.strftime("%Y-%m-%d", time.gmtime())
    try: out["n_sigma_from_1"] = f(abs(float(c["mu_common"]) - 1) / (0.5 * (float(c68[1]) - float(c68[0]))), 1)
    except Exception: out["n_sigma_from_1"] = "—"
    kp = U.get("kappa_partner") or {}
    out["kappa_lines"] = "\n".join(f"- więz κ partnera σ={k}: wspólne μ = {f(g(v,'combined','mu_common'))}, Z = {f(g(v,'combined','Z'))}σ" for k, v in kp.items()) or "—"
    info = []
    Y = jl(a.ysyst) if a.ysyst else None; SA = jl(a.salad) if a.salad else None
    if Y: info.append(f"- Y (model tła D⁺ exp vs pol2, poziom pojedynczego D): plik `{os.path.basename(a.ysyst)}` — wchodzi przez U jako „bkg_poly3”/kształty; tu informacyjnie.")
    if SA:
        ent = SA.get("entries") or {}; vals = [abs(v.get("rel_diff_salad_vs_SB", 0)) for k, v in ent.items() if "wp70" in k]
        info.append(f"- S8/B4 SALAD vs SB (ε_B pod pikiem, D⁺, WP70): max |Δ| = {f(max(vals) if vals else None)} — systematyka ε_B pojedynczego D (B4 → U).")
    out["extra_syst"] = "\n".join(info) or "—"
    figs = {}
    if not a.no_plots and os.path.exists(os.path.join(ubp, "pairs_full.npz")) and os.path.exists(os.path.join(ubp, "pair_inputs.npz")):
        try: figs = plots(ubp, a.out, a.src)
        except Exception as e: out["fig_note"] = f"(wykresy: błąd {e})"
    out["figures"] = "\n".join(f"![{NAME[ch]}]({figs[ch]})" for ch in CH if ch in figs) or out.get("fig_note", "— (brak pairs_full.npz / pair_inputs.npz w fetch)")
    tpl = open(a.template, encoding="utf-8").read()
    for k, v in out.items(): tpl = tpl.replace("{{" + k + "}}", str(v))
    open(a.out + ".md", "w", encoding="utf-8").write(tpl); json.dump(out, open(a.out + ".json", "w"), indent=1, default=str)
    print(tpl)

if __name__ == "__main__":
    main()
