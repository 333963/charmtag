"""NDS: skan wariantów N_D* na pełnej lumi (histogramy Δm z `charmtag` validate / dstar_hist) → ndstar_full.json/.md.
Użycie: python scripts/studies/ndstar_variants.py dstar_nom.npz dstar_post.json OUT_DIR  (npz: edges, norm_all, norm_<bin>, hist. MC sygnału).
Warianty: model tła (nominał próg×pol2; pełny z e^{-bu}; czysty e^{-bu}; √u×pol3), zakres, binowanie, sygnał (2G / G+CB; skala σ stała),
okno masy D0 (±36/±55 — inna selekcja, raport), per bin pT. Systematyka = suma kwadratowa max |ΔN| per kategoria (okno wyłączone,
model „czysty wykładnik” liczony osobno — Δχ²=+5 i σ×1,17 = artefakt fit_dm)."""
import sys, json, os, numpy as np
from iminuit import Minuit
from scipy.stats import norm, crystalball
from charmtag import dstar_refit as R
MPI = R.MPI

def sig_cdf_gcb(x, mu, sc, k, f, beta, mcb):
    return f * norm.cdf(x, mu, sc) + (1 - f) * (1 - crystalball.cdf(-(x - mu) / (sc * k), beta, mcb))

def bkg_cdf_gen(xe, kind, a, b, c1, c2):
    xe = np.asarray(xe, float); lo, hi = xe[:-1], xe[1:]
    x = 0.5 * (hi + lo)[:, None] + 0.5 * (hi - lo)[:, None] * R._GL_X[None, :]
    u = np.clip(x - MPI, 1e-9, None); t = u / 20.0
    if kind == "C":   fb = u ** a * np.exp(-b * u)
    elif kind == "E": fb = np.sqrt(u) * (1 + a * t + b * t * t + c1 * t ** 3)
    fb = np.clip(fb, 1e-12, None); ib = 0.5 * (hi - lo) * (fb * R._GL_W[None, :]).sum(1)
    cb = np.concatenate([[0.0], np.cumsum(ib)]); return cb / cb[-1]

def rebin(h, e, lo, hi, merge=1):
    h = np.asarray(h, float); e = np.asarray(e, float)
    if merge > 1:
        n = (len(h) // merge) * merge; h = h[:n].reshape(-1, merge).sum(1); e = e[:n + 1:merge]
    sel = (e[:-1] >= lo - 1e-9) & (e[1:] <= hi + 1e-9)
    return h[sel], np.r_[e[:-1][sel], e[1:][sel][-1]]

def fit_generic(h, e, shape, bkg, sig="2G", lo=R.LO, hi=160.0, merge=1, s_free=True):
    """Warianty poza fit_single: inne modele tła (C, E), sygnał G+CB, skala stała, zakres/binowanie (dla nominalnego modelu tła → fit_single)."""
    hh, ee = rebin(h, e, lo, hi, merge)
    if bkg in ("A", "B") and sig == "2G" and s_free:
        f = R.fit_single(hh, ee, shape, attempts=(("threshold_pol2", False, 2),) if bkg == "B" else (("expo_pol2", True, 2), ("expo_pol2_s1", True, 1)))
        return {**f, "bin_mev": float(np.diff(ee)[0])}
    def cdf_s(xe, mu, s):
        c = R._sig_cdf(xe, mu, shape["sc"] * s, shape["k"], shape["f"]) if sig == "2G" else sig_cdf_gcb(xe, mu, shape["sc"] * s, shape["k"], shape["f"], shape["beta"], shape["mcb"])
        return (c - c[0]) / (c[-1] - c[0])
    def cdf_b(xe, a, b, c1, c2):
        return R._bkg_cdf_t(xe, a, b, c1, c2, bkg == "A") if bkg in ("A", "B") else bkg_cdf_gen(xe, bkg, a, b, c1, c2)
    def nll(N, mu, s, nb, a, b, c1, c2):
        lam = np.clip(np.diff(N * cdf_s(ee, mu, s) + nb * cdf_b(ee, a, b, c1, c2)), 1e-9, None); return float(np.sum(lam - hh * np.log(lam)))
    best = None
    for a0, b0 in ((0.5, 0.03), (0.6, 0.1), (0.4, 0.0)):
        m = Minuit(nll, N=max(0.01 * hh.sum(), 50), mu=shape["mu"], s=1.0, nb=hh.sum(), a=a0, b=b0, c1=0.0, c2=0.0); m.errordef = 0.5; m.strategy = 2
        m.limits["N"] = (0, None); m.limits["nb"] = (0, None); m.limits["s"] = (0.8, 1.6); m.limits["mu"] = (shape["mu"] - 0.5, shape["mu"] + 0.5)
        m.limits["a"] = (0.01, 3); m.limits["b"] = (-0.5, 1.5); m.limits["c1"] = (-5, 20); m.limits["c2"] = (-10, 40)
        if bkg == "B": m.fixed["b"] = True; m.values["b"] = 0.0
        if bkg == "C": m.fixed["c1"] = True; m.fixed["c2"] = True
        if bkg == "E": m.limits["a"] = (-5, 20); m.limits["b"] = (-20, 40); m.limits["c1"] = (-40, 40); m.fixed["c2"] = True; m.values["a"] = 0.0; m.values["b"] = 0.0
        if not s_free: m.fixed["s"] = True
        m.simplex(); m.migrad(ncall=50000); m.hesse()
        if best is None or (m.valid and (not best.valid or m.fval < best.fval)): best = m
        if best.valid: break
    m = best; v = {k: float(m.values[k]) for k in m.parameters}
    lam = np.diff(v["N"] * cdf_s(ee, v["mu"], v["s"]) + v["nb"] * cdf_b(ee, v["a"], v["b"], v["c1"], v["c2"]))
    return {"N": v["N"], "N_err": float(m.errors["N"]), "scale": v["s"], "mu": v["mu"], "valid": bool(m.valid), "posdef": bool(m.fmin.has_posdef_covar),
            "edm": float(m.fmin.edm), "chi2": float(np.sum((hh - lam) ** 2 / np.maximum(lam, 1))), "ndf": int(len(hh) - m.nfit), "n_entries": int(hh.sum()),
            "bkg_model": {"A": "expo_pol2", "B": "threshold_pol2", "C": "pure_expo", "E": "sqrt_pol3"}[bkg], "bin_mev": float(np.diff(ee)[0])}

def gcb_shape_from_mc(hmc, e, s2g):
    hh, ee = rebin(hmc, e, 143.0, 148.5)
    def nll(n, mu, sc, k, f, beta, mcb):
        c = sig_cdf_gcb(ee, mu, sc, k, f, beta, mcb); lam = np.clip(np.diff(n * c), 1e-9, None); return float(np.sum(lam - hh * np.log(lam)))
    m = Minuit(nll, n=hh.sum(), mu=s2g["mu"], sc=s2g["sc"], k=2.0, f=0.6, beta=1.5, mcb=3.0); m.errordef = 0.5; m.strategy = 2
    m.limits["sc"] = (0.2, 2); m.limits["k"] = (1.0, 6); m.limits["f"] = (0.05, 0.98); m.limits["beta"] = (0.3, 5); m.limits["mcb"] = (1.1, 20); m.limits["mu"] = (144.8, 146.2)
    m.simplex(); m.migrad(ncall=50000); m.hesse()
    return {p: float(m.values[p]) for p in ("mu", "sc", "k", "f", "beta", "mcb")} | {"valid": bool(m.valid)}

if __name__ == "__main__":
    npz, post, out = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True)
    z = np.load(npz); e = z["edges"]; P = json.load(open(post)); shape = P["norm"]["shape_mc"]
    H = {"w40": z["norm_all"].astype(float)}
    for w, key in (("w55", "nom"), ("w36", "w36")):
        H[w] = sum(z[f"nom__H__data_none__{key}__{b}__all"].astype(float) for b in ("3-5", "5-8", "8-30"))
    hmc = sum(z[f"nom__H__signal_mc__nom__{b}__all"].astype(float) for b in ("3-5", "5-8", "8-30"))
    gcb = gcb_shape_from_mc(hmc, e, shape)
    res = {"shape_2g_mc": shape, "shape_gcb_mc": gcb, "variants": [], "per_bin": {}}
    def add(tag, cat, w="w40", **kw):
        r = fit_generic(H[w], e, gcb if kw.get("sig") == "GCB" else shape, **kw); r.update(tag=tag, cat=cat, window=w); r.pop("attempts", None)
        res["variants"].append(r)
        print(f"{tag:26s} N={r['N']:6.0f} ± {r['N_err']:4.0f}  s={r['scale']:.2f} μ={r['mu']:.2f} valid={r['valid']!s:5} posdef={r['posdef']!s:5} edm={r['edm']:.0e} χ²/ndf={r['chi2']:.0f}/{r['ndf']} {r['bkg_model']}", flush=True)
    add("nominal", "nominal", bkg="B")
    add("bkg_expo_pol2_fullmodel", "bkg", bkg="A"); add("bkg_sqrt_pol3", "bkg", bkg="E"); add("bkg_pure_expo", "bkg_rigid", bkg="C")
    add("range_hi158", "range", bkg="B", hi=158.0); add("range_hi165", "range", bkg="B", hi=165.0); add("range_lo141", "range", bkg="B", lo=141.0); add("range_lo141_hi165", "range", bkg="B", lo=141.0, hi=165.0)
    add("bin_0.10", "bin", bkg="B", merge=2); add("bin_0.20", "bin", bkg="B", merge=4)
    add("sig_gauss_cb", "sig", bkg="B", sig="GCB"); add("sig_scale_fixed", "sig", bkg="B", s_free=False)
    add("window_d0_36", "window", w="w36", bkg="B"); add("window_d0_55", "window", w="w55", bkg="B")
    for b in ("3-5", "5-8", "8-30"):
        r = R.fit_single(z[f"norm_{b}"], e, shape); r.pop("attempts", None); res["per_bin"][b] = r
        print(f"bin {b:5s} N={r['N']:6.0f} ± {r['N_err']:4.0f} s={r['scale']:.2f} valid={r['valid']} {r['bkg_model']}")
    nom = res["variants"][0]; cats = {}
    for r in res["variants"][1:]:
        if r["cat"] in ("window", "bkg_rigid"): continue
        cats.setdefault(r["cat"], []).append(abs(r["N"] - nom["N"]))
    syst_cat = {c: max(v) for c, v in cats.items()}; syst = float(np.sqrt(sum(x * x for x in syst_cat.values())))
    rigid = [r for r in res["variants"] if r["cat"] == "bkg_rigid"][0]
    vals = [r["N"] for r in res["variants"] if r["valid"] and r["window"] == "w40" and r["cat"] != "bkg_rigid"]
    res["summary"] = {"N_nominal": nom["N"], "N_stat": nom["N_err"], "N_syst": syst, "syst_per_cat": syst_cat,
                      "N_syst_incl_rigid_expo": float(np.sqrt(syst ** 2 + (rigid["N"] - nom["N"]) ** 2)), "rigid_expo": {"N": rigid["N"], "dchi2": rigid["chi2"] - nom["chi2"], "scale": rigid["scale"]},
                      "median_w40": float(np.median(vals)), "spread_w40": float(np.std(vals)), "n_variants_w40": len(vals), "all_valid": all(r["valid"] for r in res["variants"]),
                      "per_bin_sum": float(sum(r["N"] for r in res["per_bin"].values())),
                      "N_fit_dm_pf_pairs": 5587.82, "N_fit_dm_pf_pairs_err": 585.59, "N_fit_dm_V": 6607.13, "N_refit_V_invalid": 4471.67,
                      "mu_correction_factor_vs_fit_dm": 5587.82 / nom["N"],
                      "selection": "D0 z tagu D*, |m(D*)−Δm−m_D0|<40 MeV, Δm ∈ [m_π+0,9; 160] MeV, bez cięcia score, pełna lumi (V dstar_nom.npz norm_all; 656845 wejść w [139,6;165] = pair_inputs n_entries)",
                      "note": "N_exp par ∝ N_D* (direct/frozen/truth) → μ_new = μ_old × N_used/N_new; błąd względny N_D* (stat⊕syst) dodać w kwadraturze do błędu μ"}
    json.dump(res, open(os.path.join(out, "ndstar_full.json"), "w"), indent=1, ensure_ascii=False)
    with open(os.path.join(out, "ndstar_full.md"), "w") as f:
        f.write("# N_D* pełna lumi — skan wariantów (NDS, naprawiony `dstar_refit.fit_single`)\n\n| wariant | kat. | okno D0 | tło | N ± stat | skala σ | μ [MeV] | valid | χ²/ndf |\n|---|---|---|---|---|---|---|---|---|\n")
        for r in res["variants"]: f.write(f"| {r['tag']} | {r['cat']} | {r['window']} | {r['bkg_model']} | {r['N']:.0f} ± {r['N_err']:.0f} | {r['scale']:.2f} | {r['mu']:.2f} | {'✔' if r['valid'] else '✖'} | {r['chi2']:.0f}/{r['ndf']} |\n")
        f.write("\n| bin pT | N ± stat | skala σ | valid | tło |\n|---|---|---|---|---|\n")
        for b, r in res["per_bin"].items(): f.write(f"| {b} | {r['N']:.0f} ± {r['N_err']:.0f} | {r['scale']:.2f} | {'✔' if r['valid'] else '✖'} | {r['bkg_model']} |\n")
        s = res["summary"]
        f.write(f"\n**N_D* = {s['N_nominal']:.0f} ± {s['N_stat']:.0f} (stat) ± {s['N_syst']:.0f} (syst)**; syst per kat.: {json.dumps({k: round(v) for k, v in s['syst_per_cat'].items()})}; "
                f"z „czystym wykładnikiem” (Δχ²=+{s['rigid_expo']['dchi2']:.0f}, σ×{s['rigid_expo']['scale']:.2f}): ±{s['N_syst_incl_rigid_expo']:.0f}. Mediana w40 {s['median_w40']:.0f}, rozrzut {s['spread_w40']:.0f} ({s['n_variants_w40']} wariantów). Suma per bin {s['per_bin_sum']:.0f}.\n")
        f.write(f"Stare: fit_dm (pf_pairs, użyte w unblind) {s['N_fit_dm_pf_pairs']:.0f} ± {s['N_fit_dm_pf_pairs_err']:.0f}; fit_dm (V) {s['N_fit_dm_V']:.0f}; refit V (invalid) {s['N_refit_V_invalid']:.0f}. Korekta μ: ×{s['mu_correction_factor_vs_fit_dm']:.3f}.\n")
    print(json.dumps(res["summary"], indent=1, ensure_ascii=False))
