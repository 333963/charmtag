"""N2 etap 3: D⁺ — SF metodą inkluzywnych fitów masy (ε_data(WP)=N_sig(WP)/N_sig(bez), ε_MC z prawdy fd), per pT × WP,
+ wariant względny do WP90 (ε(WP)/ε(0.9)) jak wrzesień (etap5). Dane: scores Dplus_H ½ lumi (33 M wierszy, chunkami)."""
import os, sys, json, numpy as np, h5py
from iminuit import Minuit
from scipy.stats import norm
W = os.path.expanduser("~/mnt/charm/charm_work"); OUT = os.path.join(W, "N2_sf/out"); os.makedirs(OUT, exist_ok=True)
SC = os.path.join(W, "L_charmtag/runs/v11_full/scores/Dplus_H")
META = json.load(open(os.path.join(W, "L_charmtag/runs/v11_full/models/Dplus_H/meta.json")))["thresholds"]
PAIR = json.load(open(os.path.join(W, "L_charmtag/runs/v11_full/pairs/pair_inputs.json")))["meta"]["thresholds"]["4"]
PT = [(3, 5), (5, 8), (8, 30)]; M0 = 1.86966; WPS = ["none", "0.9", "0.65", "0.80", "0.70", "0.50"]
E = np.arange(1.700, 2.0501, 0.002)
def thr(wp, pt):
    if wp == "none": return np.full(len(pt), -np.inf)
    if wp in ("0.9", "0.65"): return np.full(len(pt), PAIR[wp])
    t = np.zeros(len(pt))
    for lo, hi in PT:
        s = (pt >= lo) & (pt < hi); t[s] = META[f"{lo}-{hi}"][f"{float(wp):.2f}"]
    return t
def hists(sample, chunk=4_000_000):
    f = h5py.File(os.path.join(SC, f"{sample}.h5"), "r"); n = f["m"].shape[0]; H = {}; T = {}
    for a in range(0, n, chunk):
        b = min(n, a + chunk); m = f["m"][a:b]; pt = f["pt"][a:b]; sc = f["score"][a:b]; vo = f["veto_ok"][a:b]; tr = f["truth"][a:b]
        sel0 = vo & (m > E[0]) & (m < E[-1])
        for lo, hi in PT:
            bb = sel0 & (pt >= lo) & (pt < hi); name = f"{lo}-{hi}"
            for wp in WPS:
                p = bb & (sc >= thr(wp, pt)); h, _ = np.histogram(m[p], bins=E); H[f"{name}_{wp}"] = H.get(f"{name}_{wp}", 0) + h
                if sample != "data_none":
                    tsel = p & np.isin(tr, [0, 1]) & (np.abs(m - M0) < 0.045); T[f"{name}_{wp}"] = T.get(f"{name}_{wp}", 0) + int(tsel.sum())
                    ht, _ = np.histogram(m[p & np.isin(tr, [0, 1])], bins=E); H[f"t_{name}_{wp}"] = H.get(f"t_{name}_{wp}", 0) + ht
        print(sample, a, b, flush=True)
    return H, T
def shape(ht):
    c = 0.5 * (E[1:] + E[:-1]); sel = np.abs(c - M0) < 0.06; x, h = c[sel], ht[sel].astype(float)
    def cdf(xe, n, mu, s, k, f): return n * (f * norm.cdf(xe, mu, s) + (1 - f) * norm.cdf(xe, mu, s * k))
    from iminuit.cost import ExtendedBinnedNLL
    ee = np.r_[E[:-1][sel], E[1:][sel][-1]]
    m = Minuit(ExtendedBinnedNLL(h, ee, cdf), n=h.sum(), mu=M0, s=0.008, k=2.5, f=0.7); m.limits["s"] = (0.003, 0.03); m.limits["k"] = (1.2, 6); m.limits["f"] = (0.1, 1); m.migrad()
    return {p: float(m.values[p]) for p in ("mu", "s", "k", "f")}
def fit(h, sh, fix_shape=True):
    def cdf_s(xe, mu, sc): a = sh["f"] * norm.cdf(xe, mu, sc) + (1 - sh["f"]) * norm.cdf(xe, mu, sc * sh["k"]); return (a - a[0]) / (a[-1] - a[0])
    def cdf_b(xe, b1, b2):
        t = (xe - M0) / 0.1; x = np.linspace(xe[0], xe[-1], 2001); tt = (x - M0) / 0.1; fb = np.exp(b1 * tt + b2 * tt * tt)
        cb = np.concatenate([[0], np.cumsum(0.5 * (fb[1:] + fb[:-1]) * np.diff(x))]); return np.interp(xe, x, cb / cb[-1])
    def nll(N, mu, s, nb, b1, b2):
        lam = np.clip(np.diff(N * cdf_s(E, mu, s) + nb * cdf_b(E, b1, b2)), 1e-9, None); return float(np.sum(lam - h * np.log(lam)))
    m = Minuit(nll, N=max(50, 0.02 * h.sum()), mu=sh["mu"], s=sh["s"], nb=h.sum(), b1=-0.5, b2=0.0); m.errordef = 0.5
    m.limits["N"] = (0, None); m.limits["nb"] = (0, None); m.limits["s"] = (sh["s"] * 0.8, sh["s"] * 1.3); m.limits["mu"] = (sh["mu"] - 0.003, sh["mu"] + 0.003)
    if fix_shape: m.fixed["s"] = True; m.fixed["mu"] = True
    m.migrad(ncall=20000); m.hesse()
    lam = np.diff(m.values["N"] * cdf_s(E, m.values["mu"], m.values["s"]) + m.values["nb"] * cdf_b(E, m.values["b1"], m.values["b2"]))
    return {"N": float(m.values["N"]), "N_err": float(m.errors["N"]), "valid": bool(m.valid), "chi2_ndf": float(np.sum((h - lam) ** 2 / np.maximum(lam, 1)) / (len(h) - m.nfit)), "n": int(h.sum())}
if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "hist":
        for s in sys.argv[2:]:
            H, T = hists(s); np.savez_compressed(os.path.join(OUT, f"dplus_hist_{s}.npz"), **H); json.dump(T, open(os.path.join(OUT, f"dplus_truth_{s}.json"), "w"))
    else:
        Hd = np.load(os.path.join(OUT, "dplus_hist_data_none.npz")); Hm = np.load(os.path.join(OUT, "dplus_hist_signal_mc.npz")); T = json.load(open(os.path.join(OUT, "dplus_truth_signal_mc.json")))
        res = {}
        for lo, hi in PT:
            name = f"{lo}-{hi}"; sh = shape(Hm[f"t_{name}_none"]); res[name] = {"shape": sh}
            fits = {wp: fit(Hd[f"{name}_{wp}"].astype(float), sh, fix_shape=(wp != "none")) for wp in WPS}
            fits_free = {wp: fit(Hd[f"{name}_{wp}"].astype(float), sh, fix_shape=False) for wp in WPS}
            for wp in WPS:
                eps_mc = T[f"{name}_{wp}"] / T[f"{name}_none"]
                ed = fits[wp]["N"] / fits["none"]["N"]; ed_err = ed * np.hypot(fits[wp]["N_err"] / fits[wp]["N"], fits["none"]["N_err"] / fits["none"]["N"])
                ed_free = fits_free[wp]["N"] / fits_free["none"]["N"]
                rel = (fits[wp]["N"] / fits["0.9"]["N"]) / (T[f"{name}_{wp}"] / T[f"{name}_0.9"])
                res[name][wp] = {"N_data": fits[wp]["N"], "N_err": fits[wp]["N_err"], "chi2_ndf": fits[wp]["chi2_ndf"], "valid": fits[wp]["valid"],
                                 "eps_data": ed, "eps_data_err": ed_err, "eps_data_freeshape": ed_free, "eps_mc": eps_mc,
                                 "SF_abs": ed / eps_mc, "SF_abs_stat": ed_err / eps_mc, "SF_abs_syst_shape": abs(ed_free - ed) / eps_mc, "SF_rel_to_wp90": rel}
                print(f"D+ {name:5s} {wp:5s}: N {fits[wp]['N']:.0f}±{fits[wp]['N_err']:.0f} chi2/ndf {fits[wp]['chi2_ndf']:.2f}  eps_data {ed:.3f}±{ed_err:.3f} (free {ed_free:.3f}) eps_mc {eps_mc:.3f}  SF_abs {ed/eps_mc:.3f}±{ed_err/eps_mc:.3f}  SF_rel90 {rel:.3f}")
        json.dump(res, open(os.path.join(OUT, "sf_Dplus.json"), "w"), indent=1, default=float)
