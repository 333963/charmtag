"""N2 etap 2: ε(WP) D0 z tagu D* — dane (odporny fit pass/fail, V) vs MC (prawda; + ten sam fit jako closure) → SF per pT × WP.
Systematyka: wariant fitu (sb_count bez kształtu), okno 143–148.5 vs ±1 MeV, mu/skala stała vs wolna; stat z fitu."""
import os, sys, json, numpy as np
W = os.path.expanduser("~/mnt/charm/charm_work"); OUT = os.path.join(W, "N2_sf/out")
sys.path.insert(0, os.path.join(W, "Q_lxplus/v_fix/py")); import dstar_refit as R
PT = ["3-5", "5-8", "8-30"]; WPS = ["0.80", "0.70", "0.50", "0.65", "0.9"]
D = {s: np.load(os.path.join(OUT, f"dstar_hist_{s}.npz")) for s in ("data_none", "signal_mc", "incl_mc")}
E = D["data_none"]["edges"]
def h(s, b, wp, lab, t=False): return D[s][("t_" if t else "h_") + f"{b}_{wp}_{lab}"].astype(float)
# kształt sygnału z MC truth (fd, wszystkie pT, pass+fail WP80)
ht = sum(h("signal_mc", b, "0.80", lab, t=True) for b in PT for lab in ("pass", "fail"))
shape = R.shape_from_mc(ht, E); print("shape MC:", shape)
res = {"shape": shape, "bins": {}}
for b in PT:
    res["bins"][b] = {}
    # kształt per bin (σ rośnie z pT? użyj per-bin jeśli dość statystyki)
    htb = h("signal_mc", b, "0.80", "pass", True) + h("signal_mc", b, "0.80", "fail", True)
    shb = R.shape_from_mc(htb, E) if htb.sum() > 2000 else shape
    for wp in WPS:
        hp, hf = h("data_none", b, wp, "pass"), h("data_none", b, wp, "fail")
        fd = R.fit_passfail(hp, hf, E, shb)
        fd_fix = R.fit_passfail(hp, hf, E, shb, float_scale=False, mu_free=False)
        sbc = R.sb_count(hp, hf, E, mu=shb["mu"])
        # MC prawda: zliczenie truth w oknie 143–148.5
        c = 0.5 * (E[1:] + E[:-1]); win = (c > 143) & (c < 148.5)
        tp, tf = h("signal_mc", b, wp, "pass", True)[win].sum(), h("signal_mc", b, wp, "fail", True)[win].sum()
        eps_mc = tp / (tp + tf); eps_mc_err = np.sqrt(eps_mc * (1 - eps_mc) / (tp + tf))
        # MC closure: ten sam fit na MC fd (z tłem) i na mb (inkluzywne, małe S)
        fm = R.fit_passfail(h("signal_mc", b, wp, "pass"), h("signal_mc", b, wp, "fail"), E, shb)
        sf = fd["eps"] / eps_mc; sf_stat = sf * np.hypot(fd["eps_err"] / max(fd["eps"], 1e-6), eps_mc_err / eps_mc)
        syst = {"fit_fixed": fd_fix["eps"] / eps_mc - sf, "sb_count": (sbc["eps"] / eps_mc - sf) if np.isfinite(sbc["eps"]) else None,
                "mc_closure": fm["eps"] / eps_mc - 1.0}
        sf_syst = float(np.sqrt(sum(v ** 2 for v in syst.values() if v is not None)))
        res["bins"][b][wp] = {"eps_data": fd["eps"], "eps_data_err": fd["eps_err"], "N_tot": fd["N_tot"], "N_tot_err": fd["N_tot_err"], "scale": fd["scale"], "mu": fd["mu"],
                              "chi2_ndf": fd["chi2"] / max(fd["ndf"], 1), "valid": fd["valid"], "n_pass": fd["n_pass"], "n_fail": fd["n_fail"],
                              "eps_mc_truth": float(eps_mc), "eps_mc_truth_err": float(eps_mc_err), "eps_mc_fit": fm["eps"], "eps_mc_fit_err": fm["eps_err"],
                              "eps_data_fixed": fd_fix["eps"], "eps_data_sbcount": sbc["eps"],
                              "SF": float(sf), "SF_stat": float(sf_stat), "SF_syst": sf_syst, "syst": syst}
        print(f"{b:5s} WP{wp}: eps_data {fd['eps']:.3f}±{fd['eps_err']:.3f} (fix {fd_fix['eps']:.3f}, sb {sbc['eps']:.3f}) eps_mc {eps_mc:.3f} (fit {fm['eps']:.3f})  SF {sf:.3f}±{sf_stat:.3f}±{sf_syst:.3f}  N_tot {fd['N_tot']:.0f}±{fd['N_tot_err']:.0f} s={fd['scale']:.2f} chi2/ndf {fd['chi2']/max(fd['ndf'],1):.2f} valid={fd['valid']}")
# zintegrowane (suma binów)
res["integrated"] = {}
for wp in WPS:
    hp = sum(h("data_none", b, wp, "pass") for b in PT); hf = sum(h("data_none", b, wp, "fail") for b in PT)
    fd = R.fit_passfail(hp, hf, E, shape); c = 0.5 * (E[1:] + E[:-1]); win = (c > 143) & (c < 148.5)
    tp = sum(h("signal_mc", b, wp, "pass", True)[win].sum() for b in PT); tf = sum(h("signal_mc", b, wp, "fail", True)[win].sum() for b in PT)
    eps_mc = tp / (tp + tf); res["integrated"][wp] = {"eps_data": fd["eps"], "eps_data_err": fd["eps_err"], "eps_mc_truth": float(eps_mc), "SF": fd["eps"] / eps_mc,
                                                     "SF_stat": fd["eps"] / eps_mc * fd["eps_err"] / max(fd["eps"], 1e-6), "N_tot": fd["N_tot"], "valid": fd["valid"]}
    print(f"ALL   WP{wp}: eps_data {fd['eps']:.3f}±{fd['eps_err']:.3f} eps_mc {eps_mc:.3f} SF {fd['eps']/eps_mc:.3f}±{res['integrated'][wp]['SF_stat']:.3f} N_tot {fd['N_tot']:.0f}")
json.dump(res, open(os.path.join(OUT, "sf_D0.json"), "w"), indent=1, default=float)
