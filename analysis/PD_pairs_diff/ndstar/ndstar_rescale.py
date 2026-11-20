"""ndstar: przeliczenie μ (nominał unblind + bieg PD diff_v12) na wspólną bazę N_D*. Tylko odczyt JSON-ów, bez fitów.
μ ∝ 1/N_D* dokładnie (SS ∝ n_pair ∝ N_D*; SB/BS ∝ n_real ∝ N_D* mnożone przez wolne κ≥0 → skala pochłonięta; zob. ndstar.md §2).
Użycie: .venv_h5/bin/python ndstar_rescale.py [N_ref]   (domyślnie N_ref = 4443 z NDS)."""
import sys, os, json, numpy as np
from scipy.interpolate import CubicSpline
H = os.path.dirname(os.path.abspath(__file__)); PD = os.path.dirname(H); W = os.path.dirname(PD)
nom = json.load(open(os.path.join(PD, "fit2d_unblind.json")))
pin_nom = json.load(open(os.path.join(PD, "pair_inputs.json")))["meta"]
dif = json.load(open(os.path.join(PD, "runs/lxplus_v12/pairs/fit2d_diff.json")))
pin_pd = json.load(open(os.path.join(PD, "runs/lxplus_v12/pairs/pair_inputs.json")))["meta"]
pdi = json.load(open(os.path.join(PD, "runs/lxplus_v12/pairs/pair_inputs_diff.json")))
nds = json.load(open(os.path.join(W, "NDS_ndstar/ndstar_full.json")))["summary"]
N_nom, N_pd = pin_nom["dstar_fit_data"]["N"], pin_pd["dstar_fit_data"]["N"]
N_ref = float(sys.argv[1]) if len(sys.argv) > 1 else nds["N_nominal"]
k_nom, k_pd = N_nom / N_ref, N_pd / N_ref
CH = ["D0D0", "DpD0", "DpDp"]; LAB = dif["labels"]
out = {"N_used_nominal": N_nom, "N_used_pd": N_pd, "N_ref": N_ref, "k_nom": k_nom, "k_pd": k_pd,
       "N_ref_rel_err": float(np.hypot(nds["N_stat"], nds["N_syst"]) / nds["N_nominal"])}

def prof(d):  # profil ΔNLL(μ) → (μ, ΔNLL)
    p = d["profile_dNLL"]; x = np.array(sorted(float(k) for k in p)); return x, np.array([p[f"{v:g}"] if f"{v:g}" in p else p[str(v)] for v in x])
def keyed(p):
    return {round(float(k), 4): v for k, v in p.items()}
def summed(curves):  # suma profili na wspólnej siatce, spline → min, 68 % (ΔNLL=0,5), Z=√(2ΔNLL(0))
    grids = [set(keyed(c).keys()) for c in curves]; g = np.array(sorted(set.intersection(*grids)))
    y = np.sum([[keyed(c)[v] for v in g] for c in curves], axis=0)
    cs = CubicSpline(g, y); xf = np.linspace(g[0], g[-1], 20001); yf = cs(xf); i = int(np.argmin(yf)); y0 = yf[i]
    lo = xf[:i][np.where(yf[:i] - y0 >= 0.5)[0][-1]] if np.any(yf[:i] - y0 >= 0.5) else np.nan
    hi = xf[i:][np.where(yf[i:] - y0 >= 0.5)[0][0]] if np.any(yf[i:] - y0 >= 0.5) else np.nan
    return {"mu": float(xf[i]), "lo": float(lo), "hi": float(hi), "Z": float(np.sqrt(max(2 * (cs(0.0) - y0), 0))), "grid_step": float(np.diff(g).min())}

R = {"nominal": {}, "pd": {}}
for ch in CH:
    c = nom["channels"][ch]
    R["nominal"][ch] = {"mu_used": c["mu_hat"], "err_used": c["mu_err_hesse"], "lo_used": c["mu_profile_68"][0], "hi_used": c["mu_profile_68"][1],
                        "S_used": c["S_expected_full_plane"], "Z": c["Z"]}
    for lab in LAB:
        b = dif["channels"][ch]["bins"][lab]
        R["pd"][f"{ch}|{lab}"] = {"mu_used": b["mu_hat"], "err_used": b["mu_err_hesse"], "lo_used": b["mu_profile_68"][0], "hi_used": b["mu_profile_68"][1],
                                  "S_used": b["S_expected_full_plane"], "Z": b["Z"], "f": b["f_fit"]}
for blk, k in (("nominal", k_nom), ("pd", k_pd)):
    for v in R[blk].values():
        for q in ("mu", "err", "lo", "hi"):
            v[q + "_ref"] = v[q + "_used"] * k
        v["S_ref"] = v["S_used"] / k
# spójność per kanał Σ μ_i f_i (PD) vs μ_nom — obie strony na bazie N_ref
cons = {}
for ch in CH:
    s = sum(R["pd"][f"{ch}|{l}"]["mu_ref"] * R["pd"][f"{ch}|{l}"]["f"] for l in LAB)
    cons[ch] = {"sum_mu_f_ref": s, "mu_nom_ref": R["nominal"][ch]["mu_ref"], "ratio": s / R["nominal"][ch]["mu_ref"],
                "ratio_used_bases": (s / k_pd) / R["nominal"][ch]["mu_used"]}
out["consistency"] = cons
# μ wspólne: z JSON (siatka 0,25 w profilu → kwantyzacja) i z sumy profili (spline)
comb_nom_spline = summed([nom["channels"][ch]["profile_dNLL"] for ch in CH])
comb_pd_spline = summed([dif["channels"][ch]["bins"][l]["profile_dNLL"] for ch in CH for l in LAB])
like_spline = summed([dif["channels"][ch]["bins"]["like"]["profile_dNLL"] for ch in CH])
opp_spline = summed([dif["channels"][ch]["bins"]["opposite"]["profile_dNLL"] for ch in CH])
G = {"nominal_json": dict(mu=nom["combined"]["mu_common"], lo=nom["combined"]["mu_common_68"][0], hi=nom["combined"]["mu_common_68"][1], Z=nom["combined"]["Z"], k=k_nom),
     "nominal_syst_json": dict(mu=nom["syst_profiled"]["combined"]["mu_common"], lo=nom["syst_profiled"]["combined"]["mu_common_68"][0],
                               hi=nom["syst_profiled"]["combined"]["mu_common_68"][1], Z=nom["syst_profiled"]["combined"]["Z"], k=k_nom),
     "nominal_spline": {**comb_nom_spline, "k": k_nom},
     "pd_pythia_mc_json": dict(mu=dif["tests"]["pythia_mc"]["mu_global"], lo=dif["tests"]["pythia_mc"]["mu_global_68"][0],
                               hi=dif["tests"]["pythia_mc"]["mu_global_68"][1], Z=dif["tests"]["pythia_mc"]["Z_global"], k=k_pd),
     "pd_dps_flat_json": dict(mu=dif["tests"]["dps_flat"]["mu_global"], lo=dif["tests"]["dps_flat"]["mu_global_68"][0],
                              hi=dif["tests"]["dps_flat"]["mu_global_68"][1], Z=dif["tests"]["dps_flat"]["Z_global"], k=k_pd),
     "pd_all_bins_spline": {**comb_pd_spline, "k": k_pd},
     "pd_like_json": dict(mu=dif["bins"]["like"]["mu"], lo=dif["bins"]["like"]["mu_68"][0], hi=dif["bins"]["like"]["mu_68"][1], Z=dif["bins"]["like"]["Z"], k=k_pd),
     "pd_like_spline": {**like_spline, "k": k_pd},
     "pd_opp_json": dict(mu=dif["bins"]["opposite"]["mu"], lo=dif["bins"]["opposite"]["mu_68"][0], hi=dif["bins"]["opposite"]["mu_68"][1], Z=dif["bins"]["opposite"]["Z"], k=k_pd),
     "pd_opp_spline": {**opp_spline, "k": k_pd}}
for v in G.values():
    for q in ("mu", "lo", "hi"):
        v[q + "_ref"] = v[q] * v["k"]
out["global"] = G
# f_LS: niezmiennik — przelicz z μ_ref i n_pair_ref (n_pair ∝ N_D*) i porównaj z raportem
fls = {}
for ch in CH:
    n_like = pdi["bins"][0]["n_pair"][ch] / k_pd; n_opp = pdi["bins"][1]["n_pair"][ch] / k_pd
    NL = R["pd"][f"{ch}|like"]["mu_ref"] * n_like; NO = R["pd"][f"{ch}|opposite"]["mu_ref"] * n_opp
    fls[ch] = {"N_LS_ref": NL, "N_OS_ref": NO, "f_LS_ref": NL / (NL + NO), "f_LS_report": dif["like_sign"]["channels"][ch]["f_LS"][0]
               if "like_sign" in dif and "channels" in dif["like_sign"] else None}
out["f_LS"] = fls
out["channels"] = R
json.dump(out, open(os.path.join(H, f"ndstar_rescale_{int(round(N_ref))}.json"), "w"), indent=1)

# --- wydruk
print(f"N_ref={N_ref:.1f}  N_nom(użyte)={N_nom:.1f} k_nom={k_nom:.4f}  N_PD(użyte)={N_pd:.1f} k_pd={k_pd:.4f}  δN_ref/N={out['N_ref_rel_err']:.3f}")
print("\nkanał/bin            μ(użyte)        μ(N_ref) ±Hesse   68% (N_ref)       Z     S_oczek(N_ref)")
for blk in ("nominal", "pd"):
    for kk, v in R[blk].items():
        print(f"{blk:7s} {kk:16s} {v['mu_used']:6.3f}  →  {v['mu_ref']:6.3f} ± {v['err_ref']:5.3f}  [{v['lo_ref']:6.3f}; {v['hi_ref']:6.3f}]  {v['Z']:5.2f}  {v['S_ref']:6.1f}")
print("\nspójność Σμ_i f_i / μ_nom (wspólna baza):", {c: round(v["ratio"], 3) for c, v in cons.items()},
      " (bazy użyte:", {c: round(v["ratio_used_bases"], 3) for c, v in cons.items()}, ")")
print("\nglobalne:")
for kk, v in G.items():
    print(f"  {kk:20s} μ={v['mu']:.3f} [{v['lo']:.3f}; {v['hi']:.3f}] → N_ref: {v['mu_ref']:.3f} [{v['lo_ref']:.3f}; {v['hi_ref']:.3f}]  Z={v['Z']:.2f}"
          + (f"  (siatka {v['grid_step']})" if "grid_step" in v else ""))
print("\nf_LS (z μ_ref·n_pair_ref):", {c: (round(v["f_LS_ref"], 4), v["f_LS_report"] and round(v["f_LS_report"], 4)) for c, v in fls.items()})
