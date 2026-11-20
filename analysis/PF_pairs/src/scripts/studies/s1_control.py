"""Studium S1: walidacja klasycznych próbek kontrolnych na MC (mb_nd, tło truth≥8) + przygotowanie klasyfikatora H."""
import json, pickle, sys, time
import numpy as np, pandas as pd
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # AA: pakiet z AA_merge/src
from charmtag import config as C
from charmtag.bkggen import data as D, control as K, validate as V, study as ST

SUB = os.environ.get("CHARM_SUB") or (sys.argv[1] if len(sys.argv) > 1 else "../S_bkggen/sub")   # AA: parametr
OUT = os.environ.get("CHARM_OUT") or (sys.argv[2] if len(sys.argv) > 2 else "runs/s1_control")
os.makedirs(OUT, exist_ok=True)
cfg = C.load(); sp = cfg["species"]["D0"]; t0 = time.time()
A = ST.load_all(cfg, SUB)
mc, sig, data = A["mc"], A["sig"], A["data"]
print({k: len(v) for k, v in A.items()}, "%.0fs" % (time.time() - t0))
feats = ST.feats_H(cfg, sig)
S = sig[D.is_sig_truth(sig, cfg)].reset_index(drop=True)
B = K.sideband(data, sp, "D0", "both").reset_index(drop=True)
wB = ST.bkg_weights_to_signal(B, S, cfg)
models, oS, oB = ST.train_quick(S, B, feats)
thr = ST.wp_thresholds(oS, S["pt"].to_numpy())
print("train", len(S), len(B), "AUC", V.auc(np.r_[oS, oB], np.r_[np.ones(len(oS)), np.zeros(len(oB))]), "%.0fs" % (time.time() - t0))
mc["score"] = ST.score(models, mc, feats)
data["score"] = ST.score(models, data, feats)
# OOF dla sidebandów danych (trenowane)
sbm = K.region_mask(data, sp, "D0", "near") | K.region_mask(data, sp, "D0", "far")
data.loc[sbm, "score"] = oB  # ta sama kolejność co B (sideband 'both' zachowuje kolejność)
pickle.dump({"models": models, "feats": feats, "thr": thr}, open(f"{OUT}/clfH.pkl", "wb"))
mc.to_pickle(f"{OUT}/mc.pkl"); data.to_pickle(f"{OUT}/data.pkl"); S.assign(score=oS).to_pickle(f"{OUT}/sig.pkl")

res = {"n": {k: int(len(v)) for k, v in A.items()}, "feats": feats}
bk = D.is_bkg_truth(mc)
ref = mc[bk & K.region_mask(mc, sp, "D0", "sr")]
cols = feats + ["score"]
summ = {}
for nm, ctrl in {"near": K.sideband(mc, sp, "D0", "near"), "far": K.sideband(mc, sp, "D0", "far"),
                 "near_truthbkg": K.sideband(mc[bk], sp, "D0", "near")}.items():
    cf = V.compare_features(ref, ctrl, cols)
    for wp in (None, 0.7):
        r_ = ref if wp is None else ref[ST.pass_wp(ref.score.to_numpy(), ref.pt.to_numpy(), thr, wp)]
        c_ = ctrl if wp is None else ctrl[ST.pass_wp(ctrl.score.to_numpy(), ctrl.pt.to_numpy(), thr, wp)]
        cf2 = V.compare_features(r_, c_, cols)
        summ[f"{nm}_{'nocut' if wp is None else 'wp70'}"] = {
            "n_ctrl": int(len(c_)), "n_ref": int(len(r_)),
            "chi2ndf_med": float(np.nanmedian([v["chi2ndf"] for v in cf2.values()])),
            "chi2ndf_max": float(np.nanmax([v["chi2ndf"] for v in cf2.values()])),
            "ks_max": float(np.nanmax([v["ks"] for v in cf2.values()])),
            "worst": max(cf2, key=lambda k: cf2[k]["chi2ndf"] if np.isfinite(cf2[k]["chi2ndf"]) else -1),
            "score": cf2["score"]}
    summ[f"{nm}_nocut"]["per_feature"] = cf
    # efektywność WP70 tła: kontrola vs prawda
    for wp in (0.8, 0.7, 0.5):
        e_ref = ST.pass_wp(ref.score.to_numpy(), ref.pt.to_numpy(), thr, wp).mean()
        e_c = ST.pass_wp(ctrl.score.to_numpy(), ctrl.pt.to_numpy(), thr, wp).mean()
        summ[f"{nm}_nocut"][f"effB_ratio_{int(wp*100)}"] = float(e_c / max(e_ref, 1e-9))
res["sideband"] = summ
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "per_feature"} for k, v in summ.items()}, indent=0, default=float)[:3000])

# --- refleksja: szablon z „danych” MC (wszyscy kandydaci) vs prawda MC (truth 100/101 = swap)
lo, hi = sp["mass_range"]
R = K.reflection(mc, sp)
h_c, v_c, e = K.reflection_template_var(R, lo, hi, 50)   # AA: bez obcinania ujemnych binów + Σw² (S7 §2.2)
h_c_clip = np.clip(h_c, 0, None)                            # tylko do rysowania / porównania ze starą liczbą
true_sw = mc[mc["truth"].isin([100, 101])]
h_t, _ = np.histogram(true_sw["m"], e)
sw_sig = mc[mc.truth.isin([0, 1])]
res["reflection"] = {"n_true_swaps": int(len(true_sw)), "n_template_net": float(h_c.sum()), "n_template_net_clipped": float(h_c_clip.sum()),
                     "n_true_signal": int(len(sw_sig)),
                     "chi2ndf": V.chi2_ndf_hist(h_t, h_c, None, v_c),           # AA: błędy Poissona (stare: σ=h → ~4× zaniżone)
                     "chi2ndf_old_method": V.chi2_ndf(0.5 * (e[1:] + e[:-1]), 0.5 * (e[1:] + e[:-1]), h_t, h_c_clip, bins=25, rng=(lo, hi)),
                     "n_template_net_err": float(np.sqrt(v_c.sum())),
                     "mean_true": float(true_sw.m.mean()), "mean_tmpl": float(np.sum(h_c * 0.5 * (e[1:] + e[:-1])) / max(h_c.sum(), 1)),
                     "rms_true": float(true_sw.m.std()),
                     "hist_true": h_t.tolist(), "hist_tmpl": h_c.tolist(), "var_tmpl": v_c.tolist(), "edges": e.tolist()}
Rd = K.reflection(data, sp); hd, _ = K.reflection_template(Rd, lo, hi, 50)
res["reflection"]["hist_data"] = hd.tolist()
print("refl", {k: v for k, v in res["reflection"].items() if not k.startswith(("hist", "edges"))})

# --- D*: losowy π_s z tego samego zdarzenia vs tło D* truth≥8 (Δm)
Ds_mc = D.load_npz(cfg, f"{SUB}/mb_p*.npz", "D0", decay_type=1, nd=3, tag="mc_dstar")
Ds_bkg = Ds_mc[D.is_bkg_truth(Ds_mc)]
mc_all = D.load_npz(cfg, f"{SUB}/mb_p*.npz", "D0", tag="mc")  # pełny zakres masy → pula torów
Rp = K.dstar_random_pis(mc_all[np.abs(mc_all.m - sp["m0"]) < 0.03], sp, n_per=3)
pis_min = float(np.nanpercentile(Ds_mc["d_pt2"], 1))
Rp = Rp[Rp.pis_pt > pis_min]
rng_dm = (0.1396, 0.165)
res["dstar"] = {"pis_pt_min_used": pis_min, "n_random": int(((Rp.dm_rand > rng_dm[0]) & (Rp.dm_rand < rng_dm[1])).sum()),
                "n_dstar_bkg": int(((Ds_bkg.dm > rng_dm[0]) & (Ds_bkg.dm < rng_dm[1])).sum()),
                "chi2ndf_dm": V.chi2_ndf(Ds_bkg.dm[(Ds_bkg.dm > rng_dm[0]) & (Ds_bkg.dm < rng_dm[1])],
                                         Rp.dm_rand[(Rp.dm_rand > rng_dm[0]) & (Rp.dm_rand < rng_dm[1])], bins=25, rng=rng_dm)}
dsb = K.dstar_wrong_dm(Ds_mc)
pk = Ds_bkg[(Ds_bkg.dm > 0.1435) & (Ds_bkg.dm < 0.1475)]
cfD = V.compare_features(pk, dsb, feats)
res["dstar"]["wrong_dm_vs_peak_bkg_chi2_med"] = float(np.nanmedian([v["chi2ndf"] for v in cfD.values()]))
res["dstar"]["wrong_dm_vs_peak_bkg_ks_max"] = float(np.nanmax([v["ks"] for v in cfD.values()]))
res["dstar"]["n_peak_bkg"] = int(len(pk)); res["dstar"]["n_wrong_dm"] = int(len(dsb))
hb, e2 = np.histogram(Ds_bkg.dm, 26, rng_dm); hr, _ = np.histogram(Rp.dm_rand, e2)
res["dstar"].update({"hist_bkg": hb.tolist(), "hist_rand": hr.tolist(), "edges": e2.tolist()})
print("dstar", {k: v for k, v in res["dstar"].items() if not k.startswith(("hist", "edges"))})

# --- mixing par D0–D0: tło–tło MC
pr = mc_all[D.in_range(mc_all, sp)].reset_index(drop=True)
I, J = K.same_event_pairs(pr)
Im, Jm = K.mixed_pairs(pr, I, J, factor=2)
SE, MX = K.pair_frame(pr, I, J, "same"), K.pair_frame(pr, Im, Jm, "mixed")
bb = lambda P: P[(P.t1 >= 8) & (P.t1 < 100) & (P.t2 >= 8) & (P.t2 < 100)]  # noqa: E731
SEb, MXb = bb(SE), bb(MX)
def abcd(P):
    a = (np.abs(P.m1 - sp["m0"]) < 0.055); b = (np.abs(P.m2 - sp["m0"]) < 0.055)
    A_, B_, C_, D_ = (a & b).sum(), (a & ~b).sum(), (~a & b).sum(), (~a & ~b).sum()
    return {"A": int(A_), "pred": float(B_ * C_ / max(D_, 1)), "closure": float(A_ / max(B_ * C_ / max(D_, 1), 1e-9)),
            "corr_m": float(np.corrcoef(P.m1, P.m2)[0, 1])}
res["mix"] = {"n_same": int(len(SEb)), "n_mixed": int(len(MXb)), "same": abcd(SEb), "mixed": abcd(MXb),
              "chi2ndf_m2D": V.chi2_ndf(SEb.m1 * 10 + SEb.m2, MXb.m1 * 10 + MXb.m2, bins=40),
              "chi2ndf_dphi": V.chi2_ndf(SEb.dphi, MXb.dphi, bins=20, rng=(0, np.pi)),
              "chi2ndf_m1": V.chi2_ndf(SEb.m1, MXb.m1, bins=25)}
# 2D: SR×SR udział
print("mix", res["mix"])
json.dump(res, open(f"{OUT}/s1_control.json", "w"), indent=1, default=float)
print("done %.0fs" % (time.time() - t0))
