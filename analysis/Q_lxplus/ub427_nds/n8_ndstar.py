#!/usr/bin/env python3
"""UB427_NDS (19.11.2026, projekt A, studium W2): systematyka N_D* w 427929 — ten sam skan wariantów co NDS dla 428770
(NDS_ndstar/rerun_2026-11-12/ndstar_variants.py: tło A/E (+C sztywny osobno), zakres, binowanie, sygnał G+CB / skala stała, okno D0 36/55),
na histogramie Δm z przefiltrowanych h5 427929 (klasa MB, ta sama selekcja co u4_ndstar.py). Kształt 2G i G+CB z MC jak dla 428770
(new/ndstar_full.json; MC wspólne dla obu runów). Syst = kwadratura max |ΔN| per kategoria (bez okna i sztywnego wykładnika) — jak 428770.
  n8_ndstar.py 'h5glob' --ref ndstar_full_428770.json --out DIR      |   n8_ndstar.py --selftest --ref ... --out DIR"""
import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import nds_var as V  # noqa: E402
import dstar_refit_nds as R  # noqa: E402

EDGES = np.round(np.arange(139.6, 170.0 + 1e-9, 0.05), 6)
M0_D0 = 1.86484
PT_BINS = {"3-5": (3, 5), "5-8": (5, 8), "8-30": (8, 30)}
WIN = {"w40": 0.040, "w36": 0.036, "w55": 0.055}


def hists_from_h5(fs):
    import h5py
    dm, pt, dmd0 = [], [], []
    for f in fs:
        with h5py.File(f, "r") as g:
            e = g["events"]
            q = e["Dmeson_decayType"][:] == 1
            if not q.any():
                continue
            x = e.fields(["Dmeson_delta_mass", "Dmeson_mass", "Dmeson_pt"])[:][q]
        s = (x["Dmeson_delta_mass"] > 0) & (x["Dmeson_delta_mass"] < 0.165)
        dm.append(x["Dmeson_delta_mass"][s] * 1e3); pt.append(x["Dmeson_pt"][s])
        dmd0.append(np.abs(x["Dmeson_mass"][s] - x["Dmeson_delta_mass"][s] - M0_D0))
    dm, pt, dmd0 = (np.concatenate(v) if v else np.zeros(0) for v in (dm, pt, dmd0))
    H = {w: np.histogram(dm[dmd0 < c], EDGES)[0].astype(float) for w, c in WIN.items()}
    for b, (lo, hi) in PT_BINS.items():
        H[f"w40_{b}"] = np.histogram(dm[(dmd0 < WIN["w40"]) & (pt >= lo) & (pt < hi)], EDGES)[0].astype(float)
    return H, {"files": len(fs), "n_w40": int(H["w40"].sum())}


def synth_hists(shape, n_sig=2000, n_bkg=110000, seed=3):
    """selftest: sygnał 2G (kształt MC, skala 1,4) + tło progowe u^0,5(1+…) — znany N."""
    rng = np.random.default_rng(seed)
    c = R._sig_cdf(EDGES, shape["mu"], shape["sc"] * 1.4, shape["k"], shape["f"]); ps = np.diff((c - c[0]) / (c[-1] - c[0]))
    u = np.clip(0.5 * (EDGES[1:] + EDGES[:-1]) - R.MPI, 1e-6, None); pb = u ** 0.5 * (1 + 0.3 * u / 20); pb /= pb.sum()
    H = {}
    for w, f in (("w40", 1.0), ("w36", 0.93), ("w55", 1.12)):
        H[w] = rng.poisson(n_sig * f * ps + n_bkg * f * 1.3 ** (f - 1) * pb).astype(float)
    for b, fr in (("3-5", 0.55), ("5-8", 0.33), ("8-30", 0.12)):
        H[f"w40_{b}"] = rng.poisson(n_sig * fr * ps + n_bkg * fr ** 1.5 * pb).astype(float)
    return H, {"files": 0, "n_w40": int(H["w40"].sum()), "true_N_w40": n_sig}


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("h5", nargs="?"); ap.add_argument("--ref", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    ref = json.load(open(a.ref)); s2g, gcb = ref["shape_2g_mc"], ref["shape_gcb_mc"]
    H, info = synth_hists(s2g) if a.selftest else hists_from_h5(sorted(glob.glob(a.h5)))
    np.savez(os.path.join(a.out, "dstar_hist_427929.npz"), edges=EDGES, **H)
    res = {"info": info, "dstar_refit_file": R.__file__, "shape_2g_mc": s2g, "shape_gcb_mc": gcb, "variants": [], "per_bin": {}}

    def add(tag, cat, w="w40", **kw):
        r = V.fit_generic(H[w], EDGES, gcb if kw.get("sig") == "GCB" else s2g, **kw); r.update(tag=tag, cat=cat, window=w); r.pop("attempts", None)
        res["variants"].append(r)
        print(f"{tag:26s} N={r['N']:6.0f} ± {r['N_err']:4.0f} s={r['scale']:.2f} valid={r['valid']} χ²/ndf={r['chi2']:.0f}/{r['ndf']} {r['bkg_model']}", flush=True)
    add("nominal", "nominal", bkg="B")
    add("bkg_expo_pol2_fullmodel", "bkg", bkg="A"); add("bkg_sqrt_pol3", "bkg", bkg="E"); add("bkg_pure_expo", "bkg_rigid", bkg="C")
    add("range_hi158", "range", bkg="B", hi=158.0); add("range_hi165", "range", bkg="B", hi=165.0); add("range_lo141", "range", bkg="B", lo=141.0)
    add("range_lo141_hi165", "range", bkg="B", lo=141.0, hi=165.0)
    add("bin_0.10", "bin", bkg="B", merge=2); add("bin_0.20", "bin", bkg="B", merge=4)
    add("sig_gauss_cb", "sig", bkg="B", sig="GCB"); add("sig_scale_fixed", "sig", bkg="B", s_free=False)
    add("window_d0_36", "window", w="w36", bkg="B"); add("window_d0_55", "window", w="w55", bkg="B")
    for b in PT_BINS:
        try:
            r = R.fit_single(H[f"w40_{b}"], EDGES, s2g); r.pop("attempts", None)
        except Exception as ex:  # noqa: BLE001
            r = {"error": repr(ex)}
        res["per_bin"][b] = r
        try:   # skala σ stała (= MC; w 428770 nominał ma 1,02) — w 427929 skala wolna wychodzi 1,1–1,5 (podejrzenie: tło wciągnięte w sygnał)
            rf = V.fit_generic(H[f"w40_{b}"], EDGES, s2g, bkg="B", s_free=False); rf.pop("attempts", None)
        except Exception as ex:  # noqa: BLE001
            rf = {"error": repr(ex)}
        res["per_bin"][f"{b}_scale_fixed"] = rf
    nom = res["variants"][0]; cats = {}
    for r in res["variants"][1:]:
        if r["cat"] in ("window", "bkg_rigid") or not r["valid"]:
            continue
        cats.setdefault(r["cat"], []).append(abs(r["N"] - nom["N"]))
    sc = {c: max(v) for c, v in cats.items()}; syst = float(np.sqrt(sum(x * x for x in sc.values())))
    shape_mc_rel = ref["summary"].get("shape_mc_rel", 336.0 / 4436.4)       # „kształt MC” 428770 (decomp), względnie — MC wspólne
    rel = {"stat": nom["N_err"] / nom["N"], "nds_syst": syst / nom["N"], "shape_mc": shape_mc_rel}
    rel["syst"] = float(np.hypot(rel["nds_syst"], rel["shape_mc"])); rel["total"] = float(np.hypot(rel["stat"], rel["syst"]))
    ref_rel = {"stat": 670.86 / 4436.39, "nds_syst": 605.0 / 4436.39, "shape_mc": 336.0 / 4436.39, "total": 0.219}
    res["summary"] = {"N": nom["N"], "N_stat": nom["N_err"], "N_syst_nds": syst, "syst_per_cat": sc, "rel": rel, "rel_428770": ref_rel,
                      "all_valid": all(r["valid"] for r in res["variants"]),
                      "per_bin_sum": float(sum(r.get("N", 0) for k, r in res["per_bin"].items() if not k.endswith("_scale_fixed"))),
                      "per_bin_sum_scale_fixed": float(sum(r.get("N", 0) for k, r in res["per_bin"].items() if k.endswith("_scale_fixed"))),
                      "scale_nominal": nom["scale"], "scale_428770": 1.02}
    if a.selftest:
        t = info["true_N_w40"]; res["summary"]["selftest"] = {"true_N": t, "pull": (nom["N"] - t) / nom["N_err"], "ok": abs(nom["N"] - t) < 3 * nom["N_err"]}
    json.dump(res, open(os.path.join(a.out, "ndstar_variants_427929.json"), "w"), indent=1, default=float)
    md = ["# N_D* 427929 — skan wariantów NDS (jak 428770)", "", "| wariant | kat. | okno | tło | N ± stat | skala | valid | χ²/ndf |", "|---|---|---|---|---|---|---|---|"]
    md += [f"| {r['tag']} | {r['cat']} | {r['window']} | {r['bkg_model']} | {r['N']:.0f} ± {r['N_err']:.0f} | {r['scale']:.2f} | {r['valid']} | {r['chi2']:.0f}/{r['ndf']} |"
           for r in res["variants"]]
    s = res["summary"]
    md += ["", f"**N_D*(427929) = {s['N']:.0f} ± {s['N_stat']:.0f} (stat) ± {s['N_syst_nds']:.0f} (syst NDS)**; per kat. {json.dumps({k: round(v) for k, v in sc.items()})}",
           f"Względnie: stat {rel['stat']:.3f}, syst NDS {rel['nds_syst']:.3f}, kształt MC (z 428770) {rel['shape_mc']:.3f} → syst {rel['syst']:.3f}, razem {rel['total']:.3f} "
           f"(428770: razem {ref_rel['total']:.3f}). Suma per bin {s['per_bin_sum']:.0f}; ze stałą skalą σ {s['per_bin_sum_scale_fixed']:.0f}.",
           f"Skala σ (wolna) nominał {s['scale_nominal']:.2f} vs 428770 1,02; per bin: "
           + ", ".join(f"{b} {res['per_bin'][b].get('scale', float('nan')):.2f} (N {res['per_bin'][b].get('N', float('nan')):.0f} → stała skala "
                       f"{res['per_bin'][b + '_scale_fixed'].get('N', float('nan')):.0f})" for b in PT_BINS)]
    if a.selftest:
        md.append(f"SELFTEST: prawdziwe N = {s['selftest']['true_N']}, pull {s['selftest']['pull']:+.2f} → {'OK' if s['selftest']['ok'] else 'FAIL'}")
    open(os.path.join(a.out, "NDS_427929.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md[-3:]))
    return 0 if (not a.selftest or s["selftest"]["ok"]) else 1


if __name__ == "__main__":
    sys.exit(main())
