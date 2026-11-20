#!/usr/bin/env python3
"""GEN_mcstat: inwentarz prawdy pythia_dd (lokalne A3 + warianty), pary D-D w akceptancji (pT>3, |eta|<2.2, naturalne BR = mode==1),
proxy wag krotności na poziomie generatora (n_trk ~ TRK*n_ch, lead_pt ~ lead gen), n_eff i c_w(w1) vs rekonstrukcja (t4_rates).
Tylko odczyt wejść; wyjście: out/gen_inventory.json.
  python gen_inventory.py [--var A3] [--max N] [--trk 0.85]"""
import argparse, glob, gzip, io, json, os, sys
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = "/home/k/Documents/charm"
GEN = f"{ROOT}/charm_work/AFP_ddpred/gen"
EFF = f"{ROOT}/runs/lxplus/truth/combined_A3/eff_map.npz"
W1 = f"{ROOT}/etap1_definicje_masy_v20260925_01/weights_ntrk.csv"
SP = {421: "D0", 411: "Dplus", 413: "Dstar"}
CH = {"D0D0": ("D0", "D0"), "DpD0": ("Dplus", "D0"), "DpDp": ("Dplus", "Dplus")}
RHO = {"D0D0": 1.0595, "DpD0": 1.0315, "DpDp": 1.0203}
PTMIN, ETAMAX = 3.0, 2.2


def read(path):
    with gzip.open(path, "rt") as f:
        L = f.read().splitlines()
    ev = pd.read_csv(io.StringIO("\n".join(l[2:] for l in L if l[0] == "E")), sep="\t", names=["entry", "proc", "nch", "lead"])
    D = pd.read_csv(io.StringIO("\n".join(l[2:] for l in L if l[0] == "D")), sep="\t",
                    names=["entry", "pdg", "pt", "eta", "phi", "fromB", "mode"])
    D = D[(D["mode"] == 1) & (D.pt > PTMIN) & (D.eta.abs() < ETAMAX)].merge(ev, on="entry")
    D["sp"] = D.pdg.abs().map(SP)
    D["file"] = int(os.path.basename(path).split("_")[1].split(".")[0])
    meta = json.load(open(path + ".json"))
    return D, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var", default="A3"); ap.add_argument("--max", type=int, default=0); ap.add_argument("--trk", type=float, default=0.85)
    a = ap.parse_args()
    fs = sorted(glob.glob(f"{GEN}/{a.var}/gen_*.tsv.gz"))
    fs = [f for f in fs if os.path.exists(f + ".json")][: a.max or None]
    with Pool(8) as p:
        res = p.map(read, fs)
    D = pd.concat([r[0] for r in res], ignore_index=True)
    metas = [r[1] for r in res]
    nev = {"all": sum(m["n_ok"] for m in metas)}
    for c in ("101", "103", "104", "105", "106"):
        nev[c] = sum(m["proc"][c]["n"] for m in metas)
    nchh = np.sum([m["nch_hist"] for m in metas], axis=0)

    # proxy reco: n_trk ~ TRK * n_ch (jak truth_combine), lead_pt ~ lead gen
    D["ntrk"] = a.trk * D.nch
    E = np.load(EFF)
    pe, ee, ne = E["pt_edges"], E["eta_edges"], E["ntrk_edges"]
    def eps(sp, d):
        i = np.clip(np.digitize(d.pt, pe) - 1, 0, len(pe) - 2); j = np.clip(np.digitize(d.eta.abs(), ee) - 1, 0, len(ee) - 2)
        k = np.clip(np.digitize(d.ntrk, ne) - 1, 0, len(ne) - 2)
        return E[f"fd_nd__{sp}"][i, j, k]
    D["eps"] = 0.0
    for sp in SP.values():
        m = D.sp == sp; D.loc[m, "eps"] = eps(sp, D[m])
    W = pd.read_csv(W1); e1 = np.r_[W.lo.values, W.hi.values[-1]]
    D["w1"] = W.w.values[np.clip(np.digitize(D.ntrk, e1) - 1, 0, len(W) - 1)]

    d = D[D.sp.isin(["D0", "Dplus"])].copy(); d["idx"] = np.arange(len(d))
    P = d.merge(d[["file", "entry", "sp", "idx", "pt", "eps", "pdg"]], on=["file", "entry"], suffixes=("_a", "_b"))
    P = P[P.idx_a < P.idx_b]
    sa, sb = P.sp_a.to_numpy(), P.sp_b.to_numpy()
    P["ch"] = np.where((sa == "D0") & (sb == "D0"), "D0D0", np.where((sa == "Dplus") & (sb == "Dplus"), "DpDp", "DpD0"))
    P["ls"] = np.sign(P.pdg_a) == np.sign(P.pdg_b)

    def neff(w):
        w = np.asarray(w, float); return float(w.sum() ** 2 / (w ** 2).sum()) if len(w) else 0.0
    out = {"var": a.var, "files": len(fs), "n_events": nev, "trk_proxy": a.trk, "acceptance": {"pt_min": PTMIN, "abs_eta_max": ETAMAX},
           "nch_hist_all_procs_quantiles": {q: int(np.searchsorted(np.cumsum(nchh) / nchh.sum(), q)) for q in (0.1, 0.5, 0.9)},
           "proc_sel": {}}
    for tag, sel in (("ND", [101]), ("all", [101, 103, 104, 105, 106])):
        dd, pp = D[D.proc.isin(sel)], P[P.proc.isin(sel)]
        r = {"singles": {}, "pairs": {}}
        for sp in SP.values():
            q = dd[dd.sp == sp]
            r["singles"][sp] = {"n_truth": int(len(q)), "sum_eps": float(q.eps.sum()), "sum_w1": float(q.w1.sum()),
                                "sum_eps_w1": float((q.eps * q.w1).sum()), "neff_eps_w1": neff(q.eps * q.w1),
                                "median_nch": float(q.nch.median()), "median_lead": float(q.lead.median())}
        ds = r["singles"]["Dstar"]
        for ch in CH:
            q = pp[pp.ch == ch]
            ew = q.eps_a * q.eps_b * RHO[ch]
            x = {"n_truth": int(len(q)), "n_LS": int(q.ls.sum()), "sum_eps": float(ew.sum()), "neff_eps": neff(ew),
                 "sum_eps_w1": float((ew * q.w1).sum()), "neff_eps_w1": neff(ew * q.w1), "neff_w1": neff(q.w1),
                 "median_nch": float(q.nch.median()), "median_ntrk_proxy": float(q.ntrk.median()), "median_lead": float(q.lead.median()),
                 "frac_lead_gt4": float((q.lead > 4).mean())}
            x["rel_stat_truth"] = 1 / np.sqrt(max(x["n_truth"], 1)); x["rel_stat_eps_w1"] = 1 / np.sqrt(max(x["neff_eps_w1"], 1e-9))
            x["pair_per_dstar_reco_unw"] = x["sum_eps"] / ds["sum_eps"]
            x["pair_per_dstar_reco_w1"] = x["sum_eps_w1"] / ds["sum_eps_w1"]
            x["c_w1"] = x["pair_per_dstar_reco_w1"] / x["pair_per_dstar_reco_unw"]
            r["pairs"][ch] = x
        for sp in ("D0", "Dplus"):
            s = r["singles"][sp]
            s["c_w1"] = (s["sum_eps_w1"] / ds["sum_eps_w1"]) / (s["sum_eps"] / ds["sum_eps"])
        out["proc_sel"][tag] = r
    os.makedirs("out", exist_ok=True)
    fn = f"out/gen_inventory_{a.var}_trk{a.trk}.json"
    json.dump(out, open(fn, "w"), indent=1, default=float)
    print(json.dumps({k: out[k] for k in ("var", "files", "n_events")}))
    for tag in ("ND", "all"):
        for ch, x in out["proc_sel"][tag]["pairs"].items():
            print(tag, ch, x["n_truth"], "LS", x["n_LS"], "neff_eps", round(x["neff_eps"]), "neff_eps_w1", round(x["neff_eps_w1"]),
                  "c_w1", round(x["c_w1"], 3), "ppd", round(x["pair_per_dstar_reco_unw"], 5), "med nch", x["median_nch"], "lead>4", round(x["frac_lead_gt4"], 2))
        print(tag, {sp: (v["n_truth"], round(v.get("c_w1", 1), 3)) for sp, v in out["proc_sel"][tag]["singles"].items()})


if __name__ == "__main__":
    main()
