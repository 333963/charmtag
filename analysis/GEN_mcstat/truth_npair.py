#!/usr/bin/env python3
"""GEN_mcstat/truth_npair.py — tryb truth dla n_pair wg PROCEDURA.md §1, §2 (lokalnie: B6, B7) i §4 (wiersze syst). 19.11.2026.

[pary/D*]_truth(ch) = ρ_ch · Σ_zd w · Σ_{pary a<b ∈ ch} μ_a μ_b / Σ_zd w · Σ_{D*} μ_D*,  μ_sp = eps_mult__sp(pT, |η|, n_trk) · k_sp,
n_trk = ε_trk·n_ch (0,85 nominalnie), lead_pt = lead (generator), w = w2d.npy[n_trk, lead]; tylko ND (proc = 101); akceptacja mapy:
mode == 1, pT > 2,5, |η| < 2,7 (gn_combine.acc). Wejścia: runs/lxplus_v20/{effmap_direct.npz, w2d.npy, w2d_alldata.npy, closure_mb.json, cw_h5.json,
truthD_mb.parquet}, AFP_ddpred/gen/{A3,syst_*}/gen_*.tsv.gz, etap1 weights_ntrk.csv (w1). BEZ danych par, bez μ (reguła 0.2).
Wyjście: out/truth_npair.json, out/TRUTH_NPAIR.md, out/truth_norm_v20.json (format charmtag.pairs.truthnorm).
  truth_npair.py [--max-files N] [--threads 8] [--boot 200]"""
import argparse
import glob
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CW, "Q_lxplus", "gen_norm"))
import gn_combine as G  # noqa: E402  (bins, eps_from, w_lookup, PT_E/ETA_E/NT_E — te same funkcje co pakiet)

RUN = os.path.join(HERE, "runs", "lxplus_v20")
GEN = os.path.join(CW, "AFP_ddpred", "gen")
W1CSV = os.path.join(os.path.dirname(CW), "etap1_definicje_masy_v20260925_01", "weights_ntrk.csv")
OUT = os.path.join(HERE, "out")
CH = ("D0D0", "DpD0", "DpDp")
WCH = {"D0D0": 0.298, "DpD0": 0.580, "DpDp": 0.122}
SCALE_MU = 1.495
RHO_FD = {"D0D0": 1.060, "DpD0": 1.032, "DpDp": 1.020}
SHAPE = (len(G.PT_E) - 1, len(G.ETA_E) - 1, len(G.NT_E) - 1)
NB = int(np.prod(SHAPE))
SPP = {421: "D0", 411: "Dplus", 413: "Dstar"}
VARIANTS = {"nom": (0.85, 1.0, "w2"), "trk080": (0.80, 1.0, "w2"), "trk090": (0.90, 1.0, "w2"), "lead095": (0.85, 0.95, "w2"),
            "w2all": (0.85, 1.0, "w2all"), "w1": (0.85, 1.0, "w1")}
T0 = time.time()
CFG = {}


def log(*a):
    print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)


def setup():
    E = np.load(os.path.join(RUN, "effmap_direct.npz"))
    c = json.load(open(os.path.join(RUN, "closure_mb.json")))["diag"]
    k = {}
    for sp in ("D0", "Dplus", "Dstar"):
        m = c[f"match_{sp}"]; k[sp] = m["n_legs"] / max(m["n_legs"] - m["unmatched"], 1)      # §1: wszystkie nogi / dopasowane w akceptacji
    W1 = pd.read_csv(W1CSV)
    return {"eps": {sp: E[f"eps_mult__{sp}"].ravel() for sp in ("D0", "Dplus", "Dstar")},
            "cnt": {sp: (E[f"{sp}__den"], E[f"{sp}__num_mult"]) for sp in ("D0", "Dplus", "Dstar")},
            "k": k, "W2": np.load(os.path.join(RUN, "w2d.npy")), "W2all": np.load(os.path.join(RUN, "w2d_alldata.npy")),
            "w1_edges": np.r_[W1.lo.values, W1.hi.values[-1]], "w1": W1.w.values}


def flat(pt, aeta, ntrk):
    i, j, k = G.bins(pt, aeta, ntrk)
    return (i * SHAPE[1] + j) * SHAPE[2] + k


def weight(c, kind, ntrk, lead):
    if kind == "w1":
        return c["w1"][np.clip(np.digitize(ntrk, c["w1_edges"]) - 1, 0, len(c["w1"]) - 1)]
    return G.w_lookup(c["W2"] if kind == "w2" else c["W2all"], ntrk, lead)


def read_gen(path):
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=list("tabcdefg"), dtype={"t": str})
    E = df[df.t == "E"]; E = E[E.b == 101]                                # a = ev, b = proc, c = nch, d = lead
    E = pd.DataFrame({"ev": E.a.astype(np.int64).values, "nch": E.c.astype(float).values, "lead": E.d.astype(float).values})
    D = df[df.t == "D"]                                                   # a ev, b pdg, c pt, d eta, e phi, f fromB, g mode
    D = pd.DataFrame({"ev": D.a.astype(np.int64).values, "pdg": np.abs(D.b.astype(np.int64).values), "pt": D.c.astype(float).values,
                      "eta": D.d.astype(float).values, "mode": D.g.astype(np.int64).values})
    D = D[D.ev.isin(E.ev)].merge(E, on="ev")
    return E, D


def process(args):
    path, full = args
    c = CFG
    E, D = read_gen(path)
    meta = json.load(open(path + ".json"))
    n_nd = int(meta["proc"]["101"]["n"])
    b7 = int(((D.pdg == 413) & (D["mode"] == 1) & (D.pt > 3) & (D.eta.abs() < 2.5)).sum())
    A = D[(D["mode"] == 1) & (D.pt > G.PT_E[0]) & (D.eta.abs() < G.ETA_E[-1]) & D.pdg.isin(list(SPP))].reset_index(drop=True)
    sp = A.pdg.map(SPP).to_numpy()
    # pary D0/D+ w zdarzeniu (a<b po wierszu), DpD0: a = D+, b = D0
    d = A[A.pdg.isin([421, 411])].reset_index().rename(columns={"index": "i"})
    P = d[["ev", "i", "pdg"]].merge(d[["ev", "i", "pdg"]], on="ev", suffixes=("a", "b"))
    P = P[P.ia < P.ib]
    swap = (P.pdga == 421) & (P.pdgb == 411)
    ia = np.where(swap, P.ib, P.ia); ib = np.where(swap, P.ia, P.ib)
    pa, pb = np.where(swap, P.pdgb, P.pdga), np.where(swap, P.pdga, P.pdgb)
    chs = np.where((pa == 421) & (pb == 421), "D0D0", np.where((pa == 411) & (pb == 411), "DpDp", "DpD0"))
    out = {"n_nd": n_nd, "b7_dstar": b7, "n_pairs_truth": {ch: int((chs == ch).sum()) for ch in CH}, "var": {}}
    for v, (et, lf, wk) in VARIANTS.items():
        if not full and v != "nom":
            continue
        ntrk = et * A.nch.to_numpy(); lead = lf * A.lead.to_numpy()
        w = weight(c, wk, ntrk, lead)
        fb = flat(A.pt.to_numpy(), np.abs(A.eta.to_numpy()), ntrk)
        mu = np.zeros(len(A))
        for s in ("D0", "Dplus", "Dstar"):
            m = sp == s; mu[m] = c["eps"][s][fb[m]] * c["k"][s]
        r = {"ds_w": float((w * mu)[sp == "Dstar"].sum()), "ds_u": float(mu[sp == "Dstar"].sum())}
        for s in ("D0", "Dplus"):
            r[f"{s}_w"] = float((w * mu)[sp == s].sum()); r[f"{s}_u"] = float(mu[sp == s].sum())
        for ch in CH:
            q = chs == ch
            r[f"{ch}_w"] = float((w[ia[q]] * mu[ia[q]] * mu[ib[q]]).sum()); r[f"{ch}_u"] = float((mu[ia[q]] * mu[ib[q]]).sum())
        out["var"][v] = r
        if v == "nom":                                                    # tensory do jackknife/bootstrapu map (nominał)
            mD = sp == "Dstar"
            out["V"] = np.bincount(fb[mD], weights=w[mD], minlength=NB)
            out["S"] = {}
            for ch in CH:
                q = chs == ch; S = np.zeros((NB, NB))
                np.add.at(S, (fb[ia[q]], fb[ib[q]]), w[ia[q]])
                out["S"][ch] = S
    return out


def init(c):
    CFG.update(c)


def ratio(S, V, e, k):
    sa = {"D0D0": ("D0", "D0"), "DpD0": ("Dplus", "D0"), "DpDp": ("Dplus", "Dplus")}
    den = V @ (e["Dstar"] * k["Dstar"])
    return {ch: float((e[a] * k[a]) @ S[ch] @ (e[b] * k[b]) / den) for ch, (a, b) in sa.items()}


def combine(R, ref=None):
    """Σ w_ch δ_ch (konwencja wspólna)."""
    return float(sum(WCH[ch] * R[ch] for ch in CH))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-files", type=int, default=0); ap.add_argument("--threads", type=int, default=8); ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--groups", type=int, default=50); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    c = setup()
    log("k_sp:", c["k"])
    files = sorted(glob.glob(os.path.join(GEN, "A3", "gen_*.tsv.gz")))
    if a.max_files:
        files = files[:a.max_files]
    tune = {t: sorted(glob.glob(os.path.join(GEN, t, "gen_*.tsv.gz")))[: (a.max_files or 10 ** 6)] for t in ("syst_MPIup", "syst_MPIdown", "syst_Monash")}
    G_ = a.groups
    agg = {"n_nd": 0, "b7": 0, "npt": {ch: 0 for ch in CH}, "var": {}}
    Vg = np.zeros((G_, NB)); Sg = {ch: np.zeros((G_, NB, NB)) for ch in CH}
    with Pool(a.threads, initializer=init, initargs=(c,)) as P:
        for n, r in enumerate(P.imap(process, [(f, True) for f in files], chunksize=2)):
            g = n % G_
            agg["n_nd"] += r["n_nd"]; agg["b7"] += r["b7_dstar"]
            for ch in CH:
                agg["npt"][ch] += r["n_pairs_truth"][ch]; Sg[ch][g] += r["S"][ch]
            Vg[g] += r["V"]
            for v, d in r["var"].items():
                for kk, x in d.items():
                    agg["var"].setdefault(v, {}).setdefault(kk, 0.0); agg["var"][v][kk] += x
            if n % 50 == 0:
                log(f"A3 {n + 1}/{len(files)}")
        tagg = {}
        for t, fs in tune.items():
            for r in P.imap(process, [(f, False) for f in fs], chunksize=2):
                for kk, x in r["var"]["nom"].items():
                    tagg.setdefault(t, {}).setdefault(kk, 0.0); tagg[t][kk] += x
            log("tune", t, len(fs))

    def rat(d):
        return {ch: d[f"{ch}_w"] / d["ds_w"] for ch in CH}
    R = {v: rat(d) for v, d in agg["var"].items()}
    nom = R["nom"]
    # zgodność tensorów z sumami skalarnymi
    V, S = Vg.sum(0), {ch: Sg[ch].sum(0) for ch in CH}
    rt = ratio(S, V, c["eps"], c["k"])
    assert all(abs(rt[ch] / nom[ch] - 1) < 1e-9 for ch in CH), (rt, nom)
    # jackknife po grupach plików A3
    jk = np.array([[ratio({ch: S[ch] - Sg[ch][g] for ch in CH}, V - Vg[g], c["eps"], c["k"])[ch] for ch in CH] for g in range(G_)])
    jk_rel = np.sqrt((G_ - 1) / G_ * ((jk - jk.mean(0)) ** 2).sum(0)) / np.array([nom[ch] for ch in CH])
    # bootstrap map (Poisson liczników den / num_mult w binach)
    rng = np.random.default_rng(19)
    bs = []
    for _ in range(a.boot):
        e = {}
        for s, (den, num) in c["cnt"].items():
            e[s] = G.eps_from({f"{s}__den": rng.poisson(den).astype(float), f"{s}__num_mult": rng.poisson(num).astype(float)}, s, "num_mult").ravel()
        rr = ratio(S, V, e, c["k"]); bs.append([rr[ch] for ch in CH])
    bs = np.array(bs); bs_rel = bs.std(0, ddof=1) / np.array([nom[ch] for ch in CH])
    stat_rel = np.hypot(jk_rel, bs_rel)
    # korelacja kanałów przez mapy/D* (do wspólnego): wspólne z replik
    comb_nom = combine(nom)
    jk_c = np.array([sum(WCH[ch] * jk[g, i] / nom[ch] for i, ch in enumerate(CH)) for g in range(G_)])
    bs_c = np.array([sum(WCH[ch] * bs[b, i] / nom[ch] for i, ch in enumerate(CH)) for b in range(a.boot)])
    stat_comb = float(np.hypot(np.sqrt((G_ - 1) / G_ * ((jk_c - jk_c.mean()) ** 2).sum()), bs_c.std(ddof=1)))
    # wiersze §4 (δ_ch względne)
    cl = json.load(open(os.path.join(RUN, "closure_mb.json")))
    Rs = {sp: cl["split_half"]["singles"][sp]["ratio_truth_over_direct_mult"] for sp in ("D0", "Dplus")}
    sa = {"D0D0": ("D0", "D0"), "DpD0": ("Dplus", "D0"), "DpDp": ("Dplus", "Dplus")}
    kk = c["k"]
    kb = max(abs(kk["D0"] - 1), abs(kk["Dplus"] - 1)) > 0.03
    rows = {
        "stat. MC (prawda)": {ch: float(stat_rel[i]) for i, ch in enumerate(CH)},
        "model wag w2↔w1": {ch: abs(R["w1"][ch] / nom[ch] - 1) for ch in CH},
        "próbka danych w2": {ch: abs(R["w2all"][ch] / nom[ch] - 1) for ch in CH},
        "proxy n_ch→n_trk": {ch: max(abs(R["trk080"][ch] / nom[ch] - 1), abs(R["trk090"][ch] / nom[ch] - 1)) for ch in CH},
        "proxy lead_pt": {ch: abs(R["lead095"][ch] / nom[ch] - 1) for ch in CH},
        "transfer/binning map": {ch: abs(Rs[x] - 1) + abs(Rs[y] - 1) for ch, (x, y) in sa.items()},
        "faktoryzacja par": {ch: RHO_FD[ch] - 1 for ch in CH},
        "k_sp (nogi poza mapą)": {ch: (1.0 if kb else 0.5) * (abs(kk[x] - 1) + abs(kk[y] - 1)) for ch, (x, y) in sa.items()},
    }
    tune_r = {t: rat(d) for t, d in tagg.items()}
    rows_info = {"tune (info; wiersz ±0,17 bez zmian)": {ch: max(abs(tune_r[t][ch] / nom[ch] - 1) for t in tune_r) for ch in CH} if tune_r else {}}
    comb = {k: combine(v) for k, v in rows.items()}
    comb["stat. MC (prawda)"] = stat_comb                                  # wspólne z replik (korelacje kanałów przez D*/mapy)
    # B6: c_w truth vs h5; B7
    cw_h5 = json.load(open(os.path.join(RUN, "cw_h5.json")))
    vn = agg["var"]["nom"]
    cw_tr = {s: (vn[f"{s}_w"] / vn["ds_w"]) / (vn[f"{s}_u"] / vn["ds_u"]) for s in ("D0", "Dplus")}
    b6 = {s: {"c_w_truth": cw_tr[s], "c_w_h5": cw_h5[f"c_w_{s}"], "ratio": cw_tr[s] / cw_h5[f"c_w_{s}"],
              "pass": abs(cw_tr[s] / cw_h5[f"c_w_{s}"] - 1) < 0.10} for s in cw_tr}
    cw_tr_p = {ch: (vn[f"{ch}_w"] / vn["ds_w"]) / (vn[f"{ch}_u"] / vn["ds_u"]) for ch in CH}
    TD = pd.read_parquet(os.path.join(RUN, "truthD_mb.parquet"))
    mb_ds = int(((TD.pdg.abs() == 413) & (TD["mode"] == 1) & (TD.pt > 3) & (TD.eta.abs() < 2.5)).sum())
    N_MB = 19220986
    b7 = {"A3_dstar_per_ND_event": agg["b7"] / agg["n_nd"], "mb_nd_dstar_per_event": mb_ds / N_MB}
    b7["ratio_mb_over_A3"] = b7["mb_nd_dstar_per_event"] / b7["A3_dstar_per_ND_event"]; b7["pass_15pct"] = abs(b7["ratio_mb_over_A3"] - 1) < 0.15
    # direct (nogi h5, w2) dla porównania truth/direct
    dir_ = {ch: cl["split_half"]["channels"][ch]["direct_pair_per_dstar_w"] for ch in CH}
    rep = {"files_A3": len(files), "n_ND_events": agg["n_nd"], "true_pairs_in_acc": agg["npt"], "k_sp": kk,
           "pair_per_dstar_truth": nom, "variants": R, "direct_pair_per_dstar_w": dir_,
           "truth_over_direct": {ch: nom[ch] / dir_[ch] for ch in CH}, "c_w_pairs_truth": cw_tr_p,
           "stat_rel": {"jackknife_files": dict(zip(CH, jk_rel.tolist())), "bootstrap_maps": dict(zip(CH, bs_rel.tolist())),
                        "total": dict(zip(CH, stat_rel.tolist())), "combined": stat_comb},
           "syst_rows_rel": rows, "syst_rows_info": rows_info, "syst_rows_combined_rel": comb,
           "syst_rows_mu": {k: SCALE_MU * v for k, v in comb.items()}, "tune_ratios": tune_r, "B6": b6, "B7": b7,
           "groups": G_, "boot": a.boot, "sec": time.time() - T0}
    tag = a.tag
    json.dump(rep, open(os.path.join(OUT, f"truth_npair{tag}.json"), "w"), indent=1, default=float)
    tn = {"nominal_tag": "A3_v20_w2_trk0.85", "generator": {"n_events": agg["n_nd"], "procs": "ND (101)", "files": len(files)},
          "results": {"pair_per_dstar_reco": {f"{ch}|A3_v20_w2_trk0.85": [nom[ch], nom[ch] * stat_rel[i]] for i, ch in enumerate(CH)},
                      "reco_per_event": {f"{s}|A3_v20_w2_trk0.85": [vn[f"{k}_w"] / agg["n_nd"], None]
                                         for s, k in (("D0", "D0"), ("Dplus", "Dplus"), ("Dstar", "ds"))}},
          "variants_rel_to_nominal_pair_per_dstar": {ch: {v: R[v][ch] / nom[ch] - 1 for v in R if v != "nom"} for ch in CH},
          "note": "GEN_mcstat/truth_npair.py (PROCEDURA §1): wagi w2 w sumach → truth_apply_cw: false; reco_per_event nieznormalizowane Σw (tylko stosunki)"}
    json.dump(tn, open(os.path.join(OUT, f"truth_norm_v20{tag}.json"), "w"), indent=1, default=float)
    log("OK", json.dumps({"pair_per_dstar": nom, "truth/direct": rep["truth_over_direct"], "stat": rep["stat_rel"]["total"], "B6": {s: b6[s]["ratio"] for s in b6}}, default=float))


if __name__ == "__main__":
    main()
