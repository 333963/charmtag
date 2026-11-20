#!/usr/bin/env python3
"""GEN_mcstat/truth_mbnd.py (19.11.2026) — WARIANT INFORMACYJNY, ODSTĘPSTWO od PROCEDURA §1 (prerejestracja przewidywała prawdę A3):
prawda z samego mb_nd (truthD_mb.parquet: prawdziwe D w 78 plikach mb_nd v2.0, reco n_trk/lead_pt, w2) × mapy direct → [pary/D*]_truth.
Powód: A3 (pythia_dd 8.3, standalone) ≠ mc23 mb_nd na poziomie prawdy (D⁺/D* 5,15 vs 3,63; D* / zdarzenie ×1,36 — B7), więc tryb A3 zmienia
definicję odniesienia μ. Ten sam wzór co truth_npair.py; mapy: split-half (pliki parzyste ↔ nieparzyste, jak closure_mb) i self; stat: jackknife
po plikach (78 → 26 grup) ⊕ bootstrap map. Bez danych, bez μ."""
import json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import truth_npair as T
G = T.G
c = T.setup(); T.CFG.update(c)
TD = pd.read_parquet(os.path.join(T.RUN, "truthD_mb.parquet"))
TD = TD[(TD["mode"] == 1) & (TD.pt > G.PT_E[0]) & (TD.eta.abs() < G.ETA_E[-1]) & TD.pdg.abs().isin([421, 411, 413])].reset_index(drop=True)
TD["pdg"] = TD.pdg.abs(); sp = TD.pdg.map(T.SPP).to_numpy()
fb = T.flat(TD.pt.to_numpy(), np.abs(TD.eta.to_numpy()), TD.ntrk.to_numpy())
tags = sorted(TD.tag.astype(str).unique()); tix = {t: i for i, t in enumerate(tags)}
ti = TD.tag.astype(str).map(tix).to_numpy()
W = {"w2": TD.w2.to_numpy(), "w2all": TD.w2_all.to_numpy(), "w1": T.weight(c, "w1", TD.ntrk.to_numpy(), TD.lead_pt.to_numpy()), "unw": np.ones(len(TD))}
d = TD[TD.pdg.isin([421, 411])].reset_index().rename(columns={"index": "i"})
P = d[["tag", "evnum", "i", "pdg"]].merge(d[["tag", "evnum", "i", "pdg"]], on=["tag", "evnum"], suffixes=("a", "b")); P = P[P.ia < P.ib]
swap = (P.pdga == 421) & (P.pdgb == 411)
ia = np.where(swap, P.ib, P.ia); ib = np.where(swap, P.ia, P.ib); pa, pb = np.where(swap, P.pdgb, P.pdga), np.where(swap, P.pdga, P.pdgb)
chs = np.where((pa == 421) & (pb == 421), "D0D0", np.where((pa == 411) & (pb == 411), "DpDp", "DpD0"))
NG = 26; grp = ti % NG
def tensors(w, sel_ev=None):
    m = np.ones(len(TD), bool) if sel_ev is None else sel_ev
    mD = (sp == "Dstar") & m
    V = np.bincount(fb[mD], weights=w[mD], minlength=T.NB); S = {}
    for ch in T.CH:
        q = (chs == ch) & m[ia]; Sx = np.zeros((T.NB, T.NB)); np.add.at(Sx, (fb[ia[q]], fb[ib[q]]), w[ia[q]]); S[ch] = Sx
    return S, V
out = {"n_true_pairs": {ch: int((chs == ch).sum()) for ch in T.CH}, "n_dstar": int((sp == "Dstar").sum())}
# n_eff par prawdy (wagi w2 zdarzenia × μ_a μ_b)
mu = np.zeros(len(TD))
for s in ("D0", "Dplus", "Dstar"):
    m = sp == s; mu[m] = c["eps"][s][fb[m]] * c["k"][s]
out["n_eff_truth_pairs_w2"] = {}
for ch in T.CH:
    q = chs == ch; x = W["w2"][ia[q]] * mu[ia[q]] * mu[ib[q]]; out["n_eff_truth_pairs_w2"][ch] = float(x.sum() ** 2 / (x ** 2).sum())
R = {}
for k, w in W.items():
    S, V = tensors(w); R[k] = T.ratio(S, V, c["eps"], c["k"])
out["pair_per_dstar"] = R
# split-half: mapy z drugiej parzystości (jak closure_mb) — tylko nominał w2
cl = json.load(open(os.path.join(T.RUN, "closure_mb.json")))
dir_ = {ch: cl["split_half"]["channels"][ch]["direct_pair_per_dstar_w"] for ch in T.CH}
out["direct_pair_per_dstar_w"] = dir_
out["closure_split_half_truth_pair_per_dstar"] = {ch: cl["split_half"]["channels"][ch]["truth_pair_per_dstar_w_mult"] for ch in T.CH}
out["truth_self_over_direct"] = {ch: R["w2"][ch] / dir_[ch] for ch in T.CH}
# jackknife po grupach plików ⊕ bootstrap map
S, V = tensors(W["w2"]); Sg = []; nom = R["w2"]
jk = []
for g in range(NG):
    Sx, Vx = tensors(W["w2"], grp[np.arange(len(TD))] != g)
    rr = T.ratio(Sx, Vx, c["eps"], c["k"]); jk.append([rr[ch] for ch in T.CH])
jk = np.array(jk); jk_rel = np.sqrt((NG - 1) / NG * ((jk - jk.mean(0)) ** 2).sum(0)) / np.array([nom[ch] for ch in T.CH])
rng = np.random.default_rng(19); bs = []
for _ in range(200):
    e = {s: G.eps_from({f"{s}__den": rng.poisson(den).astype(float), f"{s}__num_mult": rng.poisson(num).astype(float)}, s, "num_mult").ravel()
         for s, (den, num) in c["cnt"].items()}
    rr = T.ratio(S, V, e, c["k"]); bs.append([rr[ch] for ch in T.CH])
bs = np.array(bs); bs_rel = bs.std(0, ddof=1) / np.array([nom[ch] for ch in T.CH])
st = np.hypot(jk_rel, bs_rel)
jk_c = np.array([sum(T.WCH[ch] * jk[g, i] / nom[ch] for i, ch in enumerate(T.CH)) for g in range(NG)])
bs_c = np.array([sum(T.WCH[ch] * bs[b, i] / nom[ch] for i, ch in enumerate(T.CH)) for b in range(len(bs))])
out["stat_rel"] = {"jackknife": dict(zip(T.CH, jk_rel.tolist())), "bootstrap_maps": dict(zip(T.CH, bs_rel.tolist())), "total": dict(zip(T.CH, st.tolist())),
                   "combined": float(np.hypot(np.sqrt((NG - 1) / NG * ((jk_c - jk_c.mean()) ** 2).sum()), bs_c.std(ddof=1)))}
out["model_wag_w2_w1_rel"] = {ch: abs(R["w1"][ch] / nom[ch] - 1) for ch in T.CH}
out["model_wag_combined"] = T.combine(out["model_wag_w2_w1_rel"])
out["w2all_rel"] = {ch: abs(R["w2all"][ch] / nom[ch] - 1) for ch in T.CH}
out["c_w_pairs"] = {ch: (R["w2"][ch]) / (R["unw"][ch]) for ch in T.CH}
json.dump(out, open(os.path.join(T.OUT, "truth_mbnd.json"), "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
