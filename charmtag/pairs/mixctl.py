"""Kontrola tła par mixingiem nóg (etap S5, addytywnie). Włączenie: `pairs.mix_control: true` (domyślnie wyłączone).

Pary z RÓŻNYCH zdarzeń (noga z zdarzenia e × nogi zdarzenia-partnera z tej samej puli: krotność nóg 2/3/4+;
partner = przesunięcie o k = 1..mix_factor w losowej kolejności zdarzeń) mają z definicji brak korelacji wewnątrz zdarzenia.
Porównanie z parami z tego samego zdarzenia (reguły kolejności jak hists.build_pair_hists):
  SB×SB    — kształt 2D (χ²/ndf, po normalizacji mixingu do same-event w SB×SB), corr(m1, m2) same vs mix;
  ramiona  — same-event / mix (znormalizowany w SB×SB): >1 = korelacja prawdziwego D z kombinatoryką w zdarzeniu;
  centrum  — TYLKO mixing (ślepo): przewidywanie nieskorelowanego tła BB+SB w centrum (bez danych same-event z centrum).
Wynik: pairs/mix_control.json, pairs/mix_hists.npz (mix_<ch> = (n_cat_pairs, nb, nb), pełna płaszczyzna).
"""
from __future__ import annotations

import hashlib
import json
import logging
import os

import numpy as np
import pandas as pd

from . import CHANNELS
from .hists import build_pair_hists, leg_bins
from .model import cat_pairs

log = logging.getLogger("charmtag")


def _h(keys, seed):
    return np.array([int(hashlib.sha256(f"{seed}:{k}".encode()).hexdigest()[:12], 16) for k in keys], np.int64)


def mixed_pairs(L: pd.DataFrame, I, factor=2, seed=11, min_legs=2):
    """Nogi z binami (leg_bins) → ramka par mieszanych (kolumny *1, *2)."""
    L = L.reset_index(drop=True)
    ev = L["source_file"].astype(str) + ":" + L["entry"].astype(np.int64).astype(str)
    L = L.assign(ev=ev.to_numpy())
    n = L.groupby("ev")["ev"].transform("size").to_numpy()
    L = L[n >= min_legs].reset_index(drop=True)
    if not len(L):
        return pd.DataFrame()
    E = L.groupby("ev").size().rename("nl").reset_index()
    E["pool"] = np.minimum(E["nl"].to_numpy(), 4)
    E["h"] = _h(E["ev"].to_numpy(), seed)
    E = E.sort_values(["pool", "h"]).reset_index(drop=True)
    E["pos"] = E.groupby("pool").cumcount()
    E["npool"] = E.groupby("pool")["pos"].transform("size")
    L = L.merge(E[["ev", "pool", "pos", "npool"]], on="ev")
    keep = ["decay_type", "cat", "bin", "win", "m", "pool", "pos", "npool"]
    out = []
    for k in range(1, int(factor) + 1):
        A = L[keep].copy()
        A = A[A["npool"] > k]
        A["ppos"] = (A["pos"] + k) % A["npool"]
        B = L[keep].rename(columns={c: c + "2" for c in keep})
        P = A.merge(B, left_on=["pool", "ppos"], right_on=["pool2", "pos2"]).rename(
            columns={c: c + "1" for c in keep})
        out.append(P)
    P = pd.concat(out, ignore_index=True)
    return P


def hist_pairs(P: pd.DataFrame, I, seed=13):
    """Histogramy (n_cat_pairs, nb1, nb2) z reguł kolejności build_pair_hists (D⁺ jako m1; ten sam gatunek: niższa kat. m1)."""
    out = {ch: np.zeros((len(cat_pairs(I.ncat, *CHANNELS[ch])), I.nbins[CHANNELS[ch][0]], I.nbins[CHANNELS[ch][1]]))
           for ch in CHANNELS}
    if not len(P):
        return out
    a, b = P["decay_type1"].to_numpy(), P["decay_type2"].to_numpy()
    c1, c2 = P["cat1"].to_numpy(), P["cat2"].to_numpy()
    b1, b2 = P["bin1"].to_numpy(), P["bin2"].to_numpy()
    # D⁺ (4) jako pierwsza noga
    sw = a > b
    a, b = np.where(sw, b, a), np.where(sw, a, b)
    c1, c2 = np.where(sw, c2, c1), np.where(sw, c1, c2)
    b1, b2 = np.where(sw, b2, b1), np.where(sw, b1, b2)
    same = a == b
    rnd = np.random.default_rng(seed).random(len(P)) < 0.5
    sw = same & ((c1 > c2) | ((c1 == c2) & rnd))
    c1, c2 = np.where(sw, c2, c1), np.where(sw, c1, c2)
    b1, b2 = np.where(sw, b2, b1), np.where(sw, b1, b2)
    for ch, (ta, tb) in CHANNELS.items():
        ta, tb = min(ta, tb), max(ta, tb)
        s = (a == ta) & (b == tb)
        if not s.any():
            continue
        lut = np.full((8, 8), -1)
        for k, (i, j) in enumerate(cat_pairs(I.ncat, *CHANNELS[ch])):
            lut[i, j] = k
        kk = lut[c1[s], c2[s]]
        ok = kk >= 0
        np.add.at(out[ch], (kk[ok], b1[s][ok], b2[s][ok]), 1)
    return out


def _masks(I, ch):
    ta, tb = CHANNELS[ch]
    wa, wb = ~I.SBM[ta], ~I.SBM[tb]
    centre = np.outer(wa, wb)
    arms = np.logical_xor.outer(wa, wb)
    sbsb = ~(centre | arms)
    return sbsb, arms, centre


def _corr(H, I, ch, region):
    ta, tb = CHANNELS[ch]
    h = H.sum(0) * region
    if h.sum() < 10:
        return float("nan")
    x, y = I.C[ta], I.C[tb]
    w = h / h.sum()
    mx, my = (w.sum(1) * x).sum(), (w.sum(0) * y).sum()
    cxy = (w * np.outer(x - mx, y - my)).sum()
    vx, vy = (w.sum(1) * (x - mx) ** 2).sum(), (w.sum(0) * (y - my) ** 2).sum()
    return float(cxy / np.sqrt(max(vx * vy, 1e-30)))


def run(L: pd.DataFrame, I, out_dir: str, p: dict) -> dict:
    """L = nogi danych z cmd_pairs (source_file, entry, cand, decay_type, m [MeV], cat, trk*)."""
    factor, seed = int(p.get("mix_factor", 2)), int(p.get("mix_seed", 11))
    Lb = leg_bins(L, I)
    same_sbsb, _ = build_pair_hists(L, I, "sbsb")
    same_blind, _ = build_pair_hists(L, I, "blind")       # bez centrum (ślepo)
    P = mixed_pairs(Lb, I, factor=factor, seed=seed)
    mix = hist_pairs(P, I)
    rep = {"mix_factor": factor, "seed": seed, "n_mixed_pairs": int(len(P)), "channels": {}}
    for ch in CHANNELS:
        sbsb, arms, centre = _masks(I, ch)
        S_sb, M_sb = same_sbsb[ch] * sbsb, mix[ch] * sbsb
        ns, nm = float(S_sb.sum()), float(M_sb.sum())
        r = {"same_sbsb": ns, "mix_sbsb": nm}
        if ns < 5 or nm < 5:
            r["skipped"] = "za mało par"
            rep["channels"][ch] = r
            continue
        k = ns / nm
        pred = M_sb * k
        good = pred > 0
        chi2 = float(np.sum((S_sb[good] - pred[good]) ** 2 / np.maximum(pred[good], 1)))
        same_arms = float((same_blind[ch] * arms).sum())
        mix_arms = float((mix[ch] * arms).sum()) * k
        r.update({"norm_mix_to_same": k, "sbsb_chi2": chi2, "sbsb_ndf": int(good.sum() - 1),
                  "sbsb_chi2_ndf": chi2 / max(good.sum() - 1, 1),
                  "corr_m1m2_same_sbsb": _corr(same_sbsb[ch], I, ch, sbsb), "corr_m1m2_mix_sbsb": _corr(mix[ch], I, ch, sbsb),
                  "arms_same": same_arms, "arms_mix_pred": mix_arms,
                  "arms_ratio_same_over_mix": same_arms / max(mix_arms, 1e-9),
                  "arms_ratio_err": float(np.sqrt(1 / max(same_arms, 1) + 1 / max(mix_arms / k * 1.0, 1)) * same_arms / max(mix_arms, 1e-9)),
                  "centre_mix_pred_uncorrelated": float((mix[ch] * centre).sum()) * k,
                  "per_catpair_sbsb_ratio": [float(S_sb[i].sum() / max(pred[i].sum(), 1e-9)) for i in range(S_sb.shape[0])]})
        rep["channels"][ch] = r
        log.info("pairs mix %s: SB×SB same %d / mix %.0f×%.3g, χ²/ndf %.2f, corr same %.3f mix %.3f, ramiona same/mix %.3f",
                 ch, ns, nm, k, r["sbsb_chi2_ndf"], r["corr_m1m2_same_sbsb"], r["corr_m1m2_mix_sbsb"], r["arms_ratio_same_over_mix"])
    os.makedirs(out_dir, exist_ok=True)
    np.savez(os.path.join(out_dir, "mix_hists.npz"), **{f"mix_{ch}": v for ch, v in mix.items()},
             **{f"same_sbsb_{ch}": v for ch, v in same_sbsb.items()})
    with open(os.path.join(out_dir, "mix_control.json"), "w") as f:
        json.dump(rep, f, indent=1, default=float)
    return rep
