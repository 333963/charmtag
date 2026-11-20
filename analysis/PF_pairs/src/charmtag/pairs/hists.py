"""Histogramy par (n_cat_pairs, 25, 25) z nóg; ślepota stosowana PRZY BUDOWIE (centrum nie trafia do histogramu).

Nogi: DataFrame z kolumnami source_file (str), entry (int), cand (int), decay_type (10|4), m [MeV], cat, trk0..trk2
(identyfikatory śladów, −1 = brak). Reguły jak unblind_lib.build_pair_hists: D⁺ jako m1 w D⁺D0; ten sam gatunek:
m1 = niższa kategoria, równe kategorie → kolejność z hasha klucza pary (ślepa na masę); pary ze wspólnym śladem usuwane.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from . import CHANNELS
from .model import Inputs, cat_pairs
from .frozen import shared_track


def leg_bins(L: pd.DataFrame, I: Inputs) -> pd.DataFrame:
    b = np.full(len(L), -1)
    win = np.zeros(len(L), bool)
    for dt in (10, 4):
        s = (L.decay_type == dt).to_numpy()
        bi = np.digitize(L.m.to_numpy()[s], I.edges[dt]) - 1
        ok = (bi >= 0) & (bi < len(I.edges[dt]) - 1)
        b[s] = np.where(ok, bi, -1)
        w = np.zeros(s.sum(), bool)
        w[ok] = ~I.SBM[dt][bi[ok]]
        win[s] = w
    L = L.assign(bin=b, win=win)
    return L[L.bin >= 0]


def _rnd_bit(keys):
    return np.array([int(hashlib.sha256(k.encode()).hexdigest()[:2], 16) & 1 for k in keys], bool)


def build_pair_hists(L: pd.DataFrame, I: Inputs, region="blind"):
    """region: full (odsłonięte) | blind (bez centrum) | sbsb. Zwraca ({ch: array}, liczniki)."""
    L = leg_bins(L, I)
    out = {ch: np.zeros((len(cat_pairs(I.ncat, *CHANNELS[ch])), I.nbins[CHANNELS[ch][0]], I.nbins[CHANNELS[ch][1]])) for ch in CHANNELS}
    idx = {ch: {cp: k for k, cp in enumerate(cat_pairs(I.ncat, *CHANNELS[ch]))} for ch in CHANNELS}
    counts = {"pairs_kept": 0, "region": region}
    for sf, g in L.groupby("source_file", sort=False):
        P = g.merge(g, on=["source_file", "entry"], suffixes=("1", "2"))
        P = P[(P.decay_type1 < P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.cand1 < P.cand2))]
        sh = shared_track(P)
        counts["pairs_shared_track"] = counts.get("pairs_shared_track", 0) + int(sh.sum())  # PF: diagnostyka weta
        counts["pairs_before_veto"] = counts.get("pairs_before_veto", 0) + int(len(P))
        P = P[~sh]
        centre = (P.win1 & P.win2).to_numpy()
        arms = (P.win1 ^ P.win2).to_numpy()
        if region != "full":
            P = P[~centre]
            arms = arms[~centre]
        if region == "sbsb":
            P = P[~arms]
        counts["pairs_kept"] += int(len(P))
        if not len(P):
            continue
        a, b = P.decay_type1.to_numpy(), P.decay_type2.to_numpy()
        c1, c2 = P.cat1.to_numpy(), P.cat2.to_numpy()
        b1, b2 = P.bin1.to_numpy(), P.bin2.to_numpy()
        same = a == b
        key = (P.source_file.astype(str) + ":" + P.entry.astype(str) + ":" + P.cand1.astype(str) + ":" + P.cand2.astype(str)).to_numpy()
        rnd = _rnd_bit(key) if same.any() else np.zeros(len(P), bool)
        swap = same & ((c1 > c2) | ((c1 == c2) & rnd))
        c1, c2 = np.where(swap, c2, c1), np.where(swap, c1, c2)
        b1, b2 = np.where(swap, b2, b1), np.where(swap, b1, b2)
        for ch, (ta, tb) in CHANNELS.items():
            ta, tb = min(ta, tb), max(ta, tb)
            s = (a == ta) & (b == tb)
            if not s.any():
                continue
            lut = np.full((8, 8), -1)
            for (i, j), k in idx[ch].items():
                lut[i, j] = k
            kk = lut[c1[s], c2[s]]
            assert (kk >= 0).all(), "para kategorii poza modelem"
            np.add.at(out[ch], (kk, b1[s], b2[s]), 1)
    return out, counts
