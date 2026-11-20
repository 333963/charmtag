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


def shared_track_detail(P: pd.DataFrame):
    """PF3: dla każdej pary liczba wspólnych śladów (0–3) i typ: córka 0 = kaon (features: kaon_index 0), córki 1,2 = piony.
    Zwraca (n_shared, kk, kpi, pipi) — flagi: wspólny ślad jest K w obu nogach / K w jednej, π w drugiej / π w obu."""
    t1 = P[["trk01", "trk11", "trk21"]].to_numpy()
    t2 = P[["trk02", "trk12", "trk22"]].to_numpy()
    n = np.zeros(len(P), int); kk = np.zeros(len(P), bool); kpi = np.zeros(len(P), bool); pipi = np.zeros(len(P), bool)
    for i in range(3):
        for j in range(3):
            m = (t1[:, i] == t2[:, j]) & (t1[:, i] >= 0)
            n += m
            if i == 0 and j == 0:
                kk |= m
            elif i == 0 or j == 0:
                kpi |= m
            else:
                pipi |= m
    return n, kk, kpi, pipi


def _shared_detail(P: pd.DataFrame, sh: np.ndarray, counts: dict):
    """Akumuluje w counts['shared_detail'][kanał][region] liczniki par ze wspólnym śladem (przed wetem), region wg okien nóg."""
    n, kk, kpi, pipi = shared_track_detail(P)
    centre = (P.win1 & P.win2).to_numpy(); arms = (P.win1 ^ P.win2).to_numpy()
    reg = np.where(centre, "centre", np.where(arms, "arms", "sbsb"))
    a, b = P.decay_type1.to_numpy(), P.decay_type2.to_numpy()
    d = counts.setdefault("shared_detail", {})
    for ch, (ta, tb) in CHANNELS.items():
        ta, tb = min(ta, tb), max(ta, tb)
        s = (a == ta) & (b == tb)
        if not s.any():
            continue
        for r in ("sbsb", "arms", "centre"):
            sr = s & (reg == r)
            if not sr.any():
                continue
            m = sr & sh
            e = d.setdefault(ch, {}).setdefault(r, {"n_pairs": 0, "n_shared": 0, "shared_1track": 0, "shared_2tracks": 0,
                                                     "shared_K_K": 0, "shared_K_pi": 0, "shared_pi_pi": 0})
            e["n_pairs"] += int(sr.sum()); e["n_shared"] += int(m.sum())
            e["shared_1track"] += int((m & (n == 1)).sum()); e["shared_2tracks"] += int((m & (n >= 2)).sum())
            e["shared_K_K"] += int((m & kk).sum()); e["shared_K_pi"] += int((m & kpi).sum()); e["shared_pi_pi"] += int((m & pipi).sum())


def build_pair_hists(L: pd.DataFrame, I: Inputs, region="blind", diff=None):
    """region: full (odsłonięte) | blind (bez centrum) | sbsb. Zwraca ({ch: array}, liczniki).
    PD: `diff` = specyfikacja binów różniczkowych (`pairs.diff.spec_from_cfg`) → zwraca trzeci element {`<ch>__<k>`: array}
    (te same pary, maska binu PO wecie śladów i ślepocie; pary poza binami/nieokreślone liczone w counts['diff_undefined']).
    diff=None → zachowanie i wynik bitowo jak dotąd."""
    L = leg_bins(L, I)
    out = {ch: np.zeros((len(cat_pairs(I.ncat, *CHANNELS[ch])), I.nbins[CHANNELS[ch][0]], I.nbins[CHANNELS[ch][1]])) for ch in CHANNELS}
    idx = {ch: {cp: k for k, cp in enumerate(cat_pairs(I.ncat, *CHANNELS[ch]))} for ch in CHANNELS}
    counts = {"pairs_kept": 0, "region": region}
    if diff is not None:
        from .diff import diff_bin, nbins, hkey
        nb = nbins(diff)
        out_d = {hkey(ch, k): np.zeros_like(out[ch]) for ch in CHANNELS for k in range(nb)}
        counts["diff_kept"] = {str(k): 0 for k in range(nb)}
        counts["diff_undefined"] = 0
    for sf, g in L.groupby("source_file", sort=False):
        P = g.merge(g, on=["source_file", "entry"], suffixes=("1", "2"))
        P = P[(P.decay_type1 < P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.cand1 < P.cand2))]
        sh = shared_track(P)
        # PF3-FIN (05.11): w trybie ślepym liczniki diagnostyczne BEZ centrum (SR×SR) — wcześniej before−shared−kept
        # ujawniało liczbę par w centrum. Pełne liczniki tylko przy region == "full" (odsłonięcie).
        vis = np.ones(len(P), bool) if region == "full" else ~(P.win1 & P.win2).to_numpy()
        counts["pairs_shared_track"] = counts.get("pairs_shared_track", 0) + int((sh & vis).sum())  # PF: diagnostyka weta
        counts["pairs_before_veto"] = counts.get("pairs_before_veto", 0) + int(vis.sum())
        _shared_detail(P[vis], sh[vis], counts)                                              # PF3: rozbicie (bez centrum, gdy ślepo)
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
        if diff is not None:
            db = diff_bin(P, diff)                                   # zmienne symetryczne — kolejność nóg (swap) bez znaczenia
            counts["diff_undefined"] += int((db < 0).sum())
            for k in range(nb):
                counts["diff_kept"][str(k)] += int((db == k).sum())
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
            if diff is not None:
                dbs = db[s]
                for k in range(nb):
                    m = dbs == k
                    if m.any():
                        np.add.at(out_d[hkey(ch, k)], (kk[m], b1[s][m], b2[s][m]), 1)
    if diff is not None:
        return out, counts, out_d
    return out, counts
