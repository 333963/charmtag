"""Wejścia zamrożone 09-26/27 (unblinding_freeze_v20260927_01/frozen_inputs) → Inputs. Do testu domknięcia."""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from .model import Inputs, M0, SBW
from . import CHANNELS


def _sf_cum(SFJ, dt, e):
    t = SFJ[str(dt)]
    xs = sorted(float(k) for k in t)
    ys = [t[f"{x:g}"] for x in xs]
    if e < xs[0]:
        return float(ys[0] + (e - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0]))
    if e > xs[-1]:
        return float(ys[-1] + (e - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    return float(np.interp(e, xs, ys))


def _sf_cat(SFJ, T, dt, c):
    if dt == 10 and c == 3:
        return _sf_cum(SFJ, 10, 0.65)
    if dt == 10 and c == 4:
        return (_sf_cum(SFJ, 10, 0.9) * 0.9 - _sf_cum(SFJ, 10, 0.65) * 0.65) / 0.25
    e0, e1 = T["cats"][c]
    return (_sf_cum(SFJ, dt, e1) * e1 - (_sf_cum(SFJ, dt, e0) * e0 if e0 > 0 else 0.0)) / (e1 - e0)


def sbsb_counts(legs: pd.DataFrame, ncat) -> dict:
    """Surowe zliczenia par sideband×sideband (nogi już tylko z sidebandów), jak m7_fit.py."""
    L = legs[np.abs(legs.m - legs.decay_type.map(M0)) > legs.decay_type.map(SBW)]
    P = L.merge(L, on=["source_file", "entry"], suffixes=("1", "2"))
    P = P[(P.decay_type1 > P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.cand1 < P.cand2))]
    P = P[~shared_track(P)]
    out = {}
    for ch, (a, b) in CHANNELS.items():
        q = P[(P.decay_type1 == max(a, b)) & (P.decay_type2 == min(a, b))]
        for i in range(ncat[a]):
            for j in range(ncat[b]):
                if a == b and j < i:
                    continue
                if a == b:
                    s = ((q.cat1 == i) & (q.cat2 == j)) | ((q.cat1 == j) & (q.cat2 == i))
                else:
                    ca = np.where(q.decay_type1 == a, q.cat1, q.cat2)
                    cb = np.where(q.decay_type1 == a, q.cat2, q.cat1)
                    s = (ca == i) & (cb == j)
                out[(ch, i, j)] = int(np.sum(s))
    return out


def shared_track(P):
    t1 = P[["trk01", "trk11", "trk21"]].to_numpy()
    t2 = P[["trk02", "trk12", "trk22"]].to_numpy()
    sh = np.zeros(len(P), bool)
    for i in range(3):
        for j in range(3):
            sh |= (t1[:, i] == t2[:, j]) & (t1[:, i] >= 0)
    return sh


def load_frozen(freeze_dir, kk_note="R"):
    F = os.path.join(freeze_dir, "frozen_inputs") if os.path.isdir(os.path.join(freeze_dir, "frozen_inputs")) else freeze_dir
    T = json.load(open(os.path.join(F, "templates_v3.json")))
    A = np.load(os.path.join(F, "templates_v3.npz"))
    MB = json.load(open(os.path.join(F, "mb_rates_2d.json")))
    SFJ = json.load(open(os.path.join(F, "sf_for_pairs_v2.json")))
    ncat = {10: len(T.get("cats10", T["cats"])), 4: len(T["cats"])}
    R = 7.65e-4 / MB["dstar_true_per_event"]["w"]
    NDATA = T.get("ndata", 7_989_496)
    sig = {(dt, c): A[f"sig_{dt}_{c}"].astype(float) for dt in (10, 4) for c in range(ncat[dt])}
    dl = {(dt, c): A[f"dataleg_{dt}_{c}"].astype(float) for dt in (10, 4) for c in range(ncat[dt])}
    Y = {(dt, c): T[f"y_{dt}_{c}"] * _sf_cat(SFJ, T, dt, c) for dt in (10, 4) for c in range(ncat[dt])}
    SINGLE = {10: MB["D0_correct_per_event"]["w"], 4: MB["Dplus_correct_per_event"]["w"]}
    PR = {"D0D0": MB["pairs_D0D0"]["per_event_w"], "DpD0": MB["pairs_DpD0"]["per_event_w"], "DpDp": MB["pairs_DpDp"]["per_event_w"]}
    legs = pd.concat([pd.read_parquet(os.path.join(F, T["legs_files"]["10"])), pd.read_parquet(os.path.join(F, T["legs_files"]["4"]))],
                     ignore_index=True)
    return Inputs(edges=T["edges"], ncat=ncat, sig=sig, sig_comp={}, Y=Y, dataleg=dl, nev_legs=float(T["legs_events"]),
                  ndata=float(NDATA), n_real={dt: NDATA * SINGLE[dt] * R for dt in SINGLE},
                  n_pair={ch: NDATA * PR[ch] * R for ch in PR}, sbsb=sbsb_counts(legs, ncat),
                  meta={"source": "frozen 2026-09-26 (M7 v3)", "R": R, "ndata": NDATA, "tag_sf": T.get("tag_sf"),
                        "mb_rates": MB, "sf": SFJ, "cats": T.get("cats10", T["cats"])})


def frozen_legs_all(freeze_dir):
    F = os.path.join(freeze_dir, "frozen_inputs") if os.path.isdir(os.path.join(freeze_dir, "frozen_inputs")) else freeze_dir
    return pd.concat([pd.read_parquet(os.path.join(F, f"legs_all_v3_{dt}.parquet")) for dt in (10, 4)], ignore_index=True)
