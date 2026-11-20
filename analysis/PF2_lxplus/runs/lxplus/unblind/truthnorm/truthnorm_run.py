#!/usr/bin/env python3
"""TRUTHNORM (19.11.2026, projekt A; decyzja z 19.11 „zgoda”: normalizacja par z prawdy mb_nd × mapy direct, GEN_mcstat/PROCEDURA.md §5).

n_pair(ch) → n_pair(ch) · f_ch,  f_ch = [pary/D*]_truth,mb_nd(w2) / [pary/D*]_direct(w2)  (GEN_mcstat/out/truth_mbnd.json: truth_self_over_direct;
direct ≡ nominał 05.11: n_pair/N_D* zgodne z direct_pair_per_dstar_w do < 0,6 %). Tło (n_real, SB×SB, dataleg) bez zmian → pełny ponowny fit
(ta sama maszyneria: ub_jk.replica_fit, charmtag z PF_pairs/src). Kontrola: f ≡ 1 odtwarza μ nominalne.
427929: f jak w 428770 (metoda stosunku UB427: c_w(427929) = c_w(nom)·R_X — R_X wspólny dla prawdy i direct; przybliżenie, zob. raport).
  truthnorm_run.py --what nom428|nom427|cls:<i> --out DIR [--f-one]"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(CW, "Q_lxplus", "ub_trig")); sys.path.insert(0, os.path.join(CW, "Q_lxplus", "ub_jk"))
SP = os.environ.get("TN_SP", "")
P428 = os.path.join(CW, "PF2_lxplus", "runs", "lxplus", "unblind", "pairs")
P427 = os.path.join(SP, "p427")
TRUTH = json.load(open(os.path.join(CW, "GEN_mcstat", "out", "truth_mbnd.json")))
F = TRUTH["truth_self_over_direct"]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def scale_npair(I, f):
    import copy
    I2 = copy.copy(I)
    I2.n_pair = {k: v * f[k] for k, v in I.n_pair.items()}
    return I2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--what", required=True); ap.add_argument("--out", required=True); ap.add_argument("--f-one", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    import ubjk_run as U
    charmtag, CHANNELS, Model, load_inputs, fit_channel, Fm = U._charmtag()
    f = {ch: 1.0 for ch in CHANNELS} if a.f_one else {ch: float(F[ch]) for ch in CHANNELS}
    res = {"what": a.what, "f": f, "charmtag": charmtag.__file__}
    if a.what in ("nom428", "nom427"):
        pd_ = P428 if a.what == "nom428" else P427
        I = load_inputs(os.path.join(pd_, "pair_inputs"))
        if a.what == "nom428":
            H = dict(np.load(os.path.join(pd_, "pairs_full.npz")))
        else:
            A = dict(np.load(os.path.join(pd_, "perfile.npz")))
            H = {k[2:]: v.sum(0) for k, v in A.items() if k.startswith("H_")}
        res["n_pair_nominal"] = dict(I.n_pair)
        I2 = scale_npair(I, f)
        res["n_pair"] = dict(I2.n_pair)
        res["fit"] = U.replica_fit(Model(I2), H, list(CHANNELS), combined=True, interval=True)
    else:
        import ubtrig_run as T
        import pandas as pd
        import btrig_common as C
        i = int(a.what.split(":")[1]); c = T.CLASSES[i]
        I = load_inputs(os.path.join(P428, "pair_inputs"))
        rt = json.load(open(os.path.join(CW, "PF2_lxplus", "runs", "lxplus", "unblind", "ubtrig", "ubtrig.json")))
        r = rt["dstar_fraction"][c]["r"]; rerr = rt["dstar_fraction"][c]["r_err"]
        tkey, tcls = T.load_trg([os.path.join(CW, "B_trig", "runs", "lxplus", "trg_data_428770.npz")])
        L = pd.read_parquet(os.path.join(SP, "..", "ubtrig", "legs_data.parquet"))
        lk = C.ev_key(pd.to_numeric(L.source_file).to_numpy(np.int64), L.entry.to_numpy(np.int64))
        pos = np.clip(np.searchsorted(tkey, lk), 0, len(tkey) - 1); hit = tkey[pos] == lk
        L = L.assign(file_idx=np.where(hit, tcls[pos], len(T.CLASSES)).astype(np.int64))
        H, DL, SB, keys = U.perfile_from_legs(L, I, len(T.CLASSES) + 1)
        fev = float((tcls == i).mean())
        res.update({"class": c, "r": r, "r_err": rerr, "fev": fev, "legs_unmatched": int((~hit).sum()),
                    "pairs_class": {ch: int(H[ch][i].sum()) for ch in CHANNELS}, "pairs_class_ubtrig": {ch: rt["pairs_per_class"][ch][c] for ch in CHANNELS}})
        Ic = T.class_inputs(scale_npair(I, f), r, fev, DL, SB, keys, i)
        res["fit"] = U.replica_fit(Model(Ic), {ch: H[ch][i] for ch in CHANNELS}, list(CHANNELS), combined=True, interval=True)
    json.dump(U._clean(res), open(os.path.join(a.out, f"tn_{a.what.replace(':', '')}{'_f1' if a.f_one else ''}.json"), "w"), indent=1)
    cm = res["fit"]["combined"]
    log(a.what, "f1" if a.f_one else "truth", {ch: round(res["fit"][ch]["mu_hat"], 4) for ch in CHANNELS}, "wspólne", round(cm["mu"], 4), cm.get("mu_68"), "Z", round(cm["Z"], 2))


if __name__ == "__main__":
    main()
