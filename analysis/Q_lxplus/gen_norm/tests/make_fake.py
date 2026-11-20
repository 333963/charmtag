#!/usr/bin/env python3
"""Test lokalny gen_norm: atrapy ROOT v1.1 (CharmAnalysis + Tracks) zbudowane z PRAWDZIWYCH lokalnych części h5 v1.2 mb_nd
(kandydaci D: decayType, truth, córki, trackId — bez zmian; n_trk, ślady, prawdziwe D bez dopasowania = losowe) + atrapy danych
428770 (tylko n_trk, ślady) + atrapy scores (charmtag apply z losowym score) i config. Nic nie czyta z EOS/lxplus.
  python3 tests/make_fake.py H5_PART [H5_PART ...] --out DIR"""
import argparse
import os
import re
import sys

import awkward as ak
import h5py
import numpy as np
import uproot

PDG = {10: 421, 4: 411, 1: 413}
ND = {10: 2, 4: 3, 1: 3}


def fake_mc(h5p, out_root, rng, extra_frac=0.3):
    e = h5py.File(h5p, "r")["events"][:]
    evn = e["event_number"].astype(np.uint64)
    brk = np.flatnonzero(np.r_[True, evn[1:] != evn[:-1]])
    counts = np.diff(np.r_[brk, len(evn)])
    ev_ids = evn[brk]
    if len(np.unique(ev_ids)) != len(ev_ids):
        sys.exit(f"{h5p}: zdarzenia nie są ciągłe w pliku")
    n_ev = len(ev_ids)
    dt = e["Dmeson_decayType"].astype(np.int32)
    tc = e["Dmeson_truth_category"].astype(np.int32)
    tpdg = e["Dmeson_truthPdgId"].astype(np.int32)
    q = np.where(e["Dmeson_charge"] == 0, 1, np.sign(e["Dmeson_charge"]))
    pdg = np.vectorize(PDG.get)(dt) * np.where(tpdg != 0, np.sign(tpdg), q)
    dm = e["Dmeson_delta_mass"].astype(np.float64)
    trk = np.stack([e[f"Dmeson_daughter{i}_trackId"] for i in range(3)], 1)
    dp = np.stack([e[f"Dmeson_daughter{i}_pdgId"] for i in range(3)], 1)
    dtp = np.stack([e[f"Dmeson_daughter{i}_truthPdgId"] for i in range(3)], 1)
    ev_of = np.repeat(np.arange(n_ev), counts); loc = np.arange(len(e)) - np.repeat(brk, counts)
    # D0Index: D* → kandydat D0 w tym samym zdarzeniu z tymi samymi trackId0/1 i hipotezą córki 0
    d0i = np.full(len(e), -999, np.int32)
    key = {}
    for r in np.flatnonzero(dt == 10):
        key.setdefault((ev_of[r], trk[r, 0], trk[r, 1], dp[r, 0]), loc[r])
    for r in np.flatnonzero(dt == 1):
        d0i[r] = key.get((ev_of[r], trk[r, 0], trk[r, 1], dp[r, 0]), -999)
    nd = np.vectorize(ND.get)(dt)
    nest = lambda X: ak.unflatten(ak.unflatten(np.concatenate([X[r, :nd[r]] for r in range(len(e))]), nd), counts)  # noqa: E731
    jag = lambda x: ak.unflatten(x, counts)  # noqa: E731
    # prawdziwe D: jedna na zestaw torów kandydatów z truth 0/1 (bliźniaki K↔π → to samo D) + losowe bez dopasowania
    T = {k: [[] for _ in range(n_ev)] for k in ("pdg", "pt", "eta", "phi", "dt", "fromB", "idx")}
    grp = {}
    for r in np.flatnonzero(np.isin(tc, [0, 1])):
        k = (ev_of[r], dt[r], *sorted(trk[r, :nd[r]]))
        grp.setdefault(k, []).append(r)
    for k, rows in grp.items():
        r = rows[0]; i = ev_of[r]
        T["pdg"][i].append(int(tpdg[r] if tpdg[r] != 0 else pdg[r])); T["pt"][i].append(float(e["Dmeson_pt"][r] * 1000 * (1 + rng.normal(0, 0.02))))
        T["eta"][i].append(float(e["Dmeson_eta"][r] + rng.normal(0, 0.005))); T["phi"][i].append(float(e["Dmeson_phi"][r]))
        T["dt"][i].append(int(dt[r] if rng.random() > 0.05 else -1)); T["fromB"][i].append(int(tc[r] == 1))
        T["idx"][i].append([int(loc[x]) for x in rows])
    for i in np.flatnonzero(rng.random(n_ev) < 0.5):
        s = rng.choice([10, 4, 1]); T["pdg"][i].append(int(PDG[s] * rng.choice([-1, 1]))); T["pt"][i].append(float(rng.uniform(2.6, 20) * 1000))
        T["eta"][i].append(float(rng.uniform(-2.6, 2.6))); T["phi"][i].append(float(rng.uniform(-3.1, 3.1)))
        T["dt"][i].append(int(s)); T["fromB"][i].append(0); T["idx"][i].append([])
    # zdarzenia bez kandydatów h5 (są w ntuplu, nie ma ich w h5)
    n_x = int(extra_frac * n_ev)
    ev_all = np.r_[ev_ids, ev_ids.max() + 1 + np.arange(n_x, dtype=np.uint64)]
    pad = lambda A: ak.concatenate([A, ak.Array([[]] * n_x)]) if n_x else A  # noqa: E731
    ncand = np.r_[counts, np.zeros(n_x, int)]
    ntrk = rng.poisson(6 + 5 * np.minimum(ncand, 6)).astype(np.int32)
    ntr = ntrk + rng.poisson(3, len(ntrk))
    tpt = rng.exponential(900, ntr.sum()) + 300; teta = rng.uniform(-3, 3, ntr.sum())
    D = {"decayType": dt, "TruthMatchCategory": tc, "pdgId": pdg.astype(np.int32), "truthPdgId": tpdg, "m": e["Dmeson_mass"].astype(np.float32) * 1000,
         "mKpi1": e["Dmeson_mKpi1"], "mKpi2": e["Dmeson_mKpi2"], "D0Index": d0i, "DeltaMass": np.where(dm > 0, dm * 1000, -999).astype(np.float32),
         "pt": e["Dmeson_pt"].astype(np.float32) * 1000, "fitOutput__Lxy": e["Dmeson_Lxy"]}
    ca = {f"DMesons_{k}_NOSYS": pad(jag(v)) for k, v in D.items()}
    ca.update({"DMesons_daughterInfo__pdgId_NOSYS": pad(nest(dp)), "DMesons_daughterInfo__truthPdgId_NOSYS": pad(nest(dtp)),
               "DMesons_daughterInfo__trackId_NOSYS": pad(nest(trk))})
    ca["EventInfo_nTracksPt500Eta25_NOSYS"] = ntrk; ca["EventInfo_eventNumber"] = ev_all
    Ta = lambda k, dtp_: ak.values_astype(ak.Array(T[k] + [[]] * n_x), dtp_)  # noqa: E731
    ca.update({"TruthDmesons_pdgId": Ta("pdg", np.int32), "TruthDmesons_pt": Ta("pt", np.float32), "TruthDmesons_eta": Ta("eta", np.float32),
               "TruthDmesons_phi": Ta("phi", np.float32), "TruthDmesons_decayType": Ta("dt", np.int32), "TruthDmesons_fromBottom": Ta("fromB", np.int8),
               "TruthDmesons_recoMatchIndices": ak.values_astype(ak.Array(T["idx"] + [[]] * n_x), np.int32)})
    with uproot.recreate(out_root) as f:
        f["CharmAnalysis"] = ca
        f["Tracks"] = {"Tracks_pt_NOSYS": ak.unflatten(tpt.astype(np.float32), ntr), "Tracks_eta_NOSYS": ak.unflatten(teta.astype(np.float32), ntr)}
    return n_ev + n_x


def fake_data(out_root, rng, n=30000):
    ntrk = rng.poisson(35, n).astype(np.int32); ntr = ntrk + rng.poisson(4, n)
    tpt = rng.exponential(2000, ntr.sum()) + 300; teta = rng.uniform(-3, 3, ntr.sum())
    with uproot.recreate(out_root) as f:
        f["CharmAnalysis"] = {"EventInfo_nTracksPt500Eta25_NOSYS": ntrk, "EventInfo_eventNumber": np.arange(n, dtype=np.uint64) + 10**9}
        f["Tracks"] = {"Tracks_pt_NOSYS": ak.unflatten(tpt.astype(np.float32), ntr), "Tracks_eta_NOSYS": ak.unflatten(teta.astype(np.float32), ntr)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5", nargs="+"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(7)
    R = os.path.join(a.out, "root"); H = os.path.join(a.out, "h5", "mb_nd")
    for d in (os.path.join(R, "mb_nd"), os.path.join(R, "data_428770"), H):
        os.makedirs(d, exist_ok=True)
    for p in a.h5:
        tg = re.search(r"_(\d{6})_none\.h5$", p).group(1)
        dst = os.path.join(H, os.path.basename(p))
        if not os.path.exists(dst):
            os.symlink(os.path.abspath(p), dst)
        n = fake_mc(p, os.path.join(R, "mb_nd", f"user.fake.51907047.EXT0._{tg}.CharmAnalysis.mc.root"), rng)
        print("MC", tg, n, "zdarzeń")
    for tg in ("000013", "000021", "000030"):
        fake_data(os.path.join(R, "data_428770", f"user.fake.data.EXT0._{tg}.CharmAnalysis.data.root"), rng)
    print("OK", a.out)


if __name__ == "__main__":
    main()
