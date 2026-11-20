#!/usr/bin/env python3
"""make_synth_v20.py — atrapa pliku ROOT w formacie v2.0 z lokalnego pliku MC v1.1 (test MULT_charm bez sieci).

Przenosi drzewa CharmAnalysis (wybrane gałęzie) i Tracks z pliku v1.1, a następnie dokłada/zmienia jak w charm-prod-v2.0
(CharmAnalysis/python/Branches.py @ 6007843):
  + EventInfo_nSelBL_NOSYS (int)           ≈ liczba śladów Tracks pT>0,5 |η|<2,5 z losową utratą 8 % (atrapa selekcji BL, |d0|<1,5 mm)
  + EventInfo_leadTrackPtBL_NOSYS (float)  = max pT tych śladów [MeV] (jak v2.0: MeV)
  + EventInfo_xi_cal_A/C_NOSYS             (zamiast xi_cal_plus/minus — zamienione strony, H8)
  − EventInfo_xi_cal_plus/minus*            (usunięte w v2.0)
  --data: usuwa gałęzie prawdy MC, dodaje EventInfo_passGRL_MinBiasAFP (char), trigPassed_<chain>_NOSYS, trigPS_<chain>, trigTBP_<L1>
          (emulacja: sptrk_ptX odpala, gdy leadBL > X GeV z ε 0,97 i PS 25/50 wg LB; mbts: PS 1 i ε 0,99; J12: leadBL>3 z ε 0,5).
  python3 make_synth_v20.py IN_v11.root OUT.root [--data] [--max-events N] [--seed S]
"""
import argparse

import awkward as ak
import numpy as np
import uproot

DROP_V20 = ("EventInfo_xi_cal_plus", "EventInfo_xi_cal_minus")
MC_PREFIX = ("TruthDmesons_", "EventInfo_mcEventWeights", "EventInfo_weight_mc", "EventInfo_mcChannelNumber",
             "EventInfo_pdfInfo", "EventInfo_weight_beamspot", "EventInfo_jvt", "EventInfo_fjvt")
MC_DMES = ("DMesons_TruthMatchCategory", "DMesons_truthPdgId", "DMesons_daughterInfo__truthPdgId", "DMesons_daughterInfo__truthParentPdgId",
           "DMesons_MostPopular", "DMesons_truth", "DMesons_daughter_truthMatchProb")
KEEP_DMES = ("decayType", "pdgId", "pt", "eta", "phi", "m", "DeltaMass", "D0Index", "mKpi1", "mKpi2", "fitOutput__Lxy",
             "fitOutput__Chi2", "fitOutput__ImpactSignificance", "cosThetaP", "daughterInfo__pt", "daughterInfo__eta",
             "daughterInfo__trackId", "daughterInfo__pdgId", "daughterInfo__truthPdgId", "TruthMatchCategory", "truthPdgId")
CHAINS = {"HLT_mb_mbts_L1MBTS_1": None, "HLT_noalg_mb_L1MBTS_2": None, "HLT_mb_sptrk_pt4_L1MBTS_2": 4.0,
          "HLT_mb_sptrk_pt2_L1MBTS_2": 2.0, "HLT_noalg_L1J12": "j12", "HLT_noalg_L1AFP_A_OR_C": "afp"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--data", action="store_true"); ap.add_argument("--max-events", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    F = uproot.open(a.inp); t = F["CharmAnalysis"]; n = min(a.max_events, t.num_entries)
    keys = []
    for k in t.keys():
        base = k[:-6] if k.endswith("_NOSYS") else k
        if base.startswith(DROP_V20):
            continue
        if base.startswith("DMesons_") and base[8:] not in KEEP_DMES:
            continue
        if base.startswith(("AfpProtons_",)):
            continue
        if a.data and (base.startswith(MC_PREFIX) or base.startswith(MC_DMES)):
            continue
        if base.startswith("TruthDmesons_") and base not in ("TruthDmesons_pdgId", "TruthDmesons_pt", "TruthDmesons_eta",
                                                              "TruthDmesons_fromBottom"):
            continue
        keys.append(k)
    A = t.arrays(keys, entry_stop=n)
    tk = F["Tracks"].arrays(entry_stop=n)
    pt, eta = tk["Tracks_pt_NOSYS"], tk["Tracks_eta_NOSYS"]
    keep = rng.random(len(ak.flatten(pt))) > 0.08
    keep = ak.unflatten(keep, ak.num(pt))
    bl = (pt > 500) & (abs(eta) < 2.5) & keep
    nbl = ak.to_numpy(ak.sum(bl, axis=1)).astype(np.int32)
    lead = ak.to_numpy(ak.fill_none(ak.max(pt[bl], axis=1), 0.0)).astype(np.float32)   # MeV
    out = {k: A[k] for k in keys}
    out["EventInfo_nSelBL_NOSYS"] = nbl
    out["EventInfo_leadTrackPtBL_NOSYS"] = lead
    out["EventInfo_xi_cal_A_NOSYS"] = rng.exponential(0.01, n).astype(np.float32)
    out["EventInfo_xi_cal_C_NOSYS"] = rng.exponential(0.01, n).astype(np.float32)
    if a.data:
        lb = np.sort(rng.integers(4, 324, n)).astype(np.uint32)
        out["EventInfo_lumiBlock"] = lb
        out["EventInfo_runNumber"] = np.full(n, 428770, np.uint32)
        out["EventInfo_passGRL_MinBiasAFP"] = (rng.random(n) > 0.01).astype(np.int8)
        g = lead / 1000.0
        for c, thr in CHAINS.items():
            if thr is None:
                fire, ps = rng.random(n) < 0.99, np.ones(n)
            elif thr == "j12":
                fire, ps = (g > 3.0) & (rng.random(n) < 0.5), np.ones(n)
            elif thr == "afp":
                fire, ps = rng.random(n) < 0.02, np.where(lb <= 16, 120.0, 60.0)
            else:
                ps = np.where(lb <= 20, 25.0, 50.0) if thr == 4.0 else np.where(lb <= 20, 625.0, 1250.0)
                fire = (g > thr) & (rng.random(n) < 0.97) & (rng.random(n) < 1.0 / ps * 25)   # prescale (×25 → więcej statystyki w teście)
                ps = ps / 25
            out[f"EventInfo_trigPassed_{c}_NOSYS"] = fire.astype(np.int8)
            out[f"EventInfo_trigPS_{c}"] = ps.astype(np.float32)
        for l1 in ("L1_MBTS_1", "L1_MBTS_2", "L1_J12"):
            out[f"EventInfo_trigTBP_{l1}"] = (rng.random(n) < 0.9).astype(np.int8)
    trk = {k: tk[k] for k in tk.fields}
    with uproot.recreate(a.out) as f:
        f["CharmAnalysis"] = out
        f["Tracks"] = trk
    print(f"OK {a.out}: {n} zdarzeń, {'dane' if a.data else 'MC'}, {len(out)} gałęzi")


if __name__ == "__main__":
    main()
