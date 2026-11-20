#!/usr/bin/env python3
"""mult_extract.py — MULT_charm: jeden plik ROOT CharmAnalysis (v2.0; zgodność v1.1) → npz z trzema tabelami.

  ev_*  per zdarzenie (WSZYSTKIE wpisy drzewa): entry, evnum, lb, npv, npu, pvz, n500, nbl, leadbl [GeV], lead_trk [GeV],
        n_trk500 (Tracks, |η|<2,5), nA/nC (Tracks pT>0,5, η>+0,5 / η<−0,5), xiA, xiC, grl, trig: tp_<chain>, ps_<chain>, tbp_<L1>;
        MC: nt_D0/nt_Dplus/nt_Dstar (prawdziwe D, pT 3–30, |η|<2,5, z D z b), nt_c/nt_cbar (D0+D⁺: znak c) → pary prawdy.
  c_*   kandydaci po wstępnych cięciach PRE (mult_common): entry, sp (10/4/1), m, dm, md0, pt, eta, lxy, chi2, impsig, cosp,
        dsveto (D⁺: min(m−mKπ)), ndau_acc (córki pT>0,5 |η|<2,5), ndau_A/ndau_C (córki w 0,5<η<2,5 / −2,5<η<−0,5),
        trk0..2 (trackId; pary / weto wspólnych torów), MC: true (logika etap1 mb_truth.py), tpdg.
  t_*   MC: prawdziwe D w akceptacji (entry, pdg, pt, eta, fromB) — dla E_opp na poziomie prawdy.

Dane: NIE liczy żadnych par. Pary danych = osobny krok za bramką zgody (PROCEDURA.md §G5).
  python3 mult_extract.py --kind mc|data IN.root OUT.npz [--max-events N] [--step 50000] [--tag 000021]
"""
import argparse
import json
import os
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mult_common as C  # noqa: E402


def tracks_block(F, start, stop):
    tk = F["Tracks"]; kk = set(tk.keys())
    bpt, beta = C.res(kk, "Tracks_pt"), C.res(kk, "Tracks_eta")
    a = tk.arrays([bpt, beta], entry_start=start, entry_stop=stop)
    pt, eta = a[bpt] / 1000.0, a[beta]
    acc = abs(eta) < 2.5
    lead = ak.to_numpy(ak.fill_none(ak.max(pt[acc], axis=1), 0.0)).astype(np.float64)
    h = pt > C.TRK_PT_MIN
    n500 = ak.to_numpy(ak.sum(h & acc, axis=1)).astype(np.int32)
    nA = ak.to_numpy(ak.sum(h & (eta > C.ETA_GAP) & (eta < 2.5), axis=1)).astype(np.int32)
    nC = ak.to_numpy(ak.sum(h & (eta < -C.ETA_GAP) & (eta > -2.5), axis=1)).astype(np.int32)
    return lead, n500, nA, nC


def cand_block(a, B, is_mc, off):
    g = lambda b: a[B[b]]  # noqa: E731
    dt = g("decayType"); pt = g("pt") / 1000.0; eta = g("eta"); m = g("m") / 1000.0
    acc = (pt > C.PT_MIN) & (pt < C.PT_MAX) & (abs(eta) < C.ETA_MAX)
    lxy, chi2, isig, cosp = g("fitOutput__Lxy"), g("fitOutput__Chi2"), g("fitOutput__ImpactSignificance"), g("cosThetaP")
    dsveto = ak.where(dt == 4, np.minimum(g("m") - g("mKpi1"), g("m") - g("mKpi2")) / 1000.0, 9.0)
    pD, p0, pS = C.PRE["Dplus"], C.PRE["D0"], C.PRE["Dstar"]
    sD = (dt == 4) & (lxy > pD["lxy_min"]) & (chi2 < pD["chi2_max"]) & (cosp > pD["cosp_min"]) & (m > pD["m_lo"]) & (m < pD["m_hi"])
    s0 = (dt == 10) & (lxy > p0["lxy_min"]) & (chi2 < p0["chi2_max"]) & (m > p0["m_lo"]) & (m < p0["m_hi"])
    d0i = g("D0Index"); isds = (dt == 1) & (d0i >= 0)
    d0i_safe = ak.where(isds, d0i, 0)
    md0 = ak.where(isds, g("m")[d0i_safe] / 1000.0, -1.0)
    dm = ak.where(isds, g("DeltaMass") / 1000.0, -1.0)
    sS = isds & (abs(md0 - C.M_D0) < pS["md0_win"]) & (dm > pS["dm_lo"]) & (dm < pS["dm_hi"])
    sel = acc & (sD | s0 | sS)
    # córki: akceptacja estymatorów
    dpt = a[B["dau_pt"]] / 1000.0; deta = a[B["dau_eta"]]
    inacc = (dpt > C.TRK_PT_MIN) & (abs(deta) < 2.5)
    ndau = ak.sum(inacc, axis=2)
    ndA = ak.sum(inacc & (deta > C.ETA_GAP), axis=2); ndC = ak.sum(inacc & (deta < -C.ETA_GAP), axis=2)
    trk = ak.fill_none(ak.pad_none(a[B["dau_trackId"]], 3, axis=2, clip=True), -1)
    ent = ak.broadcast_arrays(off + np.arange(len(dt)), dt)[0]
    fl = lambda x: ak.to_numpy(ak.flatten(x[sel]))  # noqa: E731
    out = {"entry": fl(ent).astype(np.int64), "sp": fl(dt).astype(np.int8), "m": fl(m).astype(np.float32),
           "dm": fl(dm).astype(np.float32), "md0": fl(md0).astype(np.float32), "pt": fl(pt).astype(np.float32),
           "eta": fl(eta).astype(np.float32), "lxy": fl(lxy).astype(np.float32), "chi2": fl(chi2).astype(np.float32),
           "impsig": fl(isig).astype(np.float32), "cosp": fl(cosp).astype(np.float32), "dsveto": fl(dsveto).astype(np.float32),
           "ndau_acc": fl(ndau).astype(np.int8), "ndau_A": fl(ndA).astype(np.int8), "ndau_C": fl(ndC).astype(np.int8),
           **{f"trk{i}": ak.to_numpy(ak.flatten(trk[:, :, i][sel])).astype(np.int32) for i in range(3)}}
    if is_mc:   # = etap1 mb_truth.py (poprawne kandydaty) + D* przez D0Index
        tc = g("TruthMatchCategory")
        ok = ak.all(a[B["dau_pdgId"]] == a[B["dau_truthPdgId"]], axis=2) & (np.sign(g("pdgId")) == np.sign(g("truthPdgId")))
        corr = ((tc == 0) | (tc == 1)) & ok
        true = ak.where(isds, corr[d0i_safe], corr)
        out["true"] = fl(true).astype(np.int8)
        out["tpdg"] = fl(g("truthPdgId")).astype(np.int32)
    return out


def truth_block(a, keys, off):
    K = {k: C.res(keys, v) for k, v in C.TRUTH_BR.items()}
    pdg = a[K["pdg"]]; pt = a[K["pt"]] / 1000.0; eta = a[K["eta"]]
    acc = (pt > C.PT_MIN) & (pt < C.PT_MAX) & (abs(eta) < C.ETA_MAX)
    ap = abs(pdg)
    ev = {f"nt_{sp}": ak.to_numpy(ak.sum(acc & (ap == C.PDG[sp]), axis=1)).astype(np.int16) for sp in C.SPECIES}
    dd = acc & ((ap == 421) | (ap == 411))            # D = D0 i D⁺ (D0 z D* liczony jako D0), jak AFP_ddpred
    ev["nt_c"] = ak.to_numpy(ak.sum(dd & (pdg > 0), axis=1)).astype(np.int16)
    ev["nt_cbar"] = ak.to_numpy(ak.sum(dd & (pdg < 0), axis=1)).astype(np.int16)
    s = acc & ((ap == 421) | (ap == 411) | (ap == 413))
    ent = ak.broadcast_arrays(off + np.arange(len(pdg)), pdg)[0]
    t = {"entry": ak.to_numpy(ak.flatten(ent[s])).astype(np.int64), "pdg": ak.to_numpy(ak.flatten(pdg[s])).astype(np.int32),
         "pt": ak.to_numpy(ak.flatten(pt[s])).astype(np.float32), "eta": ak.to_numpy(ak.flatten(eta[s])).astype(np.float32),
         "fromB": ak.to_numpy(ak.flatten(a[K["fromB"]][s])).astype(np.int8)}
    return ev, t


def run(path, out, kind, max_events, step, tag):
    is_mc = kind == "mc"
    F = uproot.open(path); t = F["CharmAnalysis"]; keys = set(t.keys())
    if F["Tracks"].num_entries != t.num_entries:
        sys.exit(f"Tracks ({F['Tracks'].num_entries}) ≠ CharmAnalysis ({t.num_entries}) — drzewa niezsynchronizowane")
    ver = C.detect_version(keys)
    EB = C.resolve_event_branches(keys)
    TB = C.resolve_trigger_branches(keys) if not is_mc else {}
    CB = C.resolve_cand_branches(keys, is_mc)
    compat = [k for k in (("nselbl", "leadbl") if is_mc else ("nselbl", "leadbl", "grl")) if EB[k] is None]
    n_all = t.num_entries if max_events is None else min(t.num_entries, max_events)
    br = sorted({v for v in EB.values() if v} | set(TB.values()) | set(CB.values()) |
                ({C.res(keys, v) for v in C.TRUTH_BR.values()} if is_mc else set()))
    E, Cc, T = [], [], []
    off = 0
    for a in t.iterate(br, step_size=step, entry_stop=n_all, library="ak"):
        n = len(a)
        lead, n500t, nA, nC = tracks_block(F, off, off + n)
        get = lambda k, dt, fill=0: (ak.to_numpy(a[EB[k]]).astype(dt) if EB[k] else np.full(n, fill, dt))  # noqa: E731
        ev = {"entry": off + np.arange(n, dtype=np.int64), "evnum": get("evnum", np.uint64), "lb": get("lb", np.int32),
              "npv": get("npv", np.int16), "npu": get("npu", np.int16), "pvz": get("pvz", np.float32),
              "n500": get("n500", np.int16), "lead_trk": lead.astype(np.float32), "n_trk500": n500t.astype(np.int16),
              "nA": nA.astype(np.int16), "nC": nC.astype(np.int16),
              "xiA": get("xiA", np.float32), "xiC": get("xiC", np.float32), "grl": get("grl", np.int8, 1)}
        # v2.0: nSelBL, leadTrackPtBL [MeV w ROOT]; v1.1: zastępstwa (flaga compat w meta)
        ev["nbl"] = get("nselbl", np.int16) if EB["nselbl"] else ev["n500"].copy()
        ev["leadbl"] = (get("leadbl", np.float32) / 1000.0).astype(np.float32) if EB["leadbl"] else ev["lead_trk"].copy()
        for k, b in TB.items():
            ev[k] = ak.to_numpy(a[b]).astype(np.float32 if k.startswith("ps_") else np.int8)
        if is_mc:
            tev, tt = truth_block(a, keys, off)
            ev.update(tev); T.append(tt)
        E.append(ev)
        Cc.append(cand_block(a, CB, is_mc, off))
        off += n
    cat = lambda L: {c: np.concatenate([x[c] for x in L]) for c in L[0]} if L else {}  # noqa: E731
    z = {f"ev_{k}": v for k, v in cat(E).items()}
    z.update({f"c_{k}": v for k, v in cat(Cc).items()})
    z.update({f"t_{k}": v for k, v in cat(T).items()})
    meta = {"kind": kind, "source": path, "tag": tag, "format": ver, "compat_missing": compat, "def": C.DEF_VERSION,
            "n_events": int(off), "n_entries_file": int(t.num_entries), "triggers": sorted(TB)}
    np.savez_compressed(out, meta=json.dumps(meta), **z)
    msg = f"OK {out}: {ver}{' (compat: ' + ','.join(compat) + ')' if compat else ''}, {off} zdarzeń, kandydaci "
    msg += ", ".join(f"{sp}={int((z['c_sp'] == C.SP_CODE[sp]).sum())}" for sp in C.SPECIES)
    if is_mc:
        msg += ", prawdziwe (akc.) " + ", ".join(f"{sp}={int(z[f'ev_nt_{sp}'].sum())}" for sp in C.SPECIES)
        msg += f", poprawne reko {int(z['c_true'].sum())}"
    print(msg, file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--kind", choices=["mc", "data"], required=True)
    ap.add_argument("--max-events", type=int, default=None)
    ap.add_argument("--step", type=int, default=50000)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    run(a.inp, a.out, a.kind, a.max_events, a.step, a.tag)


if __name__ == "__main__":
    main()
