#!/usr/bin/env python3
"""gn_root.py — gen_norm (16.11.2026): jeden plik ROOT v1.1 (drzewa CharmAnalysis + Tracks, zsynchronizowane) → npz.

  --kind mc    (mb_nd):
     ev_*   per zdarzenie: evnum, n_trk, lead_pt, n_true_dstar, n_corr_d0, n_corr_dp
            — logika 1:1 z etap1_definicje_masy_v20260925_01/code/mb_truth.py (dla odtworzenia mb_rates_2d.json / c_w),
     c_*    poprawne kandydaty etap1 (entry, decay_type, m, pt, truth_pdg, trk0..2) = mb_truth.py `_corr`,
     t_*    prawdziwe D (effmap_gn.process: entry, evnum, ntrk, pdg, pt, eta, phi, fromB, mode, matched),
     TruthParticles: gdy w pliku są gałęzie/drzewo `TruthParticles*` → n_ch (naładowane, pT > 0,5 GeV, |η| < 2,5) i lead_gen
            per zdarzenie (ev_nch, ev_lead_gen); inaczej meta has_truth_particles = 0 (wtedy pkt (c) pomijany).
  --kind data  (data_428770): ev_evnum, ev_n_trk, ev_lead_pt = etap1 data_events.py (w2d.py nie używa triggerów).
v2.0 (17.11): D0 w `sel` tylko w oknie 1,75 < m < 2,00 GeV (nominał v1.1; v2.0 ma 1,70–2,05); gdy gałęzie są w pliku, dodatkowo ev_nselbl,
  ev_leadbl [GeV], ev_npu, ev_npv, ev_grl, ev_run, ev_lb oraz w danych ev_tp_<łańcuch> / ev_ps_<łańcuch> (TRIG_STORE) — pod w2(nSelBL, leadBL)
  i w2 per trigger. Na v1.1 wynik kolumn etap1 bez zmian (okno D0 v1.1 = 1,75–2,00).

lead_pt = max pT śladu z drzewa Tracks przy |η| < 2,5 (0 gdy brak) [GeV, float64 jak w etap1]; n_trk = EventInfo_nTracksPt500Eta25.
  python3 gn_root.py --kind mc|data IN.root OUT.npz [--max-events N] [--tag 000021]"""
import argparse
import json
import os
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import effmap_gn  # noqa: E402

# okno m(Kπ) D0 nominału (v1.1: 1,75–2,00 GeV; v2.0: 1,70–2,05) — cięcie w run_mc daje zgodność v2.0 z etap1 (na v1.1 nic nie zmienia)
D0_WIN = (1750.0, 2000.0)
# kolumny zdarzeniowe v2.0 (opcjonalne; brak w v1.1 → pomijane): nazwa logiczna → (gałąź, dtype, skala)
EV_EXTRA = {"nselbl": ("EventInfo_nSelBL", np.int16, 1.0), "leadbl": ("EventInfo_leadTrackPtBL", np.float32, 1e-3),   # MeV → GeV
            "npu": ("EventInfo_nPileupVertices", np.int16, 1.0), "npv": ("EventInfo_nPrimaryVertices", np.int16, 1.0),
            "grl": ("EventInfo_passGRL_MinBiasAFP", np.int8, 1.0), "run": ("EventInfo_runNumber", np.uint32, 1.0),
            "lb": ("EventInfo_lumiBlock", np.uint16, 1.0)}
# łańcuchy triggera do zapisu (= MULT_charm/mult_common.TRIG_STORE, kopia — paczka nie importuje MULT): tp_ (trigPassed), ps_ (trigPS)
TRIG_STORE = ["HLT_j20_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
              "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_noalg_L1AFP_A_OR_C", "HLT_noalg_L1J12", "HLT_noalg_L1RD0_FILLED", "HLT_noalg_mb_L1MBTS_1",
              "HLT_noalg_mb_L1MBTS_2"]
B_MC = ["decayType", "TruthMatchCategory", "pdgId", "truthPdgId", "m", "mKpi1", "mKpi2", "D0Index", "DeltaMass", "pt"]


def res(keys, name):
    """Nazwa gałęzi: etap1 używał sufiksu _NOSYS, effmap.py — bez; przyjmij którąkolwiek."""
    for n in (name + "_NOSYS", name):
        if n in keys:
            return n
    raise KeyError(name)


def res_opt(keys, name):
    try:
        return res(keys, name)
    except KeyError:
        return None


def extra_branches(keys):
    """{kolumna ev_*: (gałąź, dtype, skala)} dla gałęzi obecnych w pliku (v2.0) + trigPassed/trigPS z TRIG_STORE."""
    out = {}
    for c, (b, dt, sc) in EV_EXTRA.items():
        r = res_opt(keys, b)
        if r:
            out[c] = (r, dt, sc)
    for ch in TRIG_STORE:
        r = res_opt(keys, f"EventInfo_trigPassed_{ch}")
        if r:
            out[f"tp_{ch}"] = (r, np.int8, 1.0)
        r = res_opt(keys, f"EventInfo_trigPS_{ch}")
        if r:
            out[f"ps_{ch}"] = (r, np.float32, 1.0)
    return out


def read_extra(t, X, start, stop):
    return {c: (t[b].array(library="np", entry_start=start, entry_stop=stop).astype(np.float64) * sc).astype(dt) if sc != 1.0
            else t[b].array(library="np", entry_start=start, entry_stop=stop).astype(dt) for c, (b, dt, sc) in X.items()}


def fmt_of(keys):
    return "v2.0" if res_opt(keys, "EventInfo_nSelBL") else "v1.1"


def lead_and_ntrk(F, t, n_trk_br, start, stop):
    tk = F["Tracks"]; kk = set(tk.keys())
    tr = tk.arrays([res(kk, "Tracks_pt"), res(kk, "Tracks_eta")], entry_start=start, entry_stop=stop)
    pt, eta = tr[res(kk, "Tracks_pt")], tr[res(kk, "Tracks_eta")]
    lead = ak.to_numpy(ak.fill_none(ak.max(pt[abs(eta) < 2.5], axis=1), 0.0)) / 1000.0     # jak etap1 (float32 MeV → /1000)
    ntrk = t[n_trk_br].array(library="np", entry_start=start, entry_stop=stop)
    return np.asarray(lead, np.float64), np.asarray(ntrk)


def truth_particles(F):
    """Szuka prawdy cząstek (drzewo albo gałęzie TruthParticles*). Zwraca (info, funkcja(start, stop) → (n_ch, lead_gen)) albo (info, None)."""
    info = {"trees": [k.split(";")[0] for k in F.keys()], "truth_particle_branches": []}
    cands = []
    for tn in info["trees"]:
        try:
            obj = F[tn]
            ks = list(obj.keys()) if hasattr(obj, "keys") else []
        except Exception:  # noqa: BLE001
            continue
        tp = [k for k in ks if k.startswith("TruthParticles") or k.startswith("TruthParticle_") or k.startswith("TruthCharged")]
        if tp:
            cands.append((tn, tp)); info["truth_particle_branches"] += [f"{tn}/{k}" for k in tp]
    if cands:
        tn, tp = cands[0]
        pick = lambda *subs: next((k for k in tp if all(s.lower() in k.lower() for s in subs)), None)  # noqa: E731
        bpt, beta, bq, bst = pick("_pt"), pick("_eta"), pick("charge"), pick("status")
        info.update(tree=tn, pt=bpt, eta=beta, charge=bq, status=bst)
        if not (bpt and beta and bq):
            info["note"] = "gałęzie TruthParticles są, ale brak pt/eta/charge — n_ch nie liczone"
            return info, None

        def f(start, stop):
            a = F[tn].arrays([b for b in (bpt, beta, bq, bst) if b], entry_start=start, entry_stop=stop)
            sel = (a[bq] != 0) & (a[bpt] > 500.0) & (abs(a[beta]) < 2.5)
            if bst:
                sel = sel & (a[bst] == 1)
            n = ak.to_numpy(ak.sum(sel, axis=1)).astype(np.int32)
            lg = ak.to_numpy(ak.fill_none(ak.max(a[bpt][sel], axis=1), 0.0)) / 1000.0
            return n, np.asarray(lg, np.float64)
        return info, f
    return info, None


def run_mc(path, out, max_events, step, tag):
    F = uproot.open(path); t = F["CharmAnalysis"]; keys = set(t.keys())
    n_all = t.num_entries if max_events is None else min(t.num_entries, max_events)
    if F["Tracks"].num_entries != t.num_entries:
        sys.exit(f"Tracks ({F['Tracks'].num_entries}) ≠ CharmAnalysis ({t.num_entries}) — drzewa niezsynchronizowane")
    D = {b: res(keys, f"DMesons_{b}") for b in B_MC}
    DI = {b: res(keys, f"DMesons_daughterInfo__{b}") for b in ("pdgId", "truthPdgId", "trackId")}
    NTRK, EVN = res(keys, "EventInfo_nTracksPt500Eta25"), res(keys, "EventInfo_eventNumber")
    tp_info, tp_fun = truth_particles(F)
    X = extra_branches(keys)
    # effmap (prawdziwe D) — te same gałęzie co Q_lxplus/truth/effmap.py + evnum
    ke = {"pdg": "TruthDmesons_pdgId", "pt": "TruthDmesons_pt", "eta": "TruthDmesons_eta", "phi": "TruthDmesons_phi",
          "dt": "TruthDmesons_decayType", "fromB": "TruthDmesons_fromBottom", "idx": "TruthDmesons_recoMatchIndices"}
    ke = {n: effmap_gn.resolve(keys, v) for n, v in ke.items()}
    ke.update(rdt=effmap_gn.resolve(keys, "DMesons_decayType"), lxy=effmap_gn.resolve(keys, "DMesons_fitOutput__Lxy"),
              ntrk=effmap_gn.resolve(keys, "EventInfo_nTracksPt500Eta25"), evnum=effmap_gn.resolve(keys, "EventInfo_eventNumber"))
    E, C, T = [], [], []
    off = 0
    branches = sorted(set(list(D.values()) + list(DI.values()) + [NTRK, EVN] + list(ke.values())))
    for a in t.iterate(branches, step_size=step, entry_stop=n_all, library="ak"):
        n = len(a); lead, ntrk_np = lead_and_ntrk(F, t, NTRK, off, off + n)
        g = lambda b: a[D[b]]  # noqa: E731
        # ---- mb_truth.py (etap1), bez zmian logiki
        dt = g("decayType"); tc = g("TruthMatchCategory")
        ok = ak.all(a[DI["pdgId"]] == a[DI["truthPdgId"]], axis=2) & (np.sign(g("pdgId")) == np.sign(g("truthPdgId")))
        corr = ((tc == 0) | (tc == 1)) & ok
        veto = (dt == 4) & (np.minimum(g("m") - g("mKpi1"), g("m") - g("mKpi2")) < 160)
        d0win = (g("m") > D0_WIN[0]) & (g("m") < D0_WIN[1])           # v2.0: okno D0 szersze → przywrócenie okna nominału
        sel = corr & (((dt == 10) & d0win) | ((dt == 4) & ~veto))
        ntrk = ak.to_numpy(a[NTRK])
        d0i = g("D0Index")
        isds = (dt == 1) & (d0i >= 0)
        d0i_safe = ak.where(isds, d0i, 0)
        d0corr = corr[d0i_safe]; d0m = g("m")[d0i_safe]
        tds = isds & d0corr & (abs(d0m - 1864.84) < 55) & (abs(g("DeltaMass") - 145.45) < 1.5)
        ev = {"entry": off + np.arange(n), "evnum": ak.to_numpy(a[EVN]).astype(np.uint64), "n_trk": ntrk.astype(np.int32),
              "lead_pt": lead, "n_true_dstar": ak.to_numpy(ak.sum(tds, axis=1)).astype(np.int16),
              "n_corr_d0": ak.to_numpy(ak.sum(sel & (dt == 10), axis=1)).astype(np.int16),
              "n_corr_dp": ak.to_numpy(ak.sum(sel & (dt == 4), axis=1)).astype(np.int16)}
        if not np.array_equal(ntrk, ntrk_np):
            sys.exit("n_trk z iterate ≠ n_trk z odczytu zakresu — przesunięcie wpisów")
        ev.update(read_extra(t, X, off, off + n))
        if tp_fun is not None:
            ev["nch"], ev["lead_gen"] = tp_fun(off, off + n)
        E.append(ev)
        evi = ak.to_numpy(ak.flatten(ak.broadcast_arrays(off + np.arange(n), dt)[0][sel]))
        trk = ak.sort(a[DI["trackId"]][sel], axis=2)
        trk = ak.fill_none(ak.pad_none(trk, 3, axis=2, clip=True), -1)
        C.append({"entry": evi, "decay_type": ak.to_numpy(ak.flatten(dt[sel])), "m": ak.to_numpy(ak.flatten(g("m")[sel])),
                  "pt": ak.to_numpy(ak.flatten(g("pt")[sel])) / 1000, "truth_pdg": ak.to_numpy(ak.flatten(g("truthPdgId")[sel])),
                  **{f"trk{i}": ak.to_numpy(ak.flatten(trk[:, :, i])) for i in range(3)}})
        # ---- prawdziwe D (effmap_gn = effmap.py + evnum)
        r, _ = effmap_gn.process(a, ke, off)
        T.append(r)
        off += n
    cat = lambda L, c: np.concatenate([x[c] for x in L]) if L else np.zeros(0)  # noqa: E731
    z = {f"ev_{c}": cat(E, c) for c in E[0]} if E else {}
    z.update({f"c_{c}": cat(C, c) for c in C[0]} if C else {})
    z.update({f"t_{c}": cat(T, c) for c in T[0]} if T else {})
    meta = {"kind": "mc", "source": path, "tag": tag, "n_events": int(off), "n_entries_file": int(t.num_entries), "fmt": fmt_of(keys),
            "d0_window_mev": list(D0_WIN), "extra_columns": sorted(X),
            "has_truth_particles": int(tp_fun is not None), "truth_particles": tp_info, "cutflow": hist_values(F)}
    np.savez_compressed(out, meta=json.dumps(meta), **z)
    print(f"OK {out}: {off} zdarzeń, D* prawdziwe (etap1) {int(z['ev_n_true_dstar'].sum())}, poprawne D0 {int(z['ev_n_corr_d0'].sum())}, "
          f"D+ {int(z['ev_n_corr_dp'].sum())}, prawdziwe D (effmap) {len(z.get('t_pdg', []))}, TruthParticles: {meta['has_truth_particles']}",
          file=sys.stderr)


def hist_values(F):
    out = {}
    for k in ("DSelection_Vertices", "metadata"):
        try:
            h = F[k]
            v = h.values(flow=False); lab = [str(x) for x in (h.axis().labels() or [])]
            out[k] = {"values": [float(x) for x in v], "labels": lab}
        except Exception as e:  # noqa: BLE001
            out[k] = {"error": repr(e)}
    return out


def run_data(path, out, max_events, step, tag):
    F = uproot.open(path); t = F["CharmAnalysis"]; keys = set(t.keys())
    NTRK, EVN = res(keys, "EventInfo_nTracksPt500Eta25"), res(keys, "EventInfo_eventNumber")
    n_all = t.num_entries if max_events is None else min(t.num_entries, max_events)
    X = extra_branches(keys)
    E = []
    for start in range(0, n_all, step):
        stop = min(n_all, start + step)
        lead, ntrk = lead_and_ntrk(F, t, NTRK, start, stop)
        E.append({"evnum": t[EVN].array(library="np", entry_start=start, entry_stop=stop).astype(np.uint64),
                  "n_trk": ntrk.astype(np.int32), "lead_pt": lead, **read_extra(t, X, start, stop)})
    z = {f"ev_{c}": np.concatenate([x[c] for x in E]) for c in E[0]}
    meta = {"kind": "data", "source": path, "tag": tag, "n_events": int(n_all), "n_entries_file": int(t.num_entries), "fmt": fmt_of(keys),
            "extra_columns": sorted(X)}
    np.savez_compressed(out, meta=json.dumps(meta), **z)
    print(f"OK {out}: {n_all} zdarzeń danych, mediana n_trk {np.median(z['ev_n_trk']):.0f}, lead_pt p50 {np.median(z['ev_lead_pt']):.2f}",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--kind", choices=["mc", "data"], required=True)
    ap.add_argument("--max-events", type=int, default=None)
    ap.add_argument("--step", type=int, default=50000)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    (run_mc if a.kind == "mc" else run_data)(a.inp, a.out, a.max_events, a.step, a.tag)


if __name__ == "__main__":
    main()
