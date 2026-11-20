#!/usr/bin/env python3
"""effmap_gn.py (gen_norm, 16.11.2026) = kopia Q_lxplus/truth/effmap.py + kolumna `evnum` (EventInfo_eventNumber) do joinu z nogami h5.

effmap.py — wiersze prawdy D z istniejącego MC v1.1 (pełna symulacja ATLAS) do mapy wydajności. Etap Q, 2026-11-02.
Na 1 plik ROOT: wszystkie TruthDmesons |pdg| ∈ {421, 411, 413}, pT > 2.5 GeV, |η| < 2.7, z flagą `matched`:
  matched = ≥1 kandydat reco (TruthDmesons_recoMatchIndices) o decayType kanału (D0 10, D+ 4, D* 1) i Lxy ≥ 0
            (= preselekcja charmpp Donly; bez cięć klasyfikatora — ε klasyfikatora (WP) bierze się z charmtag).
  mode    = truth decayType == kanał analizy (D0→Kπ 10, D+→Kππ 4, D*→D0(Kπ)π 1).
  ntrk    = EventInfo_nTracksPt500Eta25 (reco) — oś krotności mapy.
Wyjście: npz (entry, ntrk, pdg, pt, eta, phi, fromB, mode, matched) + n_events.
  python3 effmap.py in.root out.npz [--max-events N]"""
import argparse, sys
import awkward as ak
import numpy as np
import uproot

SPEC = {421: 10, 411: 4, 413: 1}   # |pdg| → decayType kanału


def resolve(keys, name):
    for n in (name, name + "_NOSYS"):
        if n in keys:
            return n
    raise KeyError(name)


def process(a, k, off):
    pid = a[k["pdg"]]; n_ev = len(pid)
    pt = a[k["pt"]] / 1000.0; eta = a[k["eta"]]
    sel = ((abs(pid) == 421) | (abs(pid) == 411) | (abs(pid) == 413)) & (pt > 2.5) & (abs(eta) < 2.7)
    # --- dopasowanie reco: płaskie indeksy
    idx = a[k["idx"]]                                     # (ev, truth, k)
    rdt = a[k["rdt"]]; rlxy = a[k["lxy"]]                 # (ev, reco)
    nreco = ak.to_numpy(ak.num(rdt))
    roff = np.concatenate([[0], np.cumsum(nreco)])[:-1]
    ntr = ak.to_numpy(ak.num(pid))
    nk = ak.to_numpy(ak.flatten(ak.num(idx, axis=2)))     # per truth
    flat_idx = ak.to_numpy(ak.flatten(idx, axis=None)).astype(np.int64)
    ev_of_truth = np.repeat(np.arange(n_ev), ntr)
    truth_of_k = np.repeat(np.arange(ntr.sum()), nk)
    ev_of_k = ev_of_truth[truth_of_k]
    valid = (flat_idx >= 0) & (flat_idx < nreco[ev_of_k])
    g = np.where(valid, flat_idx + roff[ev_of_k], 0)
    frdt = ak.to_numpy(ak.flatten(rdt)); flxy = ak.to_numpy(ak.flatten(rlxy))
    if frdt.size == 0:
        frdt = np.zeros(1, int); flxy = np.zeros(1)
    fpid = np.abs(ak.to_numpy(ak.flatten(pid)))
    want = np.vectorize(lambda x: SPEC.get(int(x), -99))(fpid) if fpid.size else np.zeros(0, int)
    ok_k = valid & (frdt[g] == want[truth_of_k]) & (flxy[g] >= 0)
    matched = np.bincount(truth_of_k[ok_k], minlength=ntr.sum()) > 0
    fsel = ak.to_numpy(ak.flatten(sel))
    ntrk = ak.to_numpy(a[k["ntrk"]])
    evn = ak.to_numpy(a[k["evnum"]]).astype(np.uint64) if "evnum" in k else np.zeros(n_ev, np.uint64)   # gen_norm: klucz joinu z h5
    rows = dict(entry=(off + ev_of_truth)[fsel].astype(np.int64),
                evnum=evn[ev_of_truth][fsel],
                ntrk=ntrk[ev_of_truth][fsel].astype(np.int32),
                pdg=ak.to_numpy(ak.flatten(pid))[fsel].astype(np.int32),
                pt=ak.to_numpy(ak.flatten(pt))[fsel].astype(np.float32),
                eta=ak.to_numpy(ak.flatten(eta))[fsel].astype(np.float32),
                phi=ak.to_numpy(ak.flatten(a[k["phi"]]))[fsel].astype(np.float32),
                fromB=ak.to_numpy(ak.flatten(a[k["fromB"]]))[fsel].astype(np.int8),
                mode=(ak.to_numpy(ak.flatten(a[k["dt"]])) == want)[fsel].astype(np.int8),
                matched=matched[fsel].astype(np.int8))
    return rows, n_ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--max-events", type=int, default=None)
    ap.add_argument("--step", default="200 MB")
    a = ap.parse_args()
    t = uproot.open(a.inp)["CharmAnalysis"]; keys = set(t.keys())
    k = {"pdg": "TruthDmesons_pdgId", "pt": "TruthDmesons_pt", "eta": "TruthDmesons_eta", "phi": "TruthDmesons_phi",
         "dt": "TruthDmesons_decayType", "fromB": "TruthDmesons_fromBottom", "idx": "TruthDmesons_recoMatchIndices"}
    k = {n: resolve(keys, v) for n, v in k.items()}
    k.update(rdt=resolve(keys, "DMesons_decayType"), lxy=resolve(keys, "DMesons_fitOutput__Lxy"),
             ntrk=resolve(keys, "EventInfo_nTracksPt500Eta25"), evnum=resolve(keys, "EventInfo_eventNumber"))
    parts, off = [], 0
    for arr in t.iterate(list(k.values()), step_size=a.step, entry_stop=a.max_events, library="ak"):
        r, n = process(arr, k, off); parts.append(r); off += n
    rows = {c: np.concatenate([p[c] for p in parts]) if parts else np.zeros(0) for c in (parts[0] if parts else [])}
    np.savez_compressed(a.out, n_events=off, source=a.inp, **rows)
    m = rows.get("mode", np.zeros(0)) == 1
    msg = " ".join(f"{s}: {int(((abs(rows['pdg']) == s) & m & (rows['pt'] > 3) & (abs(rows['eta']) < 2.5)).sum())} w kanale, "
                   f"ε={((abs(rows['pdg']) == s) & m & (rows['pt'] > 3) & (abs(rows['eta']) < 2.5) & (rows['matched'] == 1)).sum() / max(((abs(rows['pdg']) == s) & m & (rows['pt'] > 3) & (abs(rows['eta']) < 2.5)).sum(), 1):.3f};"
                   for s in (421, 411, 413)) if off else ""
    print(f"OK {a.out}: {off} zdarzeń, {len(rows.get('pdg', []))} prawdziwych D. {msg}", file=sys.stderr)


if __name__ == "__main__":
    main()
