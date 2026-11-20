#!/usr/bin/env python3
"""B_trig — tabela zdarzeń z ROOT v2.0 (CharmAnalysis, drzewo `CharmAnalysis`) → npz: run, event, lb, tp_/ps_<chain>,
nSelBL, lead_pt [GeV], actual_mu, grl, npu, ntrk500, + liczniki D* w oknie m(Kπ) (diagnostyka uzysku).

Ścieżka A (nie czeka na h5_v20): tabelę łączy się z istniejącymi parami EV (v1.1/v1.2, run 428770) po (run, event).
Dla 427929 potrzebne są i tak pary z score H → h5_v20 + btrig_extract (ścieżka B).
Tylko zdarzenia z ≥ 1 kandydatem (tak jak drzewo). MC: run = 0 → klucz (file_idx, entry) — zapisywany zawsze.
Użycie: trg_root_extract.py --out trg_<próbka>.npz [--max-files N] PLIK.root [...]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from btrig_common import EVENT_EXTRA, TRIG_CHAINS, find_col, find_trig_cols  # noqa: E402

ROOT_EVENT = {"run": ["EventInfo_runNumber"], "event": ["EventInfo_eventNumber"]}
ROOT_EXTRA = {k: ["EventInfo_" + x for x in v if not x.startswith("EventInfo_")] + v for k, v in EVENT_EXTRA.items()}


def resolve(keys):
    m = {}
    for k, c in {**ROOT_EVENT, **ROOT_EXTRA}.items():
        b = find_col(keys, c)
        if b:
            m[k] = b
    tr = find_trig_cols(keys, TRIG_CHAINS)
    for c, v in tr.items():
        if v["pass"]:
            m[f"tp_{c}"] = v["pass"]
        if v["ps"]:
            m[f"ps_{c}"] = v["ps"]
    return m


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--tree", default="CharmAnalysis")
    ap.add_argument("--max-files", type=int, default=0)
    ap.add_argument("--no-dstar", action="store_true", help="bez liczników D* (szybciej)")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    import awkward as ak
    import uproot
    files = sorted(a.files)[: a.max_files or None]
    parts, missing, t0 = [], set(), time.time()
    for fi, path in enumerate(files):
        t = uproot.open(path)[a.tree]
        keys = t.keys()
        m = resolve(keys)
        missing |= {c for c in TRIG_CHAINS if f"tp_{c}" not in m}
        want = list(m.values())
        dsb = {}
        if not a.no_dstar:
            for k in ("decayType", "pt", "m", "DeltaMass", "D0Index"):
                b = find_col(keys, [f"DMesons_{k}"])
                if b:
                    dsb[k] = b
            want += list(dsb.values())
        arr = t.arrays(sorted(set(want)), library="ak")
        n = len(arr)
        P = {"file_idx": np.full(n, fi, np.int32), "entry": np.arange(n, dtype=np.int64)}
        for k, b in m.items():
            v = ak.to_numpy(arr[b])
            if k.startswith("tp_") or k == "grl":
                v = v.astype(np.int8)
            elif k in ("run", "event", "lb", "nSelBL", "npu", "ntrk500"):
                v = v.astype(np.int64)
            else:
                v = v.astype(np.float32)
            P[k] = v
        if "lead_pt" in P:
            P["lead_pt"] = (P["lead_pt"] / 1000.0).astype(np.float32)   # MeV → GeV (ROOT v2.0)
        if len(dsb) == 5:
            dt, pt, dm = arr[dsb["decayType"]], arr[dsb["pt"]] / 1e3, arr[dsb["DeltaMass"]] / 1e3
            md0 = (arr[dsb["m"]] - arr[dsb["DeltaMass"]]) / 1e3
            ok = (dt == 1) & (arr[dsb["D0Index"]] >= 0) & (abs(md0 - 1.86484) < 0.055)
            for lo, hi in ((3, 5), (5, 8), (8, 30)):
                b = ok & (pt >= lo) & (pt < hi)
                P[f"nds_sr_{lo}_{hi}"] = ak.to_numpy(ak.sum(b & (dm > 0.1440) & (dm < 0.1470), axis=1)).astype(np.int16)
                P[f"nds_sb_{lo}_{hi}"] = ak.to_numpy(ak.sum(b & (dm > 0.1500) & (dm < 0.1650), axis=1)).astype(np.int16)
        parts.append(P)
        print(f"[{fi + 1}/{len(files)}] {os.path.basename(path)} {n} zd. {time.time() - t0:.0f}s", flush=True)
    keys = sorted(set.intersection(*[set(p) for p in parts]))
    out = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    out["meta"] = np.array(json.dumps({"files": [os.path.basename(f) for f in files], "n_events": int(len(out["file_idx"])),
                                       "missing_chains": sorted(missing), "columns": keys,
                                       "note": "tp_=trigPassed (po prescale), ps_=trigPS (HLT×L1), lead_pt w GeV"}))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"zapisano {a.out}: {len(out['file_idx'])} zdarzeń, brak łańcuchów: {sorted(missing)}")


if __name__ == "__main__":
    sys.exit(main())
