#!/usr/bin/env python3
"""UBPNET — zbiór zdarzeń próbki nominalnej odsłonięcia 05.11: unikalne event_number z h5 v1.1 data_428770_*_none.h5 (39 plików;
zdarzenia z ≥ 1 kandydatem dowolnego gatunku — nadzbiór zdarzeń z parami). Run jeden (428770) → klucz = event_number.
Użycie: upn_evset.py --out evset.npz FILE.h5 [...]"""
import argparse
import json
import sys

import numpy as np


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("files", nargs="+")
    a = ap.parse_args(argv)
    import h5py
    ev, info = [], []
    for p in a.files:
        with h5py.File(p, "r") as f:
            d = f["events"]
            e = np.unique(d["event_number"]).astype(np.uint64)
            runs = np.unique(d["run_number"]).tolist() if "run_number" in (d.dtype.names or ()) else None
        ev.append(e); info.append({"file": p.rsplit("/", 1)[-1], "events": int(len(e)), "runs": runs})
        print(f"{p.rsplit('/', 1)[-1]}: {len(e)} zdarzeń, runy {runs}", flush=True)
    allr = sorted({r for i in info for r in (i["runs"] or [])})
    if allr and allr != [428770]:
        raise SystemExit(f"STOP: runy {allr} ≠ [428770]")
    E = np.unique(np.concatenate(ev))
    np.savez_compressed(a.out, event=E, meta=np.array(json.dumps({"files": info, "n_files": len(info), "n_events": int(len(E))})))
    print(f"zapisano {a.out}: {len(E)} zdarzeń z {len(info)} plików")
    return 0


if __name__ == "__main__":
    sys.exit(main())
