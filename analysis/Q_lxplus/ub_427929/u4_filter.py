#!/usr/bin/env python3
"""UB427 filtr (18.11.2026, projekt A): h5 v2.0 data_427929 *_none.h5 → h5 o selekcji nominału (v1.1/v1.2) i klasie zdarzeń MB427929.

Cięcia (README Q_lxplus/h5_v20: „po Lxy ≥ 0,25 i D0 1,75–2,00 v2.0 ≡ v1.2”):
  D0 (decayType 10): Dmeson_Lxy ≥ 0,25 mm i 1,75 < m < 2,00 GeV;  D* (1): Dmeson_Lxy ≥ 0,25 mm (D0-w-D* okno bez zmian v1.1→v2.0);  D⁺ (4): bez zmian.
  Zdarzenie: passGRL_MinBiasAFP == 1 i trigPassed_HLT_mb_mbts_L1MBTS_1 == 1 (klasa MB427929, B_trig/btrig_common.classify).
Wyjście: ten sam compound `events` (podzbiór wierszy, kolejność zachowana) + `evtable` (podzbiór zdarzeń klasy; mianownik), atrybuty + u4_filter.
  u4_filter.py IN.h5 [IN.h5 ...] --outdir DIR [--stats stats.json] [--no-class]"""
import argparse
import json
import os
import sys

import h5py
import numpy as np

LXY_MIN = 0.25
D0_WIN = (1.75, 2.00)
GRL, MB = "passGRL_MinBiasAFP", "trigPassed_HLT_mb_mbts_L1MBTS_1"


def masks(e, use_class=True):
    dt = e["Dmeson_decayType"].astype(np.int32); lxy = e["Dmeson_Lxy"].astype(np.float64); m = e["Dmeson_mass"].astype(np.float64)
    cand = np.ones(len(e), bool)
    cand[(dt == 10) & ~((lxy >= LXY_MIN) & (m > D0_WIN[0]) & (m < D0_WIN[1]))] = False
    cand[(dt == 1) & ~(lxy >= LXY_MIN)] = False
    evok = np.ones(len(e), bool)
    if use_class:
        evok = (e[GRL] == 1) & (e[MB] == 1)
    return cand, evok, dt


def one(path, outdir, use_class=True):
    with h5py.File(path, "r") as f:
        e = f["events"][:]
        ev = f["evtable"][:] if "evtable" in f else None
        attrs = dict(f.attrs)
    if use_class and (GRL not in e.dtype.names or MB not in e.dtype.names):
        sys.exit(f"{path}: brak kolumn klasy ({GRL}, {MB})")
    cand, evok, dt = masks(e, use_class)
    keep = cand & evok
    st = {"file": os.path.basename(path), "rows_in": int(len(e)), "rows_out": int(keep.sum())}
    for t, nm in ((10, "D0"), (4, "Dplus"), (1, "Dstar")):
        st[f"{nm}_in"] = int((dt == t).sum()); st[f"{nm}_lxy_win_ok"] = int(((dt == t) & cand).sum()); st[f"{nm}_out"] = int(((dt == t) & keep).sum())
    if ev is not None:
        eok = ((ev[GRL] == 1) & (ev[MB] == 1)) if use_class else np.ones(len(ev), bool)
        st.update({"evtable_in": int(len(ev)), "evtable_class": int(eok.sum()), "evtable_grl": int((ev[GRL] == 1).sum())})
        ev = ev[eok]
    out = os.path.join(outdir, os.path.basename(path).replace("_none.h5", "_mb.h5"))
    tmp = out + ".part"
    with h5py.File(tmp, "w") as g:
        g.create_dataset("events", data=e[keep], chunks=True if keep.sum() else None, compression="gzip", compression_opts=4)
        if ev is not None:
            g.create_dataset("evtable", data=ev, compression="gzip", compression_opts=4)
        for k, v in attrs.items():
            try:
                g.attrs[k] = v
            except TypeError:
                g.attrs[k] = str(v)
        g.attrs["u4_filter"] = json.dumps({"lxy_min_d0_dstar": LXY_MIN, "d0_window_gev": D0_WIN, "event_class": None if not use_class else f"{GRL}==1 & {MB}==1",
                                           "source": path})
    os.replace(tmp, out)
    return out, st


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+"); ap.add_argument("--outdir", required=True); ap.add_argument("--stats"); ap.add_argument("--no-class", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    S = []
    for p in a.files:
        out, st = one(p, a.outdir, not a.no_class)
        S.append(st)
        print(f"{os.path.basename(p)} → {os.path.basename(out)}: wiersze {st['rows_in']} → {st['rows_out']}; D0 {st['D0_in']}→{st['D0_out']}, "
              f"D* {st['Dstar_in']}→{st['Dstar_out']}, D+ {st['Dplus_in']}→{st['Dplus_out']}", flush=True)
    if a.stats:
        json.dump(S, open(a.stats, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
