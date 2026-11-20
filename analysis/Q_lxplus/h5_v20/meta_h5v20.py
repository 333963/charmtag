#!/usr/bin/env python3
"""meta_h5v20.py jobs.txt PARTS METADIR — per próbka: pliki, części, wiersze per gatunek, n_events_read, evtable, rozmiar, schemat → METADIR/<próbka>.json, summary.json."""
import sys, os, json, time, collections, h5py, numpy as np
jobs, parts, md = sys.argv[1:4]
os.makedirs(md, exist_ok=True)
S = collections.OrderedDict()
for line in open(jobs):
    k, f, t, p = line.split()
    s = S.setdefault(k, dict(sample=k, n_root=0, parts={}, missing=[]))
    s["n_root"] += 1
    for x in p.replace("+", ",").split(","):
        path = os.path.join(parts, k, f"{k}_{t}_{x}.h5")
        d = s["parts"].setdefault(x, dict(n=0, rows=0, n_events_read=0, n_evtable=0, bytes=0, rows_by_decayType={}, columns=None, schema=None))
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            s["missing"].append(os.path.basename(path)); continue
        try:
            with h5py.File(path, "r") as h:
                e = h["events"]; d["n"] += 1; d["rows"] += int(e.shape[0]); d["bytes"] += os.path.getsize(path)
                d["n_events_read"] += int(h.attrs.get("n_events_read", 0)); d["n_evtable"] += int(h.attrs.get("n_evtable", 0))
                d["columns"] = d["columns"] or len(e.dtype.names); d["schema"] = d["schema"] or str(h.attrs.get("schema"))
                u, c = np.unique(e["Dmeson_decayType"][:], return_counts=True)
                for a, b in zip(u, c): d["rows_by_decayType"][str(int(a))] = d["rows_by_decayType"].get(str(int(a)), 0) + int(b)
        except Exception as ex:
            s["missing"].append(f"{os.path.basename(path)} (nieczytelny: {ex})")
for k, s in S.items():
    s["complete"] = not s["missing"]; s["time"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    json.dump(s, open(os.path.join(md, f"{k}.json"), "w"), indent=1, ensure_ascii=False)
    print(k, "kompletna" if s["complete"] else f"BRAK {len(s['missing'])}", {x: (d["n"], d["rows"], round(d["bytes"] / 1e9, 2)) for x, d in s["parts"].items()})
summ = {k: dict(complete=s["complete"], n_root=s["n_root"], **{f"{x}_rows": d["rows"] for x, d in s["parts"].items()},
                **{f"{x}_GB": round(d["bytes"] / 1e9, 3) for x, d in s["parts"].items()}, n_events_read=s["parts"].get("none", {}).get("n_events_read"))
        for k, s in S.items()}
json.dump(dict(time=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), samples=summ), open(os.path.join(md, "summary.json"), "w"), indent=1, ensure_ascii=False)
