#!/usr/bin/env python3
"""Scala części PAIR_PLANES: parts/<próbka>.p<k>.npz → <próbka>.npz (+ meta), PP_SUMMARY.json. Brak części = błąd."""
import argparse, glob, json, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--parts", required=True); ap.add_argument("--out", required=True); ap.add_argument("--samples", required=True)
a = ap.parse_args()
summ = {}
for kv in a.samples.replace(",", " ").split():
    s, n = kv.split(":"); n = int(n)
    fs = [os.path.join(a.parts, f"{s}.p{k}.npz") for k in range(n)]
    miss = [f for f in fs if not os.path.exists(f)]
    if miss:
        sys.exit(f"brak części: {miss}")
    H, meta = None, None
    for f in fs:
        z = np.load(f); m = json.load(open(f.replace(".npz", ".meta.json")))
        if H is None:
            H = {k: z[k].astype(np.float64) for k in z.files}; meta = dict(m); meta["files"] = list(m["files"]); meta["seconds"] = m["seconds"]
        else:
            for k in z.files:
                H[k] += z[k]
            meta["files"] += m["files"]; meta["n_cand_sel"] += m["n_cand_sel"]; meta["n_pairs_shared_any"] += m["n_pairs_shared_any"]
            meta["seconds"] += m["seconds"]
    meta["n_files"] = len(meta["files"]); meta["part"] = "all"; meta["files"] = sorted(meta["files"])
    np.savez_compressed(os.path.join(a.out, f"{s}.npz"), **H)
    json.dump(meta, open(os.path.join(a.out, f"{s}.meta.json"), "w"), indent=1)
    summ[s] = {"n_files": meta["n_files"], "n_cand_sel": meta["n_cand_sel"], "cpu_s": round(meta["seconds"])}
    print(s, summ[s], flush=True)
json.dump(summ, open(os.path.join(a.out, "PP_SUMMARY.json"), "w"), indent=1)
