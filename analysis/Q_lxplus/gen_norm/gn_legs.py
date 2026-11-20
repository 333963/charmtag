#!/usr/bin/env python3
"""gn_legs.py — gen_norm (16.11.2026): nogi trybu `direct` z MC minbias ND (h5 + nominalne scores), DOKŁADNIE kodem charmtag.

Selekcja = `charmtag.pairs.build.legs_for_truth` + `mc_rates` (snapshot src/, bez zmian):
  D0: truth ∈ {0,1} (po mark_swaps closest_mass z apply), m ∈ [1750, 2000) MeV, pT 3–30 (in_range z apply), bez cięcia score;
  D+: to samo, m ∈ [1700, 2200) MeV, veto_ok (weto D*);
  D*: decayType 1, 0 < Δm < 0,165 GeV (scan_source), |m − Δm − m_D0| < 0,040 GeV, truth ∈ {0,1}, Δm < fit_range[1] (= `mc_rates.dsel`).
Wyjście (--out DIR): legs_mb.parquet (D0+D+: tag, evnum, row, decay_type, m, pt, eta, truth, trk0..2, in_pair_* liczniki),
  dstar_mb.parquet (D* w oknie: tag, evnum, row, m, dm, pt, eta, phi, truth), mc_rates.json (= meta.mc_rates_unweighted nominału),
  legs_meta.json (pliki, offsety, kontrola zgodności z listą h5).
  python3 gn_legs.py --cfg BASE_CFG [--out-dir NOMINAL] --out DIR [--ref-rates 1344,9534,8665,22,41,15]"""
import argparse
import json
import logging
import os
import re
import sys

import h5py
import numpy as np
import pandas as pd

from charmtag import config as C
from charmtag.io import list_columns, n_rows, open_events, resolve_columns
from charmtag.pairs import CHANNELS
from charmtag.pairs import build as B
from charmtag.pairs.frozen import shared_track
from charmtag.pairs.model import M0

log = logging.getLogger("gen_norm")
KEY = "incl_mc"


def tag_of(path):
    m = re.search(r"_(\d{6})_none\.h5$", os.path.basename(path))
    return m.group(1) if m else os.path.basename(path)


def file_map(cfg, sdir):
    """Lista plików h5 incl_mc (jak w apply) + offsety; kontrola z attrs scores (source, file_offsets) — inna lista = STOP."""
    paths = C.sample_paths(cfg, KEY) if C.is_multi(cfg, KEY) else [C.sample_path(cfg, KEY)]
    offs = np.concatenate([[0], np.cumsum([n_rows(p, cfg["io"]["dataset"]) for p in paths])]).astype(np.int64)
    chk = {}
    for sp in ("D0", "Dplus"):
        with h5py.File(os.path.join(sdir(sp), f"{KEY}.h5"), "r") as f:
            src = str(f.attrs.get("source", "")); fo = f.attrs.get("file_offsets")
        bn = [os.path.basename(x) for x in src.split(";")]          # nazwy plików (ścieżka EOS mogła mieć inny prefiks /eos/home-k/…)
        same_src = bn == [os.path.basename(x) for x in paths]
        same_off = (fo is None) or (json.loads(fo) == offs.tolist())
        chk[sp] = {"same_source_list": bool(same_src), "same_offsets": bool(same_off), "n_files_scores": len(src.split(";"))}
        if not (same_src and same_off):
            sys.exit(f"STOP: scores {sp}/{KEY}.h5 policzone na innej liście plików h5 niż {cfg['paths']['samples'][KEY]} "
                     f"(source zgodne: {same_src}, offsety zgodne: {same_off}) — wiersze `row` nie pasują")
    return paths, offs, chk


def dstar_scan(cfg, p, paths, chunk):
    """D* jak scan_source (dsm) + okno mc_rates (dsel), z wierszem, η, φ (do joinu z prawdą). Liczba musi = mc_rates.n_dstar."""
    ds_name = cfg["io"]["dataset"]
    cm = resolve_columns(list_columns(paths[0], ds_name), cfg["columns"], 3)
    need = ["run", "event", "decay_type", "m", "dm", "truth", "pt", "eta", "phi"]
    cols = {k: cm.get(k) for k in need if cm.get(k)}
    miss = [k for k in need if k not in cols]
    if miss:
        sys.exit(f"dstar_scan: brak kolumn {miss}")
    w = p["dstar_norm"]["d0_window"]; hi = p["dstar_norm"]["fit_range"][1]
    out = []
    with open_events(paths, ds_name) as ds:
        n = ds.shape[0]
        if cfg["io"].get("max_rows"):
            n = min(n, int(cfg["io"]["max_rows"]))
        for start in range(0, n, chunk):
            stop = min(n, start + chunk)
            a = ds.fields(list(cols.values()))[start:stop]
            g = {k: np.asarray(a[v]) for k, v in cols.items()}
            dsm = (g["decay_type"] == cfg["dstar_control"]["dstar_decay_type"]) & (g["dm"] > 0) & (g["dm"] < 0.165)
            m, dm, tr = g["m"].astype(float), g["dm"].astype(float), g["truth"].astype(float)
            dsel = dsm & (np.abs(m - dm - M0[10] / 1000) < w) & np.isin(tr, [0, 1]) & (dm < hi)
            if dsel.any():
                out.append(pd.DataFrame({"row": np.flatnonzero(dsel) + start, "evnum": g["event"][dsel].astype(np.uint64),
                                         "m": m[dsel], "dm": dm[dsel], "pt": g["pt"][dsel].astype(float),
                                         "eta": g["eta"][dsel].astype(float), "phi": g["phi"][dsel].astype(float), "truth": tr[dsel]}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["row", "evnum", "m", "dm", "pt", "eta", "phi", "truth"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfg", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--out-dir", default=None, help="nadpisz paths.out_dir (nominalny bieg: scores/<sp>_<var>/incl_mc.h5)")
    ap.add_argument("--ref-rates", default="1344,9534,8665,22,41,15",
                    help="oczekiwane n_dstar, n_D0, n_D+, pary D0D0, DpD0, DpDp (PF2 pairs/pair_inputs.json meta.mc_rates_unweighted); '' = bez")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    ov = {"paths": {"out_dir": a.out_dir}} if a.out_dir else None
    cfg = C.load(a.cfg, ov)
    p = B.pcfg(cfg); E = B.edges(p); var = p["variant"]
    sdir = lambda sp: os.path.join(cfg["paths"].get("scores_dir") or os.path.join(cfg["paths"]["out_dir"], "scores"), f"{sp}_{var}")  # noqa: E731
    cols = ["row", "run", "event", "m", "pt", "eta", "score", "veto_ok", "truth"]
    os.makedirs(a.out, exist_ok=True)
    paths, offs, chk = file_map(cfg, sdir)
    legs, D, nev = B.legs_for_truth(cfg, p, KEY, E, sdir, cols)
    r = B.mc_rates(legs, D, p)
    log.info("mc_rates: %s", r)
    # pary (jak mc_rates) — liczniki per noga: ile par z daną nogą w każdym kanale
    L = pd.concat([legs[10], legs[4]], ignore_index=True)
    P = L.merge(L, on=["source_file", "entry"], suffixes=("1", "2"))
    P = P[(P.decay_type1 < P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.cand1 < P.cand2))]
    P = P[~shared_track(P)]
    fi = np.searchsorted(offs, L.row.to_numpy(np.int64), side="right") - 1
    out = pd.DataFrame({"tag": [tag_of(paths[i]) for i in fi], "evnum": L.entry.to_numpy().astype(np.uint64), "row": L.row.to_numpy(np.int64),
                        "decay_type": L.decay_type.to_numpy(np.int32), "m": L.m.to_numpy(float), "pt": L.pt.to_numpy(float),
                        "eta": L.eta.to_numpy(float), "truth": L.truth.to_numpy(float),
                        **{f"trk{i}": L[f"trk{i}"].to_numpy() for i in range(3)}})
    for ch, (x, y) in CHANNELS.items():
        q = P[(P.decay_type1 == min(x, y)) & (P.decay_type2 == max(x, y))]
        cnt = pd.concat([q.cand1, q.cand2]).value_counts()
        out[f"npair_{ch}"] = out.row.map(cnt).fillna(0).astype(np.int32).to_numpy()
    out.to_parquet(os.path.join(a.out, "legs_mb.parquet"), index=False)
    Ds = dstar_scan(cfg, p, paths, p["chunk_rows"])
    Ds.insert(0, "tag", [tag_of(paths[i]) for i in np.searchsorted(offs, Ds.row.to_numpy(np.int64), side="right") - 1])
    Ds.to_parquet(os.path.join(a.out, "dstar_mb.parquet"), index=False)
    got = [r["n_dstar"], r["n_single"][10], r["n_single"][4]] + [r["n_pairs"][ch] for ch in CHANNELS]
    ok_ds = len(Ds) == r["n_dstar"]
    ref = [int(x) for x in a.ref_rates.split(",")] if a.ref_rates else None
    meta = {"cfg": os.path.abspath(a.cfg), "out_dir": cfg["paths"]["out_dir"], "n_files": len(paths), "files": [os.path.basename(x) for x in paths],
            "file_offsets": offs.tolist(), "scores_check": chk, "n_events_with_candidates": int(nev),
            "mc_rates": {"n_dstar": r["n_dstar"], "n_single": {str(k): v for k, v in r["n_single"].items()}, "n_pairs": r["n_pairs"]},
            "dstar_scan_equals_mc_rates": bool(ok_ds), "ref_rates": ref, "got_rates": got,
            "ref_match": (got == ref) if ref else None, "pairs_track_id": p["track_id"]}
    json.dump(meta, open(os.path.join(a.out, "legs_meta.json"), "w"), indent=1, default=int)
    json.dump(meta["mc_rates"], open(os.path.join(a.out, "mc_rates.json"), "w"), indent=1, default=int)
    print(f"legs: D0 {r['n_single'][10]}, D+ {r['n_single'][4]}, D* {r['n_dstar']} (skan D* {len(Ds)}), pary {r['n_pairs']}; "
          f"zgodność z nominałem {ref}: {meta['ref_match']}", file=sys.stderr)
    if not ok_ds:
        sys.exit("STOP: własny skan D* ≠ mc_rates.n_dstar — selekcja okna D* nie odtworzona")


if __name__ == "__main__":
    main()
