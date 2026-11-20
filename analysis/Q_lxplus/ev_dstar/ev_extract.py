#!/usr/bin/env python3
"""EV: zdarzeniowa próbka par D0←D* (tag Δm) z h5 v1.1 (lxplus, EOS) albo z podpróbek npz (Mac/test).

Wejście: lista plików h5 (compound `events`) lub npz (klucze = kolumny h5). Dla każdego pliku:
  D0 (decayType 10): cechy H (`charmtag.features.compute`), score OOF modelem nominalnym D0_H (fold = fold_of(run,event),
  jak w treningu), m_swap (hipoteza K↔π), swapy MC (`mark_swaps`, truth+100);
  D* (decayType 1): m, Δm, pT, η, truth, pT(π_s);
  dopasowanie D0↔D* jak `pipeline.dstar_control`: ten sam (run,event), |m_D0 − (m_D* − Δm)| < match_tol_m,
  0.8·pT(D*) < pT(D0) ≤ 1.001·pT(D*), |Δη| < 0.1.
Selekcja par: Δm < DM_MAX (0.165 GeV), |m_D0 − m0| < M_WIN (0.075 GeV). Wyjście: npz (float32) z kolumnami
  run,event,evblk (indeks zdarzenia w pliku; klucz zdarzenia = (file_idx, evblk)),file_idx,row,row_ds,m,pt,eta,phi,dm,score,fold,truth,truth_ds,pt_ds,eta_ds,pis_pt,weight,m_swap, 18 cech H
+ 'all_*' (score,pt,m,truth,weight wszystkich D0 w oknie 75 MeV — do progów WP i ε_B), + meta (progi WP, cechy, pliki).
Użycie: ev_extract.py --models DIR --out OUT.npz --tag data [--chunk 500000] FILE [FILE ...]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import yaml

DM_MAX = 0.165
M_WIN = 0.075


def _load_npz(path):
    d = np.load(path)
    ks = [k for k in d.files if not k.startswith("Dmeson_ptcone")
          and k not in ("actual_mu", "average_mu", "pdfInfo_pdgid1", "pdfInfo_pdgid2", "Dmeson_maxdR", "Dmeson_species")]
    df = pd.DataFrame({k: d[k] for k in ks})
    if "_row" in df:
        df = df.rename(columns={"_row": "_local_row"})
    else:
        df["_local_row"] = np.arange(len(df))
    df["_evblk"] = event_blocks({c: df[c].to_numpy() for c in EV_COLS if c in df})   # podpróbka npz: wiersze w kolejności h5
    return [df]


EV_COLS = ("run_number", "event_number", "event_weight", "pdfInfo_pdgid1", "pdfInfo_pdgid2", "actual_mu")


def event_blocks(cols):
    """Indeks zdarzenia w pliku = numer ciągłego bloku wierszy o stałych kolumnach zdarzeniowych (EV_COLS obecne w `cols`).
    19.11: w fd_nd run = 0 i event_number powtarza się (także w jednym pliku, scalone joby generatora) → (run, event) skleja
    różne zdarzenia (AUDYT_evjoin_2026-11-18). Wiersze jednego zdarzenia leżą w h5 ciągiem, więc zmiana któregokolwiek pola = nowe
    zdarzenie. Zwraca int64 [n]."""
    keys = [np.asarray(cols[c]) for c in EV_COLS if c in cols]
    if not keys:
        raise ValueError("event_blocks: brak kolumn zdarzeniowych")
    n = len(keys[0])
    if n == 0:
        return np.zeros(0, np.int64)
    ch = np.zeros(n, bool); ch[0] = True
    for k in keys:
        ch[1:] |= k[1:] != k[:-1]
    return np.cumsum(ch).astype(np.int64) - 1


def chunk_bounds(blk, chunk):
    """Granice porcji [a, b) wyrównane do granic zdarzeń (porcja nie tnie zdarzenia → D0 i D* tego samego zdarzenia w jednej porcji)."""
    n = len(blk); out = []; a = 0
    while a < n:
        b = min(n, a + chunk)
        while b < n and blk[b] == blk[b - 1]:
            b += 1
        out.append((a, b)); a = b
    return out


def _load_h5_chunks(path, cols, chunk, dt_col, keep_dt):
    import h5py
    with h5py.File(path, "r") as f:
        ds = f["events"]
        names = set(ds.dtype.names)
        cols = [c for c in cols if c in names]
        evc = [c for c in EV_COLS if c in names]
        e = ds.fields(evc)[:]
        blk = event_blocks({c: e[c] for c in evc})
        for a, b in chunk_bounds(blk, chunk):
            x = ds.fields(cols)[a:b]
            df = pd.DataFrame({c: x[c] for c in cols})
            df["_local_row"] = np.arange(a, b)
            df["_evblk"] = blk[a:b]
            df = df[np.isin(df[dt_col].to_numpy(), keep_dt)]
            yield df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True, help="katalog modelu D0_H (config_used.yaml, meta.json, fold*.json)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", required=True, help="data|fd|mb")
    ap.add_argument("--chunk", type=int, default=500_000)
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    from charmtag import features as F
    from charmtag.io import resolve_columns
    from charmtag.model import Ensemble
    from charmtag.samples import mark_swaps
    from charmtag.kinematics import swapped_mass_sp

    cfg = yaml.safe_load(open(os.path.join(args.models, "config_used.yaml")))
    meta = json.load(open(os.path.join(args.models, "meta.json")))
    sp = cfg["species"]["D0"]
    dc = cfg["dstar_control"]
    m0 = float(sp["m0"])
    units = cfg.get("units")
    ens = Ensemble(args.models)
    feats = list(ens.feats)
    dt_col = None
    out_parts, all_parts = [], []
    t0 = time.time()
    for fi, path in enumerate(args.files):
        is_npz = path.endswith(".npz")
        if is_npz:
            chunks = _load_npz(path)
            cols_all = list(chunks[0].columns)
        else:
            import h5py
            with h5py.File(path, "r") as f:
                cols_all = list(f["events"].dtype.names)
        cm2 = resolve_columns(cols_all, cfg["columns"], 2)
        cm3 = resolve_columns(cols_all, cfg["columns"], 3)
        dt_col = cm3.get("decay_type")
        need = sorted({c for c in (F.raw_columns(cm3, 3) + [dt_col]) if c})
        if not is_npz:
            chunks = _load_h5_chunks(path, need, args.chunk, dt_col, [int(sp["decay_type"]), int(dc["dstar_decay_type"])])
        n_d0 = n_ds = n_pairs = 0
        for raw in chunks:
            r0 = raw[raw[dt_col] == int(sp["decay_type"])].reset_index(drop=True)
            r1 = raw[raw[dt_col] == int(dc["dstar_decay_type"])].reset_index(drop=True)
            if len(r0) == 0 or len(r1) == 0:
                continue
            L = F.to_logical(r0, cm2, 2)
            L["weight"] = L["weight"].fillna(1.0)
            L = F.compute(L, sp, units)
            L["row"] = r0["_local_row"].to_numpy()
            L["file_idx"] = fi
            L["evblk"] = r0["_evblk"].to_numpy()
            mark_swaps(L, sp, cfg)
            L["m_swap"] = swapped_mass_sp(L, sp)
            L["score"], L["fold"] = ens.score(L)
            win = np.abs(L["m"].to_numpy() - m0) < M_WIN
            all_parts.append(L.loc[win, ["score", "pt", "m", "truth", "weight"]].to_numpy(dtype=np.float32))
            n_d0 += len(L)
            ds = F.to_logical(r1, cm3, 3)
            ds = ds[(ds["dm"] > 0) & (ds["dm"] < DM_MAX)]
            n_ds += len(ds)
            ds = ds.assign(m0_from=ds["m"] - ds["dm"], row_ds=r1["_local_row"].to_numpy()[ds.index], truth_ds=ds["truth"],
                           pt_ds=ds["pt"], eta_ds=ds["eta"], pis_pt=ds["d_pt2"] * float((units or {}).get("d_pt", 1.0)),
                           evblk=r1["_evblk"].to_numpy()[ds.index])[
                ["run", "event", "evblk", "m0_from", "pt_ds", "eta_ds", "dm", "row_ds", "truth_ds", "pis_pt"]]
            Lk = L.loc[win, ["run", "event", "evblk", "m", "pt", "eta"]].assign(_i=np.nonzero(win)[0])
            j = Lk.merge(ds, on=["run", "event", "evblk"])            # 19.11: + evblk (fd_nd: (run, event) niejednoznaczne)
            ok = (np.abs(j["m"] - j["m0_from"]) < dc["match_tol_m"]) & (j["pt"] <= j["pt_ds"] * 1.001) & \
                 (j["pt"] > 0.8 * j["pt_ds"]) & (np.abs(j["eta"] - j["eta_ds"]) < 0.1)
            j = j[ok].reset_index(drop=True)
            if len(j) == 0:
                continue
            base = L.iloc[j["_i"].to_numpy()].reset_index(drop=True)
            keep = ["run", "event", "evblk", "file_idx", "row", "m", "pt", "eta", "phi", "score", "fold", "truth", "weight", "m_swap"] + feats
            part = pd.concat([base[keep], j[["row_ds", "dm", "truth_ds", "pt_ds", "eta_ds", "pis_pt"]]], axis=1)
            out_parts.append(part)
            n_pairs += len(part)
        print(f"[{args.tag}] {os.path.basename(path)}: D0 {n_d0}  D* {n_ds}  pary {n_pairs}  {time.time()-t0:.0f}s", flush=True)
    if out_parts:
        P = pd.concat(out_parts, ignore_index=True)
    else:
        P = pd.DataFrame(columns=["run", "event", "evblk", "file_idx", "row", "m", "pt", "eta", "phi", "score", "fold", "truth", "weight",
                                  "m_swap", "row_ds", "dm", "truth_ds", "pt_ds", "eta_ds", "pis_pt"] + feats)
    A = np.concatenate(all_parts) if all_parts else np.zeros((0, 5), np.float32)
    out = {c: P[c].to_numpy(dtype=np.int64 if c in ("run", "event", "evblk", "file_idx", "row", "row_ds", "fold") else np.float32) for c in P.columns}
    out.update({"all_score": A[:, 0], "all_pt": A[:, 1], "all_m": A[:, 2], "all_truth": A[:, 3], "all_w": A[:, 4]})
    out["meta"] = np.array(json.dumps({"tag": args.tag, "files": [os.path.basename(p) for p in args.files], "features": feats,
                                       "thresholds": meta.get("thresholds"), "m0": m0, "dm_max": DM_MAX, "m_win": M_WIN,
                                       "n_pairs": int(len(P)), "n_all_d0_win": int(len(A))}))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez_compressed(args.out, **out)
    print(f"[{args.tag}] zapisano {args.out}: pary {len(P)}, D0 w oknie {len(A)}, {time.time()-t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
