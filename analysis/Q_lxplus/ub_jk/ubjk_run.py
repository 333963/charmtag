#!/usr/bin/env python3
"""UBJK — jackknife μ po plikach danych na ODSŁONIĘTEJ próbce par D–D (12.11).

Importuje funkcje z `charmtag` (snapshot pf_pairs/src — BEZ modyfikacji): syst.files (file_table, assign_files,
per_file_hists, scaled_inputs, unblinded_jackknife), pairs.fit2d.fit_channel, pairs.model, pairs.frozen.sbsb_counts.
Dokłada to, czego `charmtag.syst.files.unblinded_jackknife` nie robi: wspólne μ i Z (min. sumy profili, ciągłe — nie siatka 0,25),
podział na węzły, tryb `rederive` (nogi/SB×SB danych bez pliku k), jackknife ważony (Busing 1999), korelacje kanałów.

  prep     nogi (legs_data.parquet z odsłonięcia albo przebudowa `cmd_pairs` w scratch) → histogramy/dataleg/SB×SB per plik
           → perfile.{npz,json}; kontrole: Σ plików == pairs_full.npz (bitowo), Σ dataleg/sbsb == pair_inputs, fit nominalny == fit2d_unblind.json
  fit      repliki (nominał + bez pliku k) z porcji --chunk/--nchunk, tryby fixed|rederive → fit_<chunk>.json
  ref      wprost `charmtag.syst.files.unblinded_jackknife` (odniesienie dla trybu fixed) → ref_unblinded_jackknife.json
  finish   scalenie → ubjk.json + ubjk.md
  selftest lokalnie: sztuczny podział pairs_full.npz (wielomianowo) na N „plików” → prep/fit/finish; oczekiwane σ_jk/σ_Hesse ≈ 1

Tryby wejść repliki k (f = 1 − w_k, w_k = udział zdarzeń pliku):
  fixed    = files.scaled_inputs(I, f) — dokładnie jak unblinded_jackknife (dataleg/SB×SB z pełnej próbki)
  rederive = jak fixed + dataleg i SB×SB bez pliku k, nev_legs·f (wszystkie wejścia z danych bez pliku k; N_D* i MC nadal skalowane)
Wymaga CHARM_UNBLIND=1 (decyzja o odsłonięciu; flaga UNBLIND_OK sprawdzana przez ubjk_job.sh).
"""
from __future__ import annotations

import argparse
import copy
import glob
import json
import logging
import math
import os
import sys
import time

import numpy as np

log = logging.getLogger("ubjk")
TOL_MU = 2e-3          # zgodność fitu nominalnego z fit2d_unblind.json (05.11)


# ----------------------------------------------------------------------------- narzędzia
def _setup_log(out=None):
    log.setLevel(logging.INFO)
    if not log.handlers:
        fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
        h = logging.StreamHandler(sys.stdout)
        h.setFormatter(fmt)
        log.addHandler(h)
        if out:
            os.makedirs(out, exist_ok=True)
            fh = logging.FileHandler(os.path.join(out, "ubjk.log"))
            fh.setFormatter(fmt)
            log.addHandler(fh)


def _f(x):
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def _js(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    def dflt(o):
        if hasattr(o, "tolist"):
            return o.tolist()
        if hasattr(o, "item"):
            return o.item()
        return str(o)
    with open(path + ".tmp", "w") as f:
        json.dump(obj, f, indent=1, default=dflt, allow_nan=False)
    os.replace(path + ".tmp", path)


def _clean(o):
    """NaN/inf → None (json allow_nan=False)."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return _f(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _need_unblind():
    if os.environ.get("CHARM_UNBLIND") != "1":
        raise SystemExit("UBJK: wymaga CHARM_UNBLIND=1 (odsłonięty jackknife — decyzja o odsłonięciu)")


def _charmtag():
    import charmtag
    from charmtag.pairs import CHANNELS
    from charmtag.pairs.model import Model, load_inputs
    from charmtag.pairs.fit2d import fit_channel
    from charmtag.syst import files as F
    return charmtag, CHANNELS, Model, load_inputs, fit_channel, F


def sbsb_keys(I):
    return sorted(I.sbsb, key=lambda k: (k[0], k[1], k[2]))


# ----------------------------------------------------------------------------- per plik
def perfile_from_legs(L, I, nf):
    """L z kolumną file_idx → {H_<ch>: (nf,…)}, {DL_<dt>_<c>: (nf, nb)}, SBSB (nf, nkeys)."""
    from charmtag.pairs import CHANNELS
    from charmtag.pairs.model import M0, SBW
    from charmtag.pairs.frozen import sbsb_counts
    from charmtag.syst import files as F
    Hf = F.per_file_hists(L, I, nf, "full")
    sb = L[np.abs(L.m - L.decay_type.map(M0)) > L.decay_type.map(SBW)]
    keys = sbsb_keys(I)
    DL = {k: np.zeros((nf, len(I.edges[k[0]]) - 1)) for k in I.dataleg}
    SB = np.zeros((nf, len(keys)))
    for k in range(nf):
        s = sb[sb.file_idx == k]
        for (dt, c) in I.dataleg:
            DL[(dt, c)][k] = np.histogram(s.m[(s.decay_type == dt) & (s.cat == c)], I.edges[dt])[0]
        cnt = sbsb_counts(s, I.ncat)
        SB[k] = [cnt.get(key, 0) for key in keys]
    return {ch: Hf[ch] for ch in CHANNELS}, DL, SB, keys


def save_perfile(out, H, DL, SB, meta):
    arr = {f"H_{ch}": v for ch, v in H.items()}
    arr.update({f"DL_{dt}_{c}": v for (dt, c), v in DL.items()})
    arr["SBSB"] = SB
    np.savez_compressed(os.path.join(out, "perfile.npz"), **arr)
    _js(os.path.join(out, "perfile.json"), _clean(meta))


def load_perfile(d):
    js = json.load(open(os.path.join(d, "perfile.json")))
    A = dict(np.load(os.path.join(d, "perfile.npz")))
    H = {k[2:]: v for k, v in A.items() if k.startswith("H_")}
    DL = {(int(k.split("_")[1]), int(k.split("_")[2])): v for k, v in A.items() if k.startswith("DL_")}
    keys = [tuple([s.split("_")[0], int(s.split("_")[1]), int(s.split("_")[2])]) for s in js["sbsb_keys"]]
    return js, H, DL, A["SBSB"], keys


def inputs_for(I, mode, k, w, DL=None, SB=None, keys=None):
    from charmtag.syst import files as F
    if k < 0:
        return I
    f = 1.0 - float(w[k])
    I2 = F.scaled_inputs(I, f)                    # = unblinded_jackknife (fixed)
    if mode == "fixed":
        return I2
    if mode != "rederive":
        raise ValueError(mode)
    I2.dataleg = {key: np.asarray(I.dataleg[key], float) - DL[key][k] for key in I.dataleg}
    I2.sbsb = dict(I.sbsb)
    for j, key in enumerate(keys):
        I2.sbsb[key] = int(round(I.sbsb.get(key, 0) - SB[k, j]))
    I2.nev_legs = I.nev_legs * f                  # gęstość fałszywych nóg z próbki bez pliku k (scale = ndata/nev_legs bez zmian)
    return I2


# ----------------------------------------------------------------------------- fit repliki
def _minimize(g, lo, hi, xatol=1e-3):
    try:
        from scipy.optimize import minimize_scalar
        r = minimize_scalar(g, bounds=(lo, hi), method="bounded", options={"xatol": xatol, "maxiter": 60})
        return float(r.x), float(r.fun)
    except ImportError:                            # złoty podział (bez scipy)
        gr = (math.sqrt(5) - 1) / 2
        a, b = lo, hi
        c, d = b - gr * (b - a), a + gr * (b - a)
        gc, gd = g(c), g(d)
        while b - a > xatol:
            if gc < gd:
                b, d, gd = d, c, gc
                c = b - gr * (b - a)
                gc = g(c)
            else:
                a, c, gc = c, d, gd
                d = a + gr * (b - a)
                gd = g(d)
        x = 0.5 * (a + b)
        return x, g(x)


def _root(h, a, b, tol=1e-3):
    """h(a) < 0 < h(b) lub odwrotnie → pierwiastek (bisekcja); brak zmiany znaku → None."""
    ha, hb = h(a), h(b)
    if not (np.isfinite(ha) and np.isfinite(hb)) or ha * hb > 0:
        return None
    for _ in range(40):
        m = 0.5 * (a + b)
        hm = h(m)
        if ha * hm <= 0:
            b, hb = m, hm
        else:
            a, ha = m, hm
        if abs(b - a) < tol:
            break
    return 0.5 * (a + b)


def replica_fit(M, H, channels, combined=True, interval=False):
    from charmtag.pairs.fit2d import fit_channel
    out, fits = {}, {}
    for ch in channels:
        n = H[ch]
        full = np.ones(n.shape[1:], bool)
        r, cf, m1 = fit_channel(M, ch, n, full, scan=False)
        out[ch] = {k: _f(r[k]) for k in ("mu_hat", "mu_err_hesse", "Z", "dev_ndf")}
        out[ch]["valid"] = bool(r["valid"])
        fits[ch] = (cf, m1, n, full)
    mus = np.array([out[c]["mu_hat"] for c in channels], float)
    ers = np.array([out[c]["mu_err_hesse"] for c in channels], float)
    wm = float(np.sum(mus / ers ** 2) / np.sum(1 / ers ** 2))
    out["wmean"] = {"mu": wm, "err": float(1 / math.sqrt(np.sum(1 / ers ** 2)))}
    if not combined:
        return out
    cache, nbad = {}, [0]

    def g(mu):
        mu = round(float(mu), 6)
        if mu not in cache:
            s = 0.0
            for ch in channels:
                cf, m1, n, full = fits[ch]
                ms = cf.fit(n, full, fix_mu=mu, start=np.array(m1.values))
                nbad[0] += int(not ms.valid)
                s += ms.fval - m1.fval
            cache[mu] = s
        return cache[mu]
    lo, hi = max(-1.5, float(mus.min()) - 0.3), min(5.0, float(mus.max()) + 0.3)
    muc, gmin = _minimize(g, lo, hi)
    g0 = g(0.0)
    q0 = 2 * (g0 - gmin) if muc > 0 else 0.0
    c = {"mu": muc, "q0": q0, "Z": math.sqrt(max(q0, 0.0)), "gmin": gmin, "n_eval": len(cache), "n_invalid_fits": nbad[0],
         "at_bound": bool(abs(muc - lo) < 2e-3 or abs(muc - hi) < 2e-3)}
    if interval:
        h = lambda m: g(m) - gmin - 0.5
        s0 = out["wmean"]["err"]                       # skala z kanałów; przedział szukany blisko minimum (fity daleko od μ̂ są wolne)

        def side(sgn):
            step = 1.5 * s0
            for _ in range(4):
                x = min(5.0, max(-1.5, muc + sgn * step))
                if h(x) > 0:
                    return _root(h, *sorted((muc, x)), tol=5e-3)
                if x in (5.0, -1.5):
                    return None
                step *= 2
            return None
        up, dn = side(+1), side(-1)
        c["mu_68"] = [dn, up]
        c["sigma_68"] = 0.5 * (up - dn) if (up is not None and dn is not None) else None
    out["combined"] = c
    return out


# ----------------------------------------------------------------------------- statystyka jackknife
def jk_stats(theta0, th, w):
    """Standardowy (równe grupy) i ważony (Busing i in. 1999, delete-m_j) jackknife."""
    th = np.asarray(th, float)
    ok = np.isfinite(th)
    th, w = th[ok], np.asarray(w, float)[ok]
    n = len(th)
    if n < 2 or theta0 is None:
        return None
    sd = float(np.sqrt((n - 1) / n * np.sum((th - th.mean()) ** 2)))
    h = 1.0 / (w / w.sum())                                      # h_j = N/m_j
    tJ = float(n * theta0 - np.sum((1 - 1 / h) * th))
    pseudo = h * theta0 - (h - 1) * th
    var_w = float(np.mean((pseudo - tJ) ** 2 / (h - 1)))
    return {"n": n, "sigma_jk": sd, "sigma_jk_weighted": math.sqrt(var_w), "mu_jk_biascorr": tJ,
            "mean_drop": float(th.mean()), "rel_precision_sigma_jk": 1 / math.sqrt(2 * (n - 1)),
            "max_abs_shift": float(np.max(np.abs(th - theta0))), "max_shift_idx": int(np.argmax(np.abs(th - theta0)))}


# ----------------------------------------------------------------------------- polecenia
def cmd_prep(a):
    _need_unblind()
    import pandas as pd
    charmtag, CHANNELS, Model, load_inputs, fit_channel, F = _charmtag()
    from charmtag import config as C
    t0 = time.time()
    cfg = C.load(a.cfg, None, require_paths=True)
    I = load_inputs(os.path.join(a.pairs_dir, "pair_inputs"))
    H0 = dict(np.load(os.path.join(a.pairs_dir, "pairs_full.npz")))
    meta = {"charmtag_file": charmtag.__file__, "cfg": os.path.abspath(a.cfg), "pairs_dir": a.pairs_dir, "synthetic": False}
    # --- nogi
    if a.legs and os.path.exists(a.legs):
        legs_src = a.legs
        L = pd.read_parquet(a.legs)
        meta["legs_how"] = "legs_data.parquet z odsłonięcia (05.11)"
    else:
        if not a.rebuild_dir:
            raise SystemExit(f"STOP: brak nóg ({a.legs}) i brak --rebuild-dir")
        from charmtag.pairs.build import cmd_pairs
        cfg2 = copy.deepcopy(cfg)
        cfg2["paths"]["scores_dir"] = cfg["paths"].get("scores_dir") or os.path.join(cfg["paths"]["out_dir"], "scores")
        cfg2["paths"]["out_dir"] = a.rebuild_dir
        log.info("prep: brak legs_data.parquet → przebudowa par (cmd_pairs --unblind, norm=%s) w %s, scores_dir=%s",
                 a.norm, a.rebuild_dir, cfg2["paths"]["scores_dir"])
        cmd_pairs(cfg2, unblind=True, norm_mode=a.norm)
        legs_src = os.path.join(a.rebuild_dir, "pairs", "legs_data.parquet")
        if not os.path.exists(legs_src):
            raise SystemExit("STOP: przebudowa nie zapisała legs_data.parquet (pyarrow?)")
        L = pd.read_parquet(legs_src)
        meta["legs_how"] = "przebudowa cmd_pairs w scratch (src bieżący)"
        try:
            j2 = json.load(open(os.path.join(a.rebuild_dir, "pairs", "pair_inputs.json")))
            meta["rebuild_n_dstar"] = j2.get("meta", {}).get("dstar_fit_data", {}).get("N")
        except Exception:
            pass
    meta["legs_src"] = legs_src
    meta["n_legs"] = int(len(L))
    meta["groups"] = a.groups
    if a.groups == "files":
        # --- tabela plików (jak syst files: plik ← cand + liczby wierszy h5)
        ft, how = None, None
        try:
            ft = F.file_table(cfg)
            how = "h5 (charmtag.syst.files.file_table)"
        except Exception as e:
            log.warning("prep: file_table z h5 nie działa (%s)", e)
        fj = None
        if a.files_json and os.path.exists(a.files_json):
            fj = pd.DataFrame(json.load(open(a.files_json))["files"])
        if ft is None:
            if fj is None:
                raise SystemExit("STOP: brak tabeli plików (h5 i files.json)")
            ft, how = fj, f"files.json ({a.files_json})"
        meta["file_table_how"] = how
        if fj is not None and ft is not fj:
            same = len(fj) == len(ft) and (fj.rows.to_numpy() == ft.rows.to_numpy()).all()
            meta["file_table_matches_blind_files_json"] = bool(same)
            if not same:
                log.warning("prep: tabela plików z h5 ≠ files.json ślepego biegu (liczba/wiersze)")
        nf = len(ft)
        if nf < 2:
            raise SystemExit("STOP: < 2 pliki danych")
        L = F.assign_files(L, ft)
        ev = ft.events.to_numpy(float) if ft.events.notna().all() else ft.rows.to_numpy(float)
        w = ev / ev.sum()
    elif a.groups.startswith("events:"):
        # --- grupy zdarzeń: deterministyczny hash (run, event) mod G — zdarzenie nigdy nie jest dzielone; w = 1/G (hash jednorodny)
        nf = int(a.groups.split(":")[1])
        if nf < 2:
            raise SystemExit("STOP: events:G wymaga G ≥ 2")
        run = pd.to_numeric(L.source_file, errors="coerce").fillna(0).to_numpy(np.uint64)
        ent = L.entry.to_numpy().astype(np.uint64)
        with np.errstate(over="ignore"):
            hsh = (ent * np.uint64(0x9E3779B97F4A7C15) + run * np.uint64(0xBF58476D1CE4E5B9)) >> np.uint64(33)
        L = L.assign(file_idx=(hsh % np.uint64(nf)).astype(np.int64))
        ft = pd.DataFrame([{"file": f"evgrp_{g}", "job": f"g{g:03d}", "rows": 0, "events": None} for g in range(nf)])
        w = np.full(nf, 1.0 / nf)
        meta["file_table_how"] = f"grupy zdarzeń: hash(run, event) mod {nf}"
    else:
        raise SystemExit(f"--groups: files | events:G (jest {a.groups})")
    ue = L[["file_idx", "source_file", "entry"]].drop_duplicates().file_idx.value_counts()
    wl = np.array([ue.get(k, 0) for k in range(nf)], float)                 # zdarzenia z nogami per plik (diagnostyka wag)
    meta.update({"n_files": nf, "files": ft.to_dict("records"), "w": w.tolist(), "w_events_with_legs": (wl / wl.sum()).tolist(),
                 "w_corr": float(np.corrcoef(w, wl)[0, 1])})
    # --- per plik + kontrole
    H, DL, SB, keys = perfile_from_legs(L, I, nf)
    chk = {}
    chk["hists_sum_equal_pairs_full"] = {ch: bool(np.array_equal(H[ch].sum(0), H0[ch])) for ch in CHANNELS}
    chk["dataleg_sum_equal"] = bool(all(np.allclose(DL[k].sum(0), I.dataleg[k]) for k in I.dataleg))
    chk["sbsb_sum_equal"] = bool(all(int(SB[:, j].sum()) == int(I.sbsb.get(key, 0)) for j, key in enumerate(keys)))
    meta["checks"] = chk
    meta["sbsb_keys"] = [f"{k[0]}_{k[1]}_{k[2]}" for k in keys]
    meta["rederive_ok"] = bool(chk["dataleg_sum_equal"] and chk["sbsb_sum_equal"])
    meta["pairs_per_file"] = {ch: H[ch].sum(axis=(1, 2, 3)).tolist() for ch in CHANNELS}
    if not all(chk["hists_sum_equal_pairs_full"].values()):
        save_perfile(a.out, H, DL, SB, meta)
        raise SystemExit(f"STOP: Σ histogramów per plik ≠ pairs_full.npz {chk['hists_sum_equal_pairs_full']} — inne nogi/inna lista plików")
    # --- fit nominalny == 05.11
    ref = json.load(open(os.path.join(a.pairs_dir, "fit2d_unblind.json")))
    M = Model(I)
    rep = {}
    for ch in CHANNELS:
        r, _, _ = fit_channel(M, ch, H0[ch], np.ones(H0[ch].shape[1:], bool), scan=False)
        rep[ch] = {"mu": r["mu_hat"], "ref": ref["channels"][ch]["mu_hat"], "err": r["mu_err_hesse"], "ref_err": ref["channels"][ch]["mu_err_hesse"]}
    meta["nominal_reproduced"] = rep
    meta["nominal_ok"] = bool(all(abs(v["mu"] - v["ref"]) < TOL_MU for v in rep.values()))
    save_perfile(a.out, H, DL, SB, meta)
    log.info("prep: %d plików, nogi %d (%s), Σ==pairs_full %s, rederive_ok %s, nominał %s (%.0f s)", nf, len(L), meta["legs_how"],
             chk["hists_sum_equal_pairs_full"], meta["rederive_ok"], {c: round(v["mu"], 4) for c, v in rep.items()}, time.time() - t0)
    if not meta["nominal_ok"]:
        raise SystemExit(f"STOP: fit nominalny nie odtwarza fit2d_unblind.json (tol {TOL_MU}): {rep}")
    return 0


def run_replicas(pairs_dir, perfile_dir, chunk, nchunk, modes, combined=True, channels=None, I=None):
    _, CHANNELS, Model, load_inputs, _, _ = _charmtag()
    chs = channels or list(CHANNELS)
    I = I if I is not None else load_inputs(os.path.join(pairs_dir, "pair_inputs"))
    js, Hf, DL, SB, keys = load_perfile(perfile_dir)
    if "rederive" in modes and not js.get("rederive_ok"):
        log.warning("fit: rederive_ok=False → tylko tryb fixed")
        modes = [m for m in modes if m != "rederive"]
    nf = js["n_files"]
    w = np.array(js["w"])
    Ht = {ch: v.sum(0) for ch, v in Hf.items()}
    reps = list(range(-1, nf))[chunk::nchunk]
    out = {"chunk": chunk, "nchunk": nchunk, "modes": modes, "replicas": {}}
    for k in reps:
        H = {ch: (Ht[ch] - Hf[ch][k] if k >= 0 else Ht[ch]) for ch in Ht}
        rk = {}
        for mode in (["nominal"] if k < 0 else modes):
            t = time.time()
            Ik = inputs_for(I, mode, k, w, DL, SB, keys)
            rk[mode] = replica_fit(Model(Ik), H, chs, combined=combined, interval=(k < 0 and combined))
            rk[mode]["sec"] = time.time() - t
            cmb = rk[mode].get("combined", {})
            log.info("fit: replika %s %-8s μ %s wspólne %s Z %s (%s ewaluacji profilu; %.0f s)", "nom" if k < 0 else js["files"][k].get("job", k),
                     mode, {c: round(rk[mode][c]["mu_hat"], 4) for c in chs}, round(cmb["mu"], 4) if cmb else "-",
                     round(cmb["Z"], 3) if cmb else "-", cmb.get("n_eval", 0), rk[mode]["sec"])
        out["replicas"][str(k)] = rk
    return out


def cmd_fit(a):
    _need_unblind()
    modes = [m for m in a.modes.split(",") if m]
    r = run_replicas(a.pairs_dir, a.perfile, a.chunk, a.nchunk, modes, combined=not a.no_combined)
    _js(os.path.join(a.out, f"fit_{a.chunk}.json"), _clean(r))
    return 0


def cmd_ref(a):
    _need_unblind()
    import pandas as pd
    _, CHANNELS, Model, load_inputs, _, F = _charmtag()
    js = json.load(open(os.path.join(a.perfile, "perfile.json")))
    if js.get("groups", "files") != "files":
        raise SystemExit("ref: tylko dla --groups files (unblinded_jackknife przypisuje pliki z tabeli h5)")
    I = load_inputs(os.path.join(a.pairs_dir, "pair_inputs"))
    ft = pd.DataFrame(js["files"])
    L = F.assign_files(pd.read_parquet(js["legs_src"]), ft)
    t = time.time()
    rep = F.unblinded_jackknife({}, I, L, len(ft), np.array(js["w"]), chs=list(CHANNELS), ft=ft)
    _js(os.path.join(a.out, "ref_unblinded_jackknife.json"), _clean({"unblinded": rep, "sec": time.time() - t}))
    log.info("ref: unblinded_jackknife %s (%.0f s)", {c: (round(v["sigma_jk"], 4), round(v["ratio"], 3)) for c, v in rep.items()}, time.time() - t)
    return 0


def summarize(fits, js, ref_fit=None, ref_jk=None):
    """fits: lista wyników run_replicas → raport."""
    R = {}
    for f in fits:
        R.update(f["replicas"])
    nf = js["n_files"]
    miss = [k for k in range(-1, nf) if str(k) not in R]
    if miss:
        raise SystemExit(f"finish: brak replik {miss}")
    nom = R["-1"]["nominal"]
    chs = [c for c in nom if c not in ("wmean", "combined", "sec")]
    w = np.array(js["w"])
    jobs = [str(f.get("job", i)) for i, f in enumerate(js["files"])]
    modes = sorted({m for k in range(nf) for m in R[str(k)]})
    rep = {"n_files": nf, "modes": modes, "nominal": nom, "per_mode": {}, "files": jobs, "w": w.tolist()}
    if ref_fit:
        rep["ref_fit2d_unblind"] = {ch: {"mu": ref_fit["channels"][ch]["mu_hat"], "err": ref_fit["channels"][ch]["mu_err_hesse"],
                                         "mu_68": ref_fit["channels"][ch]["mu_profile_68"]} for ch in chs}
        rep["ref_fit2d_unblind"]["combined"] = ref_fit.get("combined")
    for mode in modes:
        pm = {}
        drops = {}
        for ch in chs + ["combined", "wmean"]:
            key = "mu_hat" if ch in chs else "mu"
            th = [R[str(k)][mode][ch][key] if ch in R[str(k)][mode] else None for k in range(nf)]
            th = np.array([np.nan if v is None else v for v in th], float)
            if ch not in nom:
                continue
            t0 = nom[ch][key]
            s = jk_stats(t0, th, w)
            if s is None:
                continue
            if ch in chs:
                sig = nom[ch]["mu_err_hesse"]
                s["sigma_hesse"] = sig
                s["ratio_hesse"] = s["sigma_jk"] / sig
                s["ratio_hesse_weighted"] = s["sigma_jk_weighted"] / sig
                if ref_fit:
                    lo, hi = ref_fit["channels"][ch]["mu_profile_68"]
                    s["sigma_profile68"] = 0.5 * (hi - lo)
                    s["ratio_profile68"] = s["sigma_jk"] / s["sigma_profile68"]
                z = np.array([R[str(k)][mode][ch]["Z"] for k in range(nf)], float)
                s["Z_nominal"] = nom[ch]["Z"]
                s["Z_drop_min"], s["Z_drop_max"] = float(np.nanmin(z)), float(np.nanmax(z))
                s["n_invalid"] = int(sum(not R[str(k)][mode][ch]["valid"] for k in range(nf)))
            elif ch == "combined":
                sig = nom["combined"].get("sigma_68")
                s["sigma_68_nominal"] = sig
                s["ratio_68"] = s["sigma_jk"] / sig if sig else None
                s["ratio_68_weighted"] = s["sigma_jk_weighted"] / sig if sig else None
                z = np.array([R[str(k)][mode]["combined"]["Z"] for k in range(nf)], float)
                s["Z_nominal"] = nom["combined"]["Z"]
                s["Z_drop_min"], s["Z_drop_max"] = float(np.nanmin(z)), float(np.nanmax(z))
                s["n_at_bound"] = int(sum(R[str(k)][mode]["combined"].get("at_bound", False) for k in range(nf)))
            else:
                sig = nom["wmean"]["err"]
                s["sigma_nominal"] = sig
                s["ratio"] = s["sigma_jk"] / sig
            s["max_shift_file"] = jobs[s.pop("max_shift_idx")]
            s["mu_nominal"] = t0
            pm[ch] = s
            drops[ch] = th
        cc = [c for c in chs if c in drops]
        D = np.array([drops[c] - np.nanmean(drops[c]) for c in cc])
        if D.shape[0] > 1 and np.isfinite(D).all():
            cov = (nf - 1) / nf * D @ D.T
            sd = np.sqrt(np.diag(cov))
            pm["corr_channels"] = {"channels": cc, "corr": (cov / np.outer(sd, sd)).tolist()}
        pm["mu_drop"] = {c: drops[c].tolist() for c in drops}
        rep["per_mode"][mode] = pm
    if ref_jk:
        u = ref_jk.get("unblinded", {})
        rep["ref_unblinded_jackknife"] = {c: {k: u[c][k] for k in ("mu", "sigma_fit", "sigma_jk", "ratio")} for c in u}
        if "fixed" in rep["per_mode"]:
            rep["ref_vs_fixed_max_dsigma"] = max(abs(u[c]["sigma_jk"] - rep["per_mode"]["fixed"][c]["sigma_jk"]) for c in u if c in rep["per_mode"]["fixed"])
    rep["headline_mode"] = "rederive" if "rederive" in rep["per_mode"] else "fixed"
    return rep


def to_md(rep):
    L = ["# UBJK — jackknife μ po plikach danych (ODSŁONIĘTE, pełna lumi run 428770)", ""]
    pr = rep.get("prep", {})
    L.append(f"- grupy: {pr.get('groups', 'files')} ({pr.get('file_table_how')}); nogi: {pr.get('legs_how')}; syntetyczne: {pr.get('synthetic')}")
    L.append(f"- grup (plików): {rep['n_files']}; tryby: {', '.join(rep['modes'])}; główny: **{rep['headline_mode']}**")
    pm0 = rep["per_mode"].get(rep["headline_mode"], {})
    any_ch = next((v for v in pm0.values() if isinstance(v, dict) and "rel_precision_sigma_jk" in v), None)
    if any_ch:
        L.append(f"- precyzja σ_jk ≈ ±{100 * any_ch['rel_precision_sigma_jk']:.0f} % (1/√(2(n−1)), gauss.)")
    if "ref_vs_fixed_max_dsigma" in rep:
        L.append(f"- kontrola: `charmtag.syst.files.unblinded_jackknife` vs tryb fixed: max |Δσ_jk| = {rep['ref_vs_fixed_max_dsigma']:.2e}")
    for mode, pm in rep["per_mode"].items():
        L += ["", f"## tryb `{mode}`", "",
              "| | μ nom. | σ (Hesse/68 %) | σ_jk | σ_jk ważony | σ_jk/σ | σ_jk,w/σ | μ_J (bias-corr.) | Z nom. | Z bez pliku (min–max) | max |Δμ| (plik) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        for ch, s in pm.items():
            if not isinstance(s, dict) or "sigma_jk" not in s:
                continue
            sig = s.get("sigma_hesse") or s.get("sigma_68_nominal") or s.get("sigma_nominal")
            rat = s.get("ratio_hesse") or s.get("ratio_68") or s.get("ratio")
            ratw = s.get("ratio_hesse_weighted") or s.get("ratio_68_weighted")
            fm = lambda x, p=3: "–" if x is None else f"{x:.{p}f}"
            zr = f"{s['Z_drop_min']:.2f}–{s['Z_drop_max']:.2f}" if "Z_drop_min" in s else "–"
            L.append(f"| {ch} | {fm(s['mu_nominal'])} | {fm(sig)} | {fm(s['sigma_jk'])} | {fm(s['sigma_jk_weighted'])} | {fm(rat, 2)} | "
                     f"{fm(ratw, 2)} | {fm(s['mu_jk_biascorr'])} | {fm(s.get('Z_nominal'), 2)} | {zr} | {fm(s['max_abs_shift'])} ({s['max_shift_file']}) |")
        if "corr_channels" in pm:
            cc = pm["corr_channels"]
            L += ["", "Korelacje kanałów z jackknife (μ bez pliku k): " + ", ".join(
                f"{cc['channels'][i]}–{cc['channels'][j]} {cc['corr'][i][j]:+.2f}" for i in range(len(cc["channels"])) for j in range(i + 1, len(cc["channels"])))]
    nc = rep["nominal"].get("combined", {})
    if nc:
        L += ["", f"Wspólne μ (ciągłe minimum sumy profili, replika nominalna): {nc['mu']:.3f}, 68 % {nc.get('mu_68')}, Z = {nc['Z']:.2f}"
              + (f"; 05.11 (siatka 0,25): {rep['ref_fit2d_unblind']['combined']['mu_common']:.2f}, Z = {rep['ref_fit2d_unblind']['combined']['Z']:.2f}"
                 if rep.get("ref_fit2d_unblind", {}).get("combined") else "")]
    L += ["", "σ_jk/σ > 1 ⇒ korelacje par w zdarzeniach (lub niejednorodność plików) poszerzają błąd ponad Poissona w komórkach;",
          f"jackknife ma rozdzielczość ±1/√(2(n−1)) ≈ ±{100 / math.sqrt(2 * max(rep['n_files'] - 1, 1)):.0f} % na σ (n = {rep['n_files']}) — różnice poniżej ~2× tego nieistotne.",
          "Tryb fixed = definicja `unblinded_jackknife` (dataleg/SB×SB z pełnej próbki); rederive = wszystkie wejścia z danych bez pliku k.",
          "N_D* (normalizacja) i wejścia MC — skalowane ∝ zdarzeniom (fluktuacja N_D* w systematyce normalizacji, nie tutaj)."]
    return "\n".join(L) + "\n"


def cmd_finish(a):
    _need_unblind()
    js = json.load(open(os.path.join(a.indir, "perfile.json")))
    fits = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(a.indir, "fit_*.json")))]
    if not fits:
        raise SystemExit("finish: brak fit_*.json")
    nch = {f["nchunk"] for f in fits}
    if len(nch) != 1 or len(fits) != nch.pop():
        raise SystemExit(f"finish: niekompletne porcje ({len(fits)} plików fit_*, nchunk {[f['nchunk'] for f in fits]})")
    rf = os.path.join(a.pairs_dir, "fit2d_unblind.json")
    ref_fit = json.load(open(rf)) if os.path.exists(rf) else None
    rj = os.path.join(a.indir, "ref_unblinded_jackknife.json")
    ref_jk = json.load(open(rj)) if os.path.exists(rj) else None
    rep = summarize(fits, js, ref_fit, ref_jk)
    rep["prep"] = {k: js.get(k) for k in ("groups", "legs_how", "file_table_how", "checks", "rederive_ok", "nominal_reproduced", "w_corr", "synthetic")}
    _js(os.path.join(a.out, "ubjk.json"), _clean(rep))
    md = to_md(rep)
    with open(os.path.join(a.out, "ubjk.md"), "w") as f:
        f.write(md)
    print(md)
    return 0


def cmd_selftest(a):
    """Lokalnie: pairs_full.npz + pair_inputs → wielomianowy podział na N „plików” (H, dataleg, SB×SB) → fit/finish."""
    _need_unblind()
    _, CHANNELS, Model, load_inputs, _, _ = _charmtag()
    rng = np.random.default_rng(a.seed)
    I = load_inputs(os.path.join(a.pairs_dir, "pair_inputs"))
    H0 = dict(np.load(os.path.join(a.pairs_dir, "pairs_full.npz")))
    nf = a.nfiles
    w = rng.uniform(0.5, 1.5, nf)
    w /= w.sum()
    split = lambda x: np.moveaxis(rng.multinomial(np.asarray(x).astype(np.int64), w), -1, 0).astype(float)
    H = {ch: split(H0[ch]) for ch in CHANNELS}
    DL = {k: split(np.round(v)) for k, v in I.dataleg.items()}
    keys = sbsb_keys(I)
    SB = split(np.array([I.sbsb[k] for k in keys]))
    os.makedirs(a.out, exist_ok=True)
    meta = {"synthetic": True, "seed": a.seed, "n_files": nf, "w": w.tolist(),
            "files": [{"file": f"syn_{k}", "job": f"syn{k:02d}", "rows": 0, "events": None} for k in range(nf)],
            "sbsb_keys": [f"{k[0]}_{k[1]}_{k[2]}" for k in keys], "rederive_ok": True,
            "checks": {"hists_sum_equal_pairs_full": {ch: bool(np.array_equal(H[ch].sum(0), H0[ch])) for ch in CHANNELS},
                       "dataleg_sum_equal": bool(all(np.allclose(DL[k].sum(0), I.dataleg[k]) for k in I.dataleg)),
                       "sbsb_sum_equal": bool(all(int(SB[:, j].sum()) == I.sbsb[k] for j, k in enumerate(keys)))}}
    save_perfile(a.out, H, DL, SB, meta)
    log.info("selftest: %d plików, kontrole %s", nf, meta["checks"])
    # tożsamość: rederive z „pustym” plikiem = nominał
    Iz = inputs_for(I, "rederive", 0, np.zeros(nf), {k: np.zeros_like(v) for k, v in DL.items()}, np.zeros_like(SB), keys)
    from charmtag.pairs.fit2d import fit_channel
    ident = {}
    for ch in CHANNELS:
        full = np.ones(H0[ch].shape[1:], bool)
        r1, _, _ = fit_channel(Model(I), ch, H0[ch], full, scan=False)
        r2, _, _ = fit_channel(Model(Iz), ch, H0[ch], full, scan=False)
        ident[ch] = abs(r1["mu_hat"] - r2["mu_hat"])
    log.info("selftest: |Δμ| rederive(pusty plik) vs nominał: %s", ident)
    modes = [m for m in a.modes.split(",") if m]
    r = run_replicas(a.pairs_dir, a.out, 0, 1, modes, combined=not a.no_combined)
    r["identity_dmu"] = ident
    _js(os.path.join(a.out, "fit_0.json"), _clean(r))
    a.indir = a.out
    return cmd_finish(a)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep")
    p.add_argument("--cfg", required=True)
    p.add_argument("--pairs-dir", required=True)
    p.add_argument("--legs", default=None)
    p.add_argument("--rebuild-dir", default=None)
    p.add_argument("--norm", default="direct")
    p.add_argument("--files-json", default=None)
    p.add_argument("--groups", default="files", help="files (jackknife po plikach) | events:G (G grup zdarzeń, hash)")
    p.add_argument("--out", required=True)
    p = sp.add_parser("fit")
    p.add_argument("--pairs-dir", required=True)
    p.add_argument("--perfile", required=True)
    p.add_argument("--chunk", type=int, required=True)
    p.add_argument("--nchunk", type=int, required=True)
    p.add_argument("--modes", default="fixed,rederive")
    p.add_argument("--no-combined", action="store_true")
    p.add_argument("--out", required=True)
    p = sp.add_parser("ref")
    p.add_argument("--pairs-dir", required=True)
    p.add_argument("--perfile", required=True)
    p.add_argument("--out", required=True)
    p = sp.add_parser("finish")
    p.add_argument("--pairs-dir", required=True)
    p.add_argument("--in", dest="indir", required=True)
    p.add_argument("--out", required=True)
    p = sp.add_parser("selftest")
    p.add_argument("--pairs-dir", required=True)
    p.add_argument("--nfiles", type=int, default=10)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--modes", default="fixed,rederive")
    p.add_argument("--no-combined", action="store_true")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    _setup_log(getattr(a, "out", None))
    return {"prep": cmd_prep, "fit": cmd_fit, "ref": cmd_ref, "finish": cmd_finish, "selftest": cmd_selftest}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
