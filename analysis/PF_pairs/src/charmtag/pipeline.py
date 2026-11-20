"""Komendy: train, apply, fit, validate (wywoływane z cli.py)."""
from __future__ import annotations

import json
import logging
import os

import numpy as np
import pandas as pd

from . import config as C
from . import features as F
from . import massfit as MF
from . import model as M
from . import samples as SA
from . import store

log = logging.getLogger("charmtag")
SCORE_COLS = ["row", "run", "event", "m", "pt", "eta", "dm", "truth", "score", "fold", "veto_ok"]


def mdir(cfg, species, variant):  # S5: paths.models_dir (opcjonalnie: modele z innego biegu, np. nominalnego)
    return os.path.join(cfg["paths"].get("models_dir") or os.path.join(cfg["paths"]["out_dir"], "models"), f"{species}_{variant}")


def sdir(cfg, species, variant):  # S5: paths.scores_dir (opcjonalnie: wyniki z innego biegu)
    return os.path.join(cfg["paths"].get("scores_dir") or os.path.join(cfg["paths"]["out_dir"], "scores"), f"{species}_{variant}")


def _json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


def _auc(s, y, w):
    o = np.argsort(s)
    s, y, w = s[o], y[o], w[o]
    ws, wb = np.where(y == 1, w, 0), np.where(y == 0, w, 0)
    cb = np.cumsum(wb) - wb / 2
    return float(np.sum(ws * cb) / max(ws.sum() * wb.sum(), 1e-12))


def bkg_eff_table(D, thr, pt_bins, wps, wcol="w_train"):
    """ε_B (ważone do pT,|η| sygnału i surowe) przy WP, per bin pT — OOF."""
    out = {}
    for lo, hi in zip(pt_bins[:-1], pt_bins[1:]):
        key = f"{lo}-{hi}"
        sel = (D["pt"] >= lo) & (D["pt"] < hi)
        b = D[sel & (D["y"] == 0)]
        s = D[sel & (D["y"] == 1)]
        row = {"n_sig": int(len(s)), "n_bkg": int(len(b)),
               "auc": _auc(D.loc[sel, "score"].to_numpy(), D.loc[sel, "y"].to_numpy(), D.loc[sel, wcol].to_numpy())}
        for wp in wps:
            t = thr[key][f"{wp:.2f}"]
            pb = b["score"] >= t
            row[f"effB_{int(wp*100)}"] = float((b[wcol] * pb).sum() / max(b[wcol].sum(), 1e-12))
            row[f"effB_raw_{int(wp*100)}"] = float(pb.mean()) if len(b) else float("nan")
            row[f"effS_{int(wp*100)}"] = float((s[wcol] * (s["score"] >= t)).sum() / max(s[wcol].sum(), 1e-12))
        out[key] = row
    return out


# ---------------------------------------------------------------- train
def cmd_train(cfg, species, variant=None):
    variant = variant or cfg["training"]["variant"]
    S, B = SA.training_set(cfg, species)
    feats = F.usable(pd.concat([S.head(200000), B.head(200000)]), F.feature_list(cfg, variant))
    out = mdir(cfg, species, variant)
    D, meta = M.train(S, B, feats, cfg, out)
    tr = cfg["training"]
    s = D[D["y"] == 1]
    thr = M.wp_thresholds(s["score"].to_numpy(), s["w_train"].to_numpy(), s["pt"].to_numpy(), tr["pt_bins"],
                          tr["working_points"] + [cfg["fit"]["wp"]])
    meta.update({"species": species, "variant": variant, "thresholds": thr, "n_sig": int(len(S)), "n_bkg": int(len(B)),
                 "missing_columns": list(S.attrs.get("missing_columns", []))})
    _json(os.path.join(out, "meta.json"), meta)
    metrics = bkg_eff_table(D, thr, tr["pt_bins"], tr["working_points"])
    _json(os.path.join(out, "metrics_oof.json"), metrics)
    store.save(os.path.join(out, "oof.h5"),
               D[["score", "y", "w_train", "weight", "pt", "eta", "m", "fold", "run", "event"]].astype(
                   {"score": np.float32, "w_train": np.float32}))
    C.dump(cfg, os.path.join(out, "config_used.yaml"))
    log.info("train %s/%s: cechy=%s", species, variant, feats)
    for k, v in metrics.items():
        log.info("  pT %s: AUC=%.4f effB@80=%.4g effB@70=%.4g effB@50=%.4g", k, v["auc"], v["effB_80"], v["effB_70"], v["effB_50"])
    return meta, metrics


# ---------------------------------------------------------------- apply
def cmd_apply(cfg, species, variant=None, samples=None):
    variant = variant or cfg["training"]["variant"]
    ens = M.Ensemble(mdir(cfg, species, variant))
    sp = cfg["species"][species]
    od = sdir(cfg, species, variant)
    os.makedirs(od, exist_ok=True)
    samples = samples or list(cfg["paths"]["samples"])
    done = {}
    for key in samples:
        path = C.sample_path(cfg, key)
        srcs = C.sample_paths(cfg, key) if C.is_multi(cfg, key) else [path]
        if not srcs or not all(os.path.exists(q) for q in srcs):
            log.warning("apply: brak %s → pomijam", path if not C.is_multi(cfg, key) else [q for q in srcs if not os.path.exists(q)] or "(glob pusty)")
            continue
        attrs = {"source": path, "species": species, "variant": variant,
                 "note": "row = indeks wiersza w datasecie events pliku źródłowego"}
        if C.is_multi(cfg, key):  # etap P: wiele plików
            from .io import n_rows as _nr
            offs = np.concatenate([[0], np.cumsum([_nr(q, cfg["io"]["dataset"]) for q in srcs])]).astype(np.int64)
            attrs.update({"source": ";".join(srcs), "n_files": len(srcs), "file_offsets": json.dumps(offs.tolist()),
                          "note": "row = wiersz GLOBALNY w sklejonej liście plików (source, kolejność); "
                                  "plik i = searchsorted(file_offsets,row,'right')-1, lokalny = row-file_offsets[i]"})
        attrs["sample_meta"] = C.sample_meta_json(cfg, key)
        app = store.Appender(os.path.join(od, f"{key}.h5"), attrs)
        n = 0
        for L, cm in SA.iter_species(cfg, key, species):
            L = L[SA.in_range(L, sp).to_numpy()]
            if not len(L):
                continue
            sc, fold = ens.score(L)
            L = L.assign(score=sc, fold=fold, veto_ok=SA.dstar_veto_mask(L, sp))
            app.add(L[SCORE_COLS].astype({"m": np.float64, "pt": np.float32, "eta": np.float32, "dm": np.float32,
                                          "truth": np.float32}))
            n += len(L)
        app.close()
        done[key] = n
        log.info("apply %s/%s %s: %d kandydatów", species, variant, key, n)
    try:  # etap P: metadane próbek (pliki, joby, lumi) obok wyników
        _json(os.path.join(od, "samples_meta.json"), {k: json.loads(C.sample_meta_json(cfg, k)) for k in samples})
    except Exception as ex:
        log.warning("samples_meta.json: %s", ex)
    return done


def load_scores(cfg, species, variant, key, columns=None):
    p = os.path.join(sdir(cfg, species, variant), f"{key}.h5")
    return store.load(p, columns) if os.path.exists(p) else None


LIGHT = ["m", "pt", "score", "veto_ok", "truth"]


# ---------------------------------------------------------------- fit
def cmd_fit(cfg, species, variant=None):
    variant = variant or cfg["training"]["variant"]
    meta = json.load(open(os.path.join(mdir(cfg, species, variant), "meta.json")))
    thr = meta["thresholds"]
    sp = cfg["species"][species]
    fc = cfg["fit"]
    lo, hi = sp["mass_range"]
    bins = cfg["training"]["pt_bins"]
    mc = load_scores(cfg, species, variant, "signal_mc", LIGHT)
    if mc is None:
        raise RuntimeError("brak score dla signal_mc — najpierw `charmtag apply`")
    mc = mc[mc["veto_ok"].astype(bool)]
    res = {"species": species, "variant": variant, "wp": fc["wp"], "bins": {}}
    dcache = {}
    for ds in fc["datasets"]:
        d = load_scores(cfg, species, variant, ds, LIGHT + ["row"])
        if d is not None:
            dcache[ds] = d[d["veto_ok"].astype(bool)].reset_index(drop=True)
    refl_src = fc.get("reflection_source", "mc")  # S5: mc (domyślnie) | data — szablon refleksji K↔π z danych
    if refl_src == "data" and sp.get("reflection"):
        from . import refldata as RD
        ro = RD.opts(fc)
        dref = dcache.get(ro["sample"])
        if dref is None:
            dref = load_scores(cfg, species, variant, ro["sample"], LIGHT + ["row"])
            dref = dref[dref["veto_ok"].astype(bool)].reset_index(drop=True)
        dref = dref.assign(ms=RD.swap_mass_rows(cfg, ro["sample"], species, dref["row"].to_numpy()))
        res["refl_data"] = {"options": ro, "n_candidates": int(len(dref))}
    for wp in (None, fc["wp"], *[w for w in fc.get("extra_wps", []) if w != fc["wp"]]):  # S5: extra_wps (opcja)
        tagwp = "nocut" if wp is None else f"wp{int(wp*100)}"
        cut_mc = np.ones(len(mc), bool) if wp is None else M.pass_wp(mc["score"].to_numpy(), mc["pt"].to_numpy(), thr, wp)
        for blo, bhi in zip(bins[:-1], bins[1:]):
            key = f"{blo}-{bhi}"
            inb = (mc["pt"] >= blo) & (mc["pt"] < bhi)
            sig = mc[inb & cut_mc & mc["truth"].isin(cfg["signal"]["truth_values"])]
            shape = MF.fit_mc_shape(sig["m"], np.ones(len(sig)), lo, hi, sp["m0"], fc["signal_model"], fc["nbins"])
            rt, rr = None, 0.0
            if sp.get("reflection") and fc["reflection_mode"] != "off":
                ref = mc[inb & cut_mc & mc["truth"].isin(cfg["signal"]["reflection_values"])]
                if len(ref) > 20:
                    rt, _ = np.histogram(ref["m"], fc["nbins"], (lo, hi))
                    rr = fc["reflection_ratio"] if fc["reflection_ratio"] is not None else len(ref) / max(len(sig), 1)
                    if fc["reflection_mode"] == "gauss":
                        g = MF.fit_mc_shape(ref["m"], np.ones(len(ref)), lo, hi, sp["m0"], "gauss", fc["nbins"])
                        xc = np.linspace(lo, hi, fc["nbins"] + 1)
                        rt = np.diff(MF._gcdf(xc, g["mu"], g["s1"]))
                else:
                    log.warning("fit %s %s: brak kandydatów-refleksji w MC (truth %s) → refleksja wyłączona",
                                key, tagwp, cfg["signal"]["reflection_values"])
            entry = res["bins"].setdefault(key, {})
            rdt = None
            if refl_src == "data" and sp.get("reflection"):
                cut_r = np.ones(len(dref), bool) if wp is None else M.pass_wp(dref["score"].to_numpy(), dref["pt"].to_numpy(), thr, wp)
                inb_r = ((dref["pt"] >= blo) & (dref["pt"] < bhi)).to_numpy()
                g = dref[cut_r & inb_r]
                h_, npk, tsum = RD.template(g["m"], g["ms"], sp["m0"], lo, hi, fc["nbins"], ro)
                incl = npk < ro["min_peak"]
                if incl:  # mało statystyki w binie pT → kształt inkluzywny w pT (ta sama WP), suma z binu
                    gi = dref[cut_r]
                    hi_, _, ti = RD.template(gi["m"], gi["ms"], sp["m0"], lo, hi, fc["nbins"], ro)
                    h_ = hi_ * (tsum / ti) if ti > 0 else hi_
                hyb = {}
                if ro.get("centre", "mc") == "mc" and rt is not None:  # centrum z MC, ogony z danych (RD.hybrid)
                    xe_ = np.linspace(lo, hi, fc["nbins"] + 1)
                    h_, tails, fcm = RD.hybrid(h_, rt, xe_, sp["m0"], float(ro["win"]))
                    tsum = tails / max(1 - fcm, 1e-6)  # Σ szablonu ≈ ogony / (1 − ułamek centrum MC)
                    hyb = {"tails_data": tails, "centre_frac_mc": fcm}
                rdt = (h_, tsum)
                entry[f"refl_data_{tagwp}"] = {"n_peak": npk, "template_sum": tsum, "inclusive_shape": bool(incl),
                                               "template": np.asarray(h_).tolist(), **hyb}
            entry[f"mc_shape_{tagwp}"] = shape
            entry[f"mc_nsig_{tagwp}"] = int(len(sig))
            entry[f"refl_ratio_{tagwp}"] = rr
            for ds in fc["datasets"]:
                d = dcache.get(ds)
                if d is None:
                    continue
                d = d[(d["pt"] >= blo) & (d["pt"] < bhi)]
                if wp is not None:
                    d = d[M.pass_wp(d["score"].to_numpy(), d["pt"].to_numpy(), thr, wp)]
                if len(d) < 30:
                    entry[f"{ds}_{tagwp}"] = {"error": f"za mało kandydatów ({len(d)})"}
                    continue
                try:
                    if rdt is not None and rdt[1] > 0 and rr > 0:
                        o2 = dict(ro, norm=ro["norm"] if ds == ro["sample"] else "mc_ratio")
                        entry[f"{ds}_{tagwp}"] = RD.fit_with_template(d["m"].to_numpy(), lo, hi, sp["m0"], shape, fc, rdt[0], rr,
                                                                      rdt[1], o2)
                    else:
                        entry[f"{ds}_{tagwp}"] = MF.fit_data(d["m"].to_numpy(), lo, hi, sp["m0"], shape, fc, rt, rr)
                except Exception as e:  # noqa: BLE001
                    entry[f"{ds}_{tagwp}"] = {"error": repr(e)}
    # ε klasyfikatora z fitów: dane vs MC
    for key, e in res["bins"].items():
        for ds in fc["datasets"]:
            a, b = e.get(f"{ds}_nocut", {}), e.get(f"{ds}_wp{int(fc['wp']*100)}", {})
            if "N_sig" in a and "N_sig" in b and a["N_sig"] > 0:
                r = b["N_sig"] / a["N_sig"]
                er = r * np.hypot(b["N_sig_err"] / max(abs(b["N_sig"]), 1e-9), a["N_sig_err"] / a["N_sig"])
                e[f"{ds}_eff_fit"] = {"eff": r, "err": float(er), "eff_mc_nominal": fc["wp"]}
    out = os.path.join(cfg["paths"]["out_dir"], "fits", f"{species}_{variant}{fc.get('output_tag') or ''}.json")  # S5: output_tag
    try:  # etap P: metadane próbek danych (liczba plików, joby, lumi) w wyniku fitu
        res.setdefault("samples_meta", {ds: json.loads(C.sample_meta_json(cfg, ds)) for ds in fc["datasets"]
                                        if ds in cfg["paths"]["samples"]})
    except Exception as ex:
        log.warning("samples_meta w fit: %s", ex)
    _json(out, res)
    for key, e in res["bins"].items():
        for ds in fc["datasets"]:
            f_ = e.get(f"{ds}_wp{int(fc['wp']*100)}", {})
            if "N_sig" in f_:
                log.info("fit %s pT %s %s: N=%.1f±%.1f σ=%.1f MeV χ²/ndf=%.1f/%d", species, key, ds, f_["N_sig"],
                         f_["N_sig_err"], 1e3 * f_["sigma_eff"], f_["chi2"], f_["ndf"])
    return res


# ---------------------------------------------------------------- validate
def rwin(d, sp, thr, wps, region="inner"):
    """Rzeźbienie masy: R = [N_A/N_B]_po / [N_A/N_B]_przed (A = pasmo wewn./okno, B = zewn. sidebandy)."""
    dm = np.abs(d["m"].to_numpy() - sp["m0"])
    Bm = dm > sp["sideband_min"]
    if region == "inner":
        A = (dm > sp["inner_band"][0]) & (dm < sp["inner_band"][1])
    elif region == "lowhigh":  # kształt tła: lewy vs prawy sideband (bez sygnału)
        A = Bm & (d["m"].to_numpy() < sp["m0"])
        Bm = Bm & (d["m"].to_numpy() > sp["m0"])
    else:
        A = dm < sp["signal_window"]
    out = {}
    s, pt = d["score"].to_numpy(), d["pt"].to_numpy()
    base = A.sum() / max(Bm.sum(), 1)
    for wp in wps:
        p = M.pass_wp(s, pt, thr, wp)
        nA, nB = (A & p).sum(), (Bm & p).sum()
        r = (nA / max(nB, 1)) / max(base, 1e-12)
        err = r * np.sqrt(1 / max(nA, 1) + 1 / max(nB, 1)) if nA and nB else float("nan")
        out[f"{int(wp*100)}"] = {"R": float(r), "err": float(err), "nA": int(nA), "nB": int(nB)}
    return out


def dstar_control(cfg, species, variant, key, thr, wps):
    """D0 z D*: kandydaci D0 dopasowani do D* (Δm) w tym samym zdarzeniu; ε(WP) z odejmowaniem sidebandu Δm."""
    dc = cfg["dstar_control"]
    sc = load_scores(cfg, species, variant, key)
    if sc is None:
        return None
    ds = SA.load_species(cfg, key, species, decay_type=dc["dstar_decay_type"], nd=3, tagx="dstar",
                         select=lambda L: (L["dm"] < dc["dm_sideband"][1]).to_numpy())
    if not len(ds) or not np.isfinite(ds["dm"]).any():
        return {"error": "brak kandydatów D* / kolumny Δm"}
    ds = ds.assign(m0_from=ds["m"] - ds["dm"])[["run", "event", "m0_from", "pt", "eta", "dm"]]
    j = sc.drop(columns=["dm"]).merge(ds, on=["run", "event"], suffixes=("", "_ds"))
    ok = (np.abs(j["m"] - j["m0_from"]) < dc["match_tol_m"]) & (j["pt"] <= j["pt_ds"] * 1.001) & \
         (j["pt"] > 0.8 * j["pt_ds"]) & (np.abs(j["eta"] - j["eta_ds"]) < 0.1)
    ok &= np.abs(j["m"] - cfg["species"][species]["m0"]) < cfg["species"][species]["signal_window"]
    j = j[ok].drop_duplicates(subset=["row"])
    out = {"method": "fit Δm (Gauss + próg) bez cięcia i po WP" if key != "signal_mc" else "MC truth: zliczenie w oknie Δm"}
    for blo, bhi in zip(cfg["training"]["pt_bins"][:-1], cfg["training"]["pt_bins"][1:]):
        b = ((j["pt"] >= blo) & (j["pt"] < bhi)).to_numpy()
        row = {}
        if key == "signal_mc":
            tm = j["truth"].isin(cfg["signal"]["truth_values"]).to_numpy()
            win = ((j["dm"] > dc["dm_signal"][0]) & (j["dm"] < dc["dm_signal"][1])).to_numpy()
            sel = b & tm & win
            n0 = sel.sum()
            row["N_all"] = int(n0)
            for wp in wps:
                p = M.pass_wp(j["score"].to_numpy(), j["pt"].to_numpy(), thr, wp)
                eff = (sel & p).sum() / n0 if n0 else float("nan")
                row[f"eff_{int(wp*100)}"] = float(eff)
                row[f"err_{int(wp*100)}"] = float(np.sqrt(max(eff * (1 - eff), 0) / n0)) if n0 else float("nan")
        elif cfg.get("validate", {}).get("dstar_refit"):   # etap V (AA): jednoczesny fit pass/fail, kształt z MC
            from . import dstar_refit as DR
            row.update(DR.control_row(j[b], thr, wps, cfg, species, M))
            out[f"{blo}-{bhi}"] = row
            continue
        else:
            f0 = MF.fit_dm(j.loc[b, "dm"].to_numpy())
            row["fit_all"] = f0
            for wp in wps:
                p = M.pass_wp(j["score"].to_numpy(), j["pt"].to_numpy(), thr, wp)
                f1 = MF.fit_dm(j.loc[b & p, "dm"].to_numpy(), mu0=f0["mu"], s0=f0["sigma"], fix_shape=True)
                row[f"fit_{int(wp*100)}"] = f1
                if f0["N"] > 0:
                    eff = f1["N"] / f0["N"]
                    # próbka po cięciu ⊂ bez cięcia: cov(N1, N0) ≈ σ1²
                    v = (f1["N_err"] ** 2 * (1 - 2 * eff) + eff ** 2 * f0["N_err"] ** 2) / f0["N"] ** 2
                    err = np.sqrt(v) if v > 0 else np.sqrt(max(eff * (1 - eff), 0) / f0["N"])
                else:
                    eff, err = float("nan"), float("nan")
                row[f"eff_{int(wp*100)}"] = float(eff)
                row[f"err_{int(wp*100)}"] = float(err)
        out[f"{blo}-{bhi}"] = row
    return out


def cmd_validate(cfg, species, variants=("H", "frozen7")):
    tr = cfg["training"]
    sp = cfg["species"][species]
    wps = tr["working_points"]
    res = {"species": species, "variants": {}}
    for v in variants:
        md = mdir(cfg, species, v)
        if not os.path.exists(os.path.join(md, "meta.json")):
            continue
        meta = json.load(open(os.path.join(md, "meta.json")))
        thr = meta["thresholds"]
        r = {"metrics_oof": json.load(open(os.path.join(md, "metrics_oof.json"))), "importance": meta.get("importance"),
             "features": meta["features"]}
        d = load_scores(cfg, species, v, "data_none", LIGHT)
        if d is not None:
            d = d[d["veto_ok"].astype(bool)]
            r["R_win_data_inner"] = rwin(d, sp, thr, wps, "inner")
            r["R_sb_data_lowhigh"] = rwin(d, sp, thr, wps, "lowhigh")
            r["R_win_data_inner_bins"] = {f"{lo}-{hi}": rwin(d[(d["pt"] >= lo) & (d["pt"] < hi)], sp, thr, wps, "inner")
                                          for lo, hi in zip(tr["pt_bins"][:-1], tr["pt_bins"][1:])}
        mb = load_scores(cfg, species, v, "incl_mc", LIGHT)
        if mb is not None:
            bk = mb[mb["veto_ok"].astype(bool) & ~mb["truth"].isin(cfg["signal"]["truth_values"] + cfg["signal"]["reflection_values"])]
            r["R_win_mcbkg_window"] = rwin(bk, sp, thr, wps, "window")
            sg = mb[mb["truth"].isin(cfg["signal"]["truth_values"])]
            r["effS_incl_mc"] = {f"{int(wp*100)}": float(M.pass_wp(sg["score"].to_numpy(), sg["pt"].to_numpy(), thr, wp).mean())
                                 if len(sg) else float("nan") for wp in wps}
        if species == "D0":
            r["dstar_data"] = dstar_control(cfg, species, v, "data_none", thr, wps)
            r["dstar_mc"] = dstar_control(cfg, species, v, "signal_mc", thr, wps)
        res["variants"][v] = r
    if len(res["variants"]) == 2:
        cmp_ = {}
        for key in res["variants"]["H"]["metrics_oof"]:
            a, b = res["variants"]["H"]["metrics_oof"][key], res["variants"]["frozen7"]["metrics_oof"][key]
            cmp_[key] = {f"gain_rej_{int(wp*100)}": (1 / max(a[f"effB_{int(wp*100)}"], 1e-12)) /
                         (1 / max(b[f"effB_{int(wp*100)}"], 1e-12)) for wp in wps}
        res["H_vs_frozen7"] = cmp_
    _json(os.path.join(cfg["paths"]["out_dir"], "validation", f"{species}.json"), res)
    return res
