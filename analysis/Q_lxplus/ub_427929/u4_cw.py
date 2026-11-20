#!/usr/bin/env python3
"""UB427 c_w (18.11.2026, projekt A): czynniki wag zdarzeń c_w dla runu 427929 i plik mb_rates dla `charmtag pairs`.

Nominał (428770): c_w = (X/D*)_w2 / (X/D*)_nieważone w MC mb_nd, w2(n_trk, lead_pt) = p_dane428770 / p_MC (etap1 w2d.py, referencja
pliki 000021/000030) — zamrożone w mb_rates_2d.json. Dla 427929 trzeba w2 do rozkładu 427929 (klasa MB427929: GRL ∧ HLT_mb_mbts_L1MBTS_1).
Metoda (stosunek; znosi różnice v2.0 ↔ v1.1 w definicji poprawnych kandydatów MC i brak 2 z 80 plików mb_nd):
  te same części MC (gn_root.py --kind mc na ROOT v2.0 mb_nd) i ta sama funkcja analyze2 (kopia etap1, gn_combine.py) z mapą w2 z
  (i) 428770 — referencja po kotwicach h5 v1.1 _000021/_000030 (jak nominał; gdy brak kotwic → wszystkie dane 428770),
  (ii) 427929 klasa MB;  R_X = c_X(427929) / c_X(428770_ref);  c_X(427929) := c_X(nominał) · R_X.
Wyjście: mb_rates_427929.json = nominalny mb_rates_2d.json z polami ważonymi X przeskalowanymi przez R_X (D* bez zmian ⇒
charmtag.pairs.build.weight_factors daje c_w(427929)), cw.json (wszystkie c_w, R_X, n_eff par w MC, mapy, kontrole).
  u4_cw.py --mc DIR --d428 DIR --d427 DIR --nominal-mb mb_rates_2d.json [--anchors A1.h5,A2.h5] --out DIR [--no-class]"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gn_combine as G  # noqa: E402  (kopia z Q_lxplus/gen_norm — analyze2, w2_map, load_parts, kotwice, weight_factors)

CLS = ("ev_grl", "ev_tp_HLT_mb_mbts_L1MBTS_1")
KEYS = {"D0": ("D0_correct_per_event", "w"), "Dplus": ("Dplus_correct_per_event", "w"),
        "D0D0": ("pairs_D0D0", "per_event_w"), "DpD0": ("pairs_DpD0", "per_event_w"), "DpDp": ("pairs_DpDp", "per_event_w")}


def mc_frames(MC):
    evs, cs = [], []
    for k, (tg, meta, z) in enumerate(MC):
        e = pd.DataFrame({c[3:]: z[c] for c in z if c.startswith("ev_")}); e["file"] = np.int16(k); evs.append(e)
        c = pd.DataFrame({c[2:]: z[c] for c in z if c.startswith("c_")}); c["file"] = np.int16(k); cs.append(c)
    return pd.concat(evs, ignore_index=True), pd.concat(cs, ignore_index=True)


def data_xy(DA, use_class, need_class):
    nt, lp, info = [], [], {"files": len(DA), "events": 0, "class": 0}
    for _, _, z in DA:
        m = np.ones(len(z["ev_n_trk"]), bool)
        if use_class:
            if all(c in z for c in CLS):
                m = (z[CLS[0]] == 1) & (z[CLS[1]] == 1)
            elif need_class:
                sys.exit(f"brak kolumn klasy {CLS} w części danych 427929 (ROOT v1.1?)")
        info["events"] += int(len(m)); info["class"] += int(m.sum())
        nt.append(z["ev_n_trk"][m]); lp.append(z["ev_lead_pt"][m])
    return np.concatenate(nt), np.concatenate(lp), info


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", required=True); ap.add_argument("--d428", required=True); ap.add_argument("--d427", required=True)
    ap.add_argument("--nominal-mb", required=True); ap.add_argument("--anchors", default=""); ap.add_argument("--out", required=True)
    ap.add_argument("--no-class", action="store_true", help="test lokalny: atrapy bez kolumn triggera/GRL")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    MC = G.load_parts(a.mc, "mc"); D8 = G.load_parts(a.d428, "data"); D7 = G.load_parts(a.d427, "data")
    ev, c = mc_frames(MC)
    G.say(f"MC: {len(MC)} plików, {len(ev)} zdarzeń, {len(c)} poprawnych kand.; dane 428770 {len(D8)} plików, 427929 {len(D7)} plików")
    maps = {}
    if a.anchors:
        A = G.load_anchors(a.anchors); masks, anch = G.ref_by_anchors(D8, A)
        dref = [{k: v[m] for k, v in z.items() if k.startswith("ev_")} for (_, _, z), m in zip(D8, masks) if m.any()]
        maps["428770_ref"] = (np.concatenate([z["ev_n_trk"] for z in dref]), np.concatenate([z["ev_lead_pt"] for z in dref]),
                              {"anchors": {k: anch[k] for k in ("anchors", "anchors_found_unique", "ref_events")}})
    maps["428770_all"] = data_xy(D8, use_class=False, need_class=False)
    maps["427929_mb"] = data_xy(D7, use_class=not a.no_class, need_class=not a.no_class)
    maps["427929_all"] = data_xy(D7, use_class=False, need_class=False)
    E = ev[["file", "entry", "n_trk", "lead_pt", "n_true_dstar"]]
    rep = {"mc_files": len(MC), "mc_events": int(len(ev)), "samples": {}}
    cw = {}
    for nm, (nt, lp, info) in maps.items():
        W, hd, hm, lost = G.w2_map(nt, lp, ev.n_trk.to_numpy(), ev.lead_pt.to_numpy())
        r = G.analyze2(E, c, W)
        cw[nm] = G.weight_factors(r)
        rep["samples"][nm] = {"info": info, "data_events_in_map": int(hd.sum()), "data_mass_in_empty_mc_cells": lost,
                              "mean_ntrk": float(np.mean(nt)) if len(nt) else None, "c_w": cw[nm],
                              "n_eff_pairs": {ch: r[f"pairs_{ch}"]["n_eff"] for ch in G.CH}, "n_pairs_mc": {ch: r[f"pairs_{ch}"]["n"] for ch in G.CH},
                              "dstar_per_event_w": r["dstar_true_per_event"]["w"], "rates": r}
        np.save(os.path.join(a.out, f"w2d_{nm}.npy"), W)
        G.say(nm, "c_w", {k: round(v, 4) for k, v in cw[nm].items()}, "n_eff par", {k: round(v, 1) for k, v in rep["samples"][nm]["n_eff_pairs"].items()})
    base = "428770_ref" if "428770_ref" in cw else "428770_all"
    R = {x: cw["427929_mb"][x] / cw[base][x] for x in KEYS}
    R_all = {x: cw["427929_mb"][x] / cw["428770_all"][x] for x in KEYS}
    bad = [x for x, v in R.items() if not np.isfinite(v) or v <= 0]
    if bad and not a.no_class:
        json.dump(rep, open(os.path.join(a.out, "cw.json"), "w"), indent=1, default=float)
        sys.exit(f"STOP: R_X nieokreślone dla {bad} (za mało par MC?) — c_w dla 427929 niepoliczalne")
    nom = json.load(open(a.nominal_mb))
    cw_nom = G.weight_factors(nom)
    new = json.loads(json.dumps(nom))
    for x, (blk, fld) in KEYS.items():
        new[blk][fld] = nom[blk][fld] * R[x]
    new["_u4_cw"] = {"note": "nominalny mb_rates_2d.json z polami ważonymi × R_X = c_X(427929 MB)/c_X(" + base + ") — c_w dla 427929", "R": R}
    cw_new = G.weight_factors(new)
    rep.update({"base_428770": base, "R": R, "R_vs_alldata": R_all, "c_w_nominal": cw_nom, "c_w_427929": cw_new,
                "check_c_w_427929_eq_nom_times_R": {x: abs(cw_new[x] - cw_nom[x] * R[x]) < 1e-9 for x in KEYS},
                "repro_428770_vs_nominal": {x: cw[base][x] / cw_nom[x] for x in KEYS}})
    json.dump(new, open(os.path.join(a.out, "mb_rates_427929.json"), "w"), indent=1)
    json.dump(rep, open(os.path.join(a.out, "cw.json"), "w"), indent=1, default=float)
    G.say("R_X (427929/428770):", {k: round(v, 4) for k, v in R.items()}, "c_w 427929:", {k: round(v, 4) for k, v in cw_new.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
