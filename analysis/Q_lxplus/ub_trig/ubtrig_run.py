#!/usr/bin/env python3
"""UBTRIG (18.11.2026, projekt A) — μ par D–D per rozłączna klasa triggera w runie 428770 (ODSŁONIĘTA próbka; zatwierdzone 18.11).

Pytanie: czy nadwyżka μ = 1,50 zależy od triggera (81 % par D⁰←D* z L1_J12). Ta sama próbka, nogi, model i fit co odsłonięcie 05.11
(`ub_jk`: pair_inputs, pairs_full.npz, legs_data.parquet, charmtag.pairs.model/fit2d); zmienia się tylko podział zdarzeń:
  klasa zdarzenia z bitów triggera ROOT v2.0 (trg_root_extract.py, złączenie po (run, event) — w danych klucz jednoznaczny),
  klasy rozłączne: J12only = J12 ∧ ¬MBTS, MBTS_j20 = HLT_j20_L1MBTS_2, MBTS_noj20 = MBTS ∧ ¬j20, other = ¬J12 ∧ ¬MBTS.
Wejścia klasy c (jak tryb `rederive` jackknife'u):
  n_pair, n_real × r_c, r_c = udział N_D* w klasie (fit Δm pass/fail, D⁰←D* z B_trig bt_data_428770: |m−m0| < 55 MeV, score ≥ WP90, pT 3–30);
  dataleg, SB×SB = tylko z klasy c; nev_legs = ndata = nev · f_c (udział zdarzeń klasy w tabeli trg).
  ⇒ μ_c mierzy [pary/D*]_dane,c / [pary/D*]_MC (MC minbias bez triggera, wagi w2 jak nominał). μ_c ≠ μ_c' ⇒ obciążenie triggerem.
Błąd μ_c: stat. fitu (profil) ⊕ δr_c/r_c (osobno). Kontrole: Σ klas = pairs_full; „all” przez tę samą ścieżkę = fit2d_unblind.
  ubtrig_run.py --pairs-dir P --legs L.parquet --trg 'trg_*.npz' --bt bt_data_428770.npz --bt-mc bt_fd_nd.npz --out DIR [--selftest]"""
import argparse
import copy
import glob
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "ub_jk"))
CLASSES = ["J12only", "MBTS_j20", "MBTS_noj20", "other"]
WP90 = 0.38729336857795715
TOL_MU = 2e-3


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def classes_from_tp(tp):
    """tp: {chain: bool array} → int klasa 0..3 (CLASSES) — rozłączne i wyczerpujące."""
    import btrig_common as C
    n = len(next(iter(tp.values())))
    g = lambda c: np.asarray(tp[c]).astype(bool) if c in tp else np.zeros(n, bool)  # noqa: E731
    mbts = np.zeros(n, bool)
    for c in C.MBTS_CHAINS:
        mbts |= g(c)
    j12, j20 = g(C.J12), g("HLT_j20_L1MBTS_2")
    k = np.full(n, 3, np.int8)
    k[mbts & ~j20] = 2
    k[j20] = 1
    k[j12 & ~mbts] = 0
    return k


def load_trg(paths):
    import btrig_common as C
    run, ev, tp = [], [], {c: [] for c in C.TRIG_CHAINS}
    for p in paths:
        z = np.load(p, allow_pickle=False)
        run.append(z["run"].astype(np.int64)); ev.append(z["event"].astype(np.int64))
        for c in tp:
            k = f"tp_{c}"
            tp[c].append((z[k] == 1) if k in z.files else np.zeros(len(z["run"]), bool))
    run, ev = np.concatenate(run), np.concatenate(ev)
    tp = {c: np.concatenate(v) for c, v in tp.items()}
    key = C.ev_key(run, ev)
    key, first = np.unique(key, return_index=True)
    return key, classes_from_tp({c: v[first] for c, v in tp.items()})


def lookup_cls(tkey, tcls, keys, fallback):
    """klasa zdarzenia po kluczu (tkey posortowane); brak klucza → fallback."""
    pos = np.clip(np.searchsorted(tkey, keys), 0, len(tkey) - 1)
    hit = tkey[pos] == keys
    return np.where(hit, tcls[pos], fallback).astype(np.int8), int((~hit).sum())


def permute_within_strata(tkey, tcls, strata_keys, seed):
    """Test zerowy: klasy przetasowane losowo w obrębie warstw (zdarzenie ma nogi? ma kandydata D*?) — liczności klas w każdej warstwie
    zachowane dokładnie, związek klasa ↔ trigger zerwany."""
    rng = np.random.default_rng(seed)
    st = np.zeros(len(tkey), np.int8)
    for b, ks in enumerate(strata_keys):
        st |= (np.isin(tkey, ks).astype(np.int8) << b)
    out = tcls.copy()
    for v in np.unique(st):
        ix = np.flatnonzero(st == v)
        out[ix] = tcls[ix][rng.permutation(len(ix))]
    return out


def dstar_fractions(bt, btmc, cls_lookup=None):
    """r_c = udział N_D* w klasie c (pass = c, fail = reszta), D⁰←D* z bt; kształt Δm z prawdy fd_nd.
    cls_lookup = (tkey, tcls): klasa z tabeli zdarzeń po (run, event) zamiast z bitów bt (test zerowy)."""
    import btrig_common as C
    import dstar_refit_aa as DR                                    # kopia AA_merge, nie charmtag z PF/src
    z = np.load(bt, allow_pickle=False)
    q = (np.abs(z["m"] - C.M0_D0) < C.NOM_WIN) & (z["score"] >= WP90) & (z["pt"] >= 3.0) & (z["pt"] < 30.0)
    tp = {c: (z[f"e_tp_{c}"] == 1) for c in C.TRIG_CHAINS if f"e_tp_{c}" in z.files}
    k = classes_from_tp(tp)
    if cls_lookup is not None:
        k, nmiss = lookup_cls(cls_lookup[0], cls_lookup[1], C.ev_key(z["run"].astype(np.int64), z["event"].astype(np.int64)), k)
        log(f"klasy D* z tabeli zdarzeń (test zerowy); bez klucza: {nmiss}")
    m = np.load(btmc, allow_pickle=False)
    qm = (np.abs(m["m"] - C.M0_D0) < C.NOM_WIN) & (m["score"] >= WP90) & (m["pt"] >= 3.0) & (m["pt"] < 30.0) \
        & np.isin(m["truth"], [0, 1]) & np.isin(m["truth_ds"], [0, 1])
    shape = DR.shape_from_mc(m["dm"][qm] * 1e3)
    dm = z["dm"][q]; kq = k[q]
    out = {}
    for i, c in enumerate(CLASSES):
        f = C.fit_eps(dm[kq == i], dm[kq != i], shape, interval=True)
        out[c] = {"r": float(f["eps"]), "r_err": float(C.eps_sym_err(f)), "N_tot": float(f["N"]), "determined": bool(f["determined"]),
                  "n_cand": int((kq == i).sum())}
        log(f"N_D* {c}: r = {out[c]['r']:.4f} ± {out[c]['r_err']:.4f} (kandydaci {out[c]['n_cand']}, określony {out[c]['determined']})")
    return out


def class_inputs(I, r, fev, DL, SB, keys, c):
    from charmtag.syst import files as F
    I2 = F.scaled_inputs(I, r)                                     # n_pair, n_real × r (ndata nadpisane niżej)
    I2.dataleg = {key: np.asarray(DL[key][c], float) for key in I.dataleg}
    I2.sbsb = {key: int(round(SB[c, j])) for j, key in enumerate(keys)}
    I2.nev_legs = I.nev_legs * fev
    I2.ndata = I.ndata * fev
    return I2


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs-dir"); ap.add_argument("--legs"); ap.add_argument("--trg"); ap.add_argument("--bt"); ap.add_argument("--bt-mc")
    ap.add_argument("--out", required=True); ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--null", type=int, default=None, help="test zerowy: klasy przetasowane w warstwach (ziarno)")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    if os.environ.get("CHARM_UNBLIND") != "1":
        raise SystemExit("UBTRIG: wymaga CHARM_UNBLIND=1")
    import pandas as pd
    import ubjk_run as U
    charmtag, CHANNELS, Model, load_inputs, fit_channel, F = U._charmtag()
    t0 = time.time()
    I = load_inputs(os.path.join(a.pairs_dir, "pair_inputs"))
    H0 = dict(np.load(os.path.join(a.pairs_dir, "pairs_full.npz")))
    L = pd.read_parquet(a.legs)
    rep = {"charmtag": charmtag.__file__, "classes": CLASSES, "n_legs": int(len(L))}
    # --- klasy zdarzeń
    if a.selftest:                                                  # sztuczne klasy z hasha zdarzenia (oczekiwane: μ_c zgodne)
        run = pd.to_numeric(L.source_file).to_numpy(np.int64); ev = L.entry.to_numpy(np.int64)
        u = np.unique(np.stack([run, ev], 1), axis=0)
        import btrig_common as C
        tkey = np.sort(C.ev_key(u[:, 0], u[:, 1]))
        with np.errstate(over="ignore"):
            tcls = ((tkey.astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)) >> np.uint64(61)).astype(np.int64) % 4
        tcls = np.choose(tcls, [0, 0, 1, 2]).astype(np.int8)        # nierówne udziały 0,5/0,25/0,25/0
        tcls[:: 7] = 3
        frac = {c: float((tcls == i).mean()) for i, c in enumerate(CLASSES)}
        dsf = {c: {"r": frac[c], "r_err": 0.0, "N_tot": np.nan, "determined": True, "n_cand": 0} for c in CLASSES}
    else:
        tkey, tcls = load_trg(sorted(glob.glob(a.trg)))
        if a.null is not None:
            import btrig_common as C
            zb = np.load(a.bt, allow_pickle=False)
            bk = C.ev_key(zb["run"].astype(np.int64), zb["event"].astype(np.int64))
            lk0 = C.ev_key(pd.to_numeric(L.source_file).to_numpy(np.int64), L.entry.to_numpy(np.int64))
            tcls = permute_within_strata(tkey, tcls, [np.unique(lk0), np.unique(bk)], a.null)
            rep["null_seed"] = a.null
            log(f"TEST ZEROWY ziarno {a.null}: klasy przetasowane w warstwach (nogi × D*)")
            dsf = dstar_fractions(a.bt, a.bt_mc, cls_lookup=(tkey, tcls))
        else:
            dsf = dstar_fractions(a.bt, a.bt_mc)
        bad = [c for c in CLASSES if not dsf[c]["determined"] or not np.isfinite(dsf[c]["r"])]
        if bad:
            raise SystemExit(f"STOP: r_c nieokreślone dla {bad} — brak wyniku (bieg 18.11 dał NaN przez starą fit_passfail z PF/src)")
    import btrig_common as C
    rep["n_trg_events"] = int(len(tkey))
    fev = {c: float((tcls == i).mean()) for i, c in enumerate(CLASSES)}
    rep["event_fraction"] = fev; rep["dstar_fraction"] = dsf
    lk = C.ev_key(pd.to_numeric(L.source_file).to_numpy(np.int64), L.entry.to_numpy(np.int64))
    pos = np.searchsorted(tkey, lk); pos = np.clip(pos, 0, len(tkey) - 1)
    hit = tkey[pos] == lk
    rep["legs_unmatched"] = int((~hit).sum()); rep["legs_unmatched_frac"] = float((~hit).mean())
    log(f"tabela trg: {len(tkey)} zdarzeń; udziały zdarzeń {({c: round(v, 4) for c, v in fev.items()})}; nogi bez klasy {rep['legs_unmatched']}")
    L = L.assign(file_idx=np.where(hit, tcls[pos], len(CLASSES)).astype(np.int64))
    nk = len(CLASSES) + 1                                           # ostatnia grupa = bez klasy (powinna być pusta)
    H, DL, SB, keys = U.perfile_from_legs(L, I, nk)
    rep["sum_equals_pairs_full"] = {ch: bool(np.array_equal(H[ch].sum(0), H0[ch])) for ch in CHANNELS}
    rep["pairs_per_class"] = {ch: {c: int(H[ch][i].sum()) for i, c in enumerate(CLASSES + ["bez_klasy"])} for ch in CHANNELS}
    ulev = L[["file_idx", "source_file", "entry"]].drop_duplicates().file_idx.value_counts()
    rep["events_with_legs_fraction"] = {c: float(ulev.get(i, 0) / ulev.sum()) for i, c in enumerate(CLASSES)}
    # --- nominał (kontrola) + „all” tą samą ścieżką
    fits = {}
    fits["nominal"] = U.replica_fit(Model(I), H0, list(CHANNELS), combined=True, interval=True)
    ref = json.load(open(os.path.join(a.pairs_dir, "fit2d_unblind.json")))
    rep["nominal_ok"] = bool(all(abs(fits["nominal"][ch]["mu_hat"] - ref["channels"][ch]["mu_hat"]) < TOL_MU for ch in CHANNELS))
    log("nominał:", {ch: round(fits["nominal"][ch]["mu_hat"], 4) for ch in CHANNELS}, "wspólne", round(fits["nominal"]["combined"]["mu"], 4), "ok", rep["nominal_ok"])
    for i, c in enumerate(CLASSES):
        r = dsf[c]["r"]
        if not (np.isfinite(r) and r > 0) or fev[c] <= 0:
            fits[c] = {"skipped": "r_c lub udział zdarzeń = 0"}; continue
        Ic = class_inputs(I, r, fev[c], DL, SB, keys, i)
        try:
            fits[c] = U.replica_fit(Model(Ic), {ch: H[ch][i] for ch in CHANNELS}, list(CHANNELS), combined=True, interval=True)
        except Exception as ex:  # noqa: BLE001
            fits[c] = {"error": repr(ex)}; log("fit", c, "BŁĄD", ex); continue
        fits[c]["r_rel_err"] = dsf[c]["r_err"] / r if r > 0 else None
        cm = fits[c].get("combined", {})
        log(f"{c}: μ {({ch: round(fits[c][ch]['mu_hat'], 3) for ch in CHANNELS})} wspólne {cm.get('mu')} 68% {cm.get('mu_68')}")
    rep["fits"] = fits
    # --- zgodność klas (wspólne μ; błąd = σ_68 ⊕ δr/r·μ)
    vals, errs, used = [], [], []
    for c in CLASSES:
        cm = fits.get(c, {}).get("combined") if isinstance(fits.get(c), dict) else None
        if cm and cm.get("sigma_68") and np.isfinite(cm["mu"]):
            vals.append(cm["mu"]); errs.append(math.hypot(cm["sigma_68"], cm["mu"] * (fits[c].get("r_rel_err") or 0))); used.append(c)
    rep["consistency"] = {"classes": used, **(C.chi2_consistency(vals, errs) if len(vals) >= 2 else {})}
    if "J12only" in used and len(used) > 1:
        j = used.index("J12only")
        rest = [x for x in range(len(used)) if x != j]
        w = np.array([1 / errs[x] ** 2 for x in rest]); mr = float(np.sum(w * np.array([vals[x] for x in rest])) / w.sum()); er = float(1 / math.sqrt(w.sum()))
        rep["J12_vs_rest"] = {"mu_J12": vals[j], "err_J12": errs[j], "mu_rest": mr, "err_rest": er, "diff_sigma": (vals[j] - mr) / math.hypot(errs[j], er)}
    rep["sec"] = time.time() - t0
    json.dump(U._clean(rep), open(os.path.join(a.out, "ubtrig.json"), "w"), indent=1)
    open(os.path.join(a.out, "UBTRIG.md"), "w").write(to_md(rep, CHANNELS))
    log("OK", os.path.join(a.out, "UBTRIG.md"))
    return 0


def to_md(rep, CHANNELS):
    f = rep["fits"]
    L = ["# UBTRIG — μ par D–D per klasa triggera, run 428770 (odsłonięta próbka; zatwierdzone 18.11)"
         + (f" — TEST ZEROWY (klasy losowe, ziarno {rep['null_seed']})" if rep.get("null_seed") is not None else ""), "",
         "Ta sama próbka/nogi/model co odsłonięcie 05.11; klasa zdarzenia z bitów ROOT v2.0 po (run, event). μ_c = [pary/D*]_dane,c / [pary/D*]_MC.",
         f"Normalizacja klasy: n_pair, n_real × r_c (udział N_D*, fit Δm); tło fałszywych nóg i SB×SB tylko z klasy. Zdarzeń w tabeli trg: {rep['n_trg_events']}; "
         f"nogi bez klasy: {rep['legs_unmatched']} ({100 * rep['legs_unmatched_frac']:.3f} %). Σ klas = pairs_full: {rep['sum_equals_pairs_full']}. "
         f"Nominał odtworzony: {rep['nominal_ok']}.", "",
         "| klasa | udział zdarzeń | udział zd. z nogami | r_c = udział N_D* | pary D0D0 / DpD0 / DpDp | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | δr/r |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    nom = f["nominal"]
    L.append(f"| nominał (wszystkie) | 1 | 1 | 1 | " + " / ".join(str(sum(rep['pairs_per_class'][ch].values())) for ch in CHANNELS) + " | "
             + " | ".join(f"{nom[ch]['mu_hat']:.2f} ± {nom[ch]['mu_err_hesse']:.2f}" for ch in CHANNELS)
             + f" | **{nom['combined']['mu']:.2f}** [{nom['combined']['mu_68'][0]:.2f}; {nom['combined']['mu_68'][1]:.2f}] | – |")
    for c in rep["classes"]:
        d = rep["dstar_fraction"][c]; x = f.get(c, {})
        pc = " / ".join(str(rep["pairs_per_class"][ch][c]) for ch in CHANNELS)
        if "combined" not in x:
            L.append(f"| {c} | {rep['event_fraction'][c]:.4f} | {rep['events_with_legs_fraction'][c]:.4f} | {d['r']:.3f} ± {d['r_err']:.3f} | {pc} | – | – | – | {x.get('skipped') or x.get('error')} | – |")
            continue
        cm = x["combined"]; lo, hi = [np.nan if v is None else v for v in (cm.get("mu_68") or [np.nan, np.nan])]
        L.append(f"| {c} | {rep['event_fraction'][c]:.4f} | {rep['events_with_legs_fraction'][c]:.4f} | {d['r']:.3f} ± {d['r_err']:.3f} | {pc} | "
                 + " | ".join(f"{x[ch]['mu_hat']:.2f} ± {x[ch]['mu_err_hesse']:.2f}" for ch in CHANNELS)
                 + f" | **{cm['mu']:.2f}** [{lo:.2f}; {hi:.2f}] | {x.get('r_rel_err') or 0:.3f} |")
    cs = rep.get("consistency", {})
    L += ["", f"Zgodność klas (wspólne μ, błąd σ_68 ⊕ δr/r): klasy {cs.get('classes')}, średnia {cs.get('mean', float('nan')):.3f} ± {cs.get('err', float('nan')):.3f}, "
          f"χ²/ndf = {cs.get('chi2', float('nan')):.2f}/{cs.get('ndf', 0)}, p = {cs.get('p', float('nan')):.3f}."]
    jr = rep.get("J12_vs_rest")
    if jr:
        L.append(f"J12only vs reszta: {jr['mu_J12']:.2f} ± {jr['err_J12']:.2f} vs {jr['mu_rest']:.2f} ± {jr['err_rest']:.2f} → różnica {jr['diff_sigma']:+.2f}σ.")
    mn = nom["combined"]["mu"]
    L += ["", f"Skala: próbka odsłonięcia 05.11 ma μ_nom = {mn:.3f}; zapis bieżący (N_D* = 4436 po NDS 12.11) to 1,50 → μ_c w zapisie bieżącym = μ_c × {1.50 / mn:.3f}. "
          "Wniosek o triggerze opiera się na STOSUNKACH μ_c / μ_nom (niezależnych od N_D*)."]
    L += ["", "Uwaga: klasa MBTS w 428770 to w ~93 % HLT_j20 (też dżetowa); jedyny run bez obciążenia dżetowego to 427929 (osobny etap).", ""]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
