#!/usr/bin/env python3
"""gn_combine.py — gen_norm (16.11.2026): części z gn_root.py (mb_nd, data_428770) + nogi z gn_legs.py → wyniki.

(a) w2(n_trk, lead_pt) = p_dane/p_MC(mb_nd) jak etap1 code/w2d.py (biny NE×PE, komórki MC < 20 → 0):
    w2d.npy (+ w2d_edges.json) z danych REFERENCYJNYCH (pliki 428770 _000021, _000030 = 441 691 zdarzeń jak etap1/t4),
    w2d_alldata.npy z wszystkich plików 428770 (wariant). Kontrola = etap1 code/analyze2.py (tryb 2d) → mb_rates_2d_repro.json,
    c_w (charmtag `weight_factors`) vs ref/mb_rates_2d.json: test akceptacji → acceptance.json (PASS/FAIL).
(b) nogi `direct` (h5, selekcja charmtag) ↔ prawdziwe D (ROOT TruthDmesons) po (tag, evnum, gatunek, η, pT):
    effmap_direct.npz: licznik/mianownik ε(gatunek; pT, |η|, n_trk) — eps_any (≥1 noga), eps_mult (Σ nóg), matched (Donly, stare),
    nieważone i ważone w2; D* = okno D* trybu direct (flaga okna). Tabele: truthD_mb.parquet, legs_mb_w.parquet, dstar_mb_w.parquet,
    events_D_mb.parquet (zdarzenia z prawdziwym D lub nogą: n_trk, lead_pt, w2, w2_all). c_w dla selekcji h5 (cw_h5.json).
    Domknięcie MC na mb_nd (direct vs truth×mapy, także „split-half”) → closure_mb.json (kryteria: GEN_mcstat/PROCEDURA.md).
(c) TruthParticles: gdy w częściach MC jest ev_nch → macierz odpowiedzi (n_ch, lead_gen) → (n_trk, lead_pt): response_nch.npz;
    inaczej response_nch.json = {"skipped": powód}.
  python3 gn_combine.py --mc DIR --data DIR --legs DIR --out DIR [--ref ref/mb_rates_2d.json] [--data-ref-tags 000021,000030]"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
NE = [-1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 1e4]           # = etap1 w2d.py
PE = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 1e4]
MINMC = 20
CH = {"D0D0": (10, 10), "DpD0": (4, 10), "DpDp": (4, 4)}                        # = charmtag.pairs.CHANNELS
SPDT = {421: 10, 411: 4, 413: 1}
# mapy ε: szerzej niż akceptacja par (nogi direct mogą pochodzić z D przy pT < 3 / |η| > 2,2 prawdy)
PT_E = np.array([2.5, 3, 4, 5, 6, 8, 12, 30, 1e9]); ETA_E = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.2, 2.5, 2.7]); NT_E = np.array([0, 10, 20, 30, 45, 70, 1e9])
MATCH_DETA, MATCH_DPT = 0.05, 0.15                                               # dopasowanie noga ↔ prawdziwe D
REF_COUNTS = {"events": 19714175, "data_events": 441691, "dstar_n": 1314, "D0_n": 9448, "Dplus_n": 8671,
              "pairs_D0D0_n": 21, "pairs_DpD0_n": 41, "pairs_DpDp_n": 15}
REF_DATA_EVENTS = 441691
ANCHOR_GAP = 1000            # v2.0: przerwa (wpisy) między kotwicami, powyżej której zaczyna się nowy blok referencji
RATE_RTOL = 2e-4             # v2.0 B1b: zgodność nieznormalizowanych częstości D+ / par D+D+ z ref (te same zdarzenia i kandydaci D+;
                             #   ≤ 10 zdarzeń granicznych z 441 691 → efekt ~1e-5)
REF_EVENTS_TOL = 10          # v2.0 B1b: |N_ref − 441 691| (skan 17.11: kotwice 303 497/303 497, bloki 441 686 bez krawędzi, ~441 690 z krawędziami)
CW_TOL = 0.005                                                                   # |c_w − c_w_ref| (odtworzenie bit-w-bit oczekiwane do ~1e-6)
# §3b (18.11, wariant 78 plików mb_nd; zatwierdzone): pliki ROOT v2.0 _000037/_000070 skasowane z gridu → --mb-exclude.
# Referencje B1a dla podzbioru = REF_COUNTS etap1 (80 plików) − wkład wykluczonych tagów policzony z nóg h5 v1.1 (te same tagi joba);
# dopuszczalna różnica = znana różnica definicji etap1 (ROOT) vs nogi h5 na pełnych 80 plikach: D+ |8671 − 8665| = 6, pary D+D+ 1 (15 = 15; zapas 1).
REF_H5 = {"Dplus_n": 8665, "pairs_DpDp_n": 15}
EXCL_TOL = {"Dplus_n": 6, "pairs_DpDp_n": 1}
N_MB_FILES_REF = 80
EXCL_NSIGMA = 3.0            # B1b przy wykluczeniu: |Δ częstości|/częstość < 3·√(f/(1−f)/n_ref), f = udział wykluczonych plików
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:6.0f}s]", *a, file=sys.stderr, flush=True)


def tag_of(f):
    return os.path.basename(f).rsplit(".", 1)[0]


def load_parts(d, kind):
    fs = sorted(glob.glob(os.path.join(d, "*.npz")))
    if not fs:
        sys.exit(f"brak części {kind} w {d}")
    out = []
    for f in fs:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        if meta["kind"] != kind:
            sys.exit(f"{f}: kind {meta['kind']} ≠ {kind}")
        out.append((tag_of(f), meta, {k: z[k] for k in z.files if k != "meta"}))
    return out


# ----------------------------------------------------------------------------- (a) w2 — etap1 w2d.py / analyze2.py
# ----------------------------------------------------------------------------- v2.0: referencja danych po zdarzeniach (B1 v2.0)
def load_anchors(spec):
    """Kotwice = unikalne event_number z h5 v1.1 data_428770_{000021,000030}_none (zdarzenia z ≥1 kandydatem; 303 497 z 441 691)."""
    import h5py
    A = [np.unique(h5py.File(f, "r")["events"]["event_number"]).astype(np.uint64) for f in spec.split(",") if f]
    return np.unique(np.concatenate(A)) if A else np.zeros(0, np.uint64)


def ref_by_anchors(DA, A, gap=ANCHOR_GAP):
    """Pliki v2.0 (inne dzielenie niż v1.1) → maska referencji: w każdym pliku ciągłe bloki wpisów od pierwszej do ostatniej kotwicy
    (blok kończy się, gdy przerwa między kotwicami > gap). Zakłada, że zdarzenia jednego AOD są w ROOT ciągłe (porządek wejścia)."""
    masks, info, found = [], [], 0
    for tg, _, z in DA:
        e = z["ev_evnum"].astype(np.uint64); hit = np.isin(e, A); idx = np.flatnonzero(hit); m = np.zeros(len(e), bool)
        if len(idx):
            cut = np.flatnonzero(np.diff(idx) > gap)
            los, his = np.r_[idx[0], idx[cut + 1]], np.r_[idx[cut], idx[-1]]
            if los[0] <= gap:                       # blok przy krawędzi pliku → do krawędzi (zdarzenia bez kandydata na brzegu)
                los[0] = 0
            if his[-1] >= len(e) - 1 - gap:
                his[-1] = len(e) - 1
            for lo, hi in zip(los, his):
                m[lo:hi + 1] = True
            info.append({"tag": tg, "anchors": int(len(idx)), "blocks": int(len(cut) + 1), "ref_events": int(m.sum()), "entries": int(len(e))})
        found += int(len(idx)); masks.append(m)
    uniq_found = int(len(np.unique(np.concatenate([z["ev_evnum"][hit_m] for (_, _, z), hit_m in
                                                    zip(DA, [np.isin(z["ev_evnum"].astype(np.uint64), A) for _, _, z in DA])]))))
    return masks, {"anchors": int(len(A)), "anchors_found_rows": found, "anchors_found_unique": uniq_found,
                   "ref_events": int(sum(int(m.sum()) for m in masks)), "files": info, "gap": gap}


def w2_map(d_ntrk, d_lead, m_ntrk, m_lead):
    hd, _, _ = np.histogram2d(d_ntrk, d_lead, [NE, PE]); hm, _, _ = np.histogram2d(m_ntrk, m_lead, [NE, PE])
    pdn = hd / hd.sum(); pmn = hm / hm.sum()
    w = np.divide(pdn, pmn, out=np.zeros_like(pdn), where=hm >= MINMC)
    return w, hd, hm, float(pdn[hm < MINMC].sum())


def w_lookup(W, ntrk, lead):
    i = np.clip(np.digitize(ntrk, NE) - 1, 0, len(NE) - 2); j = np.clip(np.digitize(lead, PE) - 1, 0, len(PE) - 2)
    return W[i, j]


def analyze2(ev, c, W2):
    """etap1 code/analyze2.py w trybie `2d` (te same typy: n_trk int16, lead_pt float32, w float32)."""
    ev = ev.copy()
    ev["w"] = w_lookup(W2, ev.n_trk.astype("int16"), ev.lead_pt.astype("float32")).astype("float32")
    N = len(ev); Wsum = float(ev.w.sum())
    M0 = {10: 1864.84, 4: 1869.66}; BOX = {10: 55.0, 4: 45.0}
    c = c.copy()
    c["inbox"] = (c.m - c.decay_type.map(M0)).abs() < c.decay_type.map(BOX)
    keys = c[["file", "entry"]].drop_duplicates()
    c = c.merge(ev.merge(keys, on=["file", "entry"])[["file", "entry", "w", "n_trk"]], on=["file", "entry"])
    res = {"files": int(ev.file.nunique()), "events": int(N), "sum_w": float(Wsum)}
    res["dstar_true_per_event"] = {"unw": float(ev.n_true_dstar.sum() / N), "w": float((ev.n_true_dstar.astype("float64") * ev.w).sum() / Wsum),
                                   "n": int(ev.n_true_dstar.sum())}
    for dt, n in ((10, "D0"), (4, "Dplus")):
        x = c[c.decay_type == dt]
        res[f"{n}_correct_per_event"] = {"unw": len(x) / N, "w": float(x.w.sum() / Wsum), "n": int(len(x))}
    P = c.merge(c, on=["file", "entry"], suffixes=("1", "2"))
    P = P[(P.decay_type1 > P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.trk01 * 100000 + P.trk11 < P.trk02 * 100000 + P.trk12))]
    t1 = P[["trk01", "trk11", "trk21"]].to_numpy(); t2 = P[["trk02", "trk12", "trk22"]].to_numpy(); sh = np.zeros(len(P), bool)
    for i in range(3):
        for j in range(3):
            sh |= (t1[:, i] == t2[:, j]) & (t1[:, i] >= 0)
    P = P[~sh]
    for (a, b), ch in (((10, 10), "D0D0"), ((10, 4), "DpD0"), ((4, 4), "DpDp")):
        q = P[(P.decay_type1 == a) & (P.decay_type2 == b)]
        qb = q[q.inbox1 & q.inbox2]
        n_eff = float(q.w1.astype(float).sum() ** 2 / (q.w1.astype(float) ** 2).sum()) if len(q) else 0.0
        res[f"pairs_{ch}"] = {"n": int(len(q)), "n_box": int(len(qb)), "per_event_unw": len(q) / N, "per_event_w": float(q.w1.sum() / Wsum),
                              "box_per_event_w": float(qb.w1.sum() / Wsum), "median_ntrk": float(q.n_trk1.median()) if len(q) else None,
                              "n_eff": n_eff}
    return res


def weight_factors(mb):
    """= charmtag.pairs.build.weight_factors (kopia, żeby combine nie wymagał charmtag)."""
    ds_w, ds_u = mb["dstar_true_per_event"]["w"], mb["dstar_true_per_event"]["unw"]
    f = lambda w, u: (w / ds_w) / (u / ds_u) if (ds_w and ds_u and u) else float("nan")  # noqa: E731  (0 par → nan, test lokalny)
    cs = {"D0": f(mb["D0_correct_per_event"]["w"], mb["D0_correct_per_event"]["unw"]),
          "Dplus": f(mb["Dplus_correct_per_event"]["w"], mb["Dplus_correct_per_event"]["unw"])}
    cp = {ch: f(mb[f"pairs_{ch}"]["per_event_w"], mb[f"pairs_{ch}"]["per_event_unw"]) for ch in CH}
    return {**cs, **cp}


# ----------------------------------------------------------------------------- (b) nogi ↔ prawda, mapy
def match_legs(legs, truth, dt, sp_pdg, use_phi=False):
    """Każda noga → najbliższe prawdziwe D tego gatunku w tym samym (tag, evnum): |Δη| < 0,05 i |pT_r/pT_t − 1| < 0,15.
    Zwraca indeks wiersza `truth` (−1 = brak)."""
    T = truth[np.abs(truth.pdg) == sp_pdg][["tag", "evnum", "pt", "eta", "phi"]].reset_index().rename(columns={"index": "tidx"})
    Lg = legs[["tag", "evnum", "pt", "eta"] + (["phi"] if use_phi else [])].reset_index(drop=True).reset_index().rename(columns={"index": "lidx"})
    M = Lg.merge(T, on=["tag", "evnum"], suffixes=("", "_t"))
    deta = np.abs(M.eta - M.eta_t); dpt = np.abs(M.pt / M.pt_t - 1)
    d = deta / MATCH_DETA + dpt / MATCH_DPT
    if use_phi:
        dphi = np.abs((M.phi - M.phi_t + np.pi) % (2 * np.pi) - np.pi); d = d + dphi / MATCH_DETA
        ok = (deta < MATCH_DETA) & (dpt < MATCH_DPT) & (dphi < MATCH_DETA)
    else:
        ok = (deta < MATCH_DETA) & (dpt < MATCH_DPT)
    M = M[ok].assign(d=d[ok]).sort_values("d").drop_duplicates("lidx")
    out = np.full(len(legs), -1, np.int64); out[M.lidx.to_numpy()] = M.tidx.to_numpy()
    return out


def bins(pt, aeta, ntrk):
    i = np.clip(np.digitize(pt, PT_E) - 1, 0, len(PT_E) - 2)
    j = np.clip(np.digitize(aeta, ETA_E) - 1, 0, len(ETA_E) - 2)
    k = np.clip(np.digitize(ntrk, NT_E) - 1, 0, len(NT_E) - 2)
    return i, j, k


def acc(T):
    return (T["mode"] == 1) & (T.pt > PT_E[0]) & (T.eta.abs() < ETA_E[-1])


def build_maps(T, w=None):
    """Liczniki na prawdziwych D w akceptacji mapy (mode == 1): den, num_any, num_mult, num_matched (opcjonalnie wagi)."""
    sh = (len(PT_E) - 1, len(ETA_E) - 1, len(NT_E) - 1)
    out = {}
    a = acc(T).to_numpy()
    for sp, pdg in (("D0", 421), ("Dplus", 411), ("Dstar", 413)):
        m = a & (np.abs(T.pdg.to_numpy()) == pdg)
        i, j, k = bins(T.pt.to_numpy()[m], np.abs(T.eta.to_numpy()[m]), T.ntrk.to_numpy()[m])
        ww = np.ones(m.sum()) if w is None else w[m]
        for name, v in (("den", np.ones(m.sum())), ("num_any", (T.n_legs.to_numpy()[m] > 0).astype(float)),
                        ("num_mult", T.n_legs.to_numpy()[m].astype(float)), ("num_matched", T.matched.to_numpy()[m].astype(float))):
            h = np.zeros(sh); np.add.at(h, (i, j, k), v * ww); out[f"{sp}__{name}"] = h
    return out


def eps_from(maps, sp, num, minden=30):
    """ε z hierarchicznym cofaniem (pT,η,ntrk) → (pT,η) → (pT) → całość (jak truth_combine.EffMap), gdy mianownik < minden."""
    den, nu = maps[f"{sp}__den"], maps[f"{sp}__{num}"]
    e3 = np.where(den >= minden, nu / np.maximum(den, 1e-300), np.nan)
    d2, n2 = den.sum(2), nu.sum(2); e2 = np.where(d2 >= minden, n2 / np.maximum(d2, 1e-300), np.nan)
    d1, n1 = d2.sum(1), n2.sum(1); e1 = np.where(d1 >= minden, n1 / np.maximum(d1, 1e-300), np.nan)
    e0 = nu.sum() / max(den.sum(), 1e-300)
    e = np.where(np.isnan(e3), e2[:, :, None], e3); e = np.where(np.isnan(e), e1[:, None, None], e); e = np.where(np.isnan(e), e0, e)
    return e


def pairs_of_legs(L):
    """Pary nóg jak mc_rates: to samo zdarzenie, (dt1 < dt2) lub (dt1 == dt2 i row1 < row2), bez wspólnych śladów."""
    P = L.reset_index().merge(L.reset_index(), on=["tag", "evnum"], suffixes=("1", "2"))
    P = P[(P.decay_type1 < P.decay_type2) | ((P.decay_type1 == P.decay_type2) & (P.row1 < P.row2))]
    t1 = P[["trk01", "trk11", "trk21"]].to_numpy(); t2 = P[["trk02", "trk12", "trk22"]].to_numpy(); sh = np.zeros(len(P), bool)
    for i in range(3):
        for j in range(3):
            sh |= (t1[:, i] == t2[:, j]) & (t1[:, i] >= 0)
    return P[~sh]


def truth_pairs(T):
    """Pary prawdziwych D0/D+ (mode == 1, akceptacja mapy) w tym samym zdarzeniu; kanał wg CH."""
    d = T[acc(T) & np.abs(T.pdg).isin([421, 411])].reset_index()
    P = d.merge(d, on=["tag", "evnum"], suffixes=("_a", "_b"))
    P = P[P["index_a"] < P["index_b"]]
    da, db = np.where(np.abs(P.pdg_a) == 421, 10, 4), np.where(np.abs(P.pdg_b) == 421, 10, 4)
    P = P.assign(ch=np.where((da == 10) & (db == 10), "D0D0", np.where((da == 4) & (db == 4), "DpDp", "DpD0")))
    return P


def closure(T, legs, Ds, P_legs, maps_for, label):
    """Direct (nogi h5) vs truth×mapy na tych samych zdarzeniach mb_nd, oba ważone w2; stosunki pary/D*.
    maps_for(tag) → mapy (pozwala na split-half: mapa z drugiej połowy plików)."""
    out = {"label": label, "channels": {}}
    w_l = legs.w2.to_numpy()
    ds_dir = float(Ds.w2.sum()); ds_dir_n = len(Ds)
    # truth: ε z mapy w punktach prawdy
    T = T.copy(); T["mu_mult"] = np.nan; T["mu_any"] = np.nan
    for tg, idx in T.groupby("tag").groups.items():
        M = maps_for(tg)
        for sp, pdg in (("D0", 421), ("Dplus", 411), ("Dstar", 413)):
            sel = idx[np.abs(T.loc[idx, "pdg"].to_numpy()) == pdg]
            if not len(sel):
                continue
            i, j, k = bins(T.loc[sel, "pt"].to_numpy(), np.abs(T.loc[sel, "eta"].to_numpy()), T.loc[sel, "ntrk"].to_numpy())
            T.loc[sel, "mu_mult"] = M[sp]["num_mult"][i, j, k]; T.loc[sel, "mu_any"] = M[sp]["num_any"][i, j, k]
    a = acc(T)
    TD = T[a & (np.abs(T.pdg) == 413)]
    ds_tr = {v: float((TD[f"mu_{v}"] * TD.w2).sum()) for v in ("mult", "any")}
    TP = truth_pairs(T[a])
    for ch in CH:
        q = P_legs[P_legs.ch == ch]
        wq = q.w21.to_numpy()                                         # waga zdarzenia (noga 1 = noga 2, to samo zdarzenie)
        n_eff = float(wq.sum() ** 2 / (wq ** 2).sum()) if len(q) else 0.0
        r_dir = float(wq.sum()) / ds_dir if ds_dir else np.nan
        tp = TP[TP.ch == ch]
        res = {"direct_n": int(len(q)), "direct_n_eff": n_eff, "direct_pair_per_dstar_w": r_dir,
               "direct_rel_stat": float(1 / np.sqrt(n_eff)) if n_eff else None, "truth_pairs_n": int(len(tp))}
        for v in ("mult", "any"):
            pr = float((tp[f"mu_{v}_a"] * tp[f"mu_{v}_b"] * tp.w2_a).sum())
            r_tr = pr / ds_tr[v] if ds_tr[v] else np.nan
            res[f"truth_pair_per_dstar_w_{v}"] = r_tr
            res[f"ratio_truth_over_direct_{v}"] = r_tr / r_dir if r_dir else None
            res[f"pull_{v}"] = ((r_tr / r_dir - 1) / res["direct_rel_stat"]) if (r_dir and n_eff) else None
        out["channels"][ch] = res
    wd = Ds.w2.to_numpy(); neff_ds = float(wd.sum() ** 2 / (wd ** 2).sum()) if len(Ds) and (wd ** 2).sum() > 0 else 0.0
    out["dstar"] = {"direct_n": ds_dir_n, "direct_w": ds_dir, "direct_n_eff": neff_ds, "truth_w_mult": ds_tr["mult"], "truth_w_any": ds_tr["any"]}
    for ch in CH:                                                   # σ testu: para ⊕ D* (oba z ważonego direct)
        v = out["channels"][ch]
        v["sigma_test"] = float(np.sqrt(1 / v["direct_n_eff"] + 1 / neff_ds)) if (v["direct_n_eff"] and neff_ds) else None
    # pojedyncze (D/D*) — test transferu map na singlach
    for sp, dt, pdg in (("D0", 10, 421), ("Dplus", 4, 411)):
        wl = w_l[legs.decay_type.to_numpy() == dt]; lw = float(wl.sum())
        ne = float(wl.sum() ** 2 / (wl ** 2).sum()) if (wl ** 2).sum() > 0 else 0.0
        out.setdefault("singles", {}).setdefault(sp, {})["direct_n_eff"] = ne
        out["singles"][sp]["sigma_test"] = float(np.sqrt(1 / ne + 1 / neff_ds)) if (ne and neff_ds) else None
        TS = T[a & (np.abs(T.pdg) == pdg)]
        for v in ("mult", "any"):
            out.setdefault("singles", {}).setdefault(sp, {})[f"ratio_truth_over_direct_{v}"] = \
                (float((TS[f"mu_{v}"] * TS.w2).sum()) / ds_tr[v]) / (lw / ds_dir) if (lw and ds_tr[v]) else None
    return out


def b1_exclude(acceptance, rep, ref, legs, Ds, ev, excl):
    """§3b: wariant z wykluczonymi plikami mb_nd. Usuwa nogi/D* wykluczonych tagów (zwraca przycięte legs, Ds) i podmienia w
    acceptance['checks'] testy zależne od pełnych 80 plików: count_events → informacyjnie (implikowane zdarzenia wykluczonych plików
    vs rozrzut plików), count_Dplus_n / count_pairs_DpDp_n → ref = etap1 − wkład wykluczonych tagów z nóg h5 (tol EXCL_TOL),
    B1b_* (częstości MC) → tolerancja statystyczna EXCL_NSIGMA·√(f/(1−f)/n_ref). Progi ustalone 18.11 przed wynikiem."""
    ch = acceptance["checks"]; info = acceptance.setdefault("info_nonblocking", {})
    exm = legs.tag.isin(excl).to_numpy(); exd = Ds.tag.isin(excl).to_numpy()
    Lx = legs[exm]
    x_dp = int((Lx.decay_type == 4).sum())
    Px = pairs_of_legs(Lx) if len(Lx) else Lx.iloc[:0]
    x_pp = int(((Px.decay_type1 == 4) & (Px.decay_type2 == 4)).sum()) if len(Px) else 0
    f = len(excl) / N_MB_FILES_REF
    per_file = ev.groupby("file").size()
    implied = REF_COUNTS["events"] - int(len(ev))
    ex = {"excluded_tags": sorted(excl), "file_fraction": f, "legs_excluded": int(exm.sum()), "legs_excluded_frac": float(exm.mean()) if len(exm) else 0.0,
          "dstar_excluded": int(exd.sum()), "Dplus_legs_excluded": x_dp, "DpDp_pairs_excluded": x_pp,
          "events_78": int(len(ev)), "events_implied_excluded": implied, "events_implied_per_excluded_file": implied / max(len(excl), 1),
          "events_per_file_min_max": [int(per_file.min()), int(per_file.max())]}
    ex["events_implied_in_file_range"] = bool(per_file.min() <= ex["events_implied_per_excluded_file"] <= per_file.max())
    info["exclusion_78"] = ex
    if "count_events" in ch:
        info["count_events"] = {**ch.pop("count_events"), "blocking": False, "note": "80 plików; przy wykluczeniu zob. exclusion_78"}
    for k, x in (("Dplus_n", x_dp), ("pairs_DpDp_n", x_pp)):
        old = ch.pop(f"count_{k}", None)
        if old is not None:
            info[f"count_{k}"] = {**old, "blocking": False}
        g = rep["Dplus_correct_per_event"]["n"] if k == "Dplus_n" else rep["pairs_DpDp"]["n"]
        r = REF_COUNTS[k] - x
        ch[f"B1a78_{k}"] = {"got": g, "ref": r, "ref_how": f"etap1 {REF_COUNTS[k]} − nogi h5 wykluczonych tagów {x}", "tol": EXCL_TOL[k],
                            "pass": abs(g - r) <= EXCL_TOL[k]}
    for blk, n_ref in (("Dplus_correct_per_event", REF_COUNTS["Dplus_n"]), ("pairs_DpDp", REF_COUNTS["pairs_DpDp_n"])):
        tol = EXCL_NSIGMA * float(np.sqrt(f / (1 - f) / n_ref))
        for fld in (("unw", "w") if blk.startswith("Dplus") else ("per_event_unw", "per_event_w")):
            old = ch.pop(f"B1b_{blk}_{fld}", None)
            if old is not None:
                info[f"B1b_{blk}_{fld}_strict"] = {**old, "blocking": False}
            g, r = rep[blk][fld], ref[blk][fld]
            rel = abs(g - r) / abs(r) if r else abs(g)
            ch[f"B1b78_{blk}_{fld}"] = {"got": g, "ref": r, "rel_diff": rel, "tol": tol, "pass": bool(rel < tol)}
    return legs[~exm].reset_index(drop=True), Ds[~exd].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc", required=True); ap.add_argument("--data", required=True); ap.add_argument("--legs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref", default=os.path.join(HERE, "ref", "mb_rates_2d.json"))
    ap.add_argument("--data-ref-tags", default="000021,000030")
    ap.add_argument("--fmt", default="v11", choices=["v11", "v20"], help="v20: referencja danych po kotwicach (--ref-anchors), B1 v2.0")
    ap.add_argument("--ref-anchors", default="", help="v20: h5 v1.1 data_428770_{000021,000030}_none.h5 (po przecinku)")
    ap.add_argument("--ref-counts", default="strict", choices=["strict", "none"], help="none = test lokalny (atrapy): tylko c_w vs --ref")
    ap.add_argument("--mb-exclude", default="", help="§3b: tagi plików mb_nd nieobecnych w ROOT (np. 000037,000070) — nogi h5 tych tagów usuwane")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    MC = load_parts(a.mc, "mc"); DA = load_parts(a.data, "data")
    say(f"części: MC {len(MC)}, dane {len(DA)}")
    EXCL = sorted({t.strip() for t in a.mb_exclude.split(",") if t.strip()})
    if EXCL:
        have = {tg[-6:] for tg, _, _ in MC}
        if have & set(EXCL):
            sys.exit(f"--mb-exclude {EXCL}: tagi są wśród części MC ({sorted(have & set(EXCL))}) — wykluczenie niespójne")
        if len(have) + len(EXCL) != N_MB_FILES_REF and a.ref_counts == "strict":
            sys.exit(f"--mb-exclude: części MC {len(have)} + wykluczone {len(EXCL)} ≠ {N_MB_FILES_REF}")
        say(f"§3b: wykluczone pliki mb_nd {EXCL}")
    # ---- zdarzenia MC
    evs, cs, ts = [], [], []
    TAGS = [tg[-6:] for tg, _, _ in MC]                               # 6 cyfr joba = tag h5 (mb_nd_<job>_none.h5) w nogach
    for k, (tg, meta, z) in enumerate(MC):
        e = pd.DataFrame({c[3:]: z[c] for c in z if c.startswith("ev_")}); e["file"] = np.int16(k); evs.append(e)
        c = pd.DataFrame({c[2:]: z[c] for c in z if c.startswith("c_")}); c["file"] = np.int16(k); cs.append(c)
        t = pd.DataFrame({c[2:]: z[c] for c in z if c.startswith("t_")}); t["file"] = np.int16(k); ts.append(t)
    ev = pd.concat(evs, ignore_index=True); c = pd.concat(cs, ignore_index=True); T = pd.concat(ts, ignore_index=True)
    ev["tag"] = pd.Categorical.from_codes(ev.file.to_numpy(), categories=TAGS)    # kategorie: bez 20 M napisów w pamięci
    T["tag"] = np.asarray(TAGS, dtype=object)[T.file.to_numpy()]
    has_tp = all(m["has_truth_particles"] for _, m, _ in MC)
    say(f"MC: {len(ev)} zdarzeń, {len(c)} poprawnych kand. etap1, {len(T)} prawdziwych D; TruthParticles: {has_tp}")
    ref_tags = [t for t in a.data_ref_tags.split(",") if t]
    dall = [z for _, _, z in DA]
    anch = None
    if a.fmt == "v20":
        if not a.ref_anchors:
            sys.exit("--fmt v20 wymaga --ref-anchors (h5 v1.1 _000021/_000030)")
        A = load_anchors(a.ref_anchors); masks, anch = ref_by_anchors(DA, A)
        dref = [{k: v[m] for k, v in z.items() if k.startswith("ev_")} for (_, _, z), m in zip(DA, masks) if m.any()]
        json.dump(anch, open(os.path.join(a.out, "ref_anchors.json"), "w"), indent=1)
        say("v2.0 referencja po kotwicach:", {k: anch[k] for k in ("anchors", "anchors_found_unique", "ref_events")})
        ref_tags = [f"anchors:{len(A)}"]
    else:
        dref = [z for tg, _, z in DA if any(tg.endswith(r) for r in ref_tags)]
    if not dref:
        sys.exit(f"brak danych referencyjnych {ref_tags} (referencja etap1) wśród {[tg for tg, _, _ in DA]}")
    cat = lambda L, k: np.concatenate([z[k] for z in L])  # noqa: E731
    W, hd, hm, lost = w2_map(cat(dref, "ev_n_trk"), cat(dref, "ev_lead_pt"), ev.n_trk.to_numpy(), ev.lead_pt.to_numpy())
    Wall, hd_all, _, lost_all = w2_map(cat(dall, "ev_n_trk"), cat(dall, "ev_lead_pt"), ev.n_trk.to_numpy(), ev.lead_pt.to_numpy())
    np.save(os.path.join(a.out, "w2d.npy"), W); np.save(os.path.join(a.out, "w2d_alldata.npy"), Wall)
    json.dump({"NE": NE, "PE": PE}, open(os.path.join(a.out, "w2d_edges.json"), "w"))
    np.savez(os.path.join(a.out, "hist2d.npz"), NE=NE, PE=PE, data_ref=hd, data_all=hd_all, mc=hm)
    # w2d.py: walidacja na tempie D* (float64, jak w2d.py)
    val = {"mc_files": len(MC), "mc_events": int(len(ev)), "data_events": int(hd.sum()), "data_events_all": int(hd_all.sum()),
           "data_ref_tags": ref_tags, "data_mass_in_empty_mc_cells": lost, "data_mass_in_empty_mc_cells_all": lost_all}
    for col, WW in (("unw", None), ("w2", W), ("w2_alldata", Wall)):
        ww = np.ones(len(ev)) if WW is None else w_lookup(WW, ev.n_trk.to_numpy(), ev.lead_pt.to_numpy())
        rate = float((ev.n_true_dstar * ww).sum() / ww.sum())
        val[col] = {"dstar_per_event": rate, "data_over_mc": 7.65e-4 / rate, "n_dstar": int(ev.n_true_dstar.sum())}
    json.dump(val, open(os.path.join(a.out, "w2d_validation_repro.json"), "w"), indent=1)
    say("w2: gotowe", {k: val[k] for k in ("mc_events", "data_events", "data_events_all")})
    # ---- odtworzenie mb_rates_2d.json (analyze2 2d) i c_w
    rep = analyze2(ev[["file", "entry", "n_trk", "lead_pt", "n_true_dstar"]], c, W)
    rep_all = analyze2(ev[["file", "entry", "n_trk", "lead_pt", "n_true_dstar"]], c, Wall)
    json.dump(rep, open(os.path.join(a.out, "mb_rates_2d_repro.json"), "w"), indent=1)
    json.dump(rep_all, open(os.path.join(a.out, "mb_rates_2d_alldata.json"), "w"), indent=1)
    ref = json.load(open(a.ref))
    cw, cw_ref, cw_all = weight_factors(rep), weight_factors(ref), weight_factors(rep_all)
    got_counts = {"events": rep["events"], "data_events": int(hd.sum()), "dstar_n": rep["dstar_true_per_event"]["n"],
                  "D0_n": rep["D0_correct_per_event"]["n"], "Dplus_n": rep["Dplus_correct_per_event"]["n"],
                  **{f"pairs_{ch}_n": rep[f"pairs_{ch}"]["n"] for ch in CH}}
    checks = {f"c_w_{k}": {"got": cw[k], "ref": cw_ref[k], "abs_diff": abs(cw[k] - cw_ref[k]),
                           "pass": bool(abs(cw[k] - cw_ref[k]) < CW_TOL or (np.isnan(cw[k]) and np.isnan(cw_ref[k])))} for k in cw}
    if a.ref_counts == "strict":
        checks.update({f"count_{k}": {"got": got_counts[k], "ref": REF_COUNTS[k], "pass": got_counts[k] == REF_COUNTS[k]} for k in REF_COUNTS})
    if a.fmt == "v20":
        # B1 v2.0 (PROCEDURA §3, zmiana 17.11): kandydaci D0/D* w v2.0 ≠ v1.1 (okno D0, D*: inna selekcja) → c_w D0/D*/par z D0 i
        # normalizacja c_w (przez D*) tylko informacyjnie (B1c). Blokujące: B1a (MC: zdarzenia, D+, pary D+D+, klucze nóg — niżej),
        # B1b (dane: kotwice, liczba zdarzeń ref, nieznormalizowane częstości D+ i D+D+ = ref → w2d.npy = mapa nominału).
        info = {k: checks.pop(k) for k in list(checks) if k.startswith("c_w_") or k in ("count_D0_n", "count_dstar_n", "count_pairs_D0D0_n", "count_pairs_DpD0_n")}
        for k in info.values():
            k["blocking"] = False
        rel = lambda g, r: abs(g - r) / abs(r) if r else abs(g)  # noqa: E731
        for blk in ("Dplus_correct_per_event", "pairs_DpDp"):
            for f in (("unw", "w") if blk.startswith("Dplus") else ("per_event_unw", "per_event_w")):
                g, r = rep[blk][f], ref[blk][f]
                checks[f"B1b_{blk}_{f}"] = {"got": g, "ref": r, "rel_diff": rel(g, r), "pass": bool(rel(g, r) < RATE_RTOL)}
        checks["B1b_anchors_found"] = {"got": anch["anchors_found_unique"], "ref": anch["anchors"], "pass": anch["anchors_found_unique"] == anch["anchors"]}
        checks["B1b_ref_events"] = {"got": anch["ref_events"], "ref": REF_DATA_EVENTS, "tol": REF_EVENTS_TOL,
                                    "pass": abs(anch["ref_events"] - REF_DATA_EVENTS) <= REF_EVENTS_TOL}
        checks.pop("count_data_events", None)
    else:
        info = {}
    acc_ok = all(v["pass"] for v in checks.values())
    acceptance = {"test": "odtworzenie c_w z mb_rates_2d.json (etap1, 09-26) — wymagane przed użyciem w2d.npy", "tolerance_cw_abs": CW_TOL,
                  "verdict": "PASS" if acc_ok else "FAIL", "checks": checks, "c_w_alldata": cw_all, "fmt": a.fmt, "info_nonblocking": info,
                  "expected_cw": "D0 0,935 / D+ 0,793 / D0D0 1,374 / DpD0 1,201 / DpDp 0,981"}
    json.dump(acceptance, open(os.path.join(a.out, "acceptance.json"), "w"), indent=1, default=float)
    say("akceptacja (a):", acceptance["verdict"], {k: round(v["got"], 4) if isinstance(v["got"], float) else v["got"] for k, v in checks.items()})
    # ---- (b) nogi direct + prawda
    legs = pd.read_parquet(os.path.join(a.legs, "legs_mb.parquet")); Ds = pd.read_parquet(os.path.join(a.legs, "dstar_mb.parquet"))
    lmeta = json.load(open(os.path.join(a.legs, "legs_meta.json")))
    for X in (legs, Ds):
        X["evnum"] = X.evnum.astype(np.uint64); X["tag"] = X.tag.astype(str)
    if EXCL:   # §3b: przed złączeniem — nogi wykluczonych tagów nie mogą liczyć się jako legs_without_event ani wchodzić do direct/closure
        legs, Ds = b1_exclude(acceptance, rep, ref, legs, Ds, ev, EXCL)
        say("§3b wykluczenie:", acceptance["info_nonblocking"]["exclusion_78"])
    ev["evnum"] = ev.evnum.astype(np.uint64)
    ev["w2"] = w_lookup(W, ev.n_trk.to_numpy(), ev.lead_pt.to_numpy()); ev["w2_all"] = w_lookup(Wall, ev.n_trk.to_numpy(), ev.lead_pt.to_numpy())
    dup = int(ev.duplicated(["file", "evnum"]).sum())
    EVK = ev[["tag", "evnum", "n_trk", "lead_pt", "w2", "w2_all"]].drop_duplicates(["tag", "evnum"])
    EVK["tag"] = EVK.tag.astype(str)
    legs = legs.merge(EVK, on=["tag", "evnum"], how="left"); Ds = Ds.merge(EVK, on=["tag", "evnum"], how="left")
    T["evnum"] = T.evnum.astype(np.uint64)
    T = T.merge(EVK[["tag", "evnum", "lead_pt", "w2", "w2_all"]], on=["tag", "evnum"], how="left")
    diag = {"events_key_duplicates": dup, "legs_without_event": int(legs.w2.isna().sum()), "dstar_without_event": int(Ds.w2.isna().sum()),
            "legs_meta_rates": lmeta.get("mc_rates"), "legs_ref_match": lmeta.get("ref_match")}
    if a.fmt == "v20":   # B1a: nogi (h5 v1.1, tag joba v1.1) muszą trafić w zdarzenia ROOT v2.0 o tym samym (tag, eventNumber)
        for nm, val in (("B1a_legs_without_event", diag["legs_without_event"]), ("B1a_dstar_without_event", diag["dstar_without_event"]),
                        ("B1a_event_key_duplicates", dup)):
            acceptance["checks"][nm] = {"got": val, "ref": 0, "pass": val == 0}
        acceptance["verdict"] = "PASS" if all(v["pass"] for v in acceptance["checks"].values()) else "FAIL"
        json.dump(acceptance, open(os.path.join(a.out, "acceptance.json"), "w"), indent=1, default=float)
        say("akceptacja (a, po B1a):", acceptance["verdict"])
    for X in (legs, Ds):
        X["w2"] = X.w2.fillna(0.0); X["w2_all"] = X.w2_all.fillna(0.0)
    # dopasowanie nóg do prawdy
    legs["tidx"] = -1
    for dt, pdg in ((10, 421), (4, 411)):
        m = (legs.decay_type == dt).to_numpy()
        legs.loc[m, "tidx"] = match_legs(legs[m], T, dt, pdg)
    Ds["tidx"] = match_legs(Ds, T, 1, 413, use_phi=True)
    T["n_legs"] = 0
    for X in (legs, Ds):
        cnt = X.tidx[X.tidx >= 0].value_counts()
        T.loc[cnt.index, "n_legs"] += cnt.to_numpy()
    accT = acc(T)
    for nm, X, pdgs in (("D0", legs[legs.decay_type == 10], [421]), ("Dplus", legs[legs.decay_type == 4], [411]), ("Dstar", Ds, [413])):
        ti = X.tidx.to_numpy(); okm = ti >= 0
        tt = T.iloc[ti[okm]]
        diag[f"match_{nm}"] = {"n_legs": int(len(X)), "unmatched": int((~okm).sum()), "unmatched_frac": float((~okm).mean()) if len(X) else None,
                               "matched_offmode": int((tt["mode"] != 1).sum()), "matched_outside_map_acc": int((~accT.iloc[ti[okm]]).sum()),
                               "truth_with_2plus_legs": int(((T.n_legs >= 2) & np.isin(np.abs(T.pdg), pdgs)).sum())}
    say("dopasowanie:", {k: v for k, v in diag.items() if k.startswith("match_")})
    # c_w dla selekcji h5 (direct): (X/D*)_w / (X/D*)_unw
    P = pairs_of_legs(legs)
    P["ch"] = np.where((P.decay_type1 == 10) & (P.decay_type2 == 10), "D0D0", np.where((P.decay_type1 == 4) & (P.decay_type2 == 4), "DpDp", "DpD0"))
    dsw, dsu = float(Ds.w2.sum()), len(Ds)
    cwh = {"n_dstar": dsu, "pairs_n": {ch: int((P.ch == ch).sum()) for ch in CH}}
    for nm, dt in (("D0", 10), ("Dplus", 4)):
        x = legs[legs.decay_type == dt]
        cwh[f"c_w_{nm}"] = (x.w2.sum() / dsw) / (len(x) / dsu) if dsu and len(x) else None
    for ch in CH:
        q = P[P.ch == ch]; wq = q.w21.to_numpy()
        cwh[f"c_w_{ch}"] = (wq.sum() / dsw) / (len(q) / dsu) if dsu and len(q) else None
        cwh[f"n_eff_{ch}"] = float(wq.sum() ** 2 / (wq ** 2).sum()) if len(q) else 0.0
    cwh["pairs_equal_legs_meta"] = (cwh["pairs_n"] == (lmeta.get("mc_rates") or {}).get("n_pairs"))
    json.dump(cwh, open(os.path.join(a.out, "cw_h5.json"), "w"), indent=1, default=float)
    say("c_w (selekcja h5 direct, w2 ref):", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in cwh.items()})
    # mapy
    maps = build_maps(T); maps_w = build_maps(T, T.w2.fillna(0).to_numpy())
    np.savez(os.path.join(a.out, "effmap_direct.npz"), pt_edges=PT_E, eta_edges=ETA_E, ntrk_edges=NT_E, **maps,
             **{f"w2__{k}": v for k, v in maps_w.items()},
             **{f"eps_{v}__{sp}": eps_from(maps, sp, f"num_{v}") for sp in ("D0", "Dplus", "Dstar") for v in ("any", "mult", "matched")})
    tot = {sp: {v: float(maps[f"{sp}__{v}"].sum()) for v in ("den", "num_any", "num_mult", "num_matched")} for sp in ("D0", "Dplus", "Dstar")}
    say("mapy ε (sumy):", tot)
    # tabele
    T.to_parquet(os.path.join(a.out, "truthD_mb.parquet"), index=False)
    legs.to_parquet(os.path.join(a.out, "legs_mb_w.parquet"), index=False); Ds.to_parquet(os.path.join(a.out, "dstar_mb_w.parquet"), index=False)
    keyset = pd.concat([T[["tag", "evnum"]], legs[["tag", "evnum"]], Ds[["tag", "evnum"]]]).drop_duplicates()
    EVK.merge(keyset, on=["tag", "evnum"]).to_parquet(os.path.join(a.out, "events_D_mb.parquet"), index=False)
    # domknięcie MC (b): mapy z całości i split-half (parzyste/nieparzyste joby)
    def mk(Tsub):
        Mm = build_maps(Tsub)
        return {sp: {"num_mult": eps_from(Mm, sp, "num_mult"), "num_any": eps_from(Mm, sp, "num_any")} for sp in ("D0", "Dplus", "Dstar")}
    full = mk(T)
    par = T.tag.str[-1].astype(int) % 2
    half = {0: mk(T[par == 1]), 1: mk(T[par == 0])}                 # mapa z DRUGIEJ połowy
    P2 = P
    clo = {"self": closure(T, legs, Ds, P2, lambda tg: full, "mapy z całego mb_nd (tautologia dla singli; test faktoryzacji par)"),
           "split_half": closure(T, legs, Ds, P2, lambda tg: half[int(tg[-1]) % 2], "mapy z drugiej połowy jobów (parzyste↔nieparzyste)")}
    json.dump({"diag": diag, **clo}, open(os.path.join(a.out, "closure_mb.json"), "w"), indent=1, default=float)
    say("domknięcie (split-half):", {ch: round(v["ratio_truth_over_direct_mult"] or np.nan, 3) for ch, v in clo["split_half"]["channels"].items()})
    # (c) odpowiedź n_ch → n_trk
    if has_tp and "nch" in ev:
        NCH = np.array([0, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 1e4])
        R4, _ = np.histogramdd(np.c_[ev.nch, ev.lead_gen, ev.n_trk, ev.lead_pt], [NCH, PE, NE, PE])
        np.savez(os.path.join(a.out, "response_nch.npz"), nch_edges=NCH, lead_edges=PE, ntrk_edges=NE, R=R4,
                 ratio_mean=float(ev.n_trk.sum() / max(ev.nch.sum(), 1)))
        resp = {"done": True, "ntrk_over_nch_mean": float(ev.n_trk.sum() / max(ev.nch.sum(), 1))}
    else:
        resp = {"skipped": "ntuple mb_nd v1.1 nie mają TruthParticles (tylko TruthDmesons) — n_ch→n_trk niedostępne; "
                           "zostaje proxy n_trk = ε_trk·n_ch (0,80/0,85/0,90)",
                "truth_particles_info_first_file": MC[0][1].get("truth_particles")}
    json.dump(resp, open(os.path.join(a.out, "response_nch.json"), "w"), indent=1, default=str)
    # podsumowanie
    L = ["# gen_norm — wyniki (MC + rozkład zdarzeń danych; BEZ par danych, bez μ)", "",
         f"MC mb_nd: {len(MC)} plików, {len(ev)} zdarzeń; dane 428770: {len(DA)} plików (ref {ref_tags}: {int(hd.sum())} zdarzeń, wszystkie {int(hd_all.sum())}).", "",
         f"## (a) w2 — test akceptacji: **{acceptance['verdict']}** (tol |Δc_w| < {CW_TOL})", "", "| | odtworzone | ref | Δ | wszystkie dane |", "|---|---|---|---|---|"]
    for k in cw:
        L.append(f"| c_w {k} | {cw[k]:.4f} | {cw_ref[k]:.4f} | {cw[k] - cw_ref[k]:+.4f} | {cw_all[k]:.4f} |")
    if a.ref_counts == "strict":
        allc = {**acceptance["checks"], **acceptance.get("info_nonblocking", {})}
        L += ["", "Liczniki: " + ", ".join(f"{k} {allc['count_' + k]['got']} ({'OK' if allc['count_' + k]['pass'] else '≠ ' + str(REF_COUNTS[k])})"
                                         for k in REF_COUNTS if 'count_' + k in allc)]
        if a.fmt == "v20":
            L += ["", "B1 v2.0 (blokujące): " + ", ".join(f"{k} {'OK' if v['pass'] else 'FAIL'} ({v['got']})" for k, v in acceptance["checks"].items()
                                                       if k.startswith("B1")), "",
                  "D0/D* i c_w (v2.0 ≠ v1.1 w kandydatach D0/D*) — tylko informacyjnie (B1c)."]
    L += ["", "## (b) nogi direct ↔ prawda", "", f"Nogi (gn_legs): {lmeta.get('mc_rates')}; zgodność z nominałem: {lmeta.get('ref_match')}.",
          "Dopasowanie: " + "; ".join(f"{k[6:]} {v['n_legs']} nóg, bez prawdy {v['unmatched']}, poza mode {v['matched_offmode']}, poza akc. mapy {v['matched_outside_map_acc']}"
                                      for k, v in diag.items() if k.startswith("match_")), "",
          "c_w selekcji h5 (w2 ref): " + ", ".join(f"{k[4:]} {v:.3f}" for k, v in cwh.items() if k.startswith("c_w_") and isinstance(v, float)), "",
          "| domknięcie (split-half) | direct n (n_eff) | truth/direct (mult) | truth/direct (any) | pull |", "|---|---|---|---|---|"]
    for ch, v in clo["split_half"]["channels"].items():
        f = lambda x: "—" if x is None or not np.isfinite(x) else f"{x:.3f}"  # noqa: E731
        L.append(f"| {ch} | {v['direct_n']} ({v['direct_n_eff']:.1f}) | {f(v['ratio_truth_over_direct_mult'])} | {f(v['ratio_truth_over_direct_any'])} | {f(v['pull_mult'])} |")
    L += ["", f"## (c) n_ch → n_trk: {'policzone' if resp.get('done') else 'POMINIĘTE — ' + resp['skipped']}"]
    open(os.path.join(a.out, "summary.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    say("koniec")


if __name__ == "__main__":
    main()
