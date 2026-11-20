#!/usr/bin/env python3
"""UBCOMB (19.11.2026, projekt A): wspólne μ par D–D z próbek BEZ obciążenia triggerem dżetowym.

Próbki (osobne N_D*, c_w, tło fałszywych nóg i SB×SB; wspólne μ = minimum sumy profili NLL wszystkich (próbka, kanał), jak `ub_jk.replica_fit`):
  427929 — cały run w klasie MB (GRL ∧ HLT_mb_mbts_L1MBTS_1; pary z UB427), N_D* metodą NDS;
  428770 — podzbiór zdarzeń S (maska z bitów trigPassed po (run, event)); normalizacja jak UBTRIG: n_pair, n_real × r_S (udział N_D* w S z fitu Δm
           pass/fail), dataleg/SB×SB tylko z S, nev_legs/ndata × udział zdarzeń; N_D* całości = NDS 4436 (μ ∝ 1/N_D*: n_pair × N_NDS/N_biegu).
Warianty 428770 (S):
  strict  — zdarzenia, w których przeszedł łańcuch bez obciążenia: noalg_mb_L1MBTS_1/2, mb_mbts_L1MBTS_1, mb_sptrk_L1MBTS_1, noalg_L1RD0_FILLED
            (niezależnie od współodpalenia J12 — zdarzenie i tak byłoby zapisane);
  nonj12  — klasy UBTRIG MBTS_noj20 ⊕ other (jak w specyfikacji; UWAGA: zawiera sptrk_pt4/6/8 = wymóg toru pT > 4–8 GeV, AFP i zdarzenia bez
            łańcucha z listy → nie jest „czyste”);
  nonj12_noafp — nonj12 bez zdarzeń z łańcuchem AFP.
Wyniki: każdy wariant osobno dla 428770(S), 427929 sam, oraz wspólnie 427929 + 428770(S). Skala μ = zapis bieżący (N_D* NDS).
Systematyki (raportowane, nie profilowane): N_D* niezależne per próbka (stat. NDS + δr/r podzbioru), czynnik korelacji zdarzeń f_jk per próbka,
wspólne normalizacje MC (budżet przejściowy VER_budget: mc_stat, tune — wspólne dla obu runów), model wag (osobno per run).
  ubcomb_run.py --p428 PAIRS428 --legs428 legs_data.parquet --trg428 'trg_*.npz' --bt428 bt_data_428770.npz --bt-mc bt_fd_nd.npz
                --p427 PAIRS427 --ndstar427 ndstar_427929.json [--jk427 ubjk.json] [--variants strict,nonj12,nonj12_noafp] --out DIR"""
import argparse
import glob
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
for p in (os.path.join(HERE, "..", "ub_jk"),):
    if os.path.isdir(p):
        sys.path.insert(0, p)

N_NDS_428770 = (4436.390108391369, 670.8561890721969)    # NDS rerun 12.11 (zapis μ = 1,50)
F_JK_428770 = 1.18                                       # jk_ev200 (ubjk.md) — czynnik korelacji zdarzeń
# budżet przejściowy (VER_budget/RAPORT.md, 16.11) przy μ = 1,50 → względne
SYST_COMMON_REL = {"mc_stat (n_eff)": 0.46 / 1.50, "tune Pythii": 0.17 / 1.50}   # wspólne (te same MC i normalizacja [pary/D*]_MC)
SYST_WAG_REL = 0.24 / 1.50                                                        # model wag (c_w): osobno per run, przyjęte to samo względne
STRICT = ["HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1", "HLT_mb_sptrk_L1MBTS_1", "HLT_noalg_L1RD0_FILLED"]
AFP = ["HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C"]
WP90 = 0.38729336857795715


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ----------------------------------------------------------------------------- maski podzbiorów 428770
def classes_from_tp(tp):
    """Kopia ub_trig/ubtrig_run.classes_from_tp (19.11): 0 J12only, 1 MBTS_j20, 2 MBTS_noj20, 3 other."""
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


def masks_from_tp(tp, n):
    """tp: {chain: bool} → {wariant: bool mask}. nonj12 = klasy UBTRIG MBTS_noj20 (2) ⊕ other (3)."""
    g = lambda c: np.asarray(tp[c], bool) if c in tp else np.zeros(n, bool)  # noqa: E731
    k = classes_from_tp(tp)
    strict = np.zeros(n, bool)
    for c in STRICT:
        strict |= g(c)
    afp = np.zeros(n, bool)
    for c in AFP:
        afp |= g(c)
    nonj12 = (k == 2) | (k == 3)
    return {"strict": strict, "nonj12": nonj12, "nonj12_noafp": nonj12 & ~afp}


def load_trg_masks(paths):
    import btrig_common as C
    run, ev, tp = [], [], {c: [] for c in C.TRIG_CHAINS}
    for p in paths:
        z = np.load(p, allow_pickle=False)
        run.append(z["run"].astype(np.int64)); ev.append(z["event"].astype(np.int64))
        for c in tp:
            tp[c].append((z[f"tp_{c}"] == 1) if f"tp_{c}" in z.files else np.zeros(len(z["run"]), bool))
    run, ev = np.concatenate(run), np.concatenate(ev)
    tp = {c: np.concatenate(v) for c, v in tp.items()}
    key, first = np.unique(C.ev_key(run, ev), return_index=True)
    return key, masks_from_tp({c: v[first] for c, v in tp.items()}, len(key))


def dstar_fraction(bt, btmc, variant):
    """r_S = udział N_D* (D⁰←D*, WP90, okno 55 MeV, pT 3–30) w podzbiorze S, fit Δm pass/fail; maska z bitów e_tp_* w bt."""
    import btrig_common as C
    import dstar_refit_aa as DR                                    # jawna kopia AA_merge (PF/src na lxplus ma starą fit_passfail)
    z = np.load(bt, allow_pickle=False)
    q = (np.abs(z["m"] - C.M0_D0) < C.NOM_WIN) & (z["score"] >= WP90) & (z["pt"] >= 3.0) & (z["pt"] < 30.0)
    tp = {c: (z[f"e_tp_{c}"] == 1) for c in C.TRIG_CHAINS if f"e_tp_{c}" in z.files}
    s = masks_from_tp(tp, len(z["m"]))[variant]
    m = np.load(btmc, allow_pickle=False)
    qm = (np.abs(m["m"] - C.M0_D0) < C.NOM_WIN) & (m["score"] >= WP90) & (m["pt"] >= 3.0) & (m["pt"] < 30.0) \
        & np.isin(m["truth"], [0, 1]) & np.isin(m["truth_ds"], [0, 1])
    f = C.fit_eps(z["dm"][q & s], z["dm"][q & ~s], DR.shape_from_mc(m["dm"][qm] * 1e3), interval=True)
    out = {"r": float(f["eps"]), "r_err": float(C.eps_sym_err(f)), "determined": bool(f["determined"]), "n_cand": int((q & s).sum())}
    log(f"N_D* 428770[{variant}]: r = {out['r']:.4f} ± {out['r_err']:.4f} (kandydaci {out['n_cand']}, określony {out['determined']})")
    return out


# ----------------------------------------------------------------------------- próbki
class Part:
    """Jedna próbka: Model, histogramy par, fity kanałów (raz), profil g_s(μ) z cache."""

    def __init__(self, label, M, H, channels, nd_rel, f_jk):
        from charmtag.pairs.fit2d import fit_channel
        self.label, self.H, self.channels, self.nd_rel, self.f_jk = label, H, channels, nd_rel, f_jk
        self.fits, self.ch = {}, {}
        for ch in channels:
            n = H[ch]; full = np.ones(n.shape[1:], bool)
            r, cf, m1 = fit_channel(M, ch, n, full, scan=False)
            self.ch[ch] = {"mu_hat": float(r["mu_hat"]), "mu_err_hesse": float(r["mu_err_hesse"]), "valid": bool(r["valid"]), "n": int(n.sum())}
            self.fits[ch] = (cf, m1, n, full)
        self.cache, self.nbad = {}, 0

    def g(self, mu):
        mu = round(float(mu), 6)
        if mu not in self.cache:
            s = 0.0
            for ch in self.channels:
                cf, m1, n, full = self.fits[ch]
                ms = cf.fit(n, full, fix_mu=mu, start=np.array(m1.values))
                self.nbad += int(not ms.valid)
                s += ms.fval - m1.fval
            self.cache[mu] = s
        return self.cache[mu]


def comb_fit(parts, interval=True):
    """Wspólne μ: minimum Σ_s g_s(μ) (profil per kanał, jak replica_fit), 68 % z Δ = 0,5, Z z q0."""
    import ubjk_run as U
    mus = np.array([p.ch[c]["mu_hat"] for p in parts for c in p.channels], float)
    ers = np.array([p.ch[c]["mu_err_hesse"] for p in parts for c in p.channels], float)
    ok = np.isfinite(mus) & np.isfinite(ers) & (ers > 0)
    g = lambda mu: sum(p.g(mu) for p in parts)  # noqa: E731
    lo = max(-1.5, float(np.min(mus[ok])) - 0.3) if ok.any() else -1.5
    hi = min(5.0, float(np.max(mus[ok])) + 0.3) if ok.any() else 5.0
    muc, gmin = U._minimize(g, lo, hi)
    g0 = g(0.0)
    q0 = 2 * (g0 - gmin) if muc > 0 else 0.0
    c = {"mu": muc, "q0": q0, "Z": math.sqrt(max(q0, 0.0)), "at_bound": bool(abs(muc - lo) < 2e-3 or abs(muc - hi) < 2e-3),
         "n_invalid_fits": int(sum(p.nbad for p in parts))}
    if interval:
        h = lambda m: g(m) - gmin - 0.5  # noqa: E731
        s0 = float(1 / math.sqrt(np.sum(1 / ers[ok] ** 2))) if ok.any() else 0.5

        def side(sgn):
            step = 1.5 * s0
            for _ in range(5):
                x = min(5.0, max(-1.5, muc + sgn * step))
                if h(x) > 0:
                    return U._root(h, *sorted((muc, x)), tol=5e-3)
                if x in (5.0, -1.5):
                    return None
                step *= 2
            return None
        up, dn = side(+1), side(-1)
        c["mu_68"] = [dn, up]
        c["sigma_68"] = 0.5 * (up - dn) if (up is not None and dn is not None) else None
        c["sigma_up"] = (up - muc) if up is not None else None
        c["sigma_dn"] = (muc - dn) if dn is not None else None
    return c


def syst_table(parts, comb, per_part):
    """Systematyki wspólnego μ: wagi próbek w = 1/σ_s² (σ_s ze statystyki profilu próbki); N_D* i model wag niezależne per run,
    MC wspólne; statystyka × f_jk per próbka (ważone)."""
    mu = comb["mu"]
    sig = {p.label: per_part[p.label].get("sigma_68") for p in parts}
    w = {k: (1 / v ** 2 if v else 0.0) for k, v in sig.items()}
    W = sum(w.values()) or 1.0
    a = {k: v / W for k, v in w.items()}
    nd = math.sqrt(sum((a[p.label] * p.nd_rel) ** 2 for p in parts)) * abs(mu)
    wag = math.sqrt(sum((a[p.label] * SYST_WAG_REL) ** 2 for p in parts)) * abs(mu)
    common = {k: v * abs(mu) for k, v in SYST_COMMON_REL.items()}
    fjk = math.sqrt(sum(a[p.label] * p.f_jk ** 2 for p in parts))     # efektywny czynnik (wariancje ważone)
    stat = (comb.get("sigma_68") or float("nan")) * fjk
    tot_syst = math.sqrt(sum(v ** 2 for v in common.values()) + wag ** 2)
    return {"weights": a, "f_jk_eff": fjk, "stat_with_fjk": stat, "N_dstar": nd, "model_wag": wag, **{f"common: {k}": v for k, v in common.items()},
            "syst_total_excl_Ndstar": tot_syst, "Z_with_fjk": (comb["Z"] / fjk) if fjk else None}


# ----------------------------------------------------------------------------- główny przebieg
def build_427(a, Model, load_inputs):
    I = load_inputs(os.path.join(a.p427, "pair_inputs"))
    H = dict(np.load(os.path.join(a.p427, "pairs_full.npz")))
    meta = json.load(open(os.path.join(a.p427, "pair_inputs.json"))).get("meta", {})
    n_run = float(meta["dstar_fit_data"]["N"])
    nd = json.load(open(a.ndstar427))["fit"]
    if not nd.get("valid"):
        raise SystemExit(f"STOP: N_D* 427929 (NDS) niepoprawny: {nd}")
    I.n_pair = {k: v * nd["N"] / n_run for k, v in I.n_pair.items()}   # μ w skali N_D* NDS (jak UB427 rescale)
    f_jk = F_JK_428770
    if a.jk427 and os.path.exists(a.jk427):
        js = json.load(open(a.jk427))
        hm = js.get("headline_mode", "fixed")
        f_jk = float(js["per_mode"][hm]["combined"]["ratio_68"])   # σ_jk/σ_profil68 wspólnego μ (jak jk_ev200 w 428770: 1,18)
    return I, H, {"N_run": n_run, "N_NDS": nd["N"], "N_NDS_err": nd["N_err"], "nd_rel": nd["N_err"] / nd["N"], "f_jk": f_jk,
                  "f_jk_source": "jk427" if (a.jk427 and os.path.exists(a.jk427)) else "założone = 428770 (1,18)"}


def build_428(a, variant, I0, L, tkey, tmask, U):
    import btrig_common as C
    from charmtag.syst import files as F
    meta = json.load(open(os.path.join(a.p428, "pair_inputs.json"))).get("meta", {})
    n_run = float(meta["dstar_fit_data"]["N"])
    m = tmask[variant]
    fev = float(m.mean())
    import pandas as pd
    lk = C.ev_key(pd.to_numeric(L.source_file).to_numpy(np.int64), L.entry.to_numpy(np.int64))
    pos = np.clip(np.searchsorted(tkey, lk), 0, len(tkey) - 1)
    hit = tkey[pos] == lk
    ins = hit & m[pos]
    Lv = L.assign(file_idx=np.where(ins, 0, 1).astype(np.int64))
    H, DL, SB, keys = U.perfile_from_legs(Lv, I0, 2)
    ds = dstar_fraction(a.bt428, a.bt_mc, variant)
    info = {"variant": variant, "event_fraction": fev, "legs_unmatched": int((~hit).sum()), "legs_in": int(ins.sum()), "dstar": ds,
            "pairs": {ch: int(H[ch][0].sum()) for ch in H}, "N_run": n_run}
    if not (ds["determined"] and np.isfinite(ds["r"]) and ds["r"] > 0):
        return None, None, info
    I2 = F.scaled_inputs(I0, ds["r"])
    I2.dataleg = {key: np.asarray(DL[key][0], float) for key in I0.dataleg}
    I2.sbsb = {key: int(round(SB[0, j])) for j, key in enumerate(keys)}
    I2.nev_legs = I0.nev_legs * fev
    I2.ndata = I0.ndata * fev
    I2.n_pair = {k: v * N_NDS_428770[0] / n_run for k, v in I2.n_pair.items()}
    info["nd_rel"] = math.hypot(N_NDS_428770[1] / N_NDS_428770[0], ds["r_err"] / ds["r"])
    return I2, {ch: H[ch][0] for ch in H}, info


def run(a):
    import pandas as pd
    import ubjk_run as U
    charmtag, CHANNELS, Model, load_inputs, _, _ = U._charmtag()
    t0 = time.time()
    rep = {"charmtag": charmtag.__file__, "variants": {}, "consts": {"N_NDS_428770": N_NDS_428770, "F_JK_428770": F_JK_428770,
                                                                       "SYST_COMMON_REL": SYST_COMMON_REL, "SYST_WAG_REL": SYST_WAG_REL,
                                                                       "STRICT": STRICT, "AFP": AFP}}
    chans = list(CHANNELS)
    # 427929
    I7, H7, m7 = build_427(a, Model, load_inputs)
    P7 = Part("427929", Model(I7), H7, chans, m7["nd_rel"], m7["f_jk"])
    r7 = comb_fit([P7], interval=not a.quick)
    rep["r427929"] = {**m7, "channels": P7.ch, "combined": r7}
    log(f"427929: μ = {r7['mu']:.3f} {r7.get('mu_68')} Z = {r7['Z']:.2f}")
    # 428770 — kontrola: pełna próbka odtwarza nominał 1,50 (skala NDS)
    I8 = load_inputs(os.path.join(a.p428, "pair_inputs"))
    H8 = dict(np.load(os.path.join(a.p428, "pairs_full.npz")))
    n8 = float(json.load(open(os.path.join(a.p428, "pair_inputs.json")))["meta"]["dstar_fit_data"]["N"])
    I8n = __import__("copy").copy(I8); I8n.n_pair = {k: v * N_NDS_428770[0] / n8 for k, v in I8.n_pair.items()}
    if not a.quick:
        Pn = Part("428770_all", Model(I8n), H8, chans, N_NDS_428770[1] / N_NDS_428770[0], F_JK_428770)
        rn = comb_fit([Pn], interval=False)
        rep["r428770_all_check"] = {"mu": rn["mu"], "expected": 1.50, "ok": bool(abs(rn["mu"] - 1.50) < 0.02)}
        log(f"kontrola 428770 całość (skala NDS): μ = {rn['mu']:.3f} (oczek. 1,50) ok {rep['r428770_all_check']['ok']}")
        if not rep["r428770_all_check"]["ok"]:
            raise SystemExit("STOP: pełna próbka 428770 nie odtwarza μ = 1,50 w skali NDS")
    L = pd.read_parquet(a.legs428)
    tkey, tmask = load_trg_masks(sorted(glob.glob(a.trg428)))
    rep["trg_events"] = int(len(tkey))
    rep["subset_event_fraction"] = {k: float(v.mean()) for k, v in tmask.items()}
    for v in a.variants.split(","):
        I8v, H8v, info = build_428(a, v, I8, L, tkey, tmask, U)
        R = {"info": info}
        if I8v is None:
            R["note"] = "N_D* podzbioru nieokreślone → wspólne = sam 427929"
            R["combined"] = r7
            rep["variants"][v] = R; log(v, R["note"]); continue
        try:
            P8 = Part(f"428770[{v}]", Model(I8v), H8v, chans, info["nd_rel"], F_JK_428770)
            R["r428770_subset"] = {"channels": P8.ch, "combined": comb_fit([P8], interval=not a.quick)}
            per = {"427929": r7, P8.label: R["r428770_subset"]["combined"]}
            R["combined"] = comb_fit([P7, P8], interval=not a.quick)
            if not a.quick:
                R["syst"] = syst_table([P7, P8], R["combined"], per)
        except Exception as ex:  # noqa: BLE001
            R["error"] = repr(ex); log(v, "BŁĄD", ex)
        rep["variants"][v] = R
        c = R.get("combined", {})
        log(f"{v}: 428770[S] μ = {R.get('r428770_subset', {}).get('combined', {}).get('mu')}; wspólne μ = {c.get('mu')} {c.get('mu_68')}")
    if not a.quick:
        rep["syst_427929_alone"] = syst_table([P7], r7, {"427929": r7})
    rep["sec"] = time.time() - t0
    return rep


def fmt(x, d=2):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def to_md(rep):
    r7 = rep["r427929"]; c7 = r7["combined"]
    L = ["# UBCOMB — wspólne μ par D–D z próbek bez triggera dżetowego", "",
         "μ w skali zapisu bieżącego (N_D* metodą NDS; 428770 całość = 1,50). Wspólne μ = minimum sumy profili NLL (próbka × kanał); osobne N_D*, c_w,",
         "tło fałszywych nóg i SB×SB per próbka. Przedział 68 % z profilu (bez czynnika korelacji zdarzeń f_jk; z f_jk w tabeli systematyk).", "",
         f"Kontrola: 428770 całość w skali NDS → μ = {fmt(rep.get('r428770_all_check', {}).get('mu'), 3)} (oczek. 1,50).", "",
         "| próbka / wariant | udział zdarzeń 428770 | r = udział N_D* | pary D0D0 / DpD0 / DpDp | μ 428770[S] [68 %] | μ wspólne [68 %] | Z wspólne |",
         "|---|---|---|---|---|---|---|",
         f"| 427929 sam (MB) | – | – | – | – | **{fmt(c7['mu'])}** [{fmt((c7.get('mu_68') or [None, None])[0])}; {fmt((c7.get('mu_68') or [None, None])[1])}] | {fmt(c7['Z'])} |"]
    for v, R in rep["variants"].items():
        info = R["info"]; ds = info["dstar"]; c = R.get("combined", {}); s8 = R.get("r428770_subset", {}).get("combined", {})
        p = " / ".join(str(info["pairs"].get(ch, 0)) for ch in ("D0D0", "DpD0", "DpDp"))
        m68 = c.get("mu_68") or [None, None]; s68 = s8.get("mu_68") or [None, None]
        L.append(f"| 427929 + 428770[{v}] | {info['event_fraction']:.4f} | {fmt(ds['r'], 4)} ± {fmt(ds['r_err'], 4)} | {p} | "
                 f"{fmt(s8.get('mu'))} [{fmt(s68[0])}; {fmt(s68[1])}] | **{fmt(c.get('mu'))}** [{fmt(m68[0])}; {fmt(m68[1])}] | {fmt(c.get('Z'))} |"
                 + (f" {R.get('note') or R.get('error')}" if (R.get("note") or R.get("error")) else ""))
    L += ["", "Warianty 428770: strict = łańcuchy bez obciążenia (noalg_mb_L1MBTS_1/2, mb_mbts_L1MBTS_1, mb_sptrk_L1MBTS_1, RD0_FILLED); "
          "nonj12 = klasy UBTRIG MBTS_noj20 ⊕ other (zawiera sptrk_pt4/6/8 = wymóg twardego toru, AFP i zdarzenia spoza listy łańcuchów → NIE jest czyste); "
          "nonj12_noafp = nonj12 bez AFP.", ""]
    L += ["## Systematyki (raportowane; budżet przejściowy VER_budget, μ_ref = 1,50)", "", "| | 427929 sam | " + " | ".join(f"+ 428770[{v}]" for v in rep["variants"]) + " |",
          "|---|---|" + "---|" * len(rep["variants"])]
    cols = [rep.get("syst_427929_alone", {})] + [R.get("syst", {}) for R in rep["variants"].values()]
    keys = ["stat_with_fjk", "f_jk_eff", "N_dstar", "model_wag"] + [k for k in (cols[0] or {}) if k.startswith("common:")] + ["syst_total_excl_Ndstar", "Z_with_fjk"]
    for k in keys:
        L.append(f"| {k} | " + " | ".join(fmt(c.get(k)) for c in cols) + " |")
    L += ["", f"f_jk 427929: {r7['f_jk']:.2f} ({r7['f_jk_source']}). N_D* 427929 (NDS) = {r7['N_NDS']:.0f} ± {r7['N_NDS_err']:.0f}; 428770 = {N_NDS_428770[0]:.0f} ± {N_NDS_428770[1]:.0f}.",
          "mc_stat (n_eff) ma zejść po normalizacji z prawdy (gen_norm, PASS 19.11) — liczby w tabeli są z budżetu przejściowego.", ""]
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser()
    for k in ("--p428", "--legs428", "--trg428", "--bt428", "--bt-mc", "--p427", "--ndstar427"):
        ap.add_argument(k, required=True)
    ap.add_argument("--jk427"); ap.add_argument("--variants", default="strict,nonj12,nonj12_noafp")
    ap.add_argument("--quick", action="store_true", help="test: bez przedziałów i kontroli pełnej próbki")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if os.environ.get("CHARM_UNBLIND") != "1":
        raise SystemExit("UBCOMB: wymaga CHARM_UNBLIND=1")
    os.makedirs(a.out, exist_ok=True)
    rep = run(a)
    import ubjk_run as U
    json.dump(U._clean(rep), open(os.path.join(a.out, "ubcomb.json"), "w"), indent=1)
    open(os.path.join(a.out, "UBCOMB.md"), "w").write(to_md(rep))
    log("OK", os.path.join(a.out, "UBCOMB.md"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
