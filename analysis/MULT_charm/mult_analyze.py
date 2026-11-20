#!/usr/bin/env python3
"""mult_analyze.py — MULT_charm: ⟨N_D⟩/N_ev vs krotność, dane vs MC (mb_nd), pary tylko MC (dane: bramka G5).

Wejście: katalogi npz z mult_extract.py (--data, --mc). Wyjście: --out/<etykieta>/{result.json, TABELA.md, *.png}.
Definicje (MULT-v1, PROCEDURA.md):
  Ŷ_s(j) = [N_s(j)/N_ev(j)] / [Σ_j N_s / Σ_j N_ev]  (samonormalizowany uzysk), x_j = ⟨n⟩_j/⟨n⟩
  dane: N_s(j) z fitu masy (kształt z fitu zintegrowanego), waga zdarzeń = trigPS łańcucha schematu
  MC  : (a) reko dopasowane do prawdy (Σ w po poprawnych kandydatach, bez fitu), (b) fit jak w danych (zamknięcie G2),
        (c) prawda: prawdziwe D w akceptacji / zdarzenie vs krotność reko (bez ε) → ε_rel(j) = Ŷ_a/Ŷ_c
  R_fit(j) = Ŷ_dane,fit / Ŷ_MC,fit (GŁÓWNE: ta sama metoda); R = Ŷ_dane/Ŷ_MC,a; R_corr = (Ŷ_dane/ε_rel)/Ŷ_MC,c (ważne przy G2 PASS)
  α        : Ŷ ∝ x^α na x ∈ [0,5; 4]; Δα_fit = α_dane,fit − α_MC,fit (główne), Z = Δα/σ
  pary MC  : N_par/N_D (prawda D0+D⁺, LS/OS) i reko-prawda vs x; κ (wpływ na μ, §5 PROCEDURY)
  python3 mult_analyze.py --data DIR --mc DIR --scheme sptrk4 --est opp --w2 none|self|PLIK.npy --out OUT [--label L]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mult_common as C  # noqa: E402
from mult_fit import Fitter, refl_from_mc  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
APPROVAL = os.path.join(HERE, "APPROVAL_PAIRS_DATA.txt")
MODEL = {"D0": "d0", "Dplus": "dplus", "Dstar": "dstar"}


# ----------------------------------------------------------------------------- wczytanie
def load(spec, kind):
    """spec: KATALOG albo lista 'KATALOG[:waga],KATALOG[:waga]' (mieszanka MC ND+SD+DD+CD; waga = σ/N_gen → kolumna ev_wproc)."""
    files, wf = [], []
    for part in spec.split(","):
        d, _, w = part.partition(":")
        fs = sorted(glob.glob(os.path.join(d, "*.npz")))
        if not fs:
            sys.exit(f"brak npz w {d}")
        files += fs; wf += [float(w) if w else 1.0] * len(fs)
    E, Cc, T, metas = {}, {}, {}, []
    off = 0
    for f in files:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        if meta["kind"] != kind:
            sys.exit(f"{f}: kind={meta['kind']} ≠ {kind}")
        metas.append(meta)
        E.setdefault("wproc", []).append(np.full(meta["n_events"], wf[len(metas) - 1]))
        for pref, D in (("ev_", E), ("c_", Cc), ("t_", T)):
            for k in z.files:
                if k.startswith(pref):
                    v = z[k]
                    if k[len(pref):] == "entry":
                        v = v + off
                    D.setdefault(k[len(pref):], []).append(v)
        off += meta["n_events"]
    cat = lambda D: {k: np.concatenate(v) for k, v in D.items()}  # noqa: E731
    E, Cc, T = cat(E), cat(Cc), cat(T)
    # wspólne kolumny triggerów: brak w części plików → błąd (mieszanie produkcji)
    nev = len(E["entry"])
    for k, v in E.items():
        if len(v) != nev:
            sys.exit(f"kolumna ev_{k} nie we wszystkich plikach ({len(v)} ≠ {nev}) — mieszane produkcje?")
    fmts = sorted({m["format"] for m in metas})
    return E, Cc, T, {"files": len(files), "n_events": nev, "formats": fmts,
                      "compat_missing": sorted({c for m in metas for c in m["compat_missing"]})}


# ----------------------------------------------------------------------------- selekcja
def event_mask(E, kind, scheme):
    s = (E["npv"] == C.EVSEL["npv_eq"]) & (np.abs(E["pvz"]) < C.EVSEL["pvz_max"])
    if C.EVSEL["npu_eq"] is not None:
        s &= E["npu"] == C.EVSEL["npu_eq"]
    s &= E["nbl"] >= C.EVSEL["nbl_min"]
    sc = C.TRIG_SCHEMES[scheme]
    s &= E["leadbl"] > sc["lead_min"]
    w = np.ones(len(s))
    if kind == "data":
        s &= E["grl"] == 1
        if sc["chain"]:
            tp, ps = f"tp_{sc['chain']}", f"ps_{sc['chain']}"
            if tp not in E:
                sys.exit(f"schemat {scheme}: brak gałęzi trigPassed_{sc['chain']} w danych")
            s &= E[tp] == 1
            if ps in E:
                w = np.where(E[ps] > 0, E[ps], 0.0).astype(float)
                s &= w > 0
    return s, w


def d0_hash_keep(c):
    """Arbitraż hash D0 (19.11): ta sama para torów zapisana z obiema hipotezami K↔π (derywacja trzyma obie, gdy obie masy są w
    ~1,80–1,92 GeV → „pudełko” w tle m(Kπ), które fit brał za sygnał). Zostaje jeden wiersz na (zdarzenie, para torów), wybór
    deterministyczny i niezależny od masy: bit hasha (entry, min trk, max trk) wskazuje, czy K = tor o mniejszym czy większym id."""
    d0 = c["sp"] == C.SP_CODE["D0"]
    a = np.minimum(c["trk0"], c["trk1"]).astype(np.int64); b = np.maximum(c["trk0"], c["trk1"]).astype(np.int64)
    e = c["entry"].astype(np.int64)
    with np.errstate(over="ignore"):
        hb = (((e * 1000003 + a) * 1000033 + b).astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15) >> np.uint64(63)).astype(bool)
    want = np.where(hb, b, a)                                      # id toru, który ma być kaonem
    keep = np.ones(len(e), bool)
    idx = np.flatnonzero(d0)
    order = np.lexsort((b[idx], a[idx], e[idx])); ii = idx[order]
    same = (e[ii][1:] == e[ii][:-1]) & (a[ii][1:] == a[ii][:-1]) & (b[ii][1:] == b[ii][:-1])
    dup = np.zeros(len(ii), bool); dup[1:] |= same; dup[:-1] |= same          # wiersze w grupach ≥ 2
    keep[ii[dup]] = c["trk0"][ii[dup]].astype(np.int64) == want[ii[dup]]
    return keep


def cand_mask(c, sp):
    n = C.NOM[sp]
    s = (c["sp"] == C.SP_CODE[sp]) & (c["pt"] > C.PT_MIN) & (c["pt"] < C.PT_MAX) & (np.abs(c["eta"]) < C.ETA_MAX)
    if sp == "Dplus":
        s &= (c["lxy"] > n["lxy_min"]) & (c["chi2"] < n["chi2_max"]) & (c["impsig"] < n["impsig_max"]) & (c["cosp"] > n["cosp_min"])
        s &= c["dsveto"] >= n["dstar_veto"]
    elif sp == "D0":
        s &= (c["lxy"] > n["lxy_min"]) & (c["chi2"] < n["chi2_max"]) & (c["impsig"] < n["impsig_max"]) & (c["pt"] > n.get("pt_min", 0))
        s &= d0_hash_keep(c)
    else:
        s &= (np.abs(c["md0"] - C.M_D0) < n["md0_win"]) & (c["dm"] > n["dm_lo"]) & (c["dm"] < n["dm_hi"]) & (c["pt"] > n.get("pt_min", 0))
    return s


def mult_event(E, est):
    if est in ("bl", "blx"):
        return E["nbl"].astype(float)
    if est == "trk500":
        return E["n500"].astype(float)
    return 0.5 * (E["nA"] + E["nC"]).astype(float)        # opp: średnia półkul (do ⟨n⟩)


def mult_cand(c, E_idx, E, est):
    if est == "bl":
        return E["nbl"][E_idx].astype(float)
    if est == "blx":
        return (E["nbl"][E_idx] - c["ndau_acc"]).astype(float)
    if est == "trk500":
        return E["n500"][E_idx].astype(float)
    neg = c["eta"] < 0          # D w C (η<0) → półkula A; córki tego kandydata w obszarze A odjęte
    return np.where(neg, E["nA"][E_idx] - c["ndau_A"], E["nC"][E_idx] - c["ndau_C"]).astype(float)


def denom(E, sel, w, est, nmean, edges_x):
    """Σ wag zdarzeń w binach x (opp: ½ do bin(n_A) i ½ do bin(n_C)); ⟨x⟩ w binie."""
    nb = len(edges_x) - 1
    if est == "opp":
        xs = np.concatenate([E["nA"][sel], E["nC"][sel]]) / nmean
        ws = np.concatenate([w[sel], w[sel]]) * 0.5
    else:
        xs = mult_event(E, est)[sel] / nmean; ws = w[sel]
    b = np.clip(np.digitize(xs, edges_x) - 1, 0, nb - 1)
    N = np.bincount(b, ws, nb); X = np.bincount(b, ws * xs, nb)
    return N, np.divide(X, N, out=np.full(nb, np.nan), where=N > 0)


# ----------------------------------------------------------------------------- jedna próbka
class Sample:
    def __init__(self, kind, E, Cc, T, scheme, est, w_ev=None, nmean=None):
        self.kind, self.E, self.c, self.T, self.est = kind, E, Cc, T, est
        sel, w = event_mask(E, kind, scheme)
        w = w * E["wproc"] / E["wproc"].max()
        if w_ev is not None:
            w = w * w_ev
        self.sel, self.w = sel, w
        # nmean: dla próbki fd (ε_rel) — ⟨n⟩ z mb_nd, bo fd jest warunkowane obecnością charmu
        self.nmean = float(np.average(mult_event(E, est)[sel], weights=w[sel])) if nmean is None else float(nmean)
        self.N_ev, self.x_ev = denom(E, sel, w, est, self.nmean, C.X_EDGES)
        # kandydaci w wybranych zdarzeniach (entry = globalny indeks zdarzenia)
        idx = self.c["entry"]
        self.cw = w[idx] * sel[idx]
        self.cx = mult_cand(self.c, idx, E, est) / self.nmean
        self.cb = C.x_bin(self.cx)

    def var(self, sp):
        return self.c["dm"] if sp == "Dstar" else self.c["m"]


def truth_ref_shape(S, sp, refl=None):
    """Kształt referencyjny sygnału: fit samego sygnału na kandydatach MC z prawdą (true == 1), po cięciach NOM (19.11)."""
    m = cand_mask(S.c, sp) & (S.c["true"] == 1) & (S.cw > 0)
    return Fitter(MODEL[sp], refl).truth_shape(S.var(sp)[m])


def fit_yields(S, sp, refl=None, fix_shape=None, ref=None):
    """Fity w binach x. Zwraca (N, σN, info). Wagi: N_w = N_fit · ⟨w⟩ kandydatów w binie (trigPS niezależne od D).
    ref: kształt sygnału z prawdy MC — fit zintegrowany ma μ/σ₁ w przedziałach wokół niego (mult_fit, poprawka 19.11)."""
    m = cand_mask(S.c, sp) & (S.cw > 0)
    F = Fitter(MODEL[sp], refl)
    x = S.var(sp)
    hint = F.hist(x[m])
    r0 = F.fit(hint, ref=ref) if fix_shape is None else {"shape": fix_shape}
    nb = len(C.X_EDGES) - 1
    N, eN, info = np.zeros(nb), np.zeros(nb), []
    for j in range(nb):
        mj = m & (S.cb == j)
        if mj.sum() < 30:
            info.append({"bin": j, "n": int(mj.sum()), "skip": True}); continue
        r = F.fit(F.hist(x[mj]), shape=r0["shape"])
        wm = float(S.cw[mj].mean())
        N[j], eN[j] = r["N_S"] * wm, r["err"] * wm
        info.append({"bin": j, "n": int(mj.sum()), "N_fit": r["N_S"], "err_fit": r["err"], "w_mean": wm,
                     "low_signif": bool(r.get("err_fallback") or not (r["err"] > 0 and r["N_S"] / r["err"] > 3)),
                     "err_fallback": bool(r.get("err_fallback"))})   # G3: poza fitem α
    return N, eN, {"integrated": {k: v for k, v in r0.items() if k != "shape"} | {"shape": r0["shape"]}, "bins": info, "ref_shape": ref}


def mc_true_yields(S, sp):
    m = cand_mask(S.c, sp) & (S.c["true"] == 1) & (S.cw > 0)
    nb = len(C.X_EDGES) - 1
    N = np.bincount(S.cb[m], S.cw[m], nb); N2 = np.bincount(S.cb[m], S.cw[m] ** 2, nb)
    return N, np.sqrt(N2)


def mc_truth_level(S, sp):
    """Prawdziwe D w akceptacji / zdarzenie vs krotność reko zdarzenia (opp: półkula przeciwna do prawdziwego D)."""
    nb = len(C.X_EDGES) - 1
    if S.est == "opp":
        t = S.T; ap = np.abs(t["pdg"]); m = (ap == C.PDG[sp]) & S.sel[t["entry"]]
        e = t["entry"][m]
        n = np.where(t["eta"][m] < 0, S.E["nA"][e], S.E["nC"][e]) / S.nmean
        b = C.x_bin(n); w = S.w[e]
    else:
        e = np.where(S.sel)[0]
        n = mult_event(S.E, S.est)[e] / S.nmean
        b = C.x_bin(n); w = S.w[e] * S.E[f"nt_{sp}"][e]
    N = np.bincount(b, w, nb); N2 = np.bincount(b, w * (S.w[e] if S.est != "opp" else w), nb)
    return N, np.sqrt(np.abs(N2))


def selfnorm(N, eN, Nev):
    Y = np.divide(N, Nev, out=np.zeros_like(N), where=Nev > 0)
    eY = np.divide(eN, Nev, out=np.zeros_like(N), where=Nev > 0)
    tot = N.sum() / Nev.sum()
    return Y / tot, eY / tot


def mc_pairs(S):
    """N_par/N_D vs x — prawda (D0+D⁺, zdarzenie, wszystkie / LS / OS) i reko-prawda (≥2 poprawne D0/D⁺ po cięciach NOM)."""
    nb = len(C.X_EDGES) - 1
    e = np.where(S.sel)[0]
    x = mult_event(S.E, "bl" if S.est == "opp" else S.est)[e] / (S.nmean if S.est != "opp" else np.average(S.E["nbl"][S.sel], weights=S.w[S.sel]))
    b = C.x_bin(x); w = S.w[e]
    nc, ncb = S.E["nt_c"][e].astype(float), S.E["nt_cbar"][e].astype(float)
    nD = nc + ncb
    out = {"N_D": np.bincount(b, w * nD, nb), "pairs": np.bincount(b, w * nD * (nD - 1) / 2, nb),
           "LS": np.bincount(b, w * (nc * (nc - 1) + ncb * (ncb - 1)) / 2, nb), "OS": np.bincount(b, w * nc * ncb, nb),
           "N_Dstar": np.bincount(b, w * S.E["nt_Dstar"][e], nb)}
    m = (cand_mask(S.c, "D0") | cand_mask(S.c, "Dplus")) & (S.c["true"] == 1) & (S.cw > 0)
    k = np.bincount(S.c["entry"][m], minlength=len(S.E["entry"])).astype(float)[e]
    out["reco_N_D"] = np.bincount(b, w * k, nb); out["reco_pairs"] = np.bincount(b, w * k * (k - 1) / 2, nb)
    out["note"] = "x pary = nSelBL (bez wyłączenia córek; dla est=opp też nSelBL)"
    return {kk: (v.tolist() if isinstance(v, np.ndarray) else v) for kk, v in out.items()}


def r_powerlaw(x, P_ev, dalpha):
    """r(j) = x_j^Δα / ⟨x^Δα⟩_zdarzenia (wygładzony stosunek dane/MC; normalizacja: średnio 1 na zdarzenie)."""
    x = np.nan_to_num(np.asarray(x, float), nan=1.0); x = np.clip(x, 1e-3, None); P_ev = np.asarray(P_ev, float)
    r = x ** dalpha
    return r / ((P_ev * r).sum() / P_ev.sum())


def kappa(r_D, r_Dstar, P_pair, P_Dstar, f_dps=0.4):
    """Wpływ na μ (PROCEDURA §5): pary/D* w MC ×κ, gdy uzysk D w danych/MC = r(j). DPS: para ∝ r², SPS: para ∝ r."""
    r_D, r_Dstar, P_pair, P_Dstar = map(lambda a: np.asarray(a, float), (r_D, r_Dstar, P_pair, P_Dstar))
    den = (P_Dstar * r_Dstar).sum() / P_Dstar.sum()
    k_dps = (P_pair * r_D ** 2).sum() / P_pair.sum() / den
    k_sps = (P_pair * r_D).sum() / P_pair.sum() / den
    return {"kappa_DPS": float(k_dps), "kappa_SPS": float(k_sps), "f_DPS": f_dps,
            "kappa_mix": float(f_dps * k_dps + (1 - f_dps) * k_sps), "mu_over_kappa_mix": float(1.50 / (f_dps * k_dps + (1 - f_dps) * k_sps))}


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data"); ap.add_argument("--mc", required=True, help="KATALOG albo 'mb_nd:2.781e-6,mb_sd:2.143e-6,…' (PROCEDURA §2)")
    ap.add_argument("--scheme", default="sptrk4", choices=list(C.TRIG_SCHEMES))
    ap.add_argument("--est", default="opp", choices=list(C.ESTIMATORS))
    ap.add_argument("--w2", default="none", help="none | self (w2 z tych danych, zmienne gen_norm: n500, lead_trk) | ścieżka w2d.npy")
    ap.add_argument("--out", required=True); ap.add_argument("--label", default=None)
    ap.add_argument("--pairs-data", action="store_true", help="pary na DANYCH — wymaga APPROVAL_PAIRS_DATA.txt (bramka G5)")
    ap.add_argument("--eps-mc", default=None, help="katalog npz fd_nd (wymuszony rozpad): ε_rel(x) z dużą statystyką (zalecane)")
    ap.add_argument("--inject", type=float, default=0.0, help="TEST: mnoży wagi kandydatów danych przez x^inject (sprawdzenie czułości Δα)")
    ap.add_argument("--dstar-tight", action="store_true", help="G2 (D* FAIL): zamrożona próba — okno m(D0) 30 MeV, pT > 4 GeV dla D*")
    ap.add_argument("--d0-tight", action="store_true", help="19.11 (D0 Δα = NaN: S/B ≈ 1–2 %%, wszystkie biny N/σ ≤ 3): zamrożona próba — pT(D0) > 6 GeV "
                    "(wybór z MC mb_nd przed danymi: s/√b ×1,58, S/B 1,6 %% → 10 %%)")
    ap.add_argument("--no-pu-veto", action="store_true", help="G1c: bez weta nPileupVertices == 0")
    a = ap.parse_args()
    if a.no_pu_veto:
        C.EVSEL["npu_eq"] = None
    if a.d0_tight:
        C.NOM["D0"]["pt_min"] = 6.0
    if a.dstar_tight:
        C.NOM["Dstar"]["md0_win"] = 0.030; C.NOM["Dstar"]["pt_min"] = 4.0
    if a.pairs_data:
        if not os.path.exists(APPROVAL):
            sys.exit("STOP (G5): pary na danych wymagają zatwierdzenia — brak APPROVAL_PAIRS_DATA.txt. Nic nie policzono.")
        sys.exit("G5: zgoda jest, ale pary danych liczy się kodem charmtag (fit2d per bin krotności), nie tym skryptem — PROCEDURA §G5.")
    label = a.label or f"{a.scheme}_{a.est}_w2{os.path.basename(a.w2).replace('.npy', '')}"
    od = os.path.join(a.out, label); os.makedirs(od, exist_ok=True)

    Em, Cm, Tm, im = load(a.mc, "mc")
    res = {"def": C.DEF_VERSION, "scheme": a.scheme, "scheme_info": C.TRIG_SCHEMES[a.scheme], "est": a.est, "w2": a.w2,
           "x_edges": C.X_EDGES.tolist(), "mc_input": im, "inject": a.inject}
    Ed = None
    if a.data:
        Ed, Cd, Td, idd = load(a.data, "data")
        res["data_input"] = idd
    # wagi w2 MC
    w_ev = None
    if a.w2 == "self":
        if Ed is None:
            sys.exit("--w2 self wymaga --data")
        sd, wd = event_mask(Ed, "data", a.scheme); sm, _ = event_mask(Em, "mc", a.scheme)
        W = C.w2_map(Ed["n500"][sd], Ed["lead_trk"][sd], Em["n500"][sm], Em["lead_trk"][sm], d_w=wd[sd], m_w=Em["wproc"][sm])
        w_ev = C.w2_lookup(W, Em["n500"], Em["lead_trk"]); np.save(os.path.join(od, "w2_self.npy"), W)
    elif a.w2 != "none":
        W = np.load(a.w2); w_ev = C.w2_lookup(W, Em["n500"], Em["lead_trk"])
    Sm = Sample("mc", Em, Cm, Tm, a.scheme, a.est, w_ev)
    res["mc"] = {"n_sel": int(Sm.sel.sum()), "nmean": Sm.nmean, "N_ev": Sm.N_ev.tolist(), "x": Sm.x_ev.tolist()}
    # refleksja D0 z MC (wiersze z zamienioną hipotezą: |tpdg|=421, true=0)
    mD0 = cand_mask(Cm, "D0")
    refl = refl_from_mc(Cm["m"][mD0 & (np.abs(Cm["tpdg"]) == 421) & (Cm["true"] == 0)], Cm["m"][mD0 & (Cm["true"] == 1)])
    res["d0_reflection_from_mc"] = refl
    Sf = None
    if a.eps_mc:
        Ef, Cf, Tf, ifd = load(a.eps_mc, "mc"); res["eps_mc_input"] = ifd
        Sf = Sample("mc", Ef, Cf, Tf, a.scheme, a.est, nmean=Sm.nmean)
    if a.data:
        Sd = Sample("data", Ed, Cd, Td, a.scheme, a.est)
        if a.inject:
            Sd.cw = Sd.cw * np.clip(Sd.cx, 1e-3, None) ** a.inject
        res["data"] = {"n_sel": int(Sd.sel.sum()), "sum_w": float(Sd.w[Sd.sel].sum()), "nmean": Sd.nmean,
                       "N_ev": Sd.N_ev.tolist(), "x": Sd.x_ev.tolist()}
    lines = [f"# MULT_charm — {label}", "", f"def {C.DEF_VERSION}; schemat `{a.scheme}` ({C.TRIG_SCHEMES[a.scheme]['role']}); estymator `{a.est}`; w2 `{a.w2}`",
             f"MC: {im['files']} plików, {im['n_events']} zd. ({', '.join(im['formats'])}; compat {im['compat_missing']}), wybrane {int(Sm.sel.sum())}, ⟨n⟩={Sm.nmean:.2f}"]
    if a.data:
        lines.append(f"Dane: {idd['files']} plików, {idd['n_events']} zd., wybrane {int(Sd.sel.sum())} (Σw {Sd.w[Sd.sel].sum():.0f}), ⟨n⟩={Sd.nmean:.2f}")
    lines += ["", "| gatunek | j | x(MC) | Ŷ MC reko-prawda | Ŷ MC prawda | ε_rel | Ŷ MC fit | x(dane) | Ŷ dane | R | R_corr |", "|" + "---|" * 11]
    res["species"] = {}
    for sp in C.SPECIES:
        Nt, eNt = mc_true_yields(Sm, sp)
        Ya, eYa = selfnorm(Nt, eNt, Sm.N_ev)
        Ng, eNg = mc_truth_level(Sm, sp)
        Yc, eYc = selfnorm(Ng, eNg, Sm.N_ev)
        eps = np.divide(Ya, Yc, out=np.full_like(Ya, np.nan), where=Yc > 0)
        R_eps = {"eps_rel_mb": eps.tolist()}
        if Sf is not None:   # ε_rel z fd_nd: (reko-prawda / prawda) w binach x mb_nd, samonormalizowane
            Ntf, _ = mc_true_yields(Sf, sp); Ngf, _ = mc_truth_level(Sf, sp)
            ef = np.divide(Ntf, Ngf, out=np.full_like(Ntf, np.nan), where=Ngf > 0) / (Ntf.sum() / max(Ngf.sum(), 1e-300))
            eef = ef * np.sqrt(np.divide(1.0, Ntf, out=np.zeros_like(Ntf), where=Ntf > 0))
            R_eps.update({"eps_rel_fd": ef.tolist(), "eps_rel_fd_err": eef.tolist(), "eps_alpha_fd": C.powerlaw_fit(Sm.x_ev, ef, eef)})
            eps = ef
        ref = truth_ref_shape(Sm, sp, refl)
        Nf, eNf, finfo_m = fit_yields(Sm, sp, refl, ref=ref)
        Yf, eYf = selfnorm(Nf, eNf, Sm.N_ev)
        R = {"mc_reco_true": [Ya.tolist(), eYa.tolist()], "mc_truth": [Yc.tolist(), eYc.tolist()], "eps_rel": eps.tolist(),
             "mc_fit": [Yf.tolist(), eYf.tolist()], "closure_fit_over_true": np.divide(Nf, Nt, out=np.full_like(Nf, np.nan), where=Nt > 0).tolist(),
             "fit_mc": finfo_m, "alpha": {}, **R_eps}
        for nm, (x, y, e) in {"mc_reco_true": (Sm.x_ev, Ya, eYa), "mc_truth": (Sm.x_ev, Yc, eYc), "mc_fit": (Sm.x_ev, Yf, eYf)}.items():
            R["alpha"][nm] = C.powerlaw_fit(x, y, e)
        # G2: zamknięcie fitu w MC (χ² N_fit vs N_true w binach z N_true > 0, σ² = σ_fit² + σ_true²) i zgodność α
        ok = (Nt > 0) & (eNf > 0)
        chi2 = float((((Nf - Nt) ** 2) / (eNf ** 2 + eNt ** 2))[ok].sum()); ndf = int(ok.sum())
        af, at_ = R["alpha"]["mc_fit"], R["alpha"]["mc_reco_true"]
        R["G2"] = {"chi2": chi2, "ndf": ndf, "dalpha_fit_minus_true": af[0] - at_[0], "sigma_alpha_fit": af[1],
                   "PASS": bool(ndf > 0 and chi2 <= 2 * ndf and abs(af[0] - at_[0]) < af[1])}
        if a.data:
            Nd, eNd, finfo_d = fit_yields(Sd, sp, refl, ref=ref)
            Yd, eYd = selfnorm(Nd, eNd, Sd.N_ev)
            low = np.array([b.get("low_signif", True) for b in finfo_d["bins"]])
            eYd = np.where(low, 0.0, eYd)          # G3: biny z N/σ ≤ 3 (lub < 30 kandydatów) nie wchodzą do α
            Rr = np.divide(Yd, Ya, out=np.full_like(Yd, np.nan), where=Ya > 0)
            Ydc = np.divide(Yd, eps, out=np.full_like(Yd, np.nan), where=eps > 0)
            Rc = np.divide(Ydc, Yc, out=np.full_like(Yd, np.nan), where=Yc > 0)
            R.update({"data": [Yd.tolist(), eYd.tolist()], "data_N": [Nd.tolist(), eNd.tolist()], "fit_data": finfo_d,
                      "R": Rr.tolist(), "R_corr": Rc.tolist()})
            R["alpha"]["data"] = C.powerlaw_fit(Sd.x_ev, Yd, eYd)
            R["alpha"]["data_corr"] = C.powerlaw_fit(Sd.x_ev, Ydc, np.divide(eYd, eps, out=np.zeros_like(eYd), where=eps > 0))
            ad, am = R["alpha"]["data"], R["alpha"]["mc_reco_true"]
            R["delta_alpha"] = {"value": ad[0] - am[0], "err": float(np.hypot(ad[1], am[1])),
                                "Z": (ad[0] - am[0]) / float(np.hypot(ad[1], am[1])) if np.isfinite(ad[1]) else np.nan}
            af = R["alpha"]["mc_fit"]       # GŁÓWNY: identyczna metoda (fit) w danych i MC — obciążenie fitu się skraca
            ef_ = float(np.hypot(ad[1], af[1]))
            R["delta_alpha_fit"] = {"value": ad[0] - af[0], "err": ef_, "Z": (ad[0] - af[0]) / ef_ if ef_ > 0 else np.nan}
            Rf = np.divide(Yd, Yf, out=np.full_like(Yd, np.nan), where=Yf > 0); R["R_fit"] = Rf.tolist()
            ac, at = R["alpha"]["data_corr"], R["alpha"]["mc_truth"]
            ec = float(np.hypot(ac[1], at[1]))
            R["delta_alpha_corr"] = {"value": ac[0] - at[0], "err": ec, "Z": (ac[0] - at[0]) / ec if ec > 0 else np.nan,
                                     "eps_source": "fd_nd" if Sf is not None else "mb_nd"}
        res["species"][sp] = R
        for j in range(len(C.X_EDGES) - 1):
            row = [sp, j, f"{Sm.x_ev[j]:.2f}", f"{Ya[j]:.3f}±{eYa[j]:.3f}", f"{Yc[j]:.3f}±{eYc[j]:.3f}", f"{eps[j]:.3f}", f"{Yf[j]:.3f}±{eYf[j]:.3f}"]
            row += ([f"{Sd.x_ev[j]:.2f}", f"{R['data'][0][j]:.3f}±{R['data'][1][j]:.3f}", f"{R['R'][j]:.3f}", f"{R['R_corr'][j]:.3f}"] if a.data else ["–"] * 4)
            lines.append("| " + " | ".join(map(str, row)) + " |")
    lines += ["", "## α (Ŷ ∝ x^α, x∈[0,5;4]): (α, σ, χ², ndf)"]
    for sp, R in res["species"].items():
        lines.append(f"- {sp}: " + "; ".join(f"{k} {v[0]:.3f}±{v[1]:.3f} (χ² {v[2]:.1f}/{v[3]})" for k, v in R["alpha"].items())
                     + (f"; **Δα_fit (główne) = {R['delta_alpha_fit']['value']:.3f} ± {R['delta_alpha_fit']['err']:.3f} (Z={R['delta_alpha_fit']['Z']:.2f})**" if "delta_alpha_fit" in R else "")
                     + f"; G2 {'PASS' if R['G2']['PASS'] else 'FAIL'} (χ² {R['G2']['chi2']:.1f}/{R['G2']['ndf']}, α_fit−α_prawda {R['G2']['dalpha_fit_minus_true']:+.2f})"
                     + (f"; Δα vs reko-prawda = {R['delta_alpha']['value']:.3f} ± {R['delta_alpha']['err']:.3f} (ważne tylko przy G2 PASS)" if "delta_alpha" in R else "")
                     + (f"; Δα_corr (ε z {R['delta_alpha_corr']['eps_source']}) = {R['delta_alpha_corr']['value']:.3f} ± {R['delta_alpha_corr']['err']:.3f}" if "delta_alpha_corr" in R else ""))
    # pary (MC)
    P = mc_pairs(Sm); res["mc_pairs"] = P
    lines += ["", "## Pary — tylko MC (prawda D0+D⁺, pT 3–30, |η|<2,5; x = nSelBL/⟨n⟩)", "| j | N_par/N_D | LS/N_D | OS/N_D | reko-prawda pary/D |", "|---|---|---|---|---|"]
    for j in range(len(C.X_EDGES) - 1):
        nd = P["N_D"][j]; rnd = P["reco_N_D"][j]
        lines.append(f"| {j} | {P['pairs'][j] / nd if nd else float('nan'):.4f} | {P['LS'][j] / nd if nd else float('nan'):.4f} | "
                     f"{P['OS'][j] / nd if nd else float('nan'):.4f} | {P['reco_pairs'][j] / rnd if rnd else float('nan'):.4f} |")
    if a.data and a.est != "opp":
        # r(j) z Δα_fit (prawo potęgowe, binowanie x zdarzeń MC); D = średnia ważona Δα_fit D0 i D⁺ (gatunki z G2 PASS i skończonym σ)
        def da(sp):
            d = res["species"][sp].get("delta_alpha_fit", {}); return d.get("value", np.nan), d.get("err", np.nan)
        vals = [da(sp) for sp in ("D0", "Dplus") if res["species"][sp]["G2"]["PASS"] and np.isfinite(da(sp)[1]) and da(sp)[1] > 0]
        dD = (sum(v / e ** 2 for v, e in vals) / sum(1 / e ** 2 for v, e in vals)) if vals else np.nan
        dS = da("Dstar")[0]
        if np.isfinite(dD) and np.isfinite(dS):
            rD = r_powerlaw(Sm.x_ev, Sm.N_ev, dD); rS = r_powerlaw(Sm.x_ev, Sm.N_ev, dS)
            res["kappa"] = kappa(rD, rS, P["pairs"], P["N_Dstar"]) | {"dalpha_D": dD, "dalpha_Dstar": dS}
        else:
            res["kappa"] = {"note": "brak Δα_fit dla D (G2) lub D* — κ nieliczone", "dalpha_D": dD, "dalpha_Dstar": dS}
        lines += ["", f"## κ (wpływ na μ, PROCEDURA §5): {json.dumps(res['kappa'])}"]
    elif a.data:
        lines += ["", "κ: liczony tylko dla est ∈ {bl, blx, trk500} (binowanie par = nSelBL)"]
    json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1, default=float)
    open(os.path.join(od, "TABELA.md"), "w").write("\n".join(lines) + "\n")
    plot(res, od)
    print("\n".join(lines))


def plot(res, od):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    for ax, sp in zip(axs, C.SPECIES):
        R = res["species"][sp]
        xm = np.asarray(res["mc"]["x"], float)
        for key, lab, mk in (("mc_reco_true", "MC reko-prawda", "s"), ("mc_truth", "MC prawda", "^")):
            y, e = map(np.asarray, R[key]); ax.errorbar(xm, y, e, fmt=mk, label=lab, ms=4)
        if "data" in R:
            xd = np.asarray(res["data"]["x"], float); y, e = map(np.asarray, R["data"])
            ax.errorbar(xd, y, e, fmt="o", color="k", label="dane (fit)", ms=4)
        ax.plot([0, 4], [0, 4], ":", color="grey", lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(f"n/⟨n⟩ ({res['est']})"); ax.set_ylabel("Ŷ"); ax.set_title(sp)
        ax.legend(fontsize=7)
    fig.suptitle(f"{res['scheme']} / {res['est']} / w2 {os.path.basename(res['w2'])}")
    fig.tight_layout(); fig.savefig(os.path.join(od, "yields_vs_mult.png"), dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
