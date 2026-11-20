"""B_trig — wspólne: klasy triggera, nazwy kolumn trigger (h5 v2.0 / ROOT v2.0), złączenie po (run, event), fit pass/fail,
testy zgodności (χ²) i przybliżenie czułości. Projekt B (charmbkg) — nie mieszać z A (pary D–D).

Klasy (dane, zdarzenie z kandydatem; bity = trigPassed po prescale, czyli „zdarzenie zapisane przez ten łańcuch”):
  run 428770 (physics_Main, zdominowany L1_J12):
    J12only   = HLT_noalg_L1J12 ∧ ¬MBTS_ANY            (jetowy, bez żadnego łańcucha MBTS)
    MBTSbased = MBTS_ANY                               (j20_L1MBTS_2 ∨ sptrk* L1MBTS ∨ mbts/noalg_mb L1MBTS)
      podklasy (nakładają się): j20MBTS = HLT_j20_L1MBTS_2, sptrk4 = HLT_mb_sptrk_pt4_L1MBTS_2
    other     = ¬J12 ∧ ¬MBTS_ANY                       (AFP, RD0 …; zwykle mała)
    all       = wszystkie zdarzenia (= nominał B, odtworzenie liczb EV — bramka G2)
  run 427929 (physics_MinBias, MBTS bez prescale; wymagany passGRL):
    MB427929      = HLT_mb_mbts_L1MBTS_1 ∧ GRL         (nominał niezależnego runu)
    sptrk4_427929 = HLT_mb_sptrk_pt4_L1MBTS_2 ∧ GRL     (TEN SAM łańcuch co sptrk4 w 428770 → test „run” przy stałym triggerze)
    all427929     = GRL
Rozłączne zestawy do χ²: {J12only, MBTSbased, other} w 428770 oraz MB427929 (niezależny run).
"""
from __future__ import annotations

import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CW = os.path.dirname(HERE)                                       # charm_work
AA_SRC = os.path.join(CW, "AA_merge", "src")                      # źródło prawdy charmtag (dstar_refit)
CHARMBKG = os.path.expanduser("~/Projects/charm/charm-bkg")      # pakiet charmbkg (val_epsb, salad_v1)

M0_D0 = 1.86484
BINS = {"3-5": (3.0, 5.0), "5-8": (5.0, 8.0), "8-30": (8.0, 30.0)}
WPS = ["0.95", "0.90", "0.80", "0.70", "0.50"]
E_DM = np.arange(139.6, 170.0001, 0.1)                            # MeV, jak EV
NOM_WIN = 0.055                                                   # |m(Kπ) − m0| < 55 MeV (EV „nom”)
# nominał B (EV pełna lumi, fit pf, pdgid) — punkt odniesienia testu zgodności
NOMINAL_SF95 = {"3-5": (float("nan"), float("nan")), "5-8": (0.83, 0.10), "8-30": (1.00, 0.04)}

MBTS_CHAINS = ["HLT_j20_L1MBTS_2", "HLT_mb_sptrk_pt4_L1MBTS_2", "HLT_mb_sptrk_pt2_L1MBTS_2", "HLT_mb_mbts_L1MBTS_1",
               "HLT_noalg_mb_L1MBTS_1", "HLT_noalg_mb_L1MBTS_2", "HLT_mb_sptrk_L1MBTS_1"]
J12 = "HLT_noalg_L1J12"
TRIG_CHAINS = sorted(set(MBTS_CHAINS + [J12, "HLT_noalg_L1AFP_A_OR_C", "HLT_mb_sptrk_pt4_L1AFP_A_OR_C",
                                        "HLT_noalg_L1RD0_FILLED", "HLT_mb_sptrk_pt6_L1MBTS_2", "HLT_mb_sptrk_pt8_L1MBTS_2"]))
EVENT_EXTRA = {   # nazwa logiczna → kandydaci (h5 v2.0 może mieć dowolną z konwencji; ROOT ma EventInfo_*_NOSYS)
    "nSelBL": ["nSelBL", "EventInfo_nSelBL", "n_sel_bl", "nsel_bl"],
    "lead_pt": ["leadTrackPtBL", "EventInfo_leadTrackPtBL", "lead_track_pt_bl", "lead_pt_bl"],
    "actual_mu": ["actual_mu", "actualInteractionsPerCrossing", "EventInfo_actualInteractionsPerCrossing"],
    "grl": ["passGRL_MinBiasAFP", "EventInfo_passGRL_MinBiasAFP", "pass_grl", "grl"],
    "npu": ["nPileupVertices", "EventInfo_nPileupVertices", "n_pileup_vertices"],
    "lb": ["lumiBlock", "EventInfo_lumiBlock", "lumi_block", "lb"],
    "ntrk500": ["nTracksPt500Eta25", "EventInfo_nTracksPt500Eta25"],
}
CLASSES_428770 = ["all", "J12only", "MBTSbased", "j20MBTS", "sptrk4", "other"]
CLASSES_427929 = ["all427929", "MB427929", "sptrk4_427929"]
DISJOINT_428770 = ["J12only", "MBTSbased", "other"]


# ----------------------------------------------------------------------------------------------------------- import ścieżek
def add_paths():   # ub_comb (19.11): bez charmbkg (projekt B) — dstar_refit z jawnej kopii dstar_refit_aa.py
    for p in (AA_SRC,):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)


# ----------------------------------------------------------------------------------------------------------- nazwy kolumn
_SUF = ("_NOSYS", "")


def find_col(names, cands):
    """Pierwsza istniejąca kolumna spośród kandydatów (z sufiksem _NOSYS albo bez); None gdy brak."""
    s = set(names)
    for c in cands:
        for suf in _SUF:
            if c + suf in s:
                return c + suf
    return None


def find_trig_cols(names, chains=None):
    """{chain: {"pass": kol|None, "ps": kol|None}} dla konwencji: tp_<c> / ps_<c> (mult_common), trigPassed_<c> / trigPS_<c>
    (charmpp v2.x), EventInfo_trigPassed_<c>[_NOSYS] / EventInfo_trigPS_<c> (ROOT v2.0). chains=None → wszystkie znalezione."""
    names = list(names)
    pat = re.compile(r"^(?:EventInfo_)?(?:tp_|trigPassed_)(.+?)(?:_NOSYS)?$")
    found = {}
    for n in names:
        m = pat.match(n)
        if m:
            found.setdefault(m.group(1), {})["pass"] = n
    pat_ps = re.compile(r"^(?:EventInfo_)?(?:ps_|trigPS_)(.+?)(?:_NOSYS)?$")
    for n in names:
        m = pat_ps.match(n)
        if m:
            found.setdefault(m.group(1), {})["ps"] = n
    if chains is not None:
        found = {c: found.get(c, {}) for c in chains}
    return {c: {"pass": v.get("pass"), "ps": v.get("ps")} for c, v in found.items()}


# ----------------------------------------------------------------------------------------------------------- klasy
def classify(ev: dict, run: int | None = None) -> dict:
    """ev: dict tablic per zdarzenie: 'tp_<chain>' (0/1), opcjonalnie 'grl'. Zwraca {klasa: maska bool}.
    Brak kolumny łańcucha = False (z ostrzeżeniem w 'missing')."""
    n = len(next(iter(ev.values())))
    miss = []

    def tp(c):
        k = f"tp_{c}"
        if k not in ev:
            miss.append(c)
            return np.zeros(n, bool)
        return np.asarray(ev[k]) == 1          # h5_v20: MC / brak gałęzi = −1 → False (nie astype(bool)!)

    mbts_any = np.zeros(n, bool)
    for c in MBTS_CHAINS:
        mbts_any |= tp(c)
    j12 = tp(J12)
    grl = (np.asarray(ev["grl"]) == 1) if "grl" in ev else np.ones(n, bool)
    out = {}
    if run is None or int(run) != 427929:
        out.update({"all": np.ones(n, bool), "J12only": j12 & ~mbts_any, "MBTSbased": mbts_any,
                    "j20MBTS": tp("HLT_j20_L1MBTS_2"), "sptrk4": tp("HLT_mb_sptrk_pt4_L1MBTS_2"), "other": ~j12 & ~mbts_any})
    if run is None or int(run) == 427929:
        out.update({"all427929": grl, "MB427929": tp("HLT_mb_mbts_L1MBTS_1") & grl,
                    "sptrk4_427929": tp("HLT_mb_sptrk_pt4_L1MBTS_2") & grl})
    out["_missing"] = sorted(set(miss))
    return out


# ----------------------------------------------------------------------------------------------------------- złączenie
def ev_key(run, event):
    """Klucz zdarzenia danych (run < 2^20, event < 2^43). W MC NIE używać (event powtarza się w fd) — tam (plik, wpis)."""
    run = np.asarray(run, np.int64)
    event = np.asarray(event, np.int64)
    return (run << 43) | (event & ((1 << 43) - 1))


def join_event_table(run, event, T: dict, cols):
    """Dla wierszy (run, event) zwraca {kol: wartości z tabeli zdarzeń T} (T: 'run','event' + kolumny) i maskę 'found'.
    Brak zdarzenia w T → NaN/0, found=False (bramka G1: odsetek nieznalezionych < 0,1 %)."""
    k = ev_key(run, event)
    kt = ev_key(T["run"], T["event"])
    order = np.argsort(kt, kind="mergesort")
    ks = kt[order]
    if len(ks) > 1 and np.any(ks[1:] == ks[:-1]):
        raise ValueError("tabela zdarzeń: zduplikowane (run, event)")
    pos = np.searchsorted(ks, k)
    pos_c = np.clip(pos, 0, max(len(ks) - 1, 0))
    found = (len(ks) > 0) & (ks[pos_c] == k) if len(ks) else np.zeros(len(k), bool)
    idx = order[pos_c]
    out = {"found": found}
    for c in cols:
        v = np.asarray(T[c])
        fill = np.nan if np.issubdtype(v.dtype, np.floating) else 0
        o = np.full(len(k), fill, dtype=np.float64 if fill is np.nan else v.dtype)
        o[found] = v[idx[found]]
        out[c] = o
    return out


# ----------------------------------------------------------------------------------------------------------- fit pass/fail
def fit_eps(dm_pass_gev, dm_fail_gev, shape, interval=True, model_check=False, min_entries=50):
    """ε z jednoczesnego fitu Δm pass/fail (charmtag.dstar_refit.fit_passfail, AA_merge). Zwraca dict (eps, eps_err, eps_lo/hi,
    N, N_err, determined, valid, n_pass, n_fail)."""
    add_paths()
    import dstar_refit_aa as DR                                    # kopia AA_merge (PF/src na lxplus ma starą fit_passfail bez interval)
    hp, _ = np.histogram(np.asarray(dm_pass_gev) * 1e3, E_DM)
    hf, _ = np.histogram(np.asarray(dm_fail_gev) * 1e3, E_DM)
    base = {"n_pass": int(hp.sum()), "n_fail": int(hf.sum())}
    if hp.sum() + hf.sum() < min_entries:
        return {**base, "eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False,
                "reason": "za mało wpisów"}
    try:
        r = DR.fit_passfail(hp, hf, E_DM, shape, interval=interval, model_check=model_check)
    except Exception as ex:  # noqa: BLE001
        return {**base, "eps": np.nan, "eps_err": np.nan, "N": 0.0, "N_err": 0.0, "valid": False, "determined": False,
                "reason": str(ex)}
    g = lambda k: float(r.get(k, np.nan)) if r.get(k) is not None else np.nan  # noqa: E731
    return {**base, "eps": g("eps"), "eps_err": g("eps_err"), "eps_lo": g("eps_lo"), "eps_hi": g("eps_hi"),
            "N": g("N_tot"), "N_err": g("N_tot_err"), "valid": bool(r.get("valid", False)),
            "determined": bool(r.get("determined", False)), "reason": r.get("undetermined_reason")}


def eps_sym_err(f):
    """Symetryczny błąd ε: z profilu (średnia połówek), gdy zamknięty; inaczej HESSE."""
    lo, hi, e = f.get("eps_lo", np.nan), f.get("eps_hi", np.nan), f.get("eps", np.nan)
    if np.isfinite(lo) and np.isfinite(hi):
        return 0.5 * (hi - lo)
    return f.get("eps_err", np.nan)


# ----------------------------------------------------------------------------------------------------------- testy zgodności
def chi2_consistency(vals, errs):
    """χ² zgodności niezależnych pomiarów ze średnią ważoną. Zwraca (średnia, σ, χ², ndf, p)."""
    from scipy.stats import chi2 as C2
    v, e = np.asarray(vals, float), np.asarray(errs, float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    v, e = v[ok], e[ok]
    if len(v) < 2:
        return {"mean": float(v[0]) if len(v) else np.nan, "err": float(e[0]) if len(v) else np.nan, "chi2": np.nan,
                "ndf": 0, "p": np.nan, "n": int(len(v))}
    w = 1 / e ** 2
    mu = float((w * v).sum() / w.sum())
    x2 = float((w * (v - mu) ** 2).sum())
    ndf = len(v) - 1
    return {"mean": mu, "err": float(w.sum() ** -0.5), "chi2": x2, "ndf": ndf, "p": float(C2.sf(x2, ndf)), "n": int(len(v))}


def pull_subset(sf_sub, err_sub, sf_all, err_all):
    """Pull podzbioru względem całości (te same zdarzenia częściowo): Var(sub − all) ≈ σ_sub² − σ_all² (BLUE/podzbiór)."""
    d = sf_sub - sf_all
    v = err_sub ** 2 - err_all ** 2
    if not np.isfinite(v) or v <= 0:
        v = err_sub ** 2      # przy fluktuacji σ_sub < σ_all — zachowawczo
    return float(d / math.sqrt(v))


def pull_indep(a, ea, b, eb):
    return float((a - b) / math.hypot(ea, eb))


# ----------------------------------------------------------------------------------------------------------- czułość
def sigma_scale(S_ref, B_ref, sig_ref, S, B):
    """σ(ε) dla próbki (S, B) w oknie Δm przeskalowane z referencji: σ ∝ √(S+B)/S (fit pass/fail ~ zliczenia w oknie)."""
    if S <= 0:
        return float("inf")
    return float(sig_ref * (math.sqrt(S + B) / S) / (math.sqrt(S_ref + B_ref) / S_ref))


# ----------------------------------------------------------------------------------------------------------- ε_B (ML)
def _pol2_ml_pred(m, m0, sr_half, sb_min, lo, hi, nbins=50, fix_bc=None, return_bc=False):
    """Binned Poisson ML: N(m) ∝ exp(b·t + c·t²), t = (m − m0)/0,1, na SB |m−m0| ≥ sb_min w [lo, hi]; zwraca przewidywaną
    liczbę w SR |m−m0| < sr_half. Odporne przy małej statystyce (val_epsb.fit_exp_pol2_sb = LSQ na log N — obciążone przy
    kilku wpisach/bin; test toy 17.11)."""
    from scipy.optimize import minimize
    e = np.linspace(lo, hi, nbins + 1)
    h, _ = np.histogram(m, e)
    c = 0.5 * (e[1:] + e[:-1])
    t = (c - m0) / 0.1
    sb = np.abs(c - m0) >= sb_min
    if h[sb].sum() < 10:
        return (np.nan, int(h[sb].sum()), None) if return_bc else (np.nan, int(h[sb].sum()))
    w = np.diff(e)

    def lam(p):
        return np.exp(np.clip(p[0] + p[1] * t + p[2] * t * t, -50, 50)) * w / w.mean()

    def nll(p):
        lm = lam(p)[sb]
        return float(np.sum(lm - h[sb] * np.log(np.clip(lm, 1e-300, None))))
    if fix_bc is not None:          # kształt ustalony (np. z próbki „all”), wolna tylko normalizacja → analitycznie
        lm1 = lam(np.array([0.0, fix_bc[0], fix_bc[1]]))[sb]
        a = float(np.log(max(h[sb].sum(), 1e-9) / lm1.sum()))
        px = np.array([a, fix_bc[0], fix_bc[1]])
    else:
        p0 = np.array([np.log(max(h[sb].mean(), 0.5)), 0.0, 0.0])
        px = minimize(nll, p0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 4000}).x

    class _R:
        x = px
    r = _R
    # SR: całka na drobnej siatce
    x = np.linspace(m0 - sr_half, m0 + sr_half, 401)
    tt = (x - m0) / 0.1
    dens = np.exp(np.clip(r.x[0] + r.x[1] * tt + r.x[2] * tt * tt, -50, 50)) / w.mean()
    val = float(np.sum(0.5 * (dens[1:] + dens[:-1]) * np.diff(x)))
    return (val, int(h[sb].sum()), (float(r.x[1]), float(r.x[2]))) if return_bc else (val, int(h[sb].sum()))


def pol2_shapes(m, passed, m0=M0_D0, sr_half=0.055, sb_min=0.075, lo=1.75, hi=2.00, nbins=50):
    """Kształty (b, c) exp(pol2) z próbki odniesienia (klasa „all”): dla wszystkich i dla pass."""
    m = np.asarray(m, float)
    p = np.asarray(passed, bool)
    _, _, bc_all = _pol2_ml_pred(m, m0, sr_half, sb_min, lo, hi, nbins, return_bc=True)
    _, _, bc_pass = _pol2_ml_pred(m[p], m0, sr_half, sb_min, lo, hi, nbins, return_bc=True)
    return {"all": bc_all, "pass": bc_pass}


def eps_b_pol2_ml(m, passed, m0=M0_D0, sr_half=0.055, sb_min=0.075, lo=1.75, hi=2.00, nbins=50, shapes=None):
    """ε_B(SR) = N_pass_pred / N_all_pred (oba z ML exp(pol2) na SB). Błąd: √(1/n_sb,pass − 1/n_sb,all)·ε (binomialny
    na SB — przybliżenie ekstrapolacji przy tym samym kształcie)."""
    m = np.asarray(m, float)
    p = np.asarray(passed, bool)
    sa, sp = (shapes or {}).get("all"), (shapes or {}).get("pass")
    na, nsa = _pol2_ml_pred(m, m0, sr_half, sb_min, lo, hi, nbins, fix_bc=sa)
    npa, nsp = _pol2_ml_pred(m[p], m0, sr_half, sb_min, lo, hi, nbins, fix_bc=sp)
    if not (np.isfinite(na) and np.isfinite(npa)) or na <= 0:
        return {"eps": np.nan, "err": np.nan, "n_sb": nsa, "n_sb_pass": nsp}
    eps = npa / na
    err = eps * math.sqrt(max(1.0 / max(nsp, 1) - 1.0 / max(nsa, 1), 0.0))
    return {"eps": float(eps), "err": float(err), "n_sb": nsa, "n_sb_pass": nsp}
