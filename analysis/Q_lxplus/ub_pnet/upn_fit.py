#!/usr/bin/env python3
"""UBPNET (19.11.2026, projekt A; pary D–D w danych zatwierdzone 18.11) — μ par D–D w 428770 (ta sama odsłonięta próbka zdarzeń co 05.11)
z kategoriami par z sieci PAIRNET v3 zamiast nominalnych kategorii nóg.

Schemat główny: pnetd (sieć par Deep Sets na 18 cechach H obu nóg, tło treningowe = pary SB×SB danych 428770, OOF po foldach zdarzeń;
SBTRAIN v3: σ_μ/σ_μ(nomcat) = 0,791 ± 0,009 przy tle z danych — najlepszy z 4 schematów); kontrolnie xgbd (ten sam trening, XGB) i base_sum.
K = 6 kategorii, progi z sygnału MC przy nominalnej łącznej wydajności (każda kategoria = 1/6 sygnału MC obszaru nominalnego).

Pary danych: nogi z h5 v2.0 (pn_extract --arb, cały zakres mas fitu; części z upn_extract.py), tylko zdarzenia próbki nominalnej
(event_number z h5 v1.1 data_428770_*_none, 39 plików — upn_evset.py), te same funkcje par/foldów co trening (PT.build_pairs, ek z smp = 2).

Model (wspólny dla wszystkich schematów, także porównawczego „nomcat” liczonego na tych samych parach), kanał × kategoria k, płaszczyzna 25×25:
  λ_k = μ·s_k·S1_k⊗S2_k + n_r1,k·S1_k⊗F2_k + n_r2,k·F1_k⊗S2_k + n_cc,k·F1_k⊗F2_k
  S = brzegowe szablony mas sygnału MC (prawda 0/1, z refleksją K↔π) per kategoria; F = exp(wielomian 2. st.) per noga i kategoria (wolne);
  n_r1, n_r2, n_cc ≥ 0 wolne (cała płaszczyzna z pasmami w ficie). Noga 1/2 jak pn_fit.prepare (DpD0: 1 = D⁺; inaczej 1 = wyższe pT).
  s_k = n_pair(kanał, pair_inputs odsłonięcia) × E_k:  nomcat_fact: E_ij = Y_i·Y_j·mult (Y = 0,30/0,15/0,20 jak nominał — odpowiednik
  modelu nominalnego bez kategorii z tagiem D*), nomcat_mc: 0,65²·(udział MC komórki), schemat sieci: 0,65²/K × SF_rel(kanał),
  SF_rel = [SF schematu / SF nomcat] z par hybrydowych (noga D⁺: HYBSF v3, noga D⁰: HYBDS v3 — tag D*), średnia pT z wagami par MC.
Wynik: μ na skali odsłonięcia 05.11 (μ_nom = 1,19) i ×1,50/μ_nom (zapis bieżący, N_D* = 4436).
Kontrole: (1) nominał charmtag odtworzony z pairs_full.npz (fit2d_unblind.json, tol 2e-3) — STOP gdy nie; (2) nomcat na parach v2.0 ≈ μ_nom;
(3) Asimov: fit odtwarza μ0, σ oczekiwane; (4) bootstrap Poissona po zdarzeniach: σ_μ i σ(μ_sieć − μ_nomcat_mc); (5) syst: SF_rel (norm.), nachylenie
SF między połówkami kategorii (HYBSF kategorie / HYBDS „góra”), kształt tła poly3.
Użycie: upn_fit.py --mc 'pn_*_nd_*.npz' --data 'pnf_data_428770_*.npz' --evset evset.npz --modelsd pn_models_sbtrain.pt --xgbd DIR
        --hybsf hybsf.json --hybds hybds.json --pairs-dir DIR [--legs legs_data.parquet] --out DIR   |   --selftest --no-torch --out DIR
Tryb 427929 (19.11, UBPNET427; domyślne wartości = bieg 428770 bez zmian): --sample data_427929 (klucze SF HYBSF/HYBDS), --npair-scale f_D0D0,f_DpD0,f_DpDp
(normalizacja z prawdy mb_nd; dodatkowo nomcat_mc_f1 bez mnożników), --scale-now 1.0 (pair_inputs UB427 już przy N_D* = 2046), --ref-sb 'pn_data_428770_*.npz'
(test transferu tła SBTRAIN), bez --evset (h5 po filtrze = tylko klasa MB)."""
import argparse
import glob
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PN = os.environ.get("UPN_PN") or os.path.join(os.path.dirname(HERE), "pairnet")
for p in (PN, os.path.join(os.path.dirname(HERE), "..", "PAIRNET")):          # lxplus: Q/pairnet; lokalnie: charm_work/PAIRNET
    if os.path.isfile(os.path.join(p, "pn_train.py")):
        sys.path.insert(0, os.path.abspath(p))
        break
sys.path.insert(0, HERE)
import pn_train as PT  # noqa: E402
import pn_fit as PF  # noqa: E402
import pn_hybsf as HS  # noqa: E402

CH = PT.CH
SPEC = {0: (0, 0), 1: (1, 0), 2: (1, 1)}                 # kanał → (gatunek nogi 1, nogi 2) po pn_fit.prepare (DpD0: 1 = D⁺)
WPINC = (0.30, 0.15, 0.20)                               # nominalne Y nóg (kategorie c0/c1/c2 = WP30, WP30–45, WP45–65)
E_NOM = 0.65 ** 2
SAME = {0: True, 1: False, 2: True}
NOMCODES = {True: [0, 1, 2, 4, 5, 8], False: list(range(9))}
PTB = ("3–5", "5–8", "8–30")
TOL_MU = 2e-3
SCALE_NOW = 1.50                                         # zapis bieżący μ (N_D* = 4436); skala = SCALE_NOW / μ_nom(05.11)
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f}s]", *a, flush=True)


def need_unblind():
    if os.environ.get("CHARM_UNBLIND") != "1":
        raise SystemExit("UBPNET: wymaga CHARM_UNBLIND=1 (pary D–D w oknie sygnału danych — zatwierdzone 18.11)")


# ----------------------------------------------------------------------------------------------------------- model kanału
def ugrid(sp):
    c = 0.5 * (PF.EDGES[sp][1:] + PF.EDGES[sp][:-1])
    return 2 * (c - c[0]) / (c[-1] - c[0]) - 1


def fshape(u, th):
    z = th[0] * u + th[1] * u ** 2 + (th[2] * u ** 3 if len(th) > 2 else 0.0)
    z = np.exp(z - z.max())
    return z / z.sum()


class ChanFit:
    """Model jednego kanału (wszystkie kategorie schematu). n: (C, 25, 25); s: (C,) sygnał przy μ = 1; S1, S2: (C, 25) szablony nóg."""

    def __init__(self, n, s, S1, S2, sp1, sp2, npoly=2):
        self.n = np.asarray(n, float); self.s = np.asarray(s, float); self.S1 = np.asarray(S1, float); self.S2 = np.asarray(S2, float)
        self.C = len(self.s); self.npoly = npoly
        self.u1, self.u2 = ugrid(sp1), ugrid(sp2)
        self.np_c = 2 * npoly + 3
        self.N = self.n.sum(axis=(1, 2))
        self.SS = self.s[:, None, None] * self.S1[:, :, None] * self.S2[:, None, :]
        self.lgn = np.where(self.n > 0, self.n * np.log(np.maximum(self.n, 1e-300)), 0.0)
        self.best = None

    def unpack(self, x):
        mu = x[0]
        P = np.asarray(x[1:]).reshape(self.C, self.np_c)
        k = self.npoly
        return mu, P[:, :k], P[:, k:2 * k], P[:, 2 * k], P[:, 2 * k + 1], P[:, 2 * k + 2]

    def lam(self, x):
        mu, t1, t2, r1, r2, cc = self.unpack(x)
        F1 = np.array([fshape(self.u1, t) for t in t1]); F2 = np.array([fshape(self.u2, t) for t in t2])
        L = mu * self.SS + r1[:, None, None] * self.S1[:, :, None] * F2[:, None, :] + r2[:, None, None] * F1[:, :, None] * self.S2[:, None, :] \
            + cc[:, None, None] * F1[:, :, None] * F2[:, None, :]
        return L

    def nll(self, x):
        L = np.maximum(self.lam(x), 1e-12)
        return float(np.sum(L - self.n * np.log(L)))

    def deviance(self, x):
        L = np.maximum(self.lam(x), 1e-12)
        return float(2 * np.sum(L - self.n + np.where(self.n > 0, self.n * np.log(np.maximum(self.n, 1e-300) / L), 0.0)))

    def x0(self, mu=1.0):
        x = [mu]
        for c in range(self.C):
            N = self.N[c]
            x += [0.0] * (2 * self.npoly) + [0.15 * N, 0.15 * N, 0.7 * N]
        return np.array(x, float)

    def fit(self, fix_mu=None, start=None):
        from iminuit import Minuit
        x = np.array(self.x0() if start is None else start, float)
        if fix_mu is not None:
            x[0] = fix_mu
        m = Minuit(self.nll, x)
        m.errordef = 0.5
        m.strategy = 1
        lim = [(-20.0, 50.0)]
        for c in range(self.C):
            lim += [(-12.0, 12.0)] * (2 * self.npoly) + [(0.0, None)] * 3
        m.limits = lim
        for c in range(self.C):                                    # pusta kategoria: tło = 0, kształty zamrożone
            if self.N[c] <= 0:
                b = 1 + c * self.np_c
                for j in range(self.np_c):
                    m.values[b + j] = 0.0
                    m.fixed[b + j] = True
        if fix_mu is not None:
            m.fixed[0] = True
        m.migrad(ncall=200000)
        if not m.valid:
            m.migrad(ncall=200000)
        r = {"fval": float(m.fval), "x": np.array(m.values), "valid": bool(m.valid), "mu": float(m.values[0])}
        if fix_mu is None:
            try:
                m.hesse()
                r["mu_err"] = float(m.errors[0])
            except Exception:  # noqa: BLE001
                r["mu_err"] = float("nan")
            r["dev"] = self.deviance(r["x"]); r["ndf"] = int((self.n.sum(axis=(1, 2)) > 0).sum() * 625 - m.nfit)
            self.best = r
        return r


def minimize_1d(g, lo, hi, xatol=2e-3):
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(g, bounds=(lo, hi), method="bounded", options={"xatol": xatol, "maxiter": 60})
    return float(r.x), float(r.fun)


def root_1d(h, a, b, tol=4e-3):
    ha, hb = h(a), h(b)
    if not (np.isfinite(ha) and np.isfinite(hb)) or ha * hb > 0:
        return None
    for _ in range(40):
        m = 0.5 * (a + b); hm = h(m)
        if ha * hm <= 0:
            b, hb = m, hm
        else:
            a, ha = m, hm
        if abs(b - a) < tol:
            break
    return 0.5 * (a + b)


def combined_fit(fits, interval=True, zero=True):
    """fits: {kanał: ChanFit}. Wspólne μ = minimum sumy profili kanałów (jak ub_jk.replica_fit), 68 % z Δ(−lnL) = 0,5, Z z q0."""
    out = {}
    for ch, f in fits.items():
        r = f.fit()
        out[ch] = {"mu": r["mu"], "err": r.get("mu_err"), "valid": r["valid"], "dev": r["dev"], "ndf": r["ndf"]}
    cache = {}

    def g(mu):
        mu = round(float(mu), 5)
        if mu not in cache:
            s = 0.0
            for f in fits.values():
                st = f.best["x"].copy(); st[0] = mu
                s += f.fit(fix_mu=mu, start=st)["fval"] - f.best["fval"]
            cache[mu] = s
        return cache[mu]
    mus = np.array([out[c]["mu"] for c in fits]); ers = np.array([out[c]["err"] or 1.0 for c in fits])
    ers = np.where(np.isfinite(ers) & (ers > 0), ers, 1.0)
    s0 = float(1 / math.sqrt(np.sum(1 / ers ** 2)))
    wm = float(np.sum(mus / ers ** 2) / np.sum(1 / ers ** 2))
    lo, hi = wm - 4 * s0, wm + 4 * s0
    muc, gmin = minimize_1d(g, lo, hi)
    c = {"mu": muc, "gmin": gmin, "wmean": wm, "wmean_err": s0}
    if zero:
        q0 = 2 * (g(0.0) - gmin) if muc > 0 else 0.0
        c["q0"] = q0; c["Z"] = math.sqrt(max(q0, 0.0))
    if interval:
        h = lambda m: g(m) - gmin - 0.5  # noqa: E731

        def side(sgn):
            step = 1.5 * s0
            for _ in range(5):
                x = muc + sgn * step
                if h(x) > 0:
                    return root_1d(h, *sorted((muc, x)))
                step *= 2
            return None
        up, dn = side(+1), side(-1)
        c["mu_68"] = [dn, up]
        c["sigma_68"] = 0.5 * (up - dn) if (up is not None and dn is not None) else None
    out["combined"] = c
    return out


# ----------------------------------------------------------------------------------------------------------- dane i MC
def load_data(paths, evset):
    L = PT.load_legs(paths)
    L["smp"][:] = 2                                                # jak pn_sbtrain.load_all (DSMP data_428770 = 2) → te same foldy OOF
    L["ek"] = PT.mix64(L["uid"], 1) ^ PT.mix64(L["ftag"].astype(np.uint64) * np.uint64(4) + L["smp"].astype(np.uint64), 2)
    d = {"legs_all": int(len(L["m"])), "events_all": int(len(np.unique(L["ek"])))}
    if evset is not None:
        keep = np.isin(L["event"].astype(np.uint64), evset)
        L = {k: v[keep] for k, v in L.items()}
    d.update({"legs_nominal_events": int(len(L["m"])), "events_nominal": int(len(np.unique(L["ek"])))})
    return L, d


def cats_of(D, S, th, c, schemes, K):
    """{schemat: kategoria pary (−1 = poza)} w kanale c (tylko pary z obiema nogami w zakresie mas)."""
    base = (D["chan"] == c) & D["ok"]
    same = SAME[c]
    codes = NOMCODES[same]
    lut = np.full(9, -1, np.int16); lut[codes] = np.arange(len(codes))
    nc = np.where(base & (D["nomcat"] >= 0), lut[np.clip(D["nomcat"], 0, 8)], -1).astype(np.int16)
    out = {"nomcat_fact": (nc, len(codes)), "nomcat_mc": (nc, len(codes))}
    for nm in schemes:
        if nm not in S or nm not in th[c]:
            continue
        cc = np.full(len(base), -1, np.int16)
        for k in range(K - 1, -1, -1):
            cc[S[nm] >= th[c][nm][k]] = k
        out[nm] = (np.where(base, cc, -1).astype(np.int16), K)
    return out


def sb_transfer(D, S, Dr, Sr, th, schemes, K):
    """Rozkład kategorii schematów dla par z obiema nogami w pasmach bocznych (SB×SB): próbka vs referencja. χ² jednorodności (K+1 klas: K kategorii
    + „poza schematem”), udziały w górnej połowie kategorii. Tło fitu jest wolne per kategoria, więc różnica nie obciąża μ — mierzy tylko, czy
    sieć rozdziela tło w tej próbce tak samo jak w próbce treningowej (zysk σ_μ)."""
    out = {}
    for c in range(3):
        cs, cr = cats_of(D, S, th, c, schemes, K), cats_of(Dr, Sr, th, c, schemes, K)
        qs = (D["chan"] == c) & D["ok"] & D["sb1"] & D["sb2"]; qr = (Dr["chan"] == c) & Dr["ok"] & Dr["sb1"] & Dr["sb2"]
        res = {}
        for nm in ["nomcat_mc"] + [x for x in schemes if x in cs]:
            (a_, na), (b_, nb) = cs[nm], cr[nm]
            hs = np.bincount(a_[qs] + 1, minlength=na + 1).astype(float); hr = np.bincount(b_[qr] + 1, minlength=nb + 1).astype(float)
            ns, nr = hs.sum(), hr.sum()
            if ns <= 0 or nr <= 0:
                continue
            p = (hs + hr) / (ns + nr); es, er = ns * p, nr * p
            ok = p > 0
            chi2 = float(np.sum((hs - es)[ok] ** 2 / es[ok]) + np.sum((hr - er)[ok] ** 2 / er[ok]))
            ndf = int(ok.sum() - 1)
            top = (lambda h, n_: float(h[1:1 + (n_ // 2)].sum() / max(h[1:].sum(), 1)))
            res[nm] = {"n": int(ns), "n_ref": int(nr), "frac": (hs / ns).tolist(), "frac_ref": (hr / nr).tolist(), "chi2": chi2, "ndf": ndf,
                       "chi2_ndf": chi2 / max(ndf, 1), "top_half": top(hs, na), "top_half_ref": top(hr, nb),
                       "in_scheme": float(hs[1:].sum() / ns), "in_scheme_ref": float(hr[1:].sum() / nr)}
        out[CH[c]] = res
    return out


def hist_cats(D, cat, ncat, w=None):
    q = cat >= 0
    flat = (cat[q].astype(np.int64) * 25 + D["b1"][q]) * 25 + D["b2"][q]
    h = np.bincount(flat, weights=None if w is None else w[q], minlength=ncat * 625)
    return h.reshape(ncat, 25, 25).astype(float)


def mc_templates(Dm, cat, ncat, c, min_n=50):
    base = (Dm["chan"] == c) & Dm["ok"]
    S1i = np.bincount(Dm["b1"][base], minlength=25).astype(float); S2i = np.bincount(Dm["b2"][base], minlength=25).astype(float)
    S1, S2, nmc = np.zeros((ncat, 25)), np.zeros((ncat, 25)), np.zeros(ncat, int)
    for k in range(ncat):
        q = base & (cat == k)
        nmc[k] = int(q.sum())
        a = np.bincount(Dm["b1"][q], minlength=25).astype(float) if nmc[k] >= min_n else S1i
        b = np.bincount(Dm["b2"][q], minlength=25).astype(float) if nmc[k] >= min_n else S2i
        S1[k], S2[k] = a / a.sum(), b / b.sum()
    return S1, S2, nmc


# ----------------------------------------------------------------------------------------------------------- SF (pary hybrydowe)
SAMPLE = "data_428770"                                   # próbka danych dla SF (HYBSF/HYBDS); --sample data_427929 (19.11)


def _hybsf_leg(J, ch, leg, nm):
    k = f"{SAMPLE}|{ch}|{leg}|{nm}"
    if J is None or k not in J.get("sf", {}):
        return None
    r = J["sf"][k]
    sf = np.asarray(r["eff_data"], float) / np.maximum(np.asarray(r["eff_mc"], float), 1e-9)   # per bin pT ("SF" = [łączne, ±stat, ±syst_bkg])
    err = sf * np.sqrt((np.asarray(r["err_data"]) / np.maximum(r["eff_data"], 1e-9)) ** 2 + (np.asarray(r["err_mc"]) / np.maximum(r["eff_mc"], 1e-9)) ** 2)
    top, bot = [], []
    for b in PTB:
        cats = r.get("cats", {}).get(b, {})
        fd = {int(c): v["f_data"] for c, v in cats.items() if str(c).isdigit()}; fm = {int(c): v["f_mc"] for c, v in cats.items() if str(c).isdigit()}
        t = [c for c in fd if 1 <= c <= 3]; o = [c for c in fd if c >= 4]
        top.append(sum(fd[c] for c in t) / max(sum(fm[c] for c in t), 1e-9) if t else np.nan)
        bot.append(sum(fd[c] for c in o) / max(sum(fm[c] for c in o), 1e-9) if o else np.nan)
    return {"SF": sf, "err": err, "pi": np.asarray(r.get("pi", [1 / 3] * 3), float), "top": np.asarray(top), "bot": np.asarray(bot)}


def _hybds_leg(J, ch, nm, sel):
    if J is None:
        return None
    rows = [r for r in J.get("rows", []) if r["sample"] == SAMPLE and r["sel"] == sel and r["chan"] == ch and r["scheme"] == nm]
    if len(rows) != 3:
        return None
    rows = sorted(rows, key=lambda r: PTB.index(r["bin"]))
    return {"SF": np.array([r["SF"] for r in rows], float), "err": np.array([r["SF_err"] for r in rows], float)}


def sf_rel(hybsf, hybds, ch, nm, rho):
    """SF_rel = SF(schemat)/SF(nomcat) dla pary (iloczyn nóg), błąd (korelacja schemat–nomcat z MC: rho), nachylenie (góra/dół).
    Noga D⁺: HYBSF; noga D⁰: HYBDS (tag D*). Wagi pT: „pi” z HYBSF (pary MC)."""
    legs = {"D0D0": ["D0", "D0"], "DpD0": ["D+", "D0"], "DpDp": ["D+", "D+"]}[ch]
    out = {"legs": {}, "available": True}
    tot, var, ttop, tbot = 1.0, 0.0, 1.0, 1.0
    for lg in sorted(set(legs)):
        npow = legs.count(lg)
        if lg == "D+":
            s, n = _hybsf_leg(hybsf, ch, "D+", nm), _hybsf_leg(hybsf, ch, "D+", "nomcat")
            pi = s["pi"] if s is not None else None
        else:
            s, n = _hybds_leg(hybds, ch, nm, "schemat"), _hybds_leg(hybds, ch, "nomcat", "schemat")
            g = _hybsf_leg(hybsf, ch, "D0", nm)
            pi = g["pi"] if g is not None else np.full(3, 1 / 3)
            gs = _hybds_leg(hybds, ch, nm, "góra")
            if s is not None and gs is not None:
                s["top"] = gs["SF"]; s["bot"] = 2 * s["SF"] - gs["SF"]
        if s is None or n is None:
            out["available"] = False
            out["legs"][lg] = None
            continue
        pi = np.asarray(pi, float) / np.sum(pi)
        rb = s["SF"] / n["SF"]
        es, en = s["err"] / s["SF"], n["err"] / n["SF"]
        erb = rb * np.sqrt(np.maximum(es ** 2 + en ** 2 - 2 * rho * es * en, 0.0))
        r = float(np.sum(pi * rb)); e = float(np.sqrt(np.sum((pi * erb) ** 2)))
        st = float(np.nansum(pi * s["top"]) / np.sum(pi * s["SF"])) if np.isfinite(s["top"]).all() else 1.0
        sb = float(np.nansum(pi * s["bot"]) / np.sum(pi * s["SF"])) if np.isfinite(s["bot"]).all() else 1.0
        out["legs"][lg] = {"rel": r, "err": e, "per_bin": rb.tolist(), "pi": pi.tolist(), "tilt_top": st, "tilt_bot": sb, "power": npow}
        tot *= r ** npow; var += (npow * e / r) ** 2; ttop *= st ** npow; tbot *= sb ** npow
    out.update({"rel": tot if out["available"] else 1.0, "err": tot * math.sqrt(var) if out["available"] else 0.0,
                "tilt_top": ttop if out["available"] else 1.0, "tilt_bot": tbot if out["available"] else 1.0})
    return out


# ----------------------------------------------------------------------------------------------------------- kontrola nominału
def gate_nominal(pairs_dir):
    """Fit nominalny charmtag na pairs_full.npz (jak ub_jk prep) — μ per kanał = fit2d_unblind.json; wspólne μ (replica_fit)."""
    import charmtag
    from charmtag.pairs import CHANNELS
    from charmtag.pairs.model import Model, load_inputs
    from charmtag.pairs.fit2d import fit_channel
    I = load_inputs(os.path.join(pairs_dir, "pair_inputs"))
    H0 = dict(np.load(os.path.join(pairs_dir, "pairs_full.npz")))
    ref = json.load(open(os.path.join(pairs_dir, "fit2d_unblind.json")))
    M = Model(I)
    rep = {"charmtag_file": charmtag.__file__}
    for ch in CHANNELS:
        r, _, _ = fit_channel(M, ch, H0[ch], np.ones(H0[ch].shape[1:], bool), scan=False)
        rep[ch] = {"mu": float(r["mu_hat"]), "ref": float(ref["channels"][ch]["mu_hat"]), "err": float(r["mu_err_hesse"])}
    rep["ok"] = bool(all(abs(rep[ch]["mu"] - rep[ch]["ref"]) < TOL_MU for ch in CHANNELS))
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), "ub_jk"))
        import ubjk_run as U
        c = U.replica_fit(M, H0, list(CHANNELS), combined=True, interval=True)["combined"]
        rep["combined"] = {"mu": c["mu"], "mu_68": c.get("mu_68"), "Z": c["Z"]}
    except Exception as ex:  # noqa: BLE001
        rep["combined_error"] = str(ex)
    rep["n_pair"] = {k: float(v) for k, v in I.n_pair.items()}
    return rep


# ----------------------------------------------------------------------------------------------------------- bootstrap
_G = {}


def _boot(seed):
    G = _G
    ek_inv, nek = G["ek_inv"], G["nek"]
    w_ev = np.random.default_rng(1000 + seed).poisson(1.0, nek).astype(float)
    w = w_ev[ek_inv]
    res = {}
    for nm in G["schemes"]:
        fits = {}
        for c in range(3):
            cat, ncat = G["cats"][c][nm]
            Spec = G["spec"][(c, nm)]
            fits[CH[c]] = ChanFit(hist_cats(G["D"], cat, ncat, w), Spec["s"], Spec["S1"], Spec["S2"], *SPEC[c])
        try:
            r = combined_fit(fits, interval=False, zero=False)
            res[nm] = r["combined"]["mu"]
        except Exception:  # noqa: BLE001
            res[nm] = float("nan")
    return res


# ----------------------------------------------------------------------------------------------------------- selftest
def selftest_inputs(td):
    mcp = []
    for i, p in enumerate(PT.selftest_legs(n_ev=20000, seed=3)):
        fn = os.path.join(td, f"mc_{i}.npz"); np.savez(fn, **p); mcp.append(fn)
    d = PT.selftest_legs(n_ev=40000, seed=11)[1]                      # „dane”: zdarzenia minbias-podobne (sygnał 15 % nóg)
    rs = np.random.default_rng(12)
    for s in (0, 1):                                                  # tło: gładkie (malejące) w całym zakresie fitu gatunku
        q = (d["sp"] == s) & (d["truth"] != 0)
        lo, hi = PF.EDGES[s][0], PF.EDGES[s][-1]
        d["m"][q] = (lo + (hi - lo) * (1 - np.sqrt(1 - rs.random(int(q.sum())) * 0.75)) / 0.5).astype(np.float32)
    d["meta"] = np.array(json.dumps({"tag": "data"}))
    d["run"] = np.full(len(d["m"]), 428770, np.int64)
    fn = os.path.join(td, "pnf_data_0.npz"); np.savez(fn, **d)
    ev = np.unique(d["event"].astype(np.uint64))
    evset = ev[np.random.default_rng(2).random(len(ev)) < 0.9]        # próbka nominalna = podzbiór zdarzeń
    return mcp, [fn], evset


# ----------------------------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mc"); ap.add_argument("--data"); ap.add_argument("--evset"); ap.add_argument("--modelsd"); ap.add_argument("--xgbd")
    ap.add_argument("--hybsf"); ap.add_argument("--hybds"); ap.add_argument("--pairs-dir"); ap.add_argument("--legs")
    ap.add_argument("--out", required=True); ap.add_argument("--K", type=int, default=6); ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--nboot", type=int, default=40); ap.add_argument("--schemes", default="pnetd,xgbd,base_sum")
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--no-torch", action="store_true")
    ap.add_argument("--sample", default="data_428770", help="próbka danych: SF z HYBSF/HYBDS dla tej próbki, etykiety raportu")
    ap.add_argument("--npair-scale", default="", help="mnożniki n_pair per kanał D0D0,DpD0,DpDp (np. normalizacja z prawdy mb_nd f_ch); "
                    "dodatkowo fit nomcat_mc_f1 bez mnożników (kontrola)")
    ap.add_argument("--scale-now", type=float, default=None, help="μ_zapis = μ × scale_now (domyślnie 1,50/μ_nom — 428770)")
    ap.add_argument("--ref-sb", default="", help="części pasm bocznych innej próbki (np. pn_data_428770_*.npz) — test transferu tła SBTRAIN")
    a = ap.parse_args(argv)
    global SAMPLE
    SAMPLE = a.sample
    need_unblind()
    os.makedirs(a.out, exist_ok=True)
    rep = {"config": {"K": a.K, "nboot": a.nboot, "selftest": bool(a.selftest), "schemes": a.schemes, "pn_dir": os.path.dirname(PT.__file__),
                      "sample": a.sample, "npair_scale": a.npair_scale, "scale_now": a.scale_now, "ref_sb": a.ref_sb}}
    fsc = [float(x) for x in a.npair_scale.split(",")] if a.npair_scale else None
    if fsc is not None and len(fsc) != 3:
        raise SystemExit("--npair-scale: 3 liczby (D0D0,DpD0,DpDp)")
    # ---- (1) nominał
    if a.pairs_dir:
        g = gate_nominal(a.pairs_dir)
        rep["gate_nominal"] = g
        log("nominał charmtag:", {ch: round(g[ch]["mu"], 4) for ch in CH}, "wspólne", g.get("combined"), "ok", g["ok"], g["charmtag_file"])
        if not g["ok"]:
            json.dump(rep, open(os.path.join(a.out, "ubpnet.json"), "w"), indent=1, default=float)
            raise SystemExit(f"STOP: nominał nie odtworzony (tol {TOL_MU}): {g}")
        n_pair = g["n_pair"]; mu_nom = (g.get("combined") or {}).get("mu") or 1.189
    else:
        if not a.selftest:
            raise SystemExit("STOP: --pairs-dir wymagane (normalizacja n_pair i kontrola nominału)")
        n_pair, mu_nom = None, 1.0
    # ---- wejścia
    if a.selftest:
        import tempfile
        td = tempfile.mkdtemp(prefix="upn_self_")
        mcp, dpaths, evset = selftest_inputs(td)
        a.no_torch = True
    else:
        mcp, dpaths = sorted(glob.glob(a.mc)), sorted(glob.glob(a.data))
        evset = np.load(a.evset)["event"].astype(np.uint64) if a.evset else None
        if not mcp or not dpaths:
            raise SystemExit(f"STOP: brak części MC ({len(mcp)}) lub danych ({len(dpaths)})")
    a.models1 = None
    if a.no_torch:
        a.modelsd = None
    U = HS.setup(argparse.Namespace(no_torch=a.no_torch, models1=None, modelsd=a.modelsd, xgbd=a.xgbd, threads=a.threads, K=a.K), mcp)
    M, Lm, Ps, th, enom = U["M"], U["Lm"], U["Ps"], U["th"], U["enom"]
    want = [s for s in a.schemes.split(",") if s]
    S0 = M.score(Lm, Ps)
    schemes = [s for s in want if s in S0]
    rep["config"]["schemes_used"] = schemes
    Dm = PF.prepare(Lm, Ps, S0)
    log("MC: par sygnału", len(Ps["a"]), "schematy", schemes, "ε_nom", enom)
    hybsf = json.load(open(a.hybsf)) if a.hybsf and os.path.exists(a.hybsf) else None
    hybds = json.load(open(a.hybds)) if a.hybds and os.path.exists(a.hybds) else None
    # ---- dane
    L, dd = load_data(dpaths, evset)
    A, B, cnt = PT.build_pairs(L)
    P = PT.pair_table(L, A, B)
    S = M.score(L, P)
    D = PF.prepare(L, P, S)
    nin = {CH[c]: int(((D["chan"] == c) & D["ok"]).sum()) for c in range(3)}
    rep["data"] = {**dd, **cnt, "pairs_in_range": nin, "pairs_total": int(len(A))}
    log("dane:", rep["data"])
    if a.legs and os.path.exists(a.legs):                          # diagnostyka: pokrycie zdarzeń nóg nominału (v1.1) przez nogi v2.0
        import pandas as pd
        Ln = pd.read_parquet(a.legs, columns=["entry"])
        en = np.unique(Ln.entry.to_numpy().astype(np.uint64)); ev2 = np.unique(L["event"].astype(np.uint64))
        rep["data"]["nominal_leg_events_in_v20_legs"] = float(np.isin(en, ev2).mean())
    # ---- test transferu tła SBTRAIN: rozkład kategorii par SB×SB tej próbki vs --ref-sb (np. 428770, na którego pasmach trenowano)
    if a.ref_sb:
        rp = dpaths if a.ref_sb == "self" else sorted(glob.glob(a.ref_sb))     # "self" = selftest (próbka sama ze sobą → χ² = 0)
        if rp:
            Lr, _ = load_data(rp, None)
            Ar, Br, _ = PT.build_pairs(Lr); Pr = PT.pair_table(Lr, Ar, Br); Sr = M.score(Lr, Pr); Dr = PF.prepare(Lr, Pr, Sr)
            rep["sb_transfer"] = sb_transfer(D, S, Dr, Sr, th, schemes, a.K)
            log("transfer SB×SB:", {k: {nm: round(v[nm]["chi2_ndf"], 2) for nm in v} for k, v in rep["sb_transfer"].items()})
    # ---- specyfikacje (szablony, sygnał oczekiwany) per kanał × schemat
    names = ["nomcat_fact", "nomcat_mc"] + (["nomcat_mc_f1"] if fsc is not None else []) + schemes
    cats, catm, spec, info = {}, {}, {}, {}
    for c in range(3):
        cats[c] = cats_of(D, S, th, c, schemes, a.K)
        catm[c] = cats_of(Dm, S0, th, c, schemes, a.K)
        if fsc is not None:
            cats[c]["nomcat_mc_f1"] = cats[c]["nomcat_mc"]; catm[c]["nomcat_mc_f1"] = catm[c]["nomcat_mc"]
        base = (Dm["chan"] == c) & Dm["ok"]
        nreg = base & (Dm["nomcat"] >= 0)
        for nm in names:
            cm, ncat = catm[c][nm]
            S1, S2, nmc = mc_templates(Dm, cm, ncat, c)
            if nm == "nomcat_fact":
                codes = NOMCODES[SAME[c]]
                E = np.array([WPINC[x // 3] * WPINC[x % 3] * (2.0 if SAME[c] and x // 3 != x % 3 else 1.0) for x in codes])
                sfr = None
            elif nm in ("nomcat_mc", "nomcat_mc_f1"):
                fr = np.array([(nreg & (cm == k)).sum() for k in range(ncat)], float) / max(nreg.sum(), 1)
                E = E_NOM * fr; sfr = None
            else:
                es = float((base & (cm >= 0)).mean()); en_ = float(nreg.mean()); esn = float((base & (cm >= 0) & nreg).mean())
                rho = (esn - es * en_) / math.sqrt(max(es * (1 - es) * en_ * (1 - en_), 1e-12))
                sfr = sf_rel(hybsf, hybds, CH[c], nm, rho); sfr["rho_mc"] = rho
                frk = np.array([(base & (cm == k)).sum() for k in range(ncat)], float) / max(nreg.sum(), 1)   # = ε_nom/K z konstrukcji progów
                E = E_NOM * frk * sfr["rel"]
            if a.selftest:                                          # wyrocznia: s_k = liczba prawdziwych par sygnału w kategorii (μ_prawda = 1)
                cd, _ = cats[c][nm]
                s = np.array([((cd == k) & (D["y"] == 1)).sum() for k in range(ncat)], float)
            else:
                s = n_pair[CH[c]] * E * (fsc[c] if (fsc is not None and nm != "nomcat_mc_f1") else 1.0)
            spec[(c, nm)] = {"s": s, "S1": S1, "S2": S2, "E": E, "nmc": nmc, "sf": sfr}
            info[(c, nm)] = {"E": E.tolist(), "s": s.tolist(), "n_mc_sig": nmc.tolist(), "sf": sfr,
                             "n_data": hist_cats(D, cats[c][nm][0], ncat).sum(axis=(1, 2)).tolist()}
    rep["spec"] = {f"{CH[c]}|{nm}": v for (c, nm), v in info.items()}
    # ---- fity
    res = {}
    for nm in names:
        fits = {CH[c]: ChanFit(hist_cats(D, cats[c][nm][0], cats[c][nm][1]), spec[(c, nm)]["s"], spec[(c, nm)]["S1"], spec[(c, nm)]["S2"], *SPEC[c])
                for c in range(3)}
        r = combined_fit(fits)
        # Asimov przy μ0 = μ̂ (oczekiwane σ dla tego samego tła)
        mu0 = r["combined"]["mu"]
        fa = {}
        for c in range(3):
            f = fits[CH[c]]; st = f.best["x"].copy(); st[0] = mu0
            fa[CH[c]] = ChanFit(f.lam(st), f.s, f.S1, f.S2, *SPEC[c])
        ra = combined_fit(fa, zero=True)
        r["asimov"] = {"mu0": mu0, "mu": ra["combined"]["mu"], "sigma_68": ra["combined"]["sigma_68"], "Z": ra["combined"]["Z"]}
        # syst kształtu tła: poly3
        f3 = {CH[c]: ChanFit(fits[CH[c]].n, fits[CH[c]].s, fits[CH[c]].S1, fits[CH[c]].S2, *SPEC[c], npoly=3) for c in range(3)}
        r["poly3_mu"] = combined_fit(f3, interval=False, zero=False)["combined"]["mu"]
        # syst SF: nachylenie (górna połowa kategorii × t_top, dolna × t_bot)
        if nm not in ("nomcat_fact", "nomcat_mc", "nomcat_mc_f1"):
            ft = {}
            for c in range(3):
                sp_ = spec[(c, nm)]; sfr = sp_["sf"]; f = fits[CH[c]]
                tilt = np.where(np.arange(a.K) < a.K // 2, sfr["tilt_top"], sfr["tilt_bot"])
                ft[CH[c]] = ChanFit(f.n, f.s * tilt, f.S1, f.S2, *SPEC[c])
            r["tilt_mu"] = combined_fit(ft, interval=False, zero=False)["combined"]["mu"]
            r["sf_norm_rel_err"] = float(np.sqrt(np.mean([(spec[(c, nm)]["sf"]["err"] / spec[(c, nm)]["sf"]["rel"]) ** 2 for c in range(3)])))
        res[nm] = r
        log(f"{nm}: μ = {r['combined']['mu']:.4f} {r['combined'].get('mu_68')} Z = {r['combined'].get('Z'):.2f}; kanały",
            {ch: round(r[ch]["mu"], 3) for ch in CH}, "Asimov σ", r["asimov"]["sigma_68"])
    rep["fits"] = res
    # ---- bootstrap po zdarzeniach
    if a.nboot > 0:
        _, inv = np.unique(P["ek"], return_inverse=True)
        _G.update({"D": D, "cats": cats, "spec": spec, "schemes": names, "ek_inv": inv, "nek": int(inv.max()) + 1})
        import multiprocessing as _mp                              # fork: _G dziedziczone (py3.14 domyślnie forkserver)
        with _mp.get_context("fork").Pool(min(a.threads, a.nboot)) as pool:
            B_ = pool.map(_boot, range(a.nboot))
        bt = {nm: np.array([b[nm] for b in B_], float) for nm in names}
        rep["bootstrap"] = {"n": a.nboot, "sigma": {nm: float(np.nanstd(v, ddof=1)) for nm, v in bt.items()},
                            "delta_vs_nomcat_mc": {nm: {"mean": float(np.nanmean(bt[nm] - bt["nomcat_mc"])), "sigma": float(np.nanstd(bt[nm] - bt["nomcat_mc"], ddof=1))}
                                                   for nm in names if nm != "nomcat_mc"},
                            "corr_vs_nomcat_mc": {nm: float(np.corrcoef(bt[nm], bt["nomcat_mc"])[0, 1]) for nm in names if nm != "nomcat_mc"}}
        log("bootstrap:", rep["bootstrap"]["sigma"])
    # ---- selftest: kontrola zamknięcia
    if a.selftest:
        bad = []
        for nm in names:
            c = res[nm]["combined"]
            if c["sigma_68"] is None or abs(c["mu"] - 1.0) > 3 * c["sigma_68"]:
                bad.append(nm)
            if abs(res[nm]["asimov"]["mu"] - res[nm]["asimov"]["mu0"]) > 0.05 * max(c["sigma_68"] or 1, 1e-3) + 2e-3:
                bad.append(nm + ":asimov")
        rep["selftest_ok"] = not bad
        rep["selftest_bad"] = bad
    rep["mu_nom_charmtag"] = mu_nom
    rep["scale_now"] = a.scale_now if a.scale_now is not None else (SCALE_NOW / mu_nom if mu_nom else None)
    js = json.loads(json.dumps(rep, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    json.dump(js, open(os.path.join(a.out, "ubpnet.json"), "w"), indent=1)
    open(os.path.join(a.out, "UBPNET.md"), "w").write(to_md(js))
    log("OK", os.path.join(a.out, "UBPNET.md"))
    if a.selftest and not rep["selftest_ok"]:
        raise SystemExit(f"SELFTEST FAIL: {rep['selftest_bad']}")
    return 0


def _f(x, p=3):
    return "–" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{p}f}"


def to_md(R):
    F = R["fits"]; names = list(F)
    ref = "nomcat_mc"
    sc = R.get("scale_now") or 1.0
    s_ref = (F[ref]["combined"].get("sigma_68") or float("nan"))
    a_ref = F[ref]["asimov"]["sigma_68"] or float("nan")
    cf = R.get("config", {}); smp = cf.get("sample", "data_428770")
    md = [f"# UBPNET — μ par D–D w {smp.replace('data_', '')} z kategoriami sieci par (PAIRNET v3) — odsłonięte (zatwierdzone 18.11)", "",
          *( [f"Próbka {smp}: SF z HYBSF/HYBDS dla tej próbki; mnożniki n_pair (normalizacja z prawdy mb_nd f_ch, PROCEDURA §6): {cf.get('npair_scale') or '—'}"
              + ("; nomcat_mc_f1 = nomcat_mc bez mnożników (kontrola zgodności z nominałem charmtag)." if cf.get("npair_scale") else "."), ""]
             if smp != "data_428770" or cf.get("npair_scale") else []),
          "Ten sam model fitu 2D (sygnał MC × μ + 3 składowe tła z wolnymi normalizacjami i kształtami exp(poly2) per kategoria) i te same pary danych "
          + ("(h5 v2.0, zdarzenia próbki nominalnej v1.1) " if smp == "data_428770" else "(h5 v2.0 427929 po filtrze UB427: L_xy ≥ 0,25, okno D0, klasa MB) ") + "dla wszystkich schematów; nomcat_* = nominalne kategorie nóg (bez kategorii z tagiem D*), "
          "pnetd = sieć par z tłem z SB×SB danych (schemat główny), xgbd / base_sum = kontrola. Normalizacja n_pair z pair_inputs odsłonięcia; "
          "schemat sieci × SF_rel (pary hybrydowe: D⁺ HYBSF v3, D⁰ HYBDS v3).", ""]
    g = R.get("gate_nominal")
    if g:
        md.append(f"**Kontrola nominału (charmtag, pairs_full.npz):** " + ", ".join(f"{ch} {g[ch]['mu']:.4f} (ref {g[ch]['ref']:.4f})" for ch in CH)
                  + f"; wspólne μ_nom = {_f((g.get('combined') or {}).get('mu'), 4)}; zgodność: {'TAK' if g['ok'] else 'NIE'}. Skala do zapisu bieżącego ×{_f(sc, 3)}.")
        md.append("")
    d = R.get("data", {})
    md += [f"Dane: nogi {d.get('legs_all')} → w zdarzeniach próbki nominalnej {d.get('legs_nominal_events')} (zdarzeń z nogami {d.get('events_nominal')} z {d.get('events_all')}); "
           f"pary w zakresie mas {d.get('pairs_in_range')}; pokrycie zdarzeń nóg nominału przez nogi v2.0: {_f(d.get('nominal_leg_events_in_v20_legs'))}.", "",
           "| schemat | μ D0D0 | μ DpD0 | μ DpDp | μ wspólne [68 %] | Z | σ_68 / σ(nomcat_mc) | σ Asimov / σ Asimov(nomcat_mc) | σ bootstrap | μ × skala (zapis bieżący) |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    bs = R.get("bootstrap", {}).get("sigma", {})
    for nm in names:
        r = F[nm]; c = r["combined"]; i = c.get("mu_68") or [None, None]
        md.append(f"| {nm} | " + " | ".join(f"{_f(r[ch]['mu'], 2)} ± {_f(r[ch]['err'], 2)}" for ch in CH)
                  + f" | **{_f(c['mu'])}** [{_f(i[0])}; {_f(i[1])}] | {_f(c.get('Z'), 2)} | {_f((c.get('sigma_68') or float('nan')) / s_ref)} "
                  f"| {_f((r['asimov']['sigma_68'] or float('nan')) / a_ref)} | {_f(bs.get(nm))} | {_f(c['mu'] * sc)} |")
    md += ["", "## Różnica względem nomcat_mc (te same pary, ten sam model) i systematyki schematu", "",
           "| schemat | Δμ = μ − μ(nomcat_mc) | σ_Δ (bootstrap) | Δ/σ_Δ | korelacja (bootstrap) | syst SF norm. (δμ) | syst nachylenie SF (δμ) | syst tło poly3 (δμ) | Asimov: μ0 → μ̂ |",
           "|---|---|---|---|---|---|---|---|---|"]
    B = R.get("bootstrap", {})
    for nm in names:
        if nm == ref:
            continue
        r = F[nm]; c = r["combined"]
        dl = c["mu"] - F[ref]["combined"]["mu"]
        sd = (B.get("delta_vs_nomcat_mc", {}).get(nm) or {}).get("sigma")
        md.append(f"| {nm} | {_f(dl)} | {_f(sd)} | {_f(dl / sd if sd else None, 2)} | {_f((B.get('corr_vs_nomcat_mc') or {}).get(nm), 2)} "
                  f"| {_f(c['mu'] * r['sf_norm_rel_err'] if 'sf_norm_rel_err' in r else None)} | {_f(r['tilt_mu'] - c['mu'] if 'tilt_mu' in r else None)} "
                  f"| {_f(r['poly3_mu'] - c['mu'])} | {_f(r['asimov']['mu0'])} → {_f(r['asimov']['mu'])} |")
    md += ["", "## SF_rel schematów (SF schematu / SF nomcat; para = iloczyn nóg)", "", "| kanał | schemat | SF_rel | ± | nachylenie góra / dół | ρ(schemat, nomcat) MC |", "|---|---|---|---|---|---|"]
    for k, v in R["spec"].items():
        ch, nm = k.split("|")
        s = v.get("sf")
        if s:
            md.append(f"| {ch} | {nm} | {_f(s['rel'])} | {_f(s['err'])} | {_f(s['tilt_top'])} / {_f(s['tilt_bot'])} | {_f(s.get('rho_mc'), 2)} |"
                      + ("" if s.get("available") else " (brak SF — 1,0)"))
    T = R.get("sb_transfer")
    if T:
        md += ["", "## Transfer tła SBTRAIN (pary SB×SB: ta próbka vs referencja 428770, ten sam model sieci)", "",
               "Tło w ficie jest wolne per kategoria — różnica rozkładów nie obciąża μ, mierzy tylko, czy sieć rozdziela tło tak jak na próbce treningowej.", "",
               "| kanał | schemat | n SB×SB (ref) | χ²/ndf | w schemacie: próbka / ref | górna połowa kategorii: próbka / ref |", "|---|---|---|---|---|---|"]
        for ch, v in T.items():
            for nm, x in v.items():
                md.append(f"| {ch} | {nm} | {x['n']} ({x['n_ref']}) | {x['chi2']:.1f}/{x['ndf']} | {_f(x['in_scheme'])} / {_f(x['in_scheme_ref'])} | {_f(x['top_half'])} / {_f(x['top_half_ref'])} |")
    md += ["", "## Dobroć fitów (dewiancja / ndf, kanały)", ""]
    for nm in names:
        md.append(f"- {nm}: " + ", ".join(f"{ch} {F[nm][ch]['dev']:.0f}/{F[nm][ch]['ndf']}" for ch in CH))
    md += ["", "Uwagi: σ z profilu bez korekty korelacji zdarzeń (bootstrap Poissona po zdarzeniach ją zawiera); normalizacja i tło D* tag poza modelem "
           "(nomcat_* bez kategorii z tagiem D* → nie jest identyczny z nominałem); progi sieci z MC, więc E_k = 0,65²/K × SF_rel."]
    if R.get("config", {}).get("selftest"):
        md += ["", f"SELFTEST: {'OK' if R.get('selftest_ok') else 'FAIL ' + str(R.get('selftest_bad'))} (wyrocznia: s_k = prawdziwe pary, μ_prawda = 1)"]
    return "\n".join(md) + "\n"


if __name__ == "__main__":
    sys.exit(main())
