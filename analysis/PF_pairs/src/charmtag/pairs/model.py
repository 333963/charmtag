"""Model fitu 2D par (port M7 v3 + plan B; liczby 1:1 z unblinding_freeze_v20260927_01/code/m7_fit.py i kcat.py).

Wejścia (`Inputs`): szablony sygnału nóg (MC: poprawne + swap dla D0), histogramy fałszywych nóg (sidebandy
danych), wydajności kategorii Y, normalizacja (oczekiwane liczby prawdziwych D i par), zliczenia par
sideband×sideband. Masy w MeV, 25 binów na nogę (D0 1750–2000, D⁺ 1700–2200).

Składowe w każdej parze kategorii (i, j):
  SS = μ · E_SS(i,j)                     szablon S_a(i) ⊗ S_b(j)
  SB = κ_(a,i) · d_ij · E_SB(i,j)        S_a ⊗ F_b   (d_ij więzione N(1, DEV))
  BS = κ_(b,j) · d'_ij · E_BS(i,j)       F_a ⊗ S_b
  BB = n_ij (wolne, start z GLM sideband×sideband)
F = exp(poly2) (opcjonalnie poly3) na gatunek×kategorię, parametry wolne.
κ: 'planB' (na kategorię prawdziwej nogi), 'channel' (jedno κ na kanał i kolejność — kontrola M7),
'partner' (κ_(noga, partner luźny?) z więzem ln κ_ciasny/κ_luźny ~ N(0, σ) — T12).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

import numpy as np
from iminuit import Minuit

from . import CHANNELS

M0 = {10: 1864.84, 4: 1869.66}
SBW = {10: 75.0, 4: 70.0}
LOOSE = {10: {2, 4}, 4: {2}}          # T3/T12: luźne kategorie partnera


# ----------------------------------------------------------------------------- wejścia
@dataclass
class Inputs:
    edges: dict                         # {dt: array krawędzi [MeV]}
    ncat: dict                          # {dt: liczba kategorii}
    sig: dict                           # {(dt,c): szablon sygnału (suma 1)}
    sig_comp: dict                      # {(dt,c): (poprawne, swap)} — znormalizowane do sumy sig (opcjonalnie)
    Y: dict                             # {(dt,c): wydajność kategorii × SF}
    dataleg: dict                       # {(dt,c): histogram nóg z sidebandów danych}
    nev_legs: float                     # zdarzenia, z których pochodzą dataleg (gęstość fałszywych nóg)
    ndata: float                        # zdarzenia analizowanej próbki
    n_real: dict                        # {dt: oczekiwana liczba prawdziwych D (przed kategorią) w próbce, scen. R}
    n_pair: dict                        # {ch: oczekiwana liczba par prawdziwych D (przed kategoriami), scen. R}
    sbsb: dict                          # {(ch,i,j): surowe zliczenia par sideband×sideband}
    meta: dict = field(default_factory=dict)

    # -- pochodne
    def __post_init__(self):
        self.edges = {int(k): np.asarray(v, float) for k, v in self.edges.items()}
        self.C = {k: 0.5 * (v[1:] + v[:-1]) for k, v in self.edges.items()}
        self.U = {k: 2 * (c - c[0]) / (c[-1] - c[0]) - 1 for k, c in self.C.items()}
        self.SBM = {k: np.abs(self.C[k] - M0[k]) > SBW[k] for k in self.C}
        self.nbins = {k: len(v) - 1 for k, v in self.edges.items()}

    @property
    def scale(self):
        return self.ndata / self.nev_legs

    def to_files(self, prefix):
        js = {"edges": {str(k): v.tolist() for k, v in self.edges.items()}, "ncat": {str(k): v for k, v in self.ncat.items()},
              "Y": {f"{dt}_{c}": v for (dt, c), v in self.Y.items()}, "nev_legs": self.nev_legs, "ndata": self.ndata,
              "n_real": {str(k): v for k, v in self.n_real.items()}, "n_pair": self.n_pair,
              "sbsb": {f"{ch}_{i}_{j}": int(v) for (ch, i, j), v in self.sbsb.items()}, "meta": self.meta}
        with open(prefix + ".json", "w") as f:
            json.dump(js, f, indent=1, default=float)
        arr = {f"sig_{dt}_{c}": v for (dt, c), v in self.sig.items()}
        arr.update({f"dataleg_{dt}_{c}": v for (dt, c), v in self.dataleg.items()})
        for (dt, c), (a, b) in self.sig_comp.items():
            arr[f"sigc_{dt}_{c}"] = a
            arr[f"sigs_{dt}_{c}"] = b
        np.savez(prefix + ".npz", **arr)

    @classmethod
    def from_files(cls, prefix):
        js = json.load(open(prefix + ".json"))
        A = dict(np.load(prefix + ".npz"))
        ncat = {int(k): int(v) for k, v in js["ncat"].items()}
        sig, comp, dl = {}, {}, {}
        for dt, n in ncat.items():
            for c in range(n):
                sig[(dt, c)] = A[f"sig_{dt}_{c}"].astype(float)
                dl[(dt, c)] = A[f"dataleg_{dt}_{c}"].astype(float)
                if f"sigc_{dt}_{c}" in A:
                    comp[(dt, c)] = (A[f"sigc_{dt}_{c}"].astype(float), A[f"sigs_{dt}_{c}"].astype(float))
        Y = {(int(k.split("_")[0]), int(k.split("_")[1])): float(v) for k, v in js["Y"].items()}
        sbsb = {}
        for k, v in js["sbsb"].items():
            ch, i, j = k.split("_")
            sbsb[(ch, int(i), int(j))] = int(v)
        return cls(edges=js["edges"], ncat=ncat, sig=sig, sig_comp=comp, Y=Y, dataleg=dl, nev_legs=float(js["nev_legs"]),
                   ndata=float(js["ndata"]), n_real={int(k): float(v) for k, v in js["n_real"].items()},
                   n_pair={k: float(v) for k, v in js["n_pair"].items()}, sbsb=sbsb, meta=js.get("meta", {}))


def cat_pairs(ncat, a, b):
    return [(i, j) for i in range(ncat[a]) for j in range(ncat[b]) if not (a == b and j < i)]


# ----------------------------------------------------------------------------- kształty
class Shapes:
    def __init__(self, I: Inputs, bkg_npar=2):
        self.I = I
        self.TH, self.NF, self.TH3 = {}, {}, {}
        for (dt, c) in I.dataleg:
            self.TH[(dt, c)], nf = self.fit_leg_shape(dt, c, 2)
            self.NF[(dt, c)] = nf / I.nev_legs
            self.TH3[(dt, c)], _ = self.fit_leg_shape(dt, c, 3)

    def fshape(self, dt, th):
        u = self.I.U[dt]
        z = th[0] * u + th[1] * u ** 2 + (th[2] * u ** 3 if len(th) > 2 else 0.0)
        z = np.exp(z - z.max())
        return z / z.sum()

    def fit_leg_shape(self, dt, c, npar=2):
        h = self.I.dataleg[(dt, c)].astype(float)
        m = self.I.SBM[dt]

        def nll(th):
            f = self.fshape(dt, th)
            lam = f[m] / f[m].sum() * h[m].sum() + 1e-12
            return float(np.sum(lam - h[m] * np.log(lam)))
        mi = Minuit(nll, np.zeros(npar))
        mi.errordef = 0.5
        mi.migrad()
        th = np.array(mi.values)
        return th, h[m].sum() / self.fshape(dt, th)[m].sum()

    def stretch(self, t, dt, k=1.0, shift=0.0):
        cdf = np.r_[0, np.cumsum(t)]
        e = self.I.edges[dt]
        s = np.diff(np.interp(M0[dt] + (e - M0[dt] - shift) / k, e, cdf, left=0, right=1))
        return s / s.sum()


# ----------------------------------------------------------------------------- oczekiwania
class Expectation:
    """expected(ch, kk): oczekiwane SS/SB/BS/BB w każdej parze kategorii na pełnej płaszczyźnie (jak M7 v3)."""

    def __init__(self, I: Inputs, S: Shapes, bbmode="glm"):
        self.I, self.S, self.bbmode = I, S, bbmode
        self.FRAC = {}
        self.GLM = self._bb_glm()

    def _sbsb_frac(self, ch, a, i, b, j):
        S, I = self.S, self.I
        BB2 = np.outer(S.fshape(a, S.TH[(a, i)]), S.fshape(b, S.TH[(b, j)]))
        if a == b and i == j:
            BB2 = 0.5 * (BB2 + BB2.T)
        return BB2[np.ix_(I.SBM[a], I.SBM[b])].sum()

    def _unit(self, a, i, b, j):
        same_diag = a == b and i == j
        return self.I.ndata * self.S.NF[(a, i)] * self.S.NF[(b, j)] * (0.5 if same_diag else 1.0)

    def _bb_glm(self):
        I = self.I
        rows = []
        for ch, (a, b) in CHANNELS.items():
            for i, j in cat_pairs(I.ncat, a, b):
                fr = self._sbsb_frac(ch, a, i, b, j)
                self.FRAC[(ch, i, j)] = fr
                rows.append((ch, a, i, b, j, I.sbsb.get((ch, i, j), 0), self._unit(a, i, b, j) * fr / I.scale))
        legs_ = sorted({(r[1], r[2]) for r in rows} | {(r[3], r[4]) for r in rows})
        ref = {dt: max(c for d, c in legs_ if d == dt) for dt in {lg[0] for lg in legs_}}
        ref = {dt: (2 if (dt, 2) in legs_ else ref[dt]) for dt in ref}     # M7: kategoria 2 = referencja
        free = [lg for lg in legs_ if lg[1] != ref[lg[0]]]
        chs = list(CHANNELS)

        def nll(x):
            c = dict(zip(chs, x[:3]))
            bb = {lg: 0.0 for lg in legs_}
            bb.update(dict(zip(free, x[3:])))
            v = 0.0
            for ch, a, i, b, j, n, base in rows:
                mu = base * np.exp(c[ch] + bb[(a, i)] + bb[(b, j)])
                v += mu - n * np.log(mu + 1e-300)
            return v
        m = Minuit(nll, np.zeros(3 + len(free)))
        m.errordef = 0.5
        m.migrad()
        x = np.array(m.values)
        bb = {lg: 0.0 for lg in legs_}
        bb.update(dict(zip(free, x[3:])))
        return {"c": dict(zip(chs, x[:3])), "b": bb, "rows": rows}

    def expected(self, ch, kk=1.0):
        I, S = self.I, self.S
        a, b = CHANNELS[ch]
        out = {}
        for i, j in cat_pairs(I.ncat, a, b):
            same_diag = a == b and i == j
            mult = 1.0 if (a != b or i == j) else 2.0
            SS = I.n_pair[ch] * kk * I.Y[(a, i)] * I.Y[(b, j)] * mult
            if same_diag:
                rf = I.n_real[a] * I.Y[(a, i)] * S.NF[(b, j)]
                sb, bs = 0.5 * rf, 0.5 * rf
            else:
                sb = I.n_real[a] * I.Y[(a, i)] * S.NF[(b, j)]
                bs = I.n_real[b] * I.Y[(b, j)] * S.NF[(a, i)]
            fr = self.FRAC[(ch, i, j)]
            raw = I.sbsb.get((ch, i, j), 0)
            unit = self._unit(a, i, b, j)
            out[(i, j)] = {"SS": SS, "SB": sb, "BS": bs, "BB_count": raw * I.scale / fr, "BB_raw": raw, "BB_unit": unit}
        rho = sum(v["BB_count"] for v in out.values()) / sum(v["BB_unit"] for v in out.values())
        for (i, j), v in out.items():
            v["rho"] = rho
            if self.bbmode == "count":
                v["BB"] = v["BB_count"]
            elif self.bbmode == "fact":
                v["BB"] = rho * v["BB_unit"]
            else:
                v["BB"] = v["BB_unit"] * float(np.exp(self.GLM["c"][ch] + self.GLM["b"][(a, i)] + self.GLM["b"][(b, j)]))
        return out


# ----------------------------------------------------------------------------- model kanału
class PairFit:
    """Model jednego kanału. kappa: planB | channel | partner. syst: profilowane systematyki kształtu."""
    WID = {10: 0.06, 4: 0.03}
    SHIFT = 1.0
    TAGREL = 0.10 / 0.72

    def __init__(self, S: Shapes, ch, exp, kappa="planB", syst=False, partner_sigma=None, dev=0.10, tag_cats=(3, 4)):
        self.S, self.I, self.ch, self.exp = S, S.I, ch, exp
        self.a, self.b = CHANNELS[ch]
        self.cps = cat_pairs(self.I.ncat, self.a, self.b)
        self.kappa, self.syst, self.DEV, self.psig = kappa, syst, dev, partner_sigma
        self.tag_cats = set(tag_cats)
        self.legs = sorted({(self.a, i) for i, _ in self.cps} | {(self.b, j) for _, j in self.cps})
        names = ["mu"] + [f"{k}_{i}{j}" for i, j in self.cps for k in ("sb", "bs", "bb")] + \
                [f"th{p}_{dt}_{c}" for dt, c in self.legs for p in (0, 1)]
        self.nb = len(names)
        if kappa == "channel":
            self.kkeys = ["kSB", "kBS"]
            names += ["kSB", "kBS"]
        elif kappa == "partner":
            keys = set()
            for i, j in self.cps:
                keys.add((self.a, i, self._pk(self.b, j)))
                keys.add((self.b, j, self._pk(self.a, i)))
            self.kkeys = sorted(keys)
            names += [f"kap_{dt}_{c}_p{p}" for dt, c, p in self.kkeys]
        else:
            self.kkeys = self.legs
            names += [f"kap_{dt}_{c}" for dt, c in self.kkeys]
        self.n0 = len(names)
        self.th2, self.alphas = [], []
        if syst:
            self.th2 = [f"th2_{dt}_{c}" for dt, c in self.legs]
            sp = sorted({self.a, self.b})
            self.alphas = [f"aw_{s}" for s in sp] + [f"as_{s}" for s in sp] + (["atag"] if 10 in sp else [])
            names += self.th2 + self.alphas
            self.var = {}
            for lg in self.legs:
                t0, dt = self.I.sig[lg], lg[0]
                self.var[lg] = {"w": (S.stretch(t0, dt, k=1 - self.WID[dt]), t0, S.stretch(t0, dt, k=1 + self.WID[dt])),
                                "s": (S.stretch(t0, dt, shift=-self.SHIFT), t0, S.stretch(t0, dt, shift=self.SHIFT))}
        self.names = names
        self.diag = [k for k, nm in enumerate(names) if self.a == self.b and nm[:3] == "bs_" and nm[3] == nm[4]]
        self.dev = [k for k, nm in enumerate(names) if nm[:3] in ("sb_", "bs_") and k not in self.diag]
        self.ia = np.array([names.index(x) for x in self.alphas], int)
        if kappa == "partner":
            idx = {k: self.nb + n for n, k in enumerate(self.kkeys)}
            self.ppairs = [(idx[(dt, c, 0)], idx[(dt, c, 1)]) for (dt, c, p) in self.kkeys if p == 0 and (dt, c, 1) in idx]

    def _pk(self, dt, c):
        return int(c in LOOSE.get(dt, set()))

    # start
    def x0(self):
        x = [1.0]
        for cp in self.cps:
            e = self.exp[cp]
            x += [1.0, 1.0, e["BB"]]
        for lg in self.legs:
            x += list(self.S.TH[lg])
        x += [1.0] * (self.n0 - self.nb)
        x += [0.0] * (len(self.th2) + len(self.alphas))
        return np.array(x, float)

    @staticmethod
    def _interp(lo, t0, hi, a):
        return t0 + 0.5 * a * (hi - lo) + 0.5 * a * a * (hi + lo - 2 * t0)

    def lam(self, x, sig=None, fshapes=None, ss_scale=None):
        S, I = self.S, self.I
        x = np.asarray(x, float)
        mu = x[0]
        k = 1
        rf = {}
        for cp in self.cps:
            rf[cp] = x[k:k + 3]
            k += 3
        al, th2 = {}, {}
        if self.syst:
            th2 = dict(zip(self.legs, x[self.n0:self.n0 + len(self.th2)]))
            al = dict(zip(self.alphas, x[self.n0 + len(self.th2):]))
        f = {}
        for lg in self.legs:
            th = x[k:k + 2] if not self.syst else np.r_[x[k:k + 2], th2[lg]]
            f[lg] = S.fshape(lg[0], th) if fshapes is None else fshapes[lg]
            k += 2
        if sig is None:
            if self.syst:
                sig = {}
                for lg in self.legs:
                    dt, v = lg[0], self.var[lg]
                    t0 = v["w"][1]
                    t = t0 + (self._interp(*v["w"], al[f"aw_{dt}"]) - t0) + (self._interp(*v["s"], al[f"as_{dt}"]) - t0)
                    t = np.clip(t, 0, None)
                    sig[lg] = t / t.sum()
            else:
                sig = I.sig
        kv = dict(zip(self.kkeys, x[self.nb:self.n0]))
        tagf = 1.0 + self.TAGREL * al.get("atag", 0.0)
        out = []
        for i, j in self.cps:
            e = self.exp[(i, j)]
            sa, sb_ = sig[(self.a, i)], sig[(self.b, j)]
            fa, fb = f[(self.a, i)], f[(self.b, j)]
            dsb, dbs, nbb = rf[(i, j)]
            if self.kappa == "channel":
                ka, kb = x[self.nb], (x[self.nb + 1] if self.a != self.b else x[self.nb])
            elif self.kappa == "partner":
                ka, kb = kv[(self.a, i, self._pk(self.b, j))], kv[(self.b, j, self._pk(self.a, i))]
            else:
                ka, kb = kv[(self.a, i)], kv[(self.b, j)]
            ss = e["SS"]
            if self.syst:
                ss *= (tagf if (self.a == 10 and i in self.tag_cats) else 1.0) * (tagf if (self.b == 10 and j in self.tag_cats) else 1.0)
            if ss_scale is not None:
                ss *= ss_scale[(i, j)]
            nsb, nbs = ka * dsb * e["SB"], kb * dbs * e["BS"]
            SS = np.outer(sa, sb_)
            SB = np.outer(sa, fb)
            BS = np.outer(fa, sb_)
            BB = np.outer(fa, fb)
            if self.a == self.b and i == j:
                if self.kappa == "channel":
                    nsb = ka * dsb * e["SB"] + ka * dbs * e["BS"]
                else:
                    nsb = ka * dsb * (e["SB"] + e["BS"])
                nbs = 0.0
                SB = 0.5 * (SB + SB.T)
                BB = 0.5 * (BB + BB.T)
            out.append(mu * ss * SS + nsb * SB + nbs * BS + nbb * BB)
        return np.array(out)

    def nll_fn(self, n, mask):
        Mk = np.broadcast_to(mask, n.shape)
        d = np.array(self.dev, int)
        nn = n[Mk]

        def nll(x):
            lam = np.clip(self.lam(x), 1e-9, None)[Mk]
            v = float(np.sum(lam - nn * np.log(lam))) + 0.5 * float(np.sum(((x[d] - 1.0) / self.DEV) ** 2))
            if len(self.ia):
                v += 0.5 * float(np.sum(x[self.ia] ** 2))
            if self.kappa == "partner" and self.psig:
                for i, j in self.ppairs:
                    v += 0.5 * (np.log(max(x[i], 1e-6) / max(x[j], 1e-6)) / self.psig) ** 2
            return v
        return nll

    def fit(self, n, mask=None, fix_mu=None, start=None, fix_names=(), ncall=60000):
        mask = np.ones(n.shape[1:], bool) if mask is None else mask
        m = Minuit(self.nll_fn(n, mask), self.x0() if start is None else start, name=self.names)
        m.errordef = 0.5
        for kidx, nm in enumerate(self.names[1:], start=1):
            if kidx in self.diag:
                m.fixed[nm] = True
            elif nm[:2] in ("sb", "bs", "bb") or nm.startswith("kap") or nm in ("kSB", "kBS"):
                m.limits[nm] = ((1e-4 if self.kappa == "partner" and nm.startswith("kap") else 0), None)
            elif nm.startswith(("aw_", "as_", "atag")):
                m.limits[nm] = (-4, 4)
            else:
                m.limits[nm] = (-8, 8)
        if self.kappa == "channel" and self.a == self.b:
            m.fixed["kBS"] = True
        for nm in fix_names:
            m.fixed[nm] = True
        m.limits["mu"] = (-20, 40)
        if fix_mu is not None:
            m.values["mu"] = fix_mu
            m.fixed["mu"] = True
        m.strategy = 1
        m.migrad(ncall=ncall)
        if not m.valid and self.kappa == "partner":
            m.strategy = 2
            m.migrad(ncall=2 * ncall)
        m.hesse()
        return m

    def deviance(self, m, n, mask=None):
        mask = np.ones(n.shape[1:], bool) if mask is None else mask
        Mk = np.broadcast_to(mask, n.shape)
        lam = np.clip(self.lam(np.array(m.values)), 1e-9, None)[Mk]
        nn = n[Mk]
        dev = 2 * np.sum(lam - nn + np.where(nn > 0, nn * np.log(np.where(nn > 0, nn, 1) / lam), 0.0))
        nfree = sum(1 for nm in self.names if not m.fixed[nm])
        return float(dev), int(Mk.sum() - nfree)


def region_mask(I: Inputs, ch, region):
    """full | blind (bez centrum: obie nogi w oknie) | sbsb (obie nogi w sidebandach) | arms."""
    a, b = CHANNELS[ch]
    wa, wb = ~I.SBM[a], ~I.SBM[b]
    centre = np.outer(wa, wb)
    arms = np.outer(wa, ~wb) | np.outer(~wa, wb)
    return {"full": np.ones_like(centre), "blind": ~centre, "sbsb": ~(centre | arms), "arms": arms, "centre": centre}[region]


class Model:
    """Wygodny pojemnik: wejścia + kształty + oczekiwania."""

    def __init__(self, I: Inputs, bbmode="glm"):
        self.I = I
        self.S = Shapes(I)
        self.E = Expectation(I, self.S, bbmode)

    def channel(self, ch, kk=1.0, **kw):
        return PairFit(self.S, ch, self.E.expected(ch, kk), **kw)


def load_inputs(path_prefix):
    if os.path.isdir(path_prefix):
        path_prefix = os.path.join(path_prefix, "pair_inputs")
    return Inputs.from_files(path_prefix)
