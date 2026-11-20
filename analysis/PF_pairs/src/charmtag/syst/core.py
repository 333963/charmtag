"""Systematyki ślepe na Asimovie: kształty, SF, tag D*, κ, wydajności kategorii (etap U)."""
from __future__ import annotations

import copy
import logging

import numpy as np

from ..pairs import CHANNELS
from ..pairs.model import Model, PairFit, region_mask

log = logging.getLogger("charmtag")

DEFAULTS = {
    "width": {10: 0.06, 4: 0.03},          # względna szerokość piku (jak kcat_syst 09-28)
    "mshift": 1.0,                          # MeV
    "swap": 0.20,                           # względna zmiana udziału swapów D0
    "refl": {"refl_Ds": 0.05, "refl_Lc": 0.01},   # ułamek sygnału D⁺ w oknie, kształt płaski (pary_pt_refleksy 09-28)
    "sf_rel": {10: 0.073, 4: 0.097},        # 09-26: D0 0,88±0,05±0,04; D⁺ 0,95±0,07±0,06 (nadpisać wynikiem etapu N)
    "tag_rel": 0.14,                        # SF tagu D* ±14% (PairFit.TAGREL)
    "cat_rel": 0.05,                        # niezależna niepewność Y na kategorię (kształt SF vs WP)
    "tag_cats": (3, 4),
}
SP = {10: "D0", 4: "Dp"}


def settings(cfg):
    s = copy.deepcopy(DEFAULTS)
    u = (cfg or {}).get("syst", {}) or {}
    for k, v in u.items():
        if k in ("width", "sf_rel") and isinstance(v, dict):
            s[k].update({(10 if str(kk) in ("10", "D0") else 4): float(x) for kk, x in v.items()})
        else:
            s[k] = v
    return s


# ----------------------------------------------------------------------------- fit z globalnymi obserwablami
class AsimovFit(PairFit):
    """PairFit, w którym więzy d_ij ~ N(c_ij, DEV) mają środki c (globalne obserwable Asimov)."""
    cdev = None

    def nll_fn(self, n, mask):
        Mk = np.broadcast_to(mask, n.shape)
        d = np.array(self.dev, int)
        c = np.ones(len(d)) if self.cdev is None else np.asarray(self.cdev, float)
        nn = n[Mk]

        def nll(x):
            lam = np.clip(self.lam(x), 1e-9, None)[Mk]
            v = float(np.sum(lam - nn * np.log(lam))) + 0.5 * float(np.sum(((x[d] - c) / self.DEV) ** 2))
            if len(self.ia):
                v += 0.5 * float(np.sum(x[self.ia] ** 2))
            if self.kappa == "partner" and self.psig:
                for i, j in self.ppairs:
                    v += 0.5 * (np.log(max(x[i], 1e-6) / max(x[j], 1e-6)) / self.psig) ** 2
            return v
        return nll


def channel(M: Model, ch, **kw):
    cf = M.channel(ch, **kw)
    cf.__class__ = AsimovFit
    return cf


def transfer(x_src, names_src, cf_dst):
    """Start/środki w modelu docelowym: wspólne nazwy z x_src, reszta z x0."""
    x = cf_dst.x0()
    idx = {n: k for k, n in enumerate(names_src)}
    for k, n in enumerate(cf_dst.names):
        if n in idx:
            x[k] = x_src[idx[n]]
    return x


def set_centres(cf, x):
    cf.cdev = np.asarray(x, float)[np.array(cf.dev, int)] if len(cf.dev) else np.zeros(0)
    return cf


def asimov_fit(cf, lam, start):
    full = np.ones(lam.shape[1:], bool)
    m1 = cf.fit(lam, full, start=start)
    m0 = cf.fit(lam, full, fix_mu=0.0, start=np.array(m1.values))
    q0 = 2 * (m0.fval - m1.fval) if m1.values["mu"] > 0 else 0.0
    return {"mu": float(m1.values["mu"]), "err": float(m1.errors["mu"]), "Z": float(np.sqrt(max(q0, 0.0))),
            "valid": bool(m1.valid), "x": np.array(m1.values)}


# ----------------------------------------------------------------------------- punkt bazowy
def base_point(M: Model, H: dict, ch, kappa="planB", psig=None):
    """Nuisance'y z obszaru ślepego przy μ=1 (jak fit2d.run_blind). Gdy brak danych → x0 (nominalne oczekiwanie)."""
    cf = channel(M, ch, kappa=kappa, partner_sigma=psig)
    n = None if H is None else H.get(ch)
    if n is None or n.sum() < 20:
        x = cf.x0()
        src = "x0"
    else:
        m = cf.fit(n, region_mask(M.I, ch, "blind"), fix_mu=1.0)
        x = np.array(m.values)
        src = "blind_fit" if m.valid else "blind_fit_invalid"
        cf._blind_n = n
    set_centres(cf, x)
    return cf, x, src


# ----------------------------------------------------------------------------- warianty prawdy
def _flat(I, dt):
    return np.ones(I.nbins[dt])


def truth_variant(M: Model, cf, x, name, sign, S=None):
    """λ prawdy dla wariantu `name` (kierunek sign=±1) w parametryzacji cf (ten sam wektor x)."""
    S = S or DEFAULTS
    I, Sh = M.I, M.S
    if name.startswith("width_"):
        dt = 10 if name.endswith("D0") else 4
        k = 1.0 + sign * S["width"][dt]
        return cf.lam(x, sig={lg: (Sh.stretch(I.sig[lg], dt, k=k) if lg[0] == dt else I.sig[lg]) for lg in I.sig})
    if name.startswith("mscale_"):
        dt = 10 if name.endswith("D0") else 4
        return cf.lam(x, sig={lg: (Sh.stretch(I.sig[lg], dt, shift=sign * S["mshift"]) if lg[0] == dt else I.sig[lg]) for lg in I.sig})
    if name == "bkg_poly3":
        # prawda = model z exp(poly3) dopasowany na obszarze ślepym (μ=1; α kształtu sygnału = 0) → fit modelem poly2
        if sign < 0:
            return None
        cfs = channel(M, cf.ch, syst=True)
        xs = transfer(x, cf.names, cfs)
        fix = [nm for nm in cfs.names if nm.startswith(("aw_", "as_", "atag"))]
        if getattr(cf, "_blind_n", None) is not None:
            m = cfs.fit(cf._blind_n, region_mask(I, cf.ch, "blind"), fix_mu=1.0, start=xs, fix_names=fix)
            xs = np.array(m.values)
        else:
            k = 1 + 3 * len(cf.cps)
            for lg in cf.legs:
                xs[cfs.names.index(f"th2_{lg[0]}_{lg[1]}")] = Sh.TH3[lg][2]
        lam3 = cfs.lam(xs)
        cf._poly3 = (cfs, xs, lam3)
        return lam3
    if name == "swap_D0":
        if not I.sig_comp:
            raise ValueError("brak składowych poprawne/swap (sig_comp)")
        fsw = 1.0 + sign * S["swap"]
        sig, ratio = dict(I.sig), {}
        for c in range(I.ncat[10]):
            cc, ss = I.sig_comp[(10, c)]
            t = cc + fsw * ss
            sig[(10, c)] = t / t.sum()
            ratio[c] = t.sum() / (cc + ss).sum()
        scale = {(i, j): (ratio[i] if cf.a == 10 else 1.0) * (ratio[j] if cf.b == 10 else 1.0) for i, j in cf.cps}
        return cf.lam(x, sig=sig, ss_scale=scale)
    if name.startswith("refl_"):
        r = S["refl"][name] * (1 if sign > 0 else 0)          # jednostronne (refleksja może tylko dodać)
        sig = {}
        for lg, t in I.sig.items():
            if lg[0] == 4:
                win = ~I.SBM[4]
                u = _flat(I, 4)
                u = u * t[win].sum() / u[win].sum()
                sig[lg] = t + r * u                            # nienormowane: dodaje pary z refleksją przy prawdziwym partnerze
            else:
                sig[lg] = t
        return cf.lam(x, sig=sig)
    raise ValueError(name)


def scaled_model(M: Model, yscale: dict, n_pair_scale=None):
    """Model ze zmienionym Y (SF) — te same kształty tła; Expectation przeliczone."""
    I2 = copy.copy(M.I)
    I2.Y = {k: v * yscale.get(k, 1.0) for k, v in M.I.Y.items()}
    if n_pair_scale:
        I2.n_pair = {k: v * n_pair_scale.get(k, 1.0) for k, v in M.I.n_pair.items()}
    M2 = copy.copy(M)
    M2.I = I2
    from ..pairs.model import Expectation
    M2.E = copy.copy(M.E)
    M2.E.I = I2
    return M2


def truth_yscale(M: Model, ch, x, yscale, tagf=1.0):
    M2 = scaled_model(M, yscale)
    cf2 = M2.channel(ch)
    return cf2.lam(x)


# ----------------------------------------------------------------------------- przebieg
SHAPE_SRC = ["width_D0", "width_Dp", "mscale_D0", "mscale_Dp", "bkg_poly3", "swap_D0", "refl_Ds", "refl_Lc"]


def applies(name, ch):
    a, b = CHANNELS[ch]
    if name.endswith("D0") or name in ("sf_tag",):
        return 10 in (a, b)
    if name.endswith("Dp") or name.startswith("refl_"):
        return 4 in (a, b)
    return True


def run_channel(M: Model, H, ch, S=None, sources=None, kappa_alts=True):
    """Wszystkie warianty Asimov dla kanału ch. Zwraca dict z nominałem i Δμ/ΔZ per źródło (up/down)."""
    S = S or DEFAULTS
    cf, x, src = base_point(M, H, ch)
    lam0 = cf.lam(x)
    nom = asimov_fit(cf, lam0, x)
    out = {"base": src, "nominal": {k: v for k, v in nom.items() if k != "x"}, "sources": {}}
    log.info("syst %s: baza=%s Asimov μ=%.4f σ=%.3f Z=%.2f", ch, src, nom["mu"], nom["err"], nom["Z"])

    def rec(name, lam_up, lam_dn=None, fit_cf=None, start=None):
        fc = fit_cf or cf
        st = x if start is None else start
        r = {}
        for tag, lam in (("up", lam_up), ("down", lam_dn)):
            if lam is None:
                continue
            a = asimov_fit(fc, lam, st)
            r[tag] = {"dmu": a["mu"] - nom["mu"], "dZ": a["Z"] - nom["Z"], "Z": a["Z"], "err": a["err"], "valid": a["valid"]}
        r["dmu_sym"] = float(max(abs(v["dmu"]) for k, v in r.items() if isinstance(v, dict)))
        r["dZ_min"] = float(min(v["dZ"] for k, v in r.items() if isinstance(v, dict)))
        out["sources"][name] = r
        log.info("  %-16s Δμ %s  ΔZ %s", name, " ".join(f"{k}:{v['dmu']:+.3f}" for k, v in r.items() if isinstance(v, dict)),
                 " ".join(f"{k}:{v['dZ']:+.2f}" for k, v in r.items() if isinstance(v, dict)))

    srcs = sources or (SHAPE_SRC + ["sf_D0", "sf_Dp", "sf_tag", "cat_Y", "kappa"])
    for name in srcs:
        if name in SHAPE_SRC:
            if not applies(name, ch):
                continue
            try:
                up = truth_variant(M, cf, x, name, +1, S)
                dn = None if name.startswith("refl_") or name == "bkg_poly3" else truth_variant(M, cf, x, name, -1, S)
                rec(name, up, dn)
            except ValueError as e:
                out["sources"][name] = {"error": str(e)}
        elif name in ("sf_D0", "sf_Dp"):
            dt = 10 if name == "sf_D0" else 4
            if not applies(name, ch):
                continue
            ys = lambda s: {k: (1 + s * S["sf_rel"][dt]) for k in M.I.Y if k[0] == dt}
            rec(name, truth_yscale(M, ch, x, ys(+1)), truth_yscale(M, ch, x, ys(-1)))
        elif name == "sf_tag":
            if not applies(name, ch):
                continue
            ys = lambda s: {(10, c): (1 + s * S["tag_rel"]) for c in S["tag_cats"]}
            rec(name, truth_yscale(M, ch, x, ys(+1)), truth_yscale(M, ch, x, ys(-1)))
        elif name == "cat_Y":
            a, b = CHANNELS[ch]
            q = []
            for lg in sorted({lg for lg in M.I.Y if lg[0] in (a, b)}):
                lu = truth_yscale(M, ch, x, {lg: 1 + S["cat_rel"]})
                r = asimov_fit(cf, lu, x)
                q.append((lg, r["mu"] - nom["mu"], r["Z"] - nom["Z"]))
            out["sources"]["cat_Y"] = {"per_cat": {f"{dt}_{c}": [dm, dz] for (dt, c), dm, dz in q},
                                       "dmu_sym": float(np.sqrt(sum(dm ** 2 for _, dm, _ in q))),
                                       "dZ_min": -float(np.sqrt(sum(min(dz, 0) ** 2 for *_, dz in q)))}
            log.info("  %-16s Δμ (kwadr. %d kat.) %.3f", "cat_Y", len(q), out["sources"]["cat_Y"]["dmu_sym"])
        elif name == "kappa" and kappa_alts:
            # prawda z planu B (κ per kategoria z ramion) → fit alternatywnym modelem κ
            for alt, kw in (("kappa_channel", {"kappa": "channel"}), ("kappa_partner03", {"kappa": "partner", "partner_sigma": 0.3}),
                            ("kappa_partner05", {"kappa": "partner", "partner_sigma": 0.5})):
                cfa = channel(M, ch, **kw)
                xa = transfer(x, cf.names, cfa)
                if kw["kappa"] == "channel":
                    kv = [x[k] for k, nm in enumerate(cf.names) if nm.startswith("kap")]
                    for k, nm in enumerate(cfa.names):
                        if nm in ("kSB", "kBS"):
                            xa[k] = float(np.mean(kv)) if kv else 1.0
                elif kw["kappa"] == "partner":
                    kd = {nm: x[k] for k, nm in enumerate(cf.names) if nm.startswith("kap")}
                    for k, nm in enumerate(cfa.names):
                        if nm.startswith("kap_"):
                            xa[k] = kd.get(nm.rsplit("_p", 1)[0], 1.0)
                set_centres(cfa, xa)
                # nominał alternatywnego modelu na własnym Asimovie (μ=1) i na prawdzie planu B
                a_self = asimov_fit(cfa, cfa.lam(xa), xa)
                a_b = asimov_fit(cfa, lam0, xa)
                # kierunek odwrotny: prawda z alternatywnego κ dopasowanego na obszarze ślepym → fit planem B
                rev = None
                if getattr(cf, "_blind_n", None) is not None:
                    mr = cfa.fit(cf._blind_n, region_mask(M.I, ch, "blind"), fix_mu=1.0, start=xa)
                    lam_r = cfa.lam(np.array(mr.values))
                    a_r = asimov_fit(cf, lam_r, x)
                    rev = {"dmu": a_r["mu"] - nom["mu"], "dZ": a_r["Z"] - nom["Z"], "Z": a_r["Z"], "err": a_r["err"], "valid": a_r["valid"]}
                rec_ = {"up": {"dmu": a_b["mu"] - nom["mu"], "dZ": a_b["Z"] - nom["Z"], "Z": a_b["Z"], "err": a_b["err"], "valid": a_b["valid"]},
                        "self": {"mu": a_self["mu"], "Z": a_self["Z"], "err": a_self["err"]}}
                if rev:
                    rec_["down"] = rev            # 'down' = prawda alternatywna, fit planem B
                ds = [v["dmu"] for k, v in rec_.items() if k in ("up", "down")]
                dz = [v["dZ"] for k, v in rec_.items() if k in ("up", "down")]
                rec_.update({"dmu_sym": float(max(abs(d) for d in ds)), "dZ_min": float(min(dz))})
                out["sources"][alt] = rec_
                log.info("  %-16s Δμ planB→alt %+.3f, alt→planB %s; Z własny Asimov %.2f", alt, ds[0],
                         f"{ds[1]:+.3f}" if len(ds) > 1 else "-", a_self["Z"])
    # tło poly3 jako model nominalny: fit modelem z th2 (α kształtu sygnału więzione) do prawdy poly3
    if getattr(cf, "_poly3", None) is not None:
        cfs3, xs3, lam3 = cf._poly3
        set_centres(cfs3, xs3)
        a3 = asimov_fit(cfs3, lam3, xs3)
        out["sources"]["bkg_poly3"]["poly3_model_on_poly3_truth"] = {"mu": a3["mu"], "err": a3["err"], "Z": a3["Z"]}
        n = cf._blind_n
        mk = region_mask(M.I, ch, "blind")
        m2 = cf.fit(n, mk, fix_mu=1.0, start=x)
        fix = [nm for nm in cfs3.names if nm.startswith(("aw_", "as_", "atag"))]
        m3 = cfs3.fit(n, mk, fix_mu=1.0, start=xs3, fix_names=fix)
        d2, nd2 = cf.deviance(m2, n, mk)
        d3, nd3 = cfs3.deviance(m3, n, mk)
        out["sources"]["bkg_poly3"]["blind_dev"] = {"poly2": [d2, nd2], "poly3": [d3, nd3], "ddev": d2 - d3, "dndf": nd2 - nd3}
        log.info("  tło poly3: Δdev(ślepo) = %.1f / %d st. sw.; model poly3 na prawdzie poly3: μ=%.3f σ=%.3f Z=%.2f",
                 d2 - d3, nd2 - nd3, a3["mu"], a3["err"], a3["Z"])
    # profilowane systematyki kształtu (PairFit syst=True) na nominalnym Asimovie
    cfs = channel(M, ch, syst=True)
    xs = transfer(x, cf.names, cfs)
    set_centres(cfs, xs)
    a = asimov_fit(cfs, lam0, xs)
    out["profiled_shape_syst"] = {"mu": a["mu"], "err": a["err"], "Z": a["Z"]}
    log.info("  profil syst. kształtu: σ_μ %.3f→%.3f, Z %.2f→%.2f", nom["err"], a["err"], nom["Z"], a["Z"])
    # bez kategorii z tagiem D* (czułość)
    return out
