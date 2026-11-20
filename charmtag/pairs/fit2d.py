"""Fit 2D par: tryb ślepy (domyślnie), odsłonięty (--unblind), Asimov, toys, systematyki, domknięcie."""
from __future__ import annotations

import json
import logging
import os
import time

import numpy as np

from . import CHANNELS
from .model import Model, region_mask

log = logging.getLogger("charmtag")
GRID = np.round(np.arange(-1.5, 5.01, 0.25), 2)
FINE = np.linspace(-1.5, 5.0, 2601)


def _js(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


# ----------------------------------------------------------------------------- narzędzia
def profile(cf, n, mask, m1, grid=GRID):
    x1 = np.array(m1.values)
    scan = {}
    for mu_s in grid:
        ms = cf.fit(n, mask, fix_mu=float(mu_s), start=x1)
        scan[float(mu_s)] = float(ms.fval - m1.fval)
    mus = np.array(sorted(scan))
    d = np.array([scan[k] for k in mus])
    d = d - min(d.min(), 0.0)
    df = np.interp(FINE, mus, d)
    ok = FINE[df <= 0.5]
    return {f"{k:g}": v for k, v in zip(mus, d)}, ((float(ok.min()), float(ok.max())) if len(ok) else (np.nan, np.nan))


def combine(profiles: dict):
    """Wspólne μ: suma krzywych profilu kanałów (jak stageB_post.py)."""
    tot = np.zeros_like(FINE)
    for scan in profiles.values():
        mus = np.array(sorted(float(k) for k in scan))
        d = np.array([scan[f"{m:g}"] for m in mus])
        tot += np.interp(FINE, mus, d)
    tot -= tot.min()
    mu_c = float(FINE[np.argmin(tot)])
    ok = FINE[tot <= 0.5]
    q0 = 2 * float(np.interp(0.0, FINE, tot)) if mu_c > 0 else 0.0
    return {"mu_common": mu_c, "mu_common_68": [float(ok.min()), float(ok.max())], "q0": q0, "Z": float(np.sqrt(max(q0, 0)))}


def z_from(m1, m0):
    q0 = 2 * (m0.fval - m1.fval) if m1.values["mu"] > 0 else 0.0
    return float(np.sqrt(max(q0, 0.0)))


def fit_channel(M, ch, n, mask, kappa="planB", syst=False, psig=None, scan=True):
    cf = M.channel(ch, kappa=kappa, syst=syst, partner_sigma=psig)
    m1 = cf.fit(n, mask)
    m0 = cf.fit(n, mask, fix_mu=0.0, start=np.array(m1.values))
    r = {"mu_hat": float(m1.values["mu"]), "mu_err_hesse": float(m1.errors["mu"]), "Z": z_from(m1, m0), "valid": bool(m1.valid),
         "S_expected_full_plane": float(sum(v["SS"] for v in cf.exp.values()))}
    dev, ndf = cf.deviance(m1, n, mask)
    r.update({"dev": dev, "ndf": ndf, "dev_ndf": dev / max(ndf, 1)})
    r["kappa"] = {nm: [float(m1.values[nm]), float(m1.errors[nm])] for nm in cf.names if nm.startswith("kap") or nm in ("kSB", "kBS")}
    if syst:
        r["alphas"] = {a: [float(m1.values[a]), float(m1.errors[a])] for a in cf.alphas}
    if scan:
        r["profile_dNLL"], r["mu_profile_68"] = profile(cf, n, mask, m1)
    return r, cf, m1


# ----------------------------------------------------------------------------- Asimov / toys
def asimov(M, kappa="planB", gen=None, kk=1.0, psig=None):
    """Oczekiwana czułość: dane Asimov z modelu nominalnego (μ=1) lub ze zmienionej prawdy `gen`, fit nominalny."""
    res = {}
    for ch in CHANNELS:
        cf = M.channel(ch, kk=kk, kappa=kappa, partner_sigma=psig)
        truth = _truth(M, cf, gen)
        full = np.ones(truth.shape[1:], bool)
        m1 = cf.fit(truth, full)
        m0 = cf.fit(truth, full, fix_mu=0.0, start=np.array(m1.values))
        res[ch] = {"mu": float(m1.values["mu"]), "mu_err": float(m1.errors["mu"]),
                   "Z": float(np.sqrt(max(2 * (m0.fval - m1.fval), 0))), "q0": float(2 * (m0.fval - m1.fval))}
    res["Z_quadrature"] = float(np.sqrt(sum(v["Z"] ** 2 for k, v in res.items() if k in CHANNELS)))
    return res


def _truth(M, cf, gen, x=None):
    S = M.S
    x = cf.x0() if x is None else x
    if not gen or gen == "nom":
        return cf.lam(x)
    if gen.startswith("sigw"):
        k = 1.0 + float(gen[4:] or 10) / 100
        return cf.lam(x, sig={lg: S.stretch(M.I.sig[lg], lg[0], k=k) for lg in M.I.sig})
    if gen.startswith("mshift"):
        return cf.lam(x, sig={lg: S.stretch(M.I.sig[lg], lg[0], shift=float(gen[6:] or 1)) for lg in M.I.sig})
    if gen == "bkg3":
        return cf.lam(x, fshapes={lg: S.fshape(lg[0], S.TH3[lg]) for lg in cf.legs})
    if gen.startswith("swap"):
        if not M.I.sig_comp:
            raise ValueError("brak składowych correct/swap w wejściach (sig_comp)")
        fsw = 1.0 + float(gen[4:]) / 100
        sig, ratio = dict(M.I.sig), {}
        for c in range(M.I.ncat[10]):
            cc, ss = M.I.sig_comp[(10, c)]
            t = cc + fsw * ss
            sig[(10, c)] = t / t.sum()
            ratio[c] = t.sum() / (cc + ss).sum()
        scale = {(i, j): (ratio[i] if cf.a == 10 else 1.0) * (ratio[j] if cf.b == 10 else 1.0) for i, j in cf.cps}
        return cf.lam(x, sig=sig, ss_scale=scale)
    raise ValueError(gen)


def toys(M, ch, lam_true, ntoy, seed, mu_true=1.0, kappa="planB", mask=None):
    cf = M.channel(ch, kappa=kappa)
    rng = np.random.default_rng(seed)
    mask = np.ones(lam_true.shape[1:], bool) if mask is None else mask
    mus, errs, zs, cov = [], [], [], 0
    for _ in range(ntoy):
        n = rng.poisson(np.clip(lam_true, 0, None)).astype(float)
        m1 = cf.fit(n, mask)
        m0 = cf.fit(n, mask, fix_mu=0.0, start=np.array(m1.values))
        mt = cf.fit(n, mask, fix_mu=mu_true, start=np.array(m1.values))
        mus.append(float(m1.values["mu"]))
        errs.append(float(m1.errors["mu"]))
        zs.append(z_from(m1, m0))
        cov += int(2 * (mt.fval - m1.fval) < 1.0)
    mus, errs, zs = map(np.array, (mus, errs, zs))
    pulls = (mus - mu_true) / np.where(errs > 0, errs, np.nan)
    return {"n": ntoy, "mu_true": mu_true, "mu_mean": float(mus.mean()), "bias": float(mus.mean() - mu_true),
            "sigma_mu_median": float(np.median(errs)), "pull_mean": float(np.nanmean(pulls)), "pull_std": float(np.nanstd(pulls)),
            "Z_median": float(np.median(zs)), "Z_16_84": [float(np.quantile(zs, 0.16)), float(np.quantile(zs, 0.84))],
            "coverage68": cov / max(ntoy, 1)}


# ----------------------------------------------------------------------------- tryby
def run_blind(M, H, out, ntoy=0, seed=2026, gens=("sigw3", "sigw10", "mshift1", "swap20", "bkg3")):
    """Ślepo: histogramy bez centrum (ślepota przy budowie). μ nigdy nie jest wypisywane ani zapisywane.
    Fit nuisance'ów na sidebandach+ramionach przy μ=1 (stała, jak etap A), jakość fitu, κ, czułość Asimov z nuisance'ami z ramion."""
    rep = {"mode": "blind", "channels": {}}
    for ch in CHANNELS:
        n = H[ch]
        mask = region_mask(M.I, ch, "blind")
        if n.sum() < 20:
            rep["channels"][ch] = {"skipped": f"za mało par w obszarze ślepym ({n.sum():.0f})"}
            log.warning("fit2d BLIND %s: %s", ch, rep["channels"][ch]["skipped"])
            continue
        cf = M.channel(ch)
        m = cf.fit(n, mask, fix_mu=1.0)
        dev, ndf = cf.deviance(m, n, mask)
        da, na = cf.deviance(m, n, np.broadcast_to(region_mask(M.I, ch, "arms"), n.shape[1:]) & mask)
        x = np.array(m.values)
        lam = cf.lam(x)
        full = np.ones(n.shape[1:], bool)
        ma = cf.fit(lam, full, start=x)
        m0 = cf.fit(lam, full, fix_mu=0.0, start=x)
        r = {"blind_region_dev_ndf": dev / max(ndf, 1), "arms_dev_per_bin": da / max(na, 1), "n_pairs_blind_region": float(n[np.broadcast_to(mask, n.shape)].sum()),
             "S_expected_full_plane": float(sum(v["SS"] for v in cf.exp.values())),
             "kappa": {nm: [float(m.values[nm]), float(m.errors[nm])] for nm in cf.names if nm.startswith("kap")},
             "asimov_from_arms": {"Z": float(np.sqrt(max(2 * (m0.fval - ma.fval), 0))), "sigma_mu": float(ma.errors["mu"])}}
        if ntoy:
            r["toys"] = toys(M, ch, lam, ntoy, seed)
        rep["channels"][ch] = r
        log.info("fit2d BLIND %s: dev/ndf=%.2f ramiona dev/bin=%.2f, Z_Asimov(ramiona)=%.2f σ_μ=%.2f", ch, r["blind_region_dev_ndf"],
                 r["arms_dev_per_bin"], r["asimov_from_arms"]["Z"], r["asimov_from_arms"]["sigma_mu"])
    rep["Z_asimov_quadrature"] = float(np.sqrt(sum(v["asimov_from_arms"]["Z"] ** 2 for v in rep["channels"].values() if "skipped" not in v)))
    rep["asimov_nominal"] = asimov(M)
    rep["asimov_shifts"] = {}
    for g in gens:
        try:
            a = asimov(M, gen=g)
            rep["asimov_shifts"][g] = {ch: a[ch]["mu"] - 1.0 for ch in CHANNELS}
        except ValueError as e:
            rep["asimov_shifts"][g] = str(e)
    for sig in (0.3, 0.5):
        rep[f"asimov_partner_sigma{sig}"] = asimov(M, kappa="partner", psig=sig)
    _js(os.path.join(out, "fit2d_blind.json"), rep)
    return rep


def run_unblind(M, H, out, ntoy=100, seed=2026, extra=True):
    rep = {"mode": "unblind", "channels": {}}
    prof = {}
    for ch in CHANNELS:
        t0 = time.time()
        n = H[ch]
        full = np.ones(n.shape[1:], bool)
        r, cf, m1 = fit_channel(M, ch, n, full)
        prof[ch] = r["profile_dNLL"]
        if ntoy:
            r["toys_at_fit"] = toys(M, ch, cf.lam(np.array(m1.values)), ntoy, seed, mu_true=r["mu_hat"])
            r["toys_at_fit"]["mu_corrected"] = r["mu_hat"] - r["toys_at_fit"]["bias"]
        if extra:
            c, _, _ = fit_channel(M, ch, n, full, kappa="channel", scan=False)
            r["crosscheck_kchannel"] = {k: c[k] for k in ("mu_hat", "mu_err_hesse", "Z", "dev_ndf")}
            # bez kategorii z tagiem D*
            notag = np.array([(i < 3 if cf.a == 10 else True) and (j < 3 if cf.b == 10 else True) for i, j in cf.cps])[:, None, None] & full[None]
            tagc = [f"{i}{j}" for i, j in cf.cps if (cf.a == 10 and i >= 3) or (cf.b == 10 and j >= 3)]
            fixn = [nm for nm in cf.names if (nm[:3] in ("sb_", "bs_", "bb_") and nm[3:] in tagc) or
                    nm.startswith(("th0_10_3", "th1_10_3", "th0_10_4", "th1_10_4", "kap_10_3", "kap_10_4"))]
            mnt = cf.fit(n, notag, fix_names=fixn)
            r["mu_no_tag"] = [float(mnt.values["mu"]), float(mnt.errors["mu"])]
        rep["channels"][ch] = r
        log.info("fit2d UNBLIND %s: μ=%.3f ± %.3f [%.2f, %.2f] Z=%.2f dev/ndf=%.2f (%.0f s)", ch, r["mu_hat"], r["mu_err_hesse"],
                 *r["mu_profile_68"], r["Z"], r["dev_ndf"], time.time() - t0)
    rep["combined"] = combine(prof)
    log.info("fit2d UNBLIND wspólne μ=%.2f %s Z=%.2f", rep["combined"]["mu_common"], rep["combined"]["mu_common_68"], rep["combined"]["Z"])
    if extra:
        rep["syst_profiled"] = {}
        sp = {}
        for ch in CHANNELS:
            n = H[ch]
            r, _, _ = fit_channel(M, ch, n, np.ones(n.shape[1:], bool), syst=True)
            sp[ch] = r["profile_dNLL"]
            rep["syst_profiled"][ch] = {k: r[k] for k in ("mu_hat", "mu_profile_68", "Z", "alphas")}
        rep["syst_profiled"]["combined"] = combine(sp)
        rep["kappa_partner"] = {}
        for sig in (0.3, 0.5):
            pp = {}
            row = {}
            for ch in CHANNELS:
                n = H[ch]
                r, _, _ = fit_channel(M, ch, n, np.ones(n.shape[1:], bool), kappa="partner", psig=sig)
                pp[ch] = r["profile_dNLL"]
                row[ch] = {k: r[k] for k in ("mu_hat", "Z")}
            row["combined"] = combine(pp)
            rep["kappa_partner"][f"{sig:g}"] = row
    _js(os.path.join(out, "fit2d_unblind.json"), rep)
    return rep


REFERENCE_0927 = {"D0D0": {"mu": 1.37, "Z": 1.9}, "DpD0": {"mu": 1.62, "Z": 2.9}, "DpDp": {"mu": 1.57, "Z": 1.0}, "combined_Z": 3.6, "mu_common": 1.50}


def run_closure(freeze_dir, out, hists_npz=None, ntoy=0):
    """Domknięcie: te same wejścia co 09-27 (zamrożone), plan B, pełna płaszczyzna (wynik już odsłonięty 27.09)."""
    from .frozen import load_frozen, frozen_legs_all
    from .hists import build_pair_hists
    I = load_frozen(freeze_dir)
    M = Model(I)
    H, counts = build_pair_hists(frozen_legs_all(freeze_dir), I, "full")
    same = None
    if hists_npz and os.path.exists(hists_npz):
        old = np.load(hists_npz)
        same = {ch: bool(np.array_equal(old[ch], H[ch])) for ch in CHANNELS if ch in old}
    rep = run_unblind(M, H, out, ntoy=ntoy, extra=False)
    chk = {}
    for ch in CHANNELS:
        r = rep["channels"][ch]
        chk[ch] = {"mu": round(r["mu_hat"], 3), "ref_mu": REFERENCE_0927[ch]["mu"], "Z": round(r["Z"], 2), "ref_Z": REFERENCE_0927[ch]["Z"],
                   "ok": abs(r["mu_hat"] - REFERENCE_0927[ch]["mu"]) < 0.006 and abs(r["Z"] - REFERENCE_0927[ch]["Z"]) < 0.051}
    c = rep["combined"]
    chk["combined"] = {"Z": round(c["Z"], 2), "ref_Z": 3.6, "mu": round(c["mu_common"], 3), "ref_mu": 1.50,
                       "ok": abs(c["Z"] - 3.6) < 0.051 and abs(c["mu_common"] - 1.50) < 0.006}
    res = {"hists_equal_stageB": same, "pairs": counts, "check": chk, "PASS": all(v["ok"] for v in chk.values()) and (same is None or all(same.values()))}
    _js(os.path.join(out, "closure.json"), {**res, "fit": rep})
    log.info("DOMKNIĘCIE 09-27: %s %s", "PASS" if res["PASS"] else "FAIL", json.dumps(chk))
    return res
