"""Warianty WP / kategorii score / tagu D* / arbitrażu: przebudowa par (charmtag.pairs.build.cmd_pairs) do osobnych
katalogów <out_dir>/pairs_syst/<wariant>/ i porównanie ŚLEPE (etap U).

Ślepo: Z_Asimov i σ_μ wariantu (nuisance'y z jego obszaru ślepego) vs nominał; oczekiwany rozrzut statystyczny
Δμ_dane ≈ √|σ_w² − σ_n²| (podpróbki zagnieżdżone). Δμ_dane — dopiero po odsłonięciu (fit2d --unblind w katalogu wariantu).
"""
from __future__ import annotations

import copy
import json
import logging
import os

import numpy as np

from ..pairs import CHANNELS
from ..pairs.model import Model, load_inputs
from . import core

log = logging.getLogger("charmtag")

VARIANTS = {
    "wp_tight": {"cats": [[0.0, 0.25], [0.25, 0.40], [0.40, 0.60]]},
    "wp_loose": {"cats": [[0.0, 0.35], [0.35, 0.50], [0.50, 0.65]]},
    "arb85": {"d0_arbitration_wp": 0.85},
    "arb95": {"d0_arbitration_wp": 0.95},
    "tagdm07": {"tag_dm": [0.14543, 0.0007]},
    "tagdm15": {"tag_dm": [0.14543, 0.0015]},
}
EXTRA_EFFS = {0.25, 0.35, 0.40, 0.50, 0.60, 0.85, 0.95}


def vdir(cfg, name):
    return os.path.join(cfg["paths"]["out_dir"], "pairs_syst", name)


def build_variant(cfg, name, norm_mode="direct"):
    """Przebudowa par dla wariantu (ślepo). Nie zmienia nominalnego <out_dir>/pairs."""
    from ..pairs import build as B
    cfg2 = copy.deepcopy(cfg)
    cfg2.setdefault("pairs", {})
    cfg2["pairs"] = dict(cfg2["pairs"] or {})
    cfg2["pairs"].update(VARIANTS[name])
    out = vdir(cfg, name)
    os.makedirs(out, exist_ok=True)
    orig_pdir, orig_thr = B.pdir, B.thresholds
    B.pdir = lambda c: out
    B.thresholds = lambda sig, effs: orig_thr(sig, sorted(set(effs) | EXTRA_EFFS))
    try:
        B.cmd_pairs(cfg2, unblind=False, norm_mode=norm_mode)
    finally:
        B.pdir, B.thresholds = orig_pdir, orig_thr
    with open(os.path.join(out, "variant.json"), "w") as f:
        json.dump({"name": name, "pairs_override": VARIANTS[name]}, f, indent=1)
    return out


def sensitivity(pdir_):
    I = load_inputs(os.path.join(pdir_, "pair_inputs"))
    H = dict(np.load(os.path.join(pdir_, "pairs_blind.npz")))
    M = Model(I)
    r = {}
    for ch in CHANNELS:
        if H[ch].sum() < 20:
            continue
        cf, x, src = core.base_point(M, H, ch)
        a = core.asimov_fit(cf, cf.lam(x), x)
        dev, ndf = cf.deviance(cf.fit(H[ch], np.ones(H[ch].shape[1:], bool) & core.region_mask(I, ch, "blind"), fix_mu=1.0, start=x),
                               H[ch], core.region_mask(I, ch, "blind"))
        r[ch] = {"Z": a["Z"], "err": a["err"], "S_exp": float(sum(v["SS"] for v in cf.exp.values())), "dev_ndf": dev / max(ndf, 1),
                 "n_blind": float(H[ch].sum())}
    r["Z_quadrature"] = float(np.sqrt(sum(v["Z"] ** 2 for k, v in r.items() if k in CHANNELS)))
    return r


def compare(cfg, nominal_dir):
    nom = sensitivity(nominal_dir)
    rep = {"nominal": nom, "variants": {}}
    for name in VARIANTS:
        d = vdir(cfg, name)
        if not os.path.exists(os.path.join(d, "pairs_blind.npz")):
            rep["variants"][name] = {"missing": d}
            continue
        v = sensitivity(d)
        for ch in CHANNELS:
            if ch in v and ch in nom:
                v[ch]["dZ"] = v[ch]["Z"] - nom[ch]["Z"]
                v[ch]["exp_stat_spread_dmu"] = float(np.sqrt(abs(v[ch]["err"] ** 2 - nom[ch]["err"] ** 2)))
        rep["variants"][name] = v
        log.info("syst wp %s: Z kwadr. %.2f (nominał %.2f)", name, v["Z_quadrature"], nom["Z_quadrature"])
    return rep
