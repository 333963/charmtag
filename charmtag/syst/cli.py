"""`charmtag syst` — systematyki par D–D, ŚLEPO (etap U).

  charmtag syst -c CFG core [--channels D0D0,DpD0,DpDp]     # kształty, SF, tag, κ, kategorie → syst/core_<kanał>.json
  charmtag syst -c CFG norm [--truth-dir DIR]               # normalizacja (tune, n_trk, dyfrakcja, stat. MC) → syst/norm.json
  charmtag syst -c CFG files [--no-jackknife] [--no-splits] # jackknife po plikach, podziały (ślepo) → syst/files.json
  charmtag syst -c CFG wp-build --variant NAME              # przebudowa par dla wariantu WP/kategorii → <out>/pairs_syst/NAME
  charmtag syst -c CFG wp                                   # porównanie wariantów (ślepo) → syst/wp.json
  charmtag syst -c CFG salad                                # B4: ε_B pod pikiem (SALAD/B1) per gatunek + pary → syst/bkg_shape_salad_*.json
  charmtag syst -c CFG report                               # tabela → syst/syst_blind.{md,json}
  charmtag syst -c CFG all                                  # core + norm + files + wp (jeśli warianty są) + report
Wejścia: <out_dir>/pairs/{pair_inputs.*, pairs_blind.npz, legs_data.parquet}. Nic nie jest odsłaniane.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time

import numpy as np

log = logging.getLogger("charmtag")


def _js(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else (o.item() if hasattr(o, "item") else str(o)))
    os.replace(path + ".tmp", path)


def _clean(r):
    if isinstance(r, dict):
        return {k: _clean(v) for k, v in r.items() if k not in ("x", "cf")}
    return r


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "syst":
        argv = argv[1:]
    ap = argparse.ArgumentParser(prog="charmtag syst", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("part", choices=["core", "norm", "files", "wp-build", "wp", "salad", "report", "all"])
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--out-dir", default=None, help="nadpisuje paths.out_dir")
    ap.add_argument("--pairs-dir", default=None, help="katalog wejść par (domyślnie <out_dir>/pairs)")
    ap.add_argument("--syst-dir", default=None, help="katalog wyników (domyślnie <pairs-dir>/syst)")
    ap.add_argument("--channels", default="D0D0,DpD0,DpDp")
    ap.add_argument("--truth-dir", default=None, help="katalog truth (combined_A3/, combined_syst_<wariant>/)")
    ap.add_argument("--variant", default=None)
    ap.add_argument("--norm", default="direct", choices=["direct", "frozen", "truth"])
    ap.add_argument("--no-jackknife", action="store_true")
    ap.add_argument("--no-splits", action="store_true")
    a = ap.parse_args(argv)

    from .. import config as C
    from ..cli import _setup_log
    from ..pairs.build import pdir
    ov = {"paths": {}}
    if a.out_dir:
        ov["paths"]["out_dir"] = a.out_dir
    cfg = C.load(a.config, ov, require_paths=True)
    pd_ = os.path.expanduser(a.pairs_dir) if a.pairs_dir else pdir(cfg)
    sd = os.path.expanduser(a.syst_dir) if a.syst_dir else os.path.join(pd_, "syst")
    os.makedirs(sd, exist_ok=True)
    _setup_log(sd)
    t0 = time.time()
    log.info("charmtag syst %s pairs=%s out=%s", a.part, pd_, sd)
    from . import core, norm, report
    from ..pairs.model import Model, load_inputs
    S = core.settings(cfg)
    chs = [c for c in a.channels.split(",") if c]
    truth_dir = a.truth_dir or (cfg.get("syst", {}) or {}).get("truth_dir")

    def do_core():
        I = load_inputs(os.path.join(pd_, "pair_inputs"))
        M = Model(I)
        H = dict(np.load(os.path.join(pd_, "pairs_blind.npz")))
        for ch in chs:
            r = core.run_channel(M, H, ch, S)
            r["settings"] = {k: (v if not isinstance(v, dict) else {str(kk): vv for kk, vv in v.items()}) for k, v in S.items()}
            _js(os.path.join(sd, f"core_{ch}.json"), _clean(r))

    def do_norm():
        I = load_inputs(os.path.join(pd_, "pair_inputs"))
        r = norm.collect(I.meta, truth_dir)
        _js(os.path.join(sd, "norm.json"), r)
        for k, v in r.items():
            log.info("syst norm %-16s δ %s  (%s)", k, {c: round(x, 3) for c, x in v["delta"].items()}, v.get("how"))

    def do_files():
        import pandas as pd
        from . import files
        I = load_inputs(os.path.join(pd_, "pair_inputs"))
        L = pd.read_parquet(os.path.join(pd_, "legs_data.parquet"))
        try:
            ft = files.file_table(cfg)
        except Exception as e:  # brak h5 (np. maszyna bez danych) → pomiń
            log.warning("syst files: brak tabeli plików (%s)", e)
            ft = pd.DataFrame([{"file": "?", "job": "?", "rows": len(L), "events": None}])
        r = files.run_files(cfg, I, L, ft, jackknife=not a.no_jackknife, splits=not a.no_splits, channels=chs)
        _js(os.path.join(sd, "files.json"), _clean(r))

    def do_wp():
        from . import wp
        _js(os.path.join(sd, "wp.json"), wp.compare(cfg, pd_))

    def do_salad():                                       # B4: informacyjne (B_pred pod pikiem), nie zmienia Δμ fit2d
        from . import bkg_shape_salad as bss
        sc = cfg.get("syst", {}) or {}
        s8 = {sp: sc.get(f"s8_json_{sp}") or os.path.join(cfg["paths"]["out_dir"], "syst", f"bkg_shape_salad_{sp}.json")
              for sp in ("D0", "Dplus")}
        blk = {sp: bss.build(sp, s8_json=s8[sp], d0_truth_min2=bool(sc.get("salad_d0_truth_min2", True))) for sp in s8}
        for sp, b in blk.items():
            _js(os.path.join(sd, f"bkg_shape_salad_syst_{sp}.json"), b)
        _js(os.path.join(sd, "bkg_shape_salad_pairs.json"), bss.pairs_block(blk["D0"], blk["Dplus"]))

    if a.part == "core":
        do_core()
    elif a.part == "norm":
        do_norm()
    elif a.part == "files":
        do_files()
    elif a.part == "wp-build":
        from . import wp
        if a.variant not in wp.VARIANTS:
            raise SystemExit(f"--variant: {sorted(wp.VARIANTS)}")
        wp.build_variant(cfg, a.variant, norm_mode=a.norm)
    elif a.part == "wp":
        do_wp()
    elif a.part == "salad":
        do_salad()
    elif a.part == "all":
        do_core()
        do_norm()
        do_files()
        from . import wp
        if any(os.path.exists(os.path.join(wp.vdir(cfg, v), "pairs_blind.npz")) for v in wp.VARIANTS):
            do_wp()
        if (cfg.get("syst") or {}).get("bkg_shape_salad"):
            do_salad()
    if a.part in ("report", "all"):
        md, _ = report.build(sd)
        if a.part in ("report", "all"):
            print(md)
    log.info("koniec syst %s (%.0f s)", a.part, time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
