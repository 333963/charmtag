"""charmtag CLI: inspect | train | apply | fit | validate | report | all | synth."""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time

from . import config as C

log = logging.getLogger("charmtag")


def _setup_log(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    log.setLevel(logging.INFO)
    if not log.handlers:
        h = logging.StreamHandler(sys.stdout)
        h.setFormatter(fmt)
        log.addHandler(h)
        fh = logging.FileHandler(os.path.join(out_dir, "charmtag.log"))
        fh.setFormatter(fmt)
        log.addHandler(fh)


def cmd_inspect(cfg):
    from .io import list_columns, n_rows, resolve_columns
    rep = {}
    for key in cfg["paths"]["samples"]:
        p = C.sample_path(cfg, key)
        if not os.path.exists(p):
            rep[key] = {"exists": False, "path": p}
            continue
        src = C.source(cfg, key)  # etap P: wiele plików
        cols = list_columns(p, cfg["io"]["dataset"])
        cm = resolve_columns(cols, cfg["columns"], 3)
        rep[key] = {"exists": True, "path": p, "rows": n_rows(src, cfg["io"]["dataset"]), "n_columns": len(cols),
                    "mapping": cm.mapping, "missing": cm.missing,
                    "unmapped_columns": sorted(set(cols) - {v for v in cm.mapping.values() if v})}
        if C.is_multi(cfg, key):
            rep[key]["sample_meta"] = json.loads(C.sample_meta_json(cfg, key))
        log.info("%s: %d wierszy, %d kolumn; brak logicznych: %s", key, rep[key]["rows"], len(cols), cm.missing)
        log.info("  niezmapowane: %s", rep[key]["unmapped_columns"])
    out = os.path.join(cfg["paths"]["out_dir"], "inspect.json")
    with open(out, "w") as f:
        json.dump(rep, f, indent=1)
    return rep


def main(argv=None):
    _argv = list(sys.argv[1:] if argv is None else argv)
    if _argv and _argv[0] in ("pairs", "fit2d"):      # etap O: moduł par (charmtag/pairs), addytywnie
        from .pairs.cli import main as _pairs_main
        return _pairs_main(_argv)
    if _argv and _argv[0] == "s5compare":            # etap S5: porównanie wariantów tła (charmtag/s5cmp.py), addytywnie
        from .s5cmp import main as _s5_main
        return _s5_main(_argv)
    if _argv and _argv[0] == "syst":                 # etap U: systematyki par (charmtag/syst), addytywnie
        from .syst.cli import main as _syst_main
        return _syst_main(_argv)
    if _argv and _argv[0] == "bkggen":               # etap S (+S2/S4/S6, scalenie AA): generator tła i próbki kontrolne
        from .bkggen.cli import main as _bkggen_main
        return _bkggen_main(_argv)
    if _argv and _argv[0] == "candrot":              # etap S6: kontrole na poziomie kandydata (= bkggen --method candrot)
        from .bkggen.candrot import main as _candrot_main
        return _candrot_main(_argv[1:])
    ap = argparse.ArgumentParser(prog="charmtag", description=__doc__)
    ap.add_argument("command", choices=["inspect", "train", "apply", "fit", "validate", "report", "all", "synth"])
    ap.add_argument("-c", "--config", default=None, help="YAML (nadpisuje config/default.yaml)")
    ap.add_argument("--species", default="D0", choices=["D0", "Dplus"])
    ap.add_argument("--variant", default=None, choices=["H", "frozen7"])
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--h5-dir", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--n-estimators", type=int, default=None)
    a = ap.parse_args(argv)

    ov = {"paths": {}, "io": {}, "training": {"xgb": {}}}
    if a.h5_dir:
        ov["paths"]["h5_dir"] = a.h5_dir
    if a.out_dir:
        ov["paths"]["out_dir"] = a.out_dir
    if a.max_rows:
        ov["io"]["max_rows"] = a.max_rows
    if a.n_estimators:
        ov["training"]["xgb"]["n_estimators"] = a.n_estimators
    cfg = C.load(a.config, ov)
    if a.command == "synth":
        from .synth import make
        print(make(cfg["paths"]["h5_dir"], n=a.max_rows or 20000))
        return 0
    _setup_log(cfg["paths"]["out_dir"])
    from . import pipeline as P
    from . import report as R
    t0 = time.time()
    log.info("charmtag %s species=%s variant=%s h5=%s out=%s", a.command, a.species, a.variant,
             cfg["paths"]["h5_dir"], cfg["paths"]["out_dir"])
    if a.command == "inspect":
        cmd_inspect(cfg)
    elif a.command == "train":
        P.cmd_train(cfg, a.species, a.variant)
    elif a.command == "apply":
        P.cmd_apply(cfg, a.species, a.variant, a.samples)
    elif a.command == "fit":
        P.cmd_fit(cfg, a.species, a.variant)
    elif a.command == "validate":
        P.cmd_validate(cfg, a.species)
    elif a.command == "report":
        md, _ = R.make(cfg, a.species)
        print(md)
    elif a.command == "all":
        cmd_inspect(cfg)
        variants = [a.variant] if a.variant else ["H", "frozen7"]
        for v in variants:
            P.cmd_train(cfg, a.species, v)
            P.cmd_apply(cfg, a.species, v, a.samples)
            P.cmd_fit(cfg, a.species, v)
        P.cmd_validate(cfg, a.species)
        md, _ = R.make(cfg, a.species)
        print(md)
    log.info("koniec %s (%.0f s)", a.command, time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
