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


def _pairs():
    from .pairs.cli import main
    return main


def _syst():
    from .syst.cli import main
    return main


def _bkggen():
    from .bkggen.cli import main
    return main


def _s5compare():
    from .s5cmp import main
    return main


# Podkomendy z własnym parserem (każda dostaje pełne argv, łącznie z nazwą podkomendy).
# Projekt A (double charm): pairs, fit2d, syst, s5compare.  Projekt B (tło podobne do sygnału): bkggen (także `charmtag-bkggen`).
def _s8salad():
    from .s8salad import main
    return main


_SUBCOMMANDS = {"pairs": _pairs, "fit2d": _pairs, "syst": _syst, "s5compare": _s5compare, "bkggen": _bkggen, "s8salad": _s8salad}
_EPILOG = """podkomendy z własną pomocą:  charmtag pairs -h | fit2d -h | syst -h | s5compare -h | s8salad -h | bkggen -h  (= charmtag-bkggen -h)"""


def bkggen_main(argv=None):
    """Wejście konsolowe projektu B: `charmtag-bkggen --method ...` (= `charmtag bkggen ...`)."""
    from .bkggen.cli import main
    return main(list(sys.argv[1:] if argv is None else argv))


def main(argv=None):
    _argv = list(sys.argv[1:] if argv is None else argv)
    sub = _SUBCOMMANDS.get(_argv[0]) if _argv else None
    if sub is not None:
        return sub()(_argv)
    ap = argparse.ArgumentParser(prog="charmtag", description=__doc__, epilog=_EPILOG,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["inspect", "train", "apply", "fit", "validate", "report", "all", "synth"])
    ap.add_argument("-c", "--config", default=None, help="YAML site-config (nadpisuje domyślny; przykłady: configs/examples/)")
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
    cfg = C.load(a.config, ov, require_paths=True)
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
