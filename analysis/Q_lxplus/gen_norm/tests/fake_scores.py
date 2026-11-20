#!/usr/bin/env python3
"""Test lokalny gen_norm: config + atrapy scores incl_mc (charmtag `cmd_apply` z src/ i losowym score zamiast modelu XGB).
Selekcja (in_range, mark_swaps, veto_ok, row globalny, attrs source/file_offsets) = prawdziwy kod apply.
  PYTHONPATH=src python3 tests/fake_scores.py --dir DIR   (DIR z make_fake.py: DIR/h5/mb_nd/*.h5) → DIR/cfg.yaml, DIR/run/scores/..."""
import argparse
import os

import numpy as np
import yaml

from charmtag import config as C
from charmtag import model as M
from charmtag import pipeline as PL


class FakeEnsemble:
    def __init__(self, *a, **k):
        self.rng = np.random.default_rng(1)

    def score(self, L):
        return self.rng.random(len(L)).astype(np.float32), np.zeros(len(L), np.int8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = os.path.abspath(a.dir)
    cfg_y = {"paths": {"h5_dir": os.path.join(d, "h5"), "out_dir": os.path.join(d, "run"),
                       "samples": {"incl_mc": "mb_nd/mb_nd_*_none.h5"}},
             "io": {"chunk_rows": 100000}}
    p = os.path.join(d, "cfg.yaml")
    yaml.safe_dump(cfg_y, open(p, "w"))
    cfg = C.load(p)
    M.Ensemble = FakeEnsemble
    PL.M.Ensemble = FakeEnsemble
    for sp in ("D0", "Dplus"):
        print(sp, PL.cmd_apply(cfg, sp, "H", ["incl_mc"]))
    print("OK", p)


if __name__ == "__main__":
    main()
