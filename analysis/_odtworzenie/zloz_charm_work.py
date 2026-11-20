#!/usr/bin/env python3
"""Składa z repo drzewo w układzie `charm_work`, w którym skrypty z `analysis/` działają bez zmian.

`analysis/<ścieżka>` (kod) i `results/<ścieżka>` (wyniki) to dwie połowy tego samego drzewa
`~/Documents/charm/charm_work/<ścieżka>`. Skrypty szukają wejść ścieżkami względnymi w tym drzewie, więc żeby je
uruchomić, trzeba obie połowy nałożyć na siebie w katalogu roboczym:

  python analysis/_odtworzenie/zloz_charm_work.py /tmp/cw
  →  /tmp/cw/charm_work/...   (kopie; repo nietknięte)
     /tmp/cw/scratch/p427/    (wejścia fitu 427929 w układzie, którego oczekuje truthnorm_run.py: TN_SP=/tmp/cw/scratch)

Pliki `<nazwa>.npz.json` są przy okazji zamieniane z powrotem na `<nazwa>.npz` (zob. npzjson.py).
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import npzjson  # noqa: E402

UNB = os.path.join("PF2_lxplus", "runs", "lxplus", "unblind")
SKIP_TOP = {"README.md", ".gitignore", "_odtworzenie"}


def assemble(dst):
    cw = os.path.join(dst, "charm_work")
    for half in ("analysis", "results"):
        src = os.path.join(REPO, half)
        for name in sorted(os.listdir(src)):
            if name in SKIP_TOP:
                continue
            p = os.path.join(src, name)
            if os.path.isdir(p):
                shutil.copytree(p, os.path.join(cw, name), dirs_exist_ok=True, ignore=shutil.ignore_patterns("__pycache__"))
            else:
                os.makedirs(cw, exist_ok=True)
                shutil.copy2(p, os.path.join(cw, name))
    n = 0
    for root, _, files in os.walk(cw):
        for f in files:
            if f.endswith(".npz.json"):
                npzjson.to_npz(os.path.join(root, f))
                n += 1
    # wejścia fitu 427929 tak, jak je czyta truthnorm_run.py (P427 = $TN_SP/p427: pair_inputs.{json,npz} + perfile.npz z kluczami H_<kanał>)
    import numpy as np
    p427 = os.path.join(dst, "scratch", "p427")
    os.makedirs(p427, exist_ok=True)
    pairs = os.path.join(cw, UNB, "r427929", "pairs")
    for f in ("pair_inputs.json", "pair_inputs.npz"):
        shutil.copy2(os.path.join(pairs, f), os.path.join(p427, f))
    with np.load(os.path.join(pairs, "pairs_full.npz")) as H:
        np.savez(os.path.join(p427, "perfile.npz"), **{"H_" + k: H[k][None] for k in H.files})   # jedna „grupa” = cała próbka; skrypt sumuje po osi 0
    return cw, n


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    cw, n = assemble(os.path.abspath(sys.argv[1]))
    print(f"złożono {cw} (odtworzono {n} plików .npz); TN_SP={os.path.join(os.path.abspath(sys.argv[1]), 'scratch')}")
