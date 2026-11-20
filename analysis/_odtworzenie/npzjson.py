#!/usr/bin/env python3
"""npz <-> JSON (tekst) dla małych wejść fitu.

Repo nie przechowuje plików binarnych (.npz). Małe tablice, od których zależy wynik (histogramy par, szablony,
histogramy Δm), leżą w `results/` jako `<nazwa>.npz.json`; ten moduł odtwarza z nich `.npz` o identycznej
zawartości tablic (te same klucze, dtype, kształt i wartości bit w bit — liczby zapisane przez `repr`).

  python npzjson.py to-json PLIK.npz [PLIK.npz.json]
  python npzjson.py to-npz  PLIK.npz.json [PLIK.npz]
"""
import hashlib
import json
import os
import sys

import numpy as np

FORMAT = "npz-as-json/1"


def to_json(npz_path, json_path=None, source_label=None):
    json_path = json_path or npz_path + ".json"
    with np.load(npz_path, allow_pickle=False) as A:
        arrays = {k: {"dtype": str(A[k].dtype), "shape": list(A[k].shape), "data": A[k].ravel().tolist()} for k in A.files}
    md5 = hashlib.md5(open(npz_path, "rb").read()).hexdigest()
    head = {"format": FORMAT, "source": source_label or os.path.basename(npz_path), "md5_source_npz": md5, "n_arrays": len(arrays)}
    with open(json_path, "w") as f:                      # jedna tablica w wierszu: czytelny diff, mały plik
        f.write("{\n")
        for k, v in head.items():
            f.write(f" {json.dumps(k)}: {json.dumps(v)},\n")
        f.write(' "arrays": {\n')
        items = list(arrays.items())
        for i, (k, v) in enumerate(items):
            f.write(f"  {json.dumps(k)}: {json.dumps(v, separators=(',', ':'))}{',' if i < len(items) - 1 else ''}\n")
        f.write(" }\n}\n")
    return json_path


def load(json_path):
    J = json.load(open(json_path))
    if J.get("format") != FORMAT:
        raise ValueError(f"{json_path}: nieznany format {J.get('format')!r}")
    return {k: np.array(v["data"], dtype=v["dtype"]).reshape(v["shape"]) for k, v in J["arrays"].items()}


def to_npz(json_path, npz_path=None):
    npz_path = npz_path or (json_path[:-5] if json_path.endswith(".json") else json_path + ".npz")
    np.savez(npz_path, **load(json_path))
    return npz_path


def same(npz_path, json_path):
    with np.load(npz_path, allow_pickle=False) as A:
        B = load(json_path)
        return sorted(A.files) == sorted(B) and all(A[k].dtype == B[k].dtype and A[k].shape == B[k].shape and np.array_equal(A[k], B[k]) for k in A.files)


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("to-json", "to-npz"):
        raise SystemExit(__doc__)
    fn = to_json if sys.argv[1] == "to-json" else to_npz
    print(fn(*sys.argv[2:4]))
