#!/usr/bin/env python3
"""Scal YAML: baza + nakładki (kolejno) + key=value (ścieżki z kropkami) → plik. Użycie:
merge_cfg.py OUT.yaml BASE.yaml [OVERLAY.yaml ...] [-s paths.out_dir=/eos/... -s training.bkg_source=mc]"""
import argparse, copy, os, yaml

def merge(a, b):
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("base"); ap.add_argument("overlays", nargs="*")
ap.add_argument("-s", "--set", action="append", default=[])
a = ap.parse_args()
cfg = yaml.safe_load(open(a.base)) or {}
bdir = os.path.dirname(os.path.abspath(a.base))
for k in ("h5_dir", "out_dir"):  # względne ścieżki bazy → absolutne (plik wynikowy leży gdzie indziej)
    p = (cfg.get("paths") or {}).get(k)
    if p and not os.path.isabs(os.path.expanduser(p)):
        cfg["paths"][k] = os.path.normpath(os.path.join(bdir, p))
for o in a.overlays:
    cfg = merge(cfg, yaml.safe_load(open(o)) or {})
for kv in a.set:
    k, v = kv.split("=", 1)
    d = cfg
    ks = k.split(".")
    for x in ks[:-1]:
        d = d.setdefault(x, {})
    d[ks[-1]] = yaml.safe_load(v)
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
yaml.safe_dump(cfg, open(a.out, "w"), sort_keys=False, allow_unicode=True)
print(a.out)
