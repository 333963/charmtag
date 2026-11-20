from __future__ import annotations

import copy
import glob
import json
import os
import re

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "data", "default.yaml")  # domyślny config w pakiecie (bez ścieżek lokalnych)


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load(path: str | None = None, overrides: dict | None = None, require_paths: bool = False) -> dict:
    """Domyślny config pakietu + site-config `path` + `overrides`. `require_paths=True` (CLI): paths.h5_dir/out_dir muszą być ustawione."""
    with open(DEFAULT) as f:
        cfg = yaml.safe_load(f)
    if path and os.path.abspath(path) != DEFAULT:
        with open(path) as f:
            cfg = _merge(cfg, yaml.safe_load(f))
    if overrides:
        cfg = _merge(cfg, overrides)
    base = os.path.dirname(os.path.abspath(path or DEFAULT))
    for k in ("h5_dir", "out_dir"):
        if not cfg["paths"].get(k):
            if require_paths:
                raise SystemExit(f"config: brak `paths.{k}` — ustaw w pliku -c (patrz configs/examples/) albo przez --{k.replace('_', '-')}")
            continue
        p = os.path.expanduser(str(cfg["paths"][k]))
        cfg["paths"][k] = p if os.path.isabs(p) else os.path.normpath(os.path.join(base, p))
    for k in ("mb_rates_json", "freeze_dir", "sf_json", "truth_norm_json"):   # pliki wejściowe par: też względem configu
        v = (cfg.get("pairs") or {}).get(k)
        if v:
            v = os.path.expanduser(str(v))
            cfg["pairs"][k] = v if os.path.isabs(v) else os.path.normpath(os.path.join(base, v))
    return cfg


def sample_path(cfg: dict, key: str) -> str:
    s = cfg["paths"]["samples"][key]
    if isinstance(s, (list, tuple)) or _is_glob(s):  # etap P: wiele plików → pierwszy (kolumny/istnienie)
        ps = sample_paths(cfg, key)
        return ps[0] if ps else _abs(cfg, s[0] if isinstance(s, (list, tuple)) and s else str(s))
    s = os.path.expanduser(s)
    return s if os.path.isabs(s) else os.path.join(cfg["paths"]["h5_dir"], s)


# ----------------------------------------------------------------------------- etap P: wiele plików
def _is_glob(s) -> bool:
    return isinstance(s, str) and any(ch in s for ch in "*?[")


def _abs(cfg, s: str) -> str:
    s = os.path.expanduser(str(s))
    return s if os.path.isabs(s) else os.path.join(cfg["paths"]["h5_dir"], s)


def sample_paths(cfg: dict, key: str) -> list[str]:
    """Lista plików próbki. W YAML: nazwa pliku, glob albo lista (nazw/globów), np.
    data_none: [data_428770_d150_Donly_none.h5, ../parts_data_1_24/data_428770_*_none.h5].
    Glob: posortowane, bez *.part; duplikaty usuwane (kolejność zachowana)."""
    s = cfg["paths"]["samples"][key]
    items = list(s) if isinstance(s, (list, tuple)) else [s]
    out = []
    for it in items:
        q = os.path.normpath(_abs(cfg, it))
        out += [g for g in sorted(glob.glob(q)) if not g.endswith(".part")] if _is_glob(q) else [q]
    return list(dict.fromkeys(out))


def is_multi(cfg: dict, key: str) -> bool:
    s = cfg["paths"]["samples"][key]
    return isinstance(s, (list, tuple)) or _is_glob(s)


def source(cfg: dict, key: str):
    """Ścieżka (1 plik, zachowanie jak dawniej) albo lista ścieżek (wiele plików)."""
    return sample_paths(cfg, key) if is_multi(cfg, key) else sample_path(cfg, key)


def _jobs(spec) -> list[int]:
    out = []
    for part in str(spec or "").replace(" ", "").split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def sample_meta(cfg: dict, key: str) -> dict:
    """Metadane próbki: pliki, wiersze, joby (z nazwy *_0000NN_* albo paths.sample_info.<key>.jobs), lumi.
    lumi_fraction_jobs = n_jobs/jobs_total — PRZYBLIŻENIE (joby nie mają równej lumi);
    lumi_pb tylko gdy podane w paths.sample_info.<key>.lumi_pb."""
    import h5py
    info = ((cfg["paths"].get("sample_info") or {}).get(key) or {})
    ds = cfg.get("io", {}).get("dataset", "events")
    files, jobs, nev, nev_known = [], set(), 0, True
    for p in sample_paths(cfg, key):
        e = {"file": os.path.basename(p), "exists": os.path.exists(p)}
        m = re.findall(r"_(\d{6})(?=_)", os.path.basename(p))  # ostatnia grupa 6 cyfr (nie run 428770)
        m = [x for x in m if x.startswith("0")]
        if m:
            e["jobs"] = [int(m[-1])]
        else:
            e["jobs"] = _jobs((info.get("jobs") or {}).get(os.path.basename(p), ""))
        jobs |= set(e["jobs"])
        if e["exists"]:
            with h5py.File(p, "r") as f:
                e["rows"] = int(f[ds].shape[0]) if ds in f else None
                a = f.attrs.get("n_events_read")
            if a is not None:
                e["n_events_read"] = int(a)
                nev += int(a)
            else:
                nev_known = False
        files.append(e)
    tot = info.get("jobs_total")
    return {"key": key, "n_files": len(files), "rows": int(sum(e.get("rows") or 0 for e in files)),
            "jobs": sorted(jobs), "n_jobs": len(jobs), "jobs_total": tot,
            "lumi_fraction_jobs": (len(jobs) / tot) if tot else None, "lumi_pb": info.get("lumi_pb"),
            "n_events_read": nev if nev_known else None, "files": files}


def sample_meta_json(cfg: dict, key: str) -> str:
    try:
        return json.dumps(sample_meta(cfg, key), ensure_ascii=False)
    except Exception as ex:  # metadane nie mogą zatrzymać biegu
        return json.dumps({"key": key, "error": str(ex)})


def dump(cfg: dict, path: str):
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True)
