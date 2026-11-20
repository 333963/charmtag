"""Syntetyczne h5 w układzie v1.1 (compound `events`) — tylko testy. NIE fizyka."""
from __future__ import annotations

import os

import h5py
import numpy as np

from .constants import DM_DSTAR, M0_GEV as M0  # GeV
RANGE = {10: (1.75, 2.0), 4: (1.70, 2.20), 1: (1.89, 2.2)}


def _cands(rng, n, dt, kind, run, ev0, nev):
    lo, hi = RANGE[dt]
    nd = 2 if dt == 10 else 3
    pt = 3 + rng.exponential(2.0 if kind == "sig" else 1.6, n)
    eta = rng.uniform(-2.2, 2.2, n)
    if kind == "sig":
        m = rng.normal(M0[dt], 0.018 + 0.001 * pt, n)
        Lxy = rng.exponential(1.2, n) + 0.3
        imp = rng.normal(0, 0.02, n)
        cts = rng.uniform(-1, 1, n)
    elif kind == "refl":
        m = rng.normal(M0[dt], 0.07, n)
        Lxy = rng.exponential(1.2, n) + 0.3
        imp = rng.normal(0, 0.02, n)
        cts = rng.uniform(-1, 1, n)
    else:
        m = lo + rng.exponential(0.35, n)
        Lxy = rng.exponential(0.6, n) + 0.25
        imp = rng.normal(0, 0.06, n)
        cts = np.clip(rng.normal(0, 0.6, n), -1, 1)
    keep = (m >= lo) & (m < hi)
    d = {
        "run_number": np.full(n, run, np.int32), "event_number": ev0 + rng.integers(0, nev, n),
        "event_weight": np.ones(n, np.float32), "Dmeson_decayType": np.full(n, dt, np.int32),
        "Dmeson_mass": m, "Dmeson_pt": pt, "Dmeson_eta": eta, "Dmeson_phi": rng.uniform(-np.pi, np.pi, n),
        "Dmeson_charge": rng.choice([-1, 1], n).astype(np.int32), "Dmeson_delta_mass": np.full(n, -0.999),
        "Dmeson_Lxy": Lxy, "Dmeson_chi2": rng.chisquare(3 if kind != "bkg" else 5, n),
        "Dmeson_impact": imp, "Dmeson_impactSignificance": imp / 0.02, "Dmeson_z0SinTheta": rng.normal(0, 0.05, n),
        "Dmeson_cosThetaStar": cts, "Dmeson_cosThetaP": 1 - rng.exponential(0.001 if kind != "bkg" else 0.01, n),
        "Dmeson_ptcone40": rng.exponential(1.0, n),
        "Dmeson_truth_category": np.full(n, {"sig": 0, "refl": 0, "bkg": 12}[kind], np.int32),
        "Dmeson_mKpi1": rng.uniform(700, 1600, n).astype(np.float32), "Dmeson_mKpi2": rng.uniform(700, 1600, n).astype(np.float32),
    }
    for i in range(3):
        on = i < nd
        f = rng.uniform(0.2, 0.6, n)
        dpt = pt * f * 1000 if on else np.full(n, np.nan)  # MeV jak v1.1
        d[f"Dmeson_daughter{i}_pt"] = dpt
        d[f"Dmeson_daughter{i}_eta"] = eta + rng.normal(0, 0.15 if kind != "bkg" else 0.4, n) if on else np.full(n, np.nan)
        d[f"Dmeson_daughter{i}_phi"] = d["Dmeson_phi"] + rng.normal(0, 0.15 if kind != "bkg" else 0.4, n)
        d[f"Dmeson_daughter{i}_d0"] = rng.normal(0, 0.08 if kind != "bkg" else 0.04, n)
        d[f"Dmeson_daughter{i}_z0SinThetaPV"] = rng.normal(0, 0.1, n)
        d[f"Dmeson_daughter{i}_npts_dEdx"] = rng.integers(0, 5, n).astype(np.int32)
        is_k = i == 0
        bb = 1.2 + (0.6 if (is_k and kind != "bkg") else 0) / np.maximum(dpt * 1.2e-3, 0.3) if on else np.full(n, np.nan)
        d[f"Dmeson_daughter{i}_dEdx"] = np.where(rng.random(n) < 0.6, bb * rng.normal(1, 0.12, n), -1)
    return {k: (v[keep] if np.ndim(v) else v) for k, v in d.items()}


def _swaps(rng, d0, frac):
    """Druga hipoteza K↔π dla części D0 sygnału: te same tory (zamienione), masa rozmyta, ta sama truth."""
    sel = rng.random(len(d0["Dmeson_mass"])) < frac
    out = {k: v[sel].copy() for k, v in d0.items()}
    for q in ("pt", "eta", "phi", "d0", "z0SinThetaPV", "dEdx", "npts_dEdx"):
        a, b = out[f"Dmeson_daughter0_{q}"].copy(), out[f"Dmeson_daughter1_{q}"].copy()
        out[f"Dmeson_daughter0_{q}"], out[f"Dmeson_daughter1_{q}"] = b, a
    out["Dmeson_mass"] = np.clip(M0[10] + rng.normal(0, 0.07, sel.sum()), 1.7501, 1.9999)
    return out


def _shared_track_clones(rng, c, frac, dt_out=None):
    """Drugi kandydat w tym samym zdarzeniu dzielący córkę 0 z kandydatem `c` (test weta wspólnych śladów, PF):
    d0/z0 córki 0 identyczne bitowo, pT/η/φ córki 0 „po refit" (±0,1 %), pozostałe córki i masa losowe (tło kombinatoryczne)."""
    sel = rng.random(len(c["Dmeson_mass"])) < frac
    n = int(sel.sum())
    if n == 0:
        return None
    dt = int(dt_out or c["Dmeson_decayType"][0])
    out = _cands(rng, n, dt, "sig", int(c["run_number"][0]), 0, 1)     # wierzchołek „sygnałopodobny" (przechodzi WP) …
    lo, hi = RANGE[dt]
    out["Dmeson_mass"] = lo + rng.exponential(0.35, len(out["Dmeson_mass"]))   # … ale masa kombinatoryczna (tło)
    out["Dmeson_truth_category"][:] = 12
    keep = out["Dmeson_mass"] < hi
    out = {q: v[keep] for q, v in out.items()}
    k = min(n, len(out["Dmeson_mass"]))
    src = {q: c[q][sel][:k] for q in c}
    out = {q: v[:k] for q, v in out.items()}
    out["event_number"] = src["event_number"].copy()
    for q in ("d0", "z0SinThetaPV", "dEdx", "npts_dEdx"):
        out[f"Dmeson_daughter0_{q}"] = src[f"Dmeson_daughter0_{q}"].copy()
    for q in ("pt", "eta", "phi"):
        out[f"Dmeson_daughter0_{q}"] = src[f"Dmeson_daughter0_{q}"] * rng.normal(1, 1e-3, k)
    return out


def _cat(parts):
    keys = parts[0].keys()
    return {k: np.concatenate([p[k] for p in parts]) for k in keys}


def _dstar(rng, d0):
    """D* dla części D0 sygnału: ta sama (run,event), m(D*) = m(D0)+Δm."""
    n = len(d0["Dmeson_mass"])
    out = {k: v.copy() for k, v in d0.items()}
    dm = rng.normal(DM_DSTAR, 0.0006, n)
    out["Dmeson_decayType"][:] = 1
    out["Dmeson_delta_mass"] = dm
    out["Dmeson_mass"] = d0["Dmeson_mass"] + dm
    out["Dmeson_pt"] = d0["Dmeson_pt"] / rng.uniform(0.88, 0.97, n)
    out["Dmeson_daughter2_pt"] = out["Dmeson_pt"] * 60.0
    return out


def _dstar_bkg(rng, n, run, nev):
    b = _cands(rng, n, 10, "bkg", run, 0, nev)
    b["Dmeson_decayType"][:] = 1
    b["Dmeson_delta_mass"] = 0.1396 + rng.exponential(0.02, len(b["Dmeson_mass"]))
    b["Dmeson_mass"] = b["Dmeson_mass"] + b["Dmeson_delta_mass"]
    return b


def _group_by_event(o, run, event):
    """Porządek wierszy `o` (losowy) → ten sam zbiór wierszy, ale kandydaci jednego zdarzenia (run, event) w ciągłym bloku,
    jak w prawdziwych h5. Zdarzenia zostają w losowej kolejności (wg pierwszego wystąpienia), kandydaci w zdarzeniu też.
    Bez tego `pairs.build.check_event_key_unique` (strażnik z 16.11) słusznie odrzucał dane syntetyczne."""
    key = (run[o].astype(np.int64) << 40) | event[o].astype(np.int64)
    _, first, inv = np.unique(key, return_index=True, return_inverse=True)
    return o[np.argsort(first[inv], kind="stable")]


def write(path, d, rng):
    o = rng.permutation(len(d["Dmeson_mass"]))
    o = _group_by_event(o, d["run_number"], d["event_number"])
    d = dict(d)
    d["Dmeson_truth_channel"] = np.array([b"x"] * len(o), dtype=object)  # kolumna vlen str jak v1.1
    dt = np.dtype([(k, v.dtype if v.dtype != object else h5py.string_dtype()) for k, v in d.items()])
    arr = np.empty(len(o), dt)
    for k, v in d.items():
        arr[k] = v[o]
    with h5py.File(path, "w") as f:
        f.create_dataset("events", data=arr, compression="gzip", compression_opts=1)


def make(outdir, n=20000, seed=1, nev_factor=2, shared_frac=0.15):
    """Syntetyczne h5 (układ v1.1). `nev_factor`: liczba zdarzeń = nev_factor·n (małe → wiele kandydatów na zdarzenie →
    pary same-event, także ze wspólnym śladem: `shared_frac` klonów na kandydata tła) — ścieżka par/weta jest ćwiczona."""
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    nev = max(int(nev_factor * n), 1000)
    # MC forcedDecay: run 0, event_number nieunikalne (jak v1.1)
    parts = []
    for dt in (10, 4):
        parts += [_cands(rng, n, dt, "sig", 0, 0, nev), _cands(rng, n // 2, dt, "bkg", 0, 0, nev)]
    d0s = parts[0]
    parts.append(_swaps(rng, d0s, 0.5))
    sel = rng.random(len(d0s["Dmeson_mass"])) < 0.3
    parts.append(_dstar(rng, {k: v[sel] for k, v in d0s.items()}))
    write(os.path.join(outdir, "fd_nd_d150_Donly_none.h5"), _cat(parts), rng)
    # MC inclusive
    parts = []
    for dt in (10, 4):
        sg, bk = _cands(rng, n // 10, dt, "sig", 0, 0, nev), _cands(rng, 3 * n, dt, "bkg", 0, 0, nev)
        parts += [sg, bk]
        cl = _shared_track_clones(rng, sg, shared_frac)
        if cl is not None:
            parts.append(cl)
        if dt == 10:  # D* prawdy w MC minbias → normalizacja 'direct' (N_D* ≥ 50) zamiast cichego 'frozen'
            sel = rng.random(len(sg["Dmeson_mass"])) < 0.3
            parts.append(_dstar(rng, {k: v[sel] for k, v in sg.items()}))
    write(os.path.join(outdir, "mb_nd_d150_Donly_none.h5"), _cat(parts), rng)
    # dane: tło + mały sygnał (truth = -1 wszędzie)
    for name, scale in (("data_428770_d150_Donly_none.h5", 1.0), ("data_428770_d150_Donly_single_tag.h5", 0.1)):
        parts = []
        for dt in (10, 4):
            s = _cands(rng, int(scale * n / 4), dt, "sig", 428770, 0, nev)
            bk = _cands(rng, int(scale * 4 * n), dt, "bkg", 428770, 0, nev)
            parts += [s, bk]
            cl = _shared_track_clones(rng, s, shared_frac)
            if cl is not None:
                parts.append(cl)
            if dt == 10:
                sel = rng.random(len(s["Dmeson_mass"])) < 0.3
                parts.append(_dstar(rng, {k: v[sel] for k, v in s.items()}))
                parts.append(_dstar_bkg(rng, int(scale * n / 2), 428770, nev))
        for p in parts:
            p["Dmeson_truth_category"][:] = -1
            p["Dmeson_ptcone40"] = p["Dmeson_ptcone40"] * 1.5
        write(os.path.join(outdir, name), _cat(parts), rng)
    return outdir
